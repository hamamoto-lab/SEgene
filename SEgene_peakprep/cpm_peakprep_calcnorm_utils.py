# === cpm_peakprep_calcnorm_utils.py ===
# --- 標準ライブラリ ---
import os
import sys
import logging
import re
import json
import time
from typing import Dict, Optional, List, Tuple, Union, Any

# --- サードパーティライブラリ ---
import numpy as np
import pandas as pd

# --- 定数定義 ---
FEATURE_COUNTS_METADATA_COLS = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
FEATURE_COUNTS_METADATA_COL_COUNT = 6  # featureCounts出力のメタデータ列数
TSV_FORMAT_COLS = ['PeakID', 'Chr', 'Start', 'End']  # 標準TSV出力のメタデータカラム

# --- rpy2のインポート（条件付き） ---
try:
    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    # pandas-rpy2変換を有効化
    pandas2ri.activate()
    RPY2_AVAILABLE = True
except ImportError as e:
    RPY2_AVAILABLE = False
    print(f"Error importing rpy2: {e}", file=sys.stderr)
    print("Please install rpy2: pip install rpy2", file=sys.stderr)

# --- ロギング関連 ---
def setup_logging(log_level_str: str, script_log_file: Optional[str], log_dir: str) -> logging.Logger:
    """
    Sets up logging configuration.
    
    Args:
        log_level_str (str): Log level as a string
        script_log_file (Optional[str]): Path to log file, None to disable file logging
        log_dir (str): Directory for log file
        
    Returns:
        logging.Logger: Configured logger
    """
    # ログレベルの決定
    log_level = getattr(logging, log_level_str.upper(), logging.INFO)
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # ルートロガーのリセット
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()
    
    root_logger.setLevel(log_level)
    
    # コンソールハンドラの追加
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    formatter = logging.Formatter(log_format, datefmt=date_format)
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)
    
    # 特定のロガーを取得
    logger = logging.getLogger("NormFactors")
    
    # 要求された場合はファイルハンドラを追加
    if script_log_file:
        log_file_path = os.path.join(log_dir, script_log_file)
        
        try:
            # ログディレクトリが存在しない場合は作成
            os.makedirs(log_dir, exist_ok=True)
            logger.debug(f"Created or verified log directory: {log_dir}")
            
            # 既存のログファイルをバックアップ
            if os.path.exists(log_file_path):
                timestamp = time.strftime("%Y%m%d-%H%M%S")
                backup_log_path = f"{log_file_path}.{timestamp}.bak"
                try:
                    os.rename(log_file_path, backup_log_path)
                    logger.debug(f"Existing log file renamed to: {os.path.basename(backup_log_path)}")
                except Exception as e:
                    logger.warning(f"Failed to rename existing log file {log_file_path}: {e}")
            
            fh = logging.FileHandler(log_file_path, mode='w')
            fh.setLevel(log_level)
            fh.setFormatter(formatter)
            root_logger.addHandler(fh)
            logger.debug(f"Logging to file: {log_file_path}")
        except Exception as e:
            logger.error(f"Failed to configure file logging to {log_file_path}: {e}")
    
    return logger

# --- ディレクトリ管理 ---
def create_output_directories(
    output_file: str, 
    cleaned_output: Optional[str], 
    log_file: Optional[str], 
    log_dir: str, 
    full_metadata_output: Optional[str],
    logger: logging.Logger
) -> bool:
    """
    Creates all necessary output directories.
    
    Args:
        output_file (str): Path to output file
        cleaned_output (Optional[str]): Path to cleaned counts file
        log_file (Optional[str]): Path to log file
        log_dir (str): Directory for log file
        full_metadata_output (Optional[str]): Path to complete featureCounts data output file
        logger (logging.Logger): Logger object
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # 出力ディレクトリを作成
        output_dir = os.path.dirname(output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created or verified output directory: {output_dir}")
        
        # クリーニング出力ディレクトリを作成（必要な場合）
        if cleaned_output:
            cleaned_dir = os.path.dirname(cleaned_output)
            if cleaned_dir:
                os.makedirs(cleaned_dir, exist_ok=True)
                logger.debug(f"Created or verified cleaned output directory: {cleaned_dir}")
        
        # ログディレクトリを作成（必要な場合）
        if log_file and log_dir:
            os.makedirs(log_dir, exist_ok=True)
            logger.debug(f"Created or verified log directory: {log_dir}")
        
        # フルメタデータ出力ディレクトリを作成（必要な場合）
        if full_metadata_output:
            full_data_dir = os.path.dirname(full_metadata_output)
            if full_data_dir:
                os.makedirs(full_data_dir, exist_ok=True)
                logger.debug(f"Created or verified full metadata output directory: {full_data_dir}")
            
        return True
        
    except OSError as e:
        logger.critical(f"Failed to create directory: {e}")
        return False

def get_temp_cleaned_file_path(output_file: str, logger: logging.Logger) -> str:
    """
    Generates a path for temporary cleaned file.
    
    Args:
        output_file (str): Path to output file
        logger (logging.Logger): Logger object
        
    Returns:
        str: Path to temporary cleaned file
    """
    output_dir = os.path.dirname(output_file) or '.'
    temp_suffix = time.strftime("%Y%m%d-%H%M%S")
    cleaned_file_path = os.path.join(output_dir, f"cleaned_counts_{temp_suffix}.txt")
    logger.info(f"No cleaned output file specified, using temporary file: {cleaned_file_path}")
    
    return cleaned_file_path

def cleanup_temp_files(temp_file_path: str, is_temp: bool, logger: logging.Logger) -> None:
    """
    Cleans up temporary files.
    
    Args:
        temp_file_path (str): Path to temporary file
        is_temp (bool): Whether the file is temporary
        logger (logging.Logger): Logger object
    """
    if is_temp and os.path.exists(temp_file_path):
        try:
            os.remove(temp_file_path)
            logger.debug(f"Temporary cleaned file removed: {temp_file_path}")
        except Exception as e:
            logger.warning(f"Failed to remove temporary file {temp_file_path}: {e}")

# --- パターンルール処理 ---
def load_pattern_rules(pattern_rules_file: Optional[str], logger: logging.Logger, remove_extensions: bool = False) -> List[Dict[str, str]]:
    """
    Loads pattern rules from JSON file.
    
    Args:
        pattern_rules_file (Optional[str]): Path to JSON file with pattern rules
        logger (logging.Logger): Logger object
        remove_extensions (bool): Flag to add patterns for removing common BAM extensions
        
    Returns:
        List[Dict[str, str]]: List of pattern rules
    """
    rules = []
    
    # BAM拡張子削除パターンの追加
    if remove_extensions:
        rules.append({"pattern": r"\.bam$", "replace": ""})
        rules.append({"pattern": r"\.sorted\.bam$", "replace": ""})
        rules.append({"pattern": r"\.mLb\.clN\.sorted\.bam$", "replace": ""})
        logger.info("Added common BAM extension removal patterns")
    
    # JSONからカスタムパターンを読み込む
    if pattern_rules_file:
        try:
            with open(pattern_rules_file, 'r') as f:
                custom_rules = json.load(f)
                
            if isinstance(custom_rules, list):
                for rule in custom_rules:
                    if isinstance(rule, dict) and "pattern" in rule and "replace" in rule:
                        rules.append(rule)
                    else:
                        logger.warning(f"Skipping invalid rule: {rule}")
                
                logger.info(f"Loaded {len(custom_rules)} custom pattern rules from {pattern_rules_file}")
            else:
                logger.error(f"Invalid pattern rules format in {pattern_rules_file}. Expected a list of rules.")
                
        except Exception as e:
            logger.error(f"Error loading pattern rules from {pattern_rules_file}: {e}")
    
    return rules

# --- R環境のセットアップ ---
def setup_r_environment(logger: logging.Logger) -> Tuple[Optional[Any], Optional[Any], Optional[Any]]:
    """
    Sets up R environment and required packages.
    
    Args:
        logger (logging.Logger): Logger object
        
    Returns:
        Tuple[Optional[Any], Optional[Any], Optional[Any]]: 
            R environment, base package, edgeR package. None if setup fails.
    """
    if not RPY2_AVAILABLE:
        logger.error("rpy2 is not available. Cannot set up R environment.")
        return None, None, None
        
    try:
        # R環境の初期化
        r = robjects.r
        base = importr('base')
        
        # edgeRパッケージのインポート
        try:
            edger = importr('edgeR')
            return r, base, edger
        except Exception as e:
            logger.error(f"Could not import edgeR package: {e}")
            logger.error("Please install edgeR in R: BiocManager::install('edgeR')")
            return None, None, None
            
    except Exception as e:
        logger.error(f"Error setting up R environment: {e}", exc_info=True)
        return None, None, None

# --- ファイル検証 ---
def validate_counts_file(counts_file: str, logger: logging.Logger) -> bool:
    """
    Validates the featureCounts output file.
    
    Args:
        counts_file (str): Path to the featureCounts output file
        logger (logging.Logger): Logger object
        
    Returns:
        bool: True if file is valid, False otherwise
    """
    try:
        # ファイルの存在確認
        if not os.path.exists(counts_file):
            logger.error(f"Counts file not found: {counts_file}")
            return False
        
        # ファイルの読み込みとバリデーション
        df = pd.read_csv(counts_file, sep='\t', comment='#', nrows=5)
        
        # 必要なカラムの確認
        missing_cols = [col for col in FEATURE_COUNTS_METADATA_COLS if col not in df.columns]
        
        if missing_cols:
            logger.error(f"Required columns missing: {missing_cols}")
            return False
        
        # カウント列が存在するか確認
        if len(df.columns) <= FEATURE_COUNTS_METADATA_COL_COUNT:
            logger.error("No count columns found in file")
            return False
        
        logger.info(f"Counts file validation passed: {counts_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error validating counts file: {e}")
        return False

# --- ファイル読み書き ---
def read_counts_file_with_comments(counts_file: str, logger: logging.Logger) -> Tuple[List[str], pd.DataFrame]:
    """
    Reads featureCounts file with comments preserved.
    
    Args:
        counts_file (str): Path to featureCounts file
        logger (logging.Logger): Logger object
        
    Returns:
        Tuple[List[str], pd.DataFrame]: Comment lines and data DataFrame
    """
    try:
        # コメント行を読み込む
        header_comments = []
        with open(counts_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    header_comments.append(line)
                else:
                    break
        
        # データをDataFrameとして読み込む
        counts_df = pd.read_csv(counts_file, sep='\t', comment='#')
        
        return header_comments, counts_df
        
    except Exception as e:
        logger.error(f"Error reading counts file: {e}", exc_info=True)
        raise

def write_counts_file_with_comments(
    output_file: str,
    header_comments: List[str],
    df: pd.DataFrame,
    logger: logging.Logger
) -> bool:
    """
    Writes featureCounts file with comments preserved.
    
    Args:
        output_file (str): Path to output file
        header_comments (List[str]): Comment lines to include
        df (pd.DataFrame): Data to write
        logger (logging.Logger): Logger object
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # 出力ディレクトリの作成（必要な場合）
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            logger.debug(f"Created output directory: {output_dir}")
        
        # コメント付きでファイルを書き込む
        with open(output_file, 'w') as f:
            # コメント行を書き込む
            for comment in header_comments:
                f.write(comment)
            
            # DataFrameをTSVとして書き込む
            df.to_csv(f, sep='\t', index=False)
            
        return True
        
    except Exception as e:
        logger.error(f"Error writing counts file: {e}", exc_info=True)
        return False

# --- サンプル名クリーニング ---
def clean_sample_names(
    df: pd.DataFrame, 
    logger: logging.Logger, 
    sample_delimiter: Optional[str] = None,
    pattern_rules: Optional[List[Dict[str, str]]] = None
) -> pd.DataFrame:
    """
    Cleans sample names (column names) in the DataFrame.
    
    Args:
        df (pd.DataFrame): Original DataFrame
        logger (logging.Logger): Logger object
        sample_delimiter (Optional[str]): Delimiter to extract sample names. If provided, 
                                         the part before this delimiter will be used as the sample name.
        pattern_rules (Optional[List[Dict[str, str]]]): Additional pattern rules for sample name cleaning.
            Example: [{"pattern": r"\.bam$", "replace": ""}]
        
    Returns:
        pd.DataFrame: DataFrame with cleaned column names
    """
    logger.info("Cleaning sample names...")
    
    # 元のカラム名を保存
    original_columns = df.columns.tolist()
    
    # メタデータカラムとカウントカラムを分離
    metadata_cols = original_columns[:FEATURE_COUNTS_METADATA_COL_COUNT]
    count_cols = original_columns[FEATURE_COUNTS_METADATA_COL_COUNT:]
    
    # カウントカラムの名前をクリーニング
    clean_count_cols = []
    for col in count_cols:
        # 基本クリーニング：最後のスラッシュまでを削る
        clean_col = os.path.basename(col)
        
        # デリミターが指定されている場合、それに基づいてサンプル名を抽出
        if sample_delimiter and sample_delimiter in clean_col:
            try:
                # デリミターで区切られた最初の部分をサンプル名として使用
                clean_col = clean_col.split(sample_delimiter)[0]
                logger.debug(f"Extracted sample name '{clean_col}' using delimiter '{sample_delimiter}'")
            except Exception as e:
                logger.warning(f"Failed to extract sample name using delimiter '{sample_delimiter}' from '{clean_col}': {e}")
                # 例外発生時は元のファイル名（ディレクトリパスなし）を使用
        
        # 追加のパターンルールを適用
        if pattern_rules:
            for rule in pattern_rules:
                if "pattern" in rule and "replace" in rule:
                    clean_col = re.sub(rule["pattern"], rule["replace"], clean_col)
        
        clean_count_cols.append(clean_col)
    
    # サンプル名の重複チェック
    if len(set(clean_count_cols)) < len(clean_count_cols):
        # 重複がある場合はログに警告
        duplicates = {}
        for idx, name in enumerate(clean_count_cols):
            if name in duplicates:
                duplicates[name].append(count_cols[idx])
            else:
                duplicates[name] = [count_cols[idx]]
        
        duplicate_names = {name: paths for name, paths in duplicates.items() if len(paths) > 1}
        if duplicate_names:
            logger.warning(f"Duplicate sample names found after cleaning: {duplicate_names}")
            logger.warning("This may cause issues when processing the data.")

    # 新しいカラム名でデータフレームを更新
    df.columns = metadata_cols + clean_count_cols
    
    logger.debug(f"Sample names before cleaning: {count_cols[:3]}...")
    logger.debug(f"Sample names after cleaning: {clean_count_cols[:3]}...")
    
    return df

def clean_count_file_sample_names(
    counts_file: str, 
    output_file: str,
    logger: logging.Logger,
    sample_delimiter: Optional[str] = None,
    pattern_rules: Optional[List[Dict[str, str]]] = None
) -> bool:
    """
    Cleans sample names in featureCounts file and saves to a new file.
    
    Args:
        counts_file (str): Path to input file (featureCounts output)
        output_file (str): Path to output file
        logger (logging.Logger): Logger object
        sample_delimiter (Optional[str]): Delimiter to extract sample names. If provided, 
                                        the part before this delimiter will be used as the sample name.
        pattern_rules (Optional[List[Dict[str, str]]]): Additional pattern rules for sample name cleaning.
            Example: [{"pattern": r"\.bam$", "replace": ""}]
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        logger.info(f"Reading counts file: {counts_file}")
        
        # コメント付きでファイルを読み込む
        header_comments, counts_df = read_counts_file_with_comments(counts_file, logger)
        
        logger.info(f"Loaded counts data with {len(counts_df)} features and {len(counts_df.columns) - FEATURE_COUNTS_METADATA_COL_COUNT} samples")
        
        # サンプル名をクリーニング
        counts_df = clean_sample_names(counts_df, logger, sample_delimiter, pattern_rules)
        
        # クリーニングしたファイルを書き込む
        logger.info(f"Saving cleaned counts file to: {output_file}")
        if not write_counts_file_with_comments(output_file, header_comments, counts_df, logger):
            return False
            
        logger.info(f"Successfully saved cleaned counts file with {len(counts_df.columns) - FEATURE_COUNTS_METADATA_COL_COUNT} samples")
        return True
        
    except Exception as e:
        logger.error(f"Error cleaning counts file: {e}", exc_info=True)
        return False

# --- 正規化処理 ---
def create_dge_list(
    counts_matrix: pd.DataFrame,
    r: Any,
    edger: Any,
    logger: logging.Logger
) -> Any:
    """
    Creates an edgeR DGEList object from a count matrix.
    
    Args:
        counts_matrix (pd.DataFrame): Count data matrix
        r (Any): R environment
        edger (Any): edgeR package
        logger (logging.Logger): Logger object
        
    Returns:
        Any: edgeR DGEList object
    """
    try:
        # Rの行列に変換
        r_counts = pandas2ri.py2rpy(counts_matrix)
        
        # DGEListを作成
        dge = edger.DGEList(counts=r_counts)
        
        return dge
        
    except Exception as e:
        logger.error(f"Error creating DGEList: {e}", exc_info=True)
        raise

def filter_low_expression_features(
    dge: Any,
    min_cpm: float,
    min_samples: int,
    r: Any,
    base: Any,
    edger: Any,
    logger: logging.Logger
) -> Tuple[Any, Optional[np.ndarray]]:
    """
    Filters low expression features from a DGEList.
    
    Args:
        dge (Any): edgeR DGEList object
        min_cpm (float): Minimum CPM threshold for filtering
        min_samples (int): Minimum number of samples with CPM > threshold
        r (Any): R environment
        base (Any): base package
        edger (Any): edgeR package
        logger (logging.Logger): Logger object
        
    Returns:
        Tuple[Any, Optional[np.ndarray]]: Filtered DGEList and boolean array of kept indices
    """
    # min_samplesが0の場合はフィルタリングをスキップ
    if min_samples == 0:
        logger.info("Skipping filtering step as min_samples=0")
        return dge, None
    
    try:
        # CPMを計算
        logger.info(f"Filtering: keeping features with CPM > {min_cpm} in at least {min_samples} samples")
        r_cpm = edger.cpm(dge)
        
        # フィルタ条件を適用
        r_keep = base.rowSums(r_cpm > min_cpm) >= min_samples
        dge_filtered = dge.rx(r_keep, True)
        
        # フィルタリング結果をログに記録
        kept_indices = np.array(pandas2ri.rpy2py(r_keep), dtype=bool)
        n_kept = int(kept_indices.sum())
        n_total = len(kept_indices)
        logger.info(f"Kept {n_kept} out of {n_total} features after filtering ({n_kept/n_total:.1%})")
        
        return dge_filtered, kept_indices
        
    except Exception as e:
        logger.error(f"Error during expression filtering: {e}", exc_info=True)
        raise

def get_normalization_factors(
    dge_norm: Any,
    count_cols: List[str],
    r: Any,
    logger: logging.Logger
) -> None:
    """
    Extracts and logs normalization factors from a normalized DGEList.
    
    Args:
        dge_norm (Any): Normalized edgeR DGEList object
        count_cols (List[str]): Count data column names
        r (Any): R environment
        logger (logging.Logger): Logger object
    """
    try:
        # 正規化係数を取得
        norm_factors_r = dge_norm.rx2('samples').rx2('norm.factors')
        conversion_method = ""
        
        # オブジェクトの型をチェックして適切に処理
        if isinstance(norm_factors_r, np.ndarray):
            norm_factors = norm_factors_r
            conversion_method = "Direct use of NumPy array"
        else:
            # RPy2の変換を試みる
            try:
                norm_factors = pandas2ri.rpy2py(norm_factors_r)
                conversion_method = "Using pandas2ri.rpy2py conversion"
            except Exception:
                # Rコマンドにフォールバック
                r_command = r("function(x) as.numeric(x$samples$norm.factors)")
                norm_factors = np.array(r_command(dge_norm))
                conversion_method = "Converted via R command"
        
        # 正規化係数がサンプル数と一致するか確認
        if len(norm_factors) == len(count_cols):
            # サンプル名と正規化係数をマッピング
            norm_factors_dict = dict(zip(count_cols, norm_factors))
            
            # 統計量を計算
            min_factor = np.min(norm_factors)
            max_factor = np.max(norm_factors)
            mean_factor = np.mean(norm_factors)
            
            logger.debug(f"Normalization factors retrieval method: {conversion_method}")
            logger.debug(f"Normalization factors - Min: {min_factor:.4f}, Max: {max_factor:.4f}, Mean: {mean_factor:.4f}")
            
            # 大きな分散について警告
            if max_factor / min_factor > 3:
                logger.warning(f"Large dispersion in normalization factors detected (Max/Min ratio: {max_factor/min_factor:.2f})")
                
            # 正規化係数の上位5サンプルをログに記録
            norm_factors_str = ", ".join([f"{sample}:{factor:.4f}" for sample, factor in sorted(norm_factors_dict.items(), 
                                        key=lambda x: x[1])[:5]])
            logger.debug(f"Normalization factors example (top 5 samples): {norm_factors_str}...")
        else:
            logger.warning(f"Number of normalization factors ({len(norm_factors)}) does not match number of samples ({len(count_cols)})")
                
    except Exception as e:
        logger.warning(f"Error occurred while retrieving normalization factors: {e}")
        logger.warning("Cannot display normalization factor details, but computation will continue")

def calculate_log_cpm(
    dge_norm: Any,
    count_cols: List[str],
    edger: Any,
    logger: logging.Logger
) -> pd.DataFrame:
    """
    Calculates log-CPM values from a normalized DGEList.
    
    Args:
        dge_norm (Any): Normalized edgeR DGEList object
        count_cols (List[str]): Count data column names
        edger (Any): edgeR package
        logger (logging.Logger): Logger object
        
    Returns:
        pd.DataFrame: DataFrame of log-CPM values
    """
    try:
        # log-CPMを計算
        logcpm = edger.cpm(dge_norm, normalized_lib_sizes=True, log=True)
        
        # pandasのDataFrameに変換
        logcpm_df = pd.DataFrame(np.array(logcpm), columns=count_cols)
        
        # 統計量をログに記録
        try:
            logcpm_min = logcpm_df.values.min()
            logcpm_max = logcpm_df.values.max()
            logcpm_mean = logcpm_df.values.mean()
            logger.info(f"Normalized logCPM statistics - Min: {logcpm_min:.2f}, Max: {logcpm_max:.2f}, Mean: {logcpm_mean:.2f}")
        except Exception as e:
            logger.warning(f"Error calculating normalized logCPM statistics: {e}")
        
        return logcpm_df
        
    except Exception as e:
        logger.error(f"Error calculating log-CPM values: {e}", exc_info=True)
        raise

def apply_edger_normalization(
    counts_df: pd.DataFrame,
    min_cpm: float,
    min_samples: int,
    logger: logging.Logger,
    method: str = "upperquartile"
) -> Optional[pd.DataFrame]:
    """
    Applies edgeR normalization to count data.
    
    Args:
        counts_df (pd.DataFrame): DataFrame containing featureCounts output
        min_cpm (float): Minimum CPM threshold for filtering
        min_samples (int): Minimum number of samples with CPM > threshold. Set to 0 to disable filtering.
        logger (logging.Logger): Logger object
        method (str): Normalization method (upperquartile, TMM, RLE, none)
        
    Returns:
        Optional[pd.DataFrame]: Normalized log-CPM DataFrame, or None if error occurs
    """
    logger.info(f"Applying edgeR normalization with method: {method}")
    
    if not RPY2_AVAILABLE:
        logger.error("rpy2 is not available. Cannot perform normalization.")
        return None
    
    try:
        # R環境のセットアップ
        r, base, edger = setup_r_environment(logger)
        if not r or not base or not edger:
            return None
        
        # カウント列の特定
        if len(counts_df.columns) <= FEATURE_COUNTS_METADATA_COL_COUNT:
            logger.error(f"No count columns found in DataFrame (expected columns > {FEATURE_COUNTS_METADATA_COL_COUNT})")
            return None
        
        count_cols = counts_df.columns[FEATURE_COUNTS_METADATA_COL_COUNT:]
        logger.info(f"Found {len(count_cols)} count columns")
        
        counts_matrix = counts_df[count_cols]
        
        # DGEListの作成
        dge = create_dge_list(counts_matrix, r, edger, logger)
        
        # フィルタリングの適用
        dge_filtered, kept_indices = filter_low_expression_features(dge, min_cpm, min_samples, r, base, edger, logger)
        
        # 正規化係数の計算
        dge_norm = edger.calcNormFactors(dge_filtered, method=method)
        
        # 正規化係数をログに記録
        get_normalization_factors(dge_norm, count_cols, r, logger)
        
        # log-CPM値を計算
        logcpm_df = calculate_log_cpm(dge_norm, count_cols, edger, logger)
        
        # メタデータと結合
        if min_samples == 0 or kept_indices is None:
            result_df = pd.concat([counts_df.iloc[:, :FEATURE_COUNTS_METADATA_COL_COUNT], logcpm_df], axis=1)
        else:
            result_df = counts_df.iloc[kept_indices, :FEATURE_COUNTS_METADATA_COL_COUNT].copy()
            result_df = pd.concat([result_df, logcpm_df], axis=1)
        
        # Start列を0-basedに変換（featureCountsは1-based）
        if 'Start' in result_df.columns:
            logger.debug("Converting 'Start' column to 0-based (BED format)")
            result_df['Start'] = result_df['Start'] - 1
        
        logger.info(f"Normalization completed. Result DataFrame has {len(result_df)} rows")
        return result_df
        
    except Exception as e:
        logger.error(f"Error during edgeR normalization: {e}", exc_info=True)
        return None

# --- 出力形式変換 ---
def convert_to_tsv_format(
    normalized_df: pd.DataFrame,
    logger: logging.Logger,
    output_id_name: str = "PeakID"
) -> pd.DataFrame:
    """
    Converts normalized DataFrame to standard TSV format with PeakID, Chr, Start, End and sample columns.
    
    Args:
        normalized_df (pd.DataFrame): Normalized DataFrame with featureCounts metadata
        logger (logging.Logger): Logger object
        output_id_name (str): Name to use for the ID column in output DataFrame. Defaults to "PeakID".
        
    Returns:
        pd.DataFrame: DataFrame in TSV format (PeakID, Chr, Start, End, sample1, sample2, ...)
    """
    logger.info("Converting normalized data to standard TSV format")
    
    try:
        # 必要なカラムを確認
        required_cols = ['Geneid', 'Chr', 'Start', 'End']
        missing_cols = [col for col in required_cols if col not in normalized_df.columns]
        
        if missing_cols:
            logger.error(f"Required columns for TSV format missing: {missing_cols}")
            logger.warning("Returning original DataFrame without conversion")
            return normalized_df
        
        # サンプルカラムの特定
        metadata_cols = [col for col in FEATURE_COUNTS_METADATA_COLS if col in normalized_df.columns]
        sample_cols = [col for col in normalized_df.columns if col not in metadata_cols]
        
        # 必要なカラムのみを選択
        tsv_df = normalized_df[['Geneid', 'Chr', 'Start', 'End'] + sample_cols].copy()
        
        # 'Geneid'を指定されたIDカラム名にリネーム
        tsv_df = tsv_df.rename(columns={'Geneid': output_id_name})
        
        logger.info(f"Converted to TSV format with {len(tsv_df)} rows and {len(tsv_df.columns)} columns")
        logger.debug(f"TSV format columns: {list(tsv_df.columns[:10])}...")
        
        return tsv_df
    
    except Exception as e:
        logger.error(f"Error converting to TSV format: {e}", exc_info=True)
        logger.warning("Returning original DataFrame without conversion")
        return normalized_df

def ensure_tsv_extension(file_path: str, logger: logging.Logger) -> str:
    """
    Ensures the file path has a .tsv extension.
    
    Args:
        file_path (str): Original file path
        logger (logging.Logger): Logger object
        
    Returns:
        str: File path with .tsv extension
    """
    if not file_path.lower().endswith('.tsv'):
        base_name, ext = os.path.splitext(file_path)
        new_path = f"{base_name}.tsv"
        logger.info(f"Changed output file extension to .tsv: {os.path.basename(new_path)}")
        return new_path
    return file_path