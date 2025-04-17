# === bigwig_peakprep_bamdiff_utils.py ===

import os
import sys
import time
import logging
import subprocess
import yaml
import pandas as pd
import datetime
from typing import Optional, List, Dict, Tuple, Any, TypeVar
from natsort import natsorted

# cpm_peakprep_utils.pyから必要な関数をインポート
try:
    from cpm_peakprep_utils import get_sample_data, invert_dictionary
except ImportError as e:
    print(f"Error importing functions from utility modules: {e}", file=sys.stderr)
    print("Please ensure 'cpm_peakprep_utils.py' exists and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

K = TypeVar('K')
V = TypeVar('V')

def _save_bamcompare_log(log_path: str, stdout_content: Optional[str], stderr_content: Optional[str], logger: logging.Logger):
    """Saves stdout and stderr content from bamCompare to the specified log file."""
    # ログ保存試行のログ
    logger.info(f"Attempting to save bamCompare log to: {log_path}")
    try:
        # 保存先ディレクトリが存在しない場合は作成
        log_dir = os.path.dirname(log_path)
        if log_dir and not os.path.exists(log_dir):
            os.makedirs(log_dir, exist_ok=True)
            # ディレクトリ作成ログ
            logger.info(f"Created directory for log file: {log_dir}")

        # ファイルに書き込み (encoding指定)
        with open(log_path, 'w', encoding='utf-8') as f:
            f.write("--- bamCompare stdout ---\n")
            # None や空文字列の場合に対応
            f.write((stdout_content if stdout_content else "[No stdout captured or content was empty]") + "\n")
            f.write("\n--- bamCompare stderr ---\n")
            # None や空文字列の場合に対応
            f.write((stderr_content if stderr_content else "[No stderr captured or content was empty]") + "\n")
        # ログ保存成功ログ
        logger.info(f"Successfully saved bamCompare log to: {log_path}")
    except Exception as e:
        # ログ保存失敗時のエラーログ
        # ログ保存の失敗は run_bamcompare の成否には影響させない
        logger.error(f"Failed to save bamCompare log to {log_path}: {e}")


def _get_bamcompare_version(bamcompare_path: str, logger: logging.Logger) -> Optional[str]:
    """Gets the version of bamCompare by running 'bamCompare --version'."""
    try:
        # バージョン確認コマンド実行
        version_cmd = [bamcompare_path, "--version"]
        result = subprocess.run(version_cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        version_info = result.stdout.strip() if result.stdout else "Unknown version"
        logger.debug(f"bamCompare version: {version_info}")
        return version_info
    except Exception as e:
        logger.warning(f"Failed to get bamCompare version: {e}")
        return None


def run_bamcompare(
    bam1_file: str,
    bam2_file: str,
    output_dir: str,
    logger: logging.Logger,
    output_filename: Optional[str] = None,
    operation: str = "log2",
    pseudocount: float = 1.0,
    bin_size: int = 50,
    effective_genome_size: int = 2913022398,
    bamcompare_path: str = 'bamCompare',
    threads: int = 1,
    log_file: Optional[str] = None,
    additional_options: Optional[List[str]] = None
) -> Optional[str]:
    """
    Runs bamCompare from deeptools to generate differential coverage tracks between two BAM files.
    
    Args:
        bam1_file (str): Path to the first BAM file (typically treatment/ChIP).
        bam2_file (str): Path to the second BAM file (typically control/Input).
        output_dir (str): Directory to save the output bigWig file.
        logger (logging.Logger): Logger object for logging messages.
        output_filename (Optional[str]): Name for the output file. 
            If None, defaults to '{bam1_basename}_vs_{bam2_basename}.{operation}.bigwig'.
        operation (str): Operation to perform between the two BAM files. 
            Options include "log2", "ratio", "subtract", "add", "mean", "reciprocal_ratio", "first", "second".
            Defaults to "log2" (log2 ratio).
        pseudocount (float): Pseudocount to add when calculating log2 or ratios. Defaults to 1.0.
        bin_size (int): Size of the bins in bases. Defaults to 50.
        effective_genome_size (int): The effective genome size used for normalization. 
            Defaults to 2913022398 (human genome).
        bamcompare_path (str): Path to the bamCompare executable. Defaults to 'bamCompare'.
        threads (int): Number of threads to use. Defaults to 1.
        log_file (Optional[str]): Path to save the stdout and stderr from bamCompare.
            If None, logs are not saved to a separate file. Defaults to None.
        additional_options (Optional[List[str]]): List of additional arguments
            to pass to bamCompare. Defaults to None.
            
    Returns:
        Optional[str]: The full path to the generated bigWig file if successful,
                      None otherwise.
    """
    # --- 入力値の検証 ---
    for bam_file, bam_label in [(bam1_file, "First BAM"), (bam2_file, "Second BAM")]:
        if not bam_file:
            logger.error(f"{bam_label} file path not provided to run_bamcompare.")
            return None
        
        if not os.path.exists(bam_file):
            logger.error(f"{bam_label} file not found: {bam_file}")
            return None
        
    # --- bamCompareのバージョン確認 ---
    _get_bamcompare_version(bamcompare_path, logger)
        
    # --- 出力ディレクトリの確認 ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory {output_dir}: {e}")
        return None
        
    # --- 出力ファイル名の決定 ---
    if output_filename is None:
        bam1_basename = os.path.splitext(os.path.basename(bam1_file))[0]
        bam2_basename = os.path.splitext(os.path.basename(bam2_file))[0]
        output_filename = f"{bam1_basename}_vs_{bam2_basename}.{operation}.bigwig"
        
    output_file_path = os.path.join(output_dir, output_filename)
    
    # --- bamCompareコマンドリストの構築 ---
    cmd = [bamcompare_path]
    
    # 基本オプション
    cmd.extend(['-b1', bam1_file])
    cmd.extend(['-b2', bam2_file])
    cmd.extend(['--outFileName', output_file_path])
    cmd.extend(['--operation', operation])
    cmd.extend(['--pseudocount', str(pseudocount)])
    cmd.extend(['--binSize', str(bin_size)])
    cmd.extend(['--effectiveGenomeSize', str(effective_genome_size)])
    cmd.extend(['-p', str(threads)])
    
    # 追加オプション
    if additional_options:
        cmd.extend(additional_options)
    
    # --- サブプロセスとしてbamCompareを実行 ---
    # 実行開始ログ
    logger.info(f"Running bamCompare for {os.path.basename(bam1_file)} vs {os.path.basename(bam2_file)}...")
    # 実行するコマンド文字列をログに出力
    logger.debug(f"Command: {' '.join(cmd)}")
    
    stdout_log: Optional[str] = None
    stderr_log: Optional[str] = None
    
    # 実行時間計測の準備
    start_time = time.time()
    
    try:
        # コマンド実行
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding='utf-8')
        
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"bamCompare execution took {execution_time:.2f} seconds.")
        
        # 実行結果のstdout/stderrを取得
        stdout_log = result.stdout
        stderr_log = result.stderr
        
        # 正常終了ログ
        logger.info("bamCompare finished successfully.")
        
        # stdoutがあればデバッグログ出力
        if stdout_log:
            logger.debug(f"bamCompare stdout:\n{stdout_log}")
        
        # stderrがあればデバッグログ出力
        if stderr_log:
            logger.debug(f"bamCompare stderr:\n{stderr_log}")
            
        # ログファイルが指定されていれば保存
        if log_file:
            _save_bamcompare_log(log_file, stdout_log, stderr_log, logger)
            
        # 成功したので出力ファイルの絶対パスを返す
        return os.path.abspath(output_file_path)
        
    # --- エラーハンドリング ---
    except FileNotFoundError:
        # コマンドが見つからない場合
        logger.error(f"Error: '{bamcompare_path}' command not found.")
        logger.error("Please ensure bamCompare is installed and in your PATH, or provide the correct full path.")
        return None
    except subprocess.CalledProcessError as e:
        # コマンドがエラー終了した場合
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed bamCompare execution took {execution_time:.2f} seconds.")
        
        stdout_log = e.stdout
        stderr_log = e.stderr
        logger.error(f"bamCompare failed with exit code {e.returncode}")
        if stdout_log: logger.error(f"Captured stdout on failure:\n{stdout_log}")
        if stderr_log: logger.error(f"Captured stderr on failure:\n{stderr_log}")
        
        # エラー時もログ保存試行
        if log_file:
            _save_bamcompare_log(log_file, stdout_log, stderr_log, logger)
            
        # 不完全ファイルの削除試行
        if os.path.exists(output_file_path):
            logger.warning(f"Attempting to remove potentially incomplete output file: {output_file_path}")
            try: 
                os.remove(output_file_path)
            except OSError as rm_err: 
                logger.error(f"Failed to remove incomplete output file {output_file_path}: {rm_err}")
        return None
    except Exception as e:
        # その他の予期せぬエラー
        # 実行時間の計算と記録
        end_time = time.time()
        execution_time = end_time - start_time
        logger.debug(f"Failed bamCompare execution took {execution_time:.2f} seconds.")
        
        logger.error(f"An unexpected error occurred while running bamCompare: {e}", exc_info=True)
        
        # ログ保存試行
        if log_file:
            _save_bamcompare_log(log_file, stdout_log, stderr_log, logger)
        return None


def find_longest_common_prefix(sample_name: str, control_name: str, min_length: int = 2) -> str:
    """
    Find the longest common prefix between sample and control names.

    This function (assumed to be defined elsewhere) removes trailing hyphens
    or underscores from the common prefix and returns it only if its length
    meets or exceeds the specified minimum length.

    Args:
        sample_name (str): The sample name string.
        control_name (str): The control name string.
        min_length (int, optional): The minimum required length for the
                                    cleaned common prefix to be considered valid.
                                    Defaults to 2.

    Returns:
        str: The cleaned common prefix if it meets the minimum length,
             otherwise an empty string.
    """
    # 例: os.path.commonprefix を使う実装
    common_prefix = os.path.commonprefix([sample_name, control_name])
    cleaned_prefix = common_prefix.rstrip('-_')
    if len(cleaned_prefix) >= min_length:
        return cleaned_prefix
    else:
        return ""


def generate_sample_control_pairs(
    sample_bams: Dict[str, str], 
    control_bams: Dict[str, str], 
    logger: logging.Logger,
    sample_control_mapping: Dict[str, str]
) -> Dict[str, Tuple[str, str, str]]:
    """
    Generates sample-control pairs for differential analysis.
    
    Args:
        sample_bams (Dict[str, str]): Dictionary mapping sample names to BAM file paths
        control_bams (Dict[str, str]): Dictionary mapping control names to BAM file paths
        logger (logging.Logger): Logger object for logging messages
        sample_control_mapping (Dict[str, str]): Mapping from sample names to control names.
            
    Returns:
        Dict[str, Tuple[str, str, str]]: Dictionary mapping pair IDs to tuples of 
            (sample BAM path, control BAM path, output name suggestion)
    """
    logger.info(f"Generating sample-control pairs for {len(sample_bams)} samples and {len(control_bams)} controls")
    
    if not sample_bams:
        logger.error("No sample BAM files provided")
        return {}
        
    if not control_bams:
        logger.error("No control BAM files provided")
        return {}
    
    pairs = {}
    
    logger.info(f"Using provided sample-control mapping for {len(sample_control_mapping)} samples")
    
    for sample_name, control_name in sample_control_mapping.items():
        if sample_name not in sample_bams:
            logger.warning(f"Sample '{sample_name}' in mapping not found in sample BAMs. Skipping.")
            continue
            
        if control_name not in control_bams:
            logger.warning(f"Control '{control_name}' in mapping not found in control BAMs. Skipping.")
            continue
            
        sample_bam = sample_bams[sample_name]
        control_bam = control_bams[control_name]
        pair_id = f"{sample_name}_vs_{control_name}"
        output_name = f"{sample_name}_vs_{control_name}.log2ratio.bigwig"
        
        pairs[pair_id] = (sample_bam, control_bam, output_name)
        logger.debug(f"Created pair '{pair_id}': {sample_name} vs {control_name}")
    
    logger.info(f"Generated {len(pairs)} sample-control pairs for differential analysis")
    return pairs


def create_sample_control_mapping(
    sample_bams: Dict[str, str],
    control_bams: Dict[str, str],
    logger: logging.Logger,
    strategy: str = "natsort",
    output_mapping_file: Optional[str] = None
) -> Dict[str, str]:
    """
    Creates a mapping between sample and control names based on the specified strategy.
    
    Args:
        sample_bams (Dict[str, str]): Dictionary mapping sample names to BAM file paths
        control_bams (Dict[str, str]): Dictionary mapping control names to BAM file paths
        logger (logging.Logger): Logger object for logging messages
        strategy (str, optional): Strategy to use for mapping. Options:
            - "natsort": Natural sort both lists and pair items in order (requires equal counts)
            - "fix_input": Use a single input control for all samples
            Defaults to "natsort".
        output_mapping_file (Optional[str], optional): If provided, write the mapping to this file
            
    Returns:
        Dict[str, str]: Mapping from sample names to control names
    """
    logger.info(f"Creating sample-control mapping using strategy '{strategy}'")
    
    mapping = {}
    
    if strategy == "natsort":
        # Natural sort both sample and control names
        sorted_sample_names = natsorted(sample_bams.keys())
        sorted_control_names = natsorted(control_bams.keys())
        
        # Check if counts match
        if len(sorted_sample_names) != len(sorted_control_names):
            logger.error(f"Cannot use natsort strategy: sample count ({len(sorted_sample_names)}) does not match control count ({len(sorted_control_names)})")
            raise ValueError(f"Sample count ({len(sorted_sample_names)}) does not match control count ({len(sorted_control_names)})")
        
        # Create pairs by position
        for i, (sample_name, control_name) in enumerate(zip(sorted_sample_names, sorted_control_names)):
            mapping[sample_name] = control_name
            
            # Find common identifier for debugging
            common_id = find_common_identifier(sample_name, control_name)
            if common_id:
                logger.debug(f"Pair {i+1}: '{sample_name}' and '{control_name}' (common identifier: '{common_id}')")
            else:
                logger.debug(f"Pair {i+1}: '{sample_name}' and '{control_name}' (no common identifier)")
    
    elif strategy == "fix_input":
        # Check if there's exactly one control
        if len(control_bams) != 1:
            logger.error(f"Cannot use fix_input strategy: expected exactly 1 control, got {len(control_bams)}")
            raise ValueError(f"The fix_input strategy requires exactly 1 control, got {len(control_bams)}")
        
        # Use the single control for all samples
        control_name = next(iter(control_bams.keys()))
        
        for sample_name in sample_bams.keys():
            mapping[sample_name] = control_name
            logger.debug(f"Mapped sample '{sample_name}' to fixed control '{control_name}'")
    
    else:
        logger.error(f"Unknown mapping strategy: '{strategy}'")
        raise ValueError(f"Unknown mapping strategy: '{strategy}'")
    
    # Write mapping to file if requested
    if output_mapping_file:
        try:
            with open(output_mapping_file, 'w') as f:
                f.write("sample_name\tcontrol_name\n")
                for sample_name, control_name in mapping.items():
                    f.write(f"{sample_name}\t{control_name}\n")
            logger.info(f"Wrote mapping to file: {output_mapping_file}")
        except Exception as e:
            logger.error(f"Failed to write mapping to file {output_mapping_file}: {str(e)}")
    
    logger.info(f"Created mapping for {len(mapping)} samples")
    return mapping


def find_common_identifier(sample_name: str, control_name: str) -> str:
    """
    Find common sample identifier between sample and control names.
    For ChIP-seq data, this is typically something like 'T1' in 'T1-H3K27ac' and 'T1-Input'.
    
    Args:
        sample_name (str): Sample name (e.g., 'T1-H3K27ac_R1.mLb.clN.sorted')
        control_name (str): Control name (e.g., 'T1-Input_R1.mLb.clN.sorted')
        
    Returns:
        str: Common identifier if found, or empty string
    """
    # Extract sample ID (typically before first hyphen)
    sample_id = sample_name.split('-')[0] if '-' in sample_name else ""
    control_id = control_name.split('-')[0] if '-' in control_name else ""
    
    if sample_id and sample_id == control_id:
        return sample_id
    
    # If not found, look for any common significant token
    sample_tokens = set(sample_name.split('_')[0].split('-'))
    control_tokens = set(control_name.split('_')[0].split('-'))
    common_tokens = sample_tokens.intersection(control_tokens)
    
    if common_tokens:
        # Filter out very short tokens
        meaningful_tokens = [t for t in common_tokens if len(t) > 1]
        if meaningful_tokens:
            return ', '.join(sorted(meaningful_tokens))
    
    common_prefix = find_longest_common_prefix(sample_name, control_name)
    return common_prefix



def scan_bam_files(
    sample_dir: str,
    control_dir: str,
    logger: logging.Logger,
    sample_pattern: str = "*.bam",
    control_pattern: str = "*Input*.bam",
    sample_delimiter: Optional[str] = None
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Scans sample and control BAM files and retrieves mappings from names to file paths.

    Args:
        sample_dir (str): Directory containing sample BAM files.
        control_dir (str): Directory containing control BAM files.
        logger (logging.Logger): Logger object.
        sample_pattern (str): Pattern for sample filenames.
        control_pattern (str): Pattern for control filenames.
        sample_delimiter (Optional[str]): Delimiter for extracting sample names.

    Returns:
        Tuple[Dict[str, str], Dict[str, str]]:
            - Dictionary mapping sample names to BAM file paths.
            - Dictionary mapping control names to BAM file paths.
    """
    logger.info(f"Scanning sample BAM directory: {sample_dir}")
    logger.info(f"Scanning control BAM directory: {control_dir}")

    # サンプルBAMファイルのスキャン
    try:
        sample_to_bam, sample_bam_files = get_sample_data( # Assume get_sample_data is defined elsewhere
            bam_folder=sample_dir,
            logger=logger,
            sample_name_delimiter=sample_delimiter,
            filename_pattern=sample_pattern
        )

        if not sample_bam_files:
            logger.warning(f"No sample BAM files matching pattern '{sample_pattern}' found in {sample_dir}")
            return {}, {}

        logger.info(f"Found {len(sample_bam_files)} sample BAM files as {len(sample_to_bam)} distinct samples")

        # 最初の数サンプルをログ出力
        for i, (sample, bam_path) in enumerate(sample_to_bam.items()):
            if i < 5:  # 最初の5つのサンプルのみログ出力
                logger.debug(f"Sample: {sample} -> BAM: {os.path.basename(bam_path)}")
    except Exception as e:
        logger.error(f"Error occurred while scanning sample BAM files: {e}", exc_info=True)
        return {}, {}

    # コントロールBAMファイルのスキャン
    try:
        control_to_bam, control_bam_files = get_sample_data( # Assume get_sample_data is defined elsewhere
            bam_folder=control_dir,
            logger=logger,
            sample_name_delimiter=sample_delimiter,
            filename_pattern=control_pattern
        )

        if not control_bam_files:
            logger.warning(f"No control BAM files matching pattern '{control_pattern}' found in {control_dir}")
            # Consider returning sample_to_bam here if samples were found but controls were not
            # return sample_to_bam, {} # Depending on desired behavior
            return {}, {} # Current behavior: return empty if controls are missing

        logger.info(f"Found {len(control_bam_files)} control BAM files as {len(control_to_bam)} distinct controls")

        # 最初の数コントロールをログ出力
        for i, (control, bam_path) in enumerate(control_to_bam.items()):
            if i < 5:  # 最初の5つのコントロールのみログ出力
                logger.debug(f"Control: {control} -> BAM: {os.path.basename(bam_path)}")
    except Exception as e:
        logger.error(f"Error occurred while scanning control BAM files: {e}", exc_info=True)
        # Consider returning sample_to_bam here if samples were found but scanning controls failed
        # return sample_to_bam, {} # Depending on desired behavior
        return {}, {} # Current behavior: return empty if control scanning fails

    return sample_to_bam, control_to_bam


def prepare_yaml(
    sample_bams: Dict[str, str],  # サンプル名 -> BAMパス
    control_bams: Dict[str, str],  # コントロール名 -> BAMパス  
    sample_control_pairs: Dict[str, Tuple[str, str, str]],  # ペアID -> (サンプルBAM, コントロールBAM, 出力名)
    output_dir: str,
    bamcompare_params: Dict[str, Any], 
    output_yaml_path: str,
    logger: logging.Logger
) -> Optional[str]:
    """
    Generates a YAML format conversion plan from sample-control pairs and processing parameters.
    
    Args:
        sample_bams (Dict[str, str]): Dictionary mapping sample names to BAM file paths
        control_bams (Dict[str, str]): Dictionary mapping control names to BAM file paths
        sample_control_pairs (Dict[str, Tuple[str, str, str]]): Dictionary mapping pair IDs to tuples of 
                                                              (sample BAM, control BAM, output name)
        output_dir (str): Output directory (where the generated bigWig files will be saved)
        bamcompare_params (Dict[str, Any]): Dictionary of parameters to pass to bamCompare
                                  (operation, bin_size, pseudocount, effective_genome_size, threads, etc.)
        output_yaml_path (str): Path for the output YAML file
        logger (logging.Logger): Logger object for logging
        
    Returns:
        Optional[str]: Path to the generated YAML file, or None if an error occurred
    """
    logger.info(f"Preparing YAML conversion plan file: {output_yaml_path}")
    
    try:
        # BAMサンプル名とコントロール名の辞書を作成（パス->名前）
        bam_path_to_sample = {path: name for name, path in sample_bams.items()}
        bam_path_to_control = {path: name for name, path in control_bams.items()}
        
        # 出力用のデータ構造を初期化
        yaml_data = {
            "meta_info": {
                "generated_date": time.strftime("%Y-%m-%d %H:%M:%S"),
                "pipeline": "bigwig_peakprep_bamdiff pipeline"
            },
            "parameters": bamcompare_params.copy(),  # bamCompareパラメータのコピー
            "samples": {},         # サンプル情報
            "controls": {},        # コントロール情報
            "conversion_plan": []  # 変換予定リスト
        }
        
        # サンプル情報を追加
        for sample_name, bam_path in sample_bams.items():
            yaml_data["samples"][sample_name] = {
                "bam_file": os.path.basename(bam_path),
                "bam_path": bam_path
            }
            
        # コントロール情報を追加
        for control_name, bam_path in control_bams.items():
            yaml_data["controls"][control_name] = {
                "bam_file": os.path.basename(bam_path),
                "bam_path": bam_path
            }
        
        # ペア変換情報を追加
        conversion_items = []
        for pair_id, (sample_bam, control_bam, output_name) in sample_control_pairs.items():
            sample_name = bam_path_to_sample.get(sample_bam, os.path.basename(sample_bam))
            control_name = bam_path_to_control.get(control_bam, os.path.basename(control_bam))
            
            # 出力ファイルのフルパスを生成
            output_path = os.path.join(output_dir, output_name)
            
            # 変換アイテムの作成
            conversion_item = {
                "pair_id": pair_id,
                "sample_name": sample_name,
                "control_name": control_name,
                "sample_bam": os.path.basename(sample_bam),
                "control_bam": os.path.basename(control_bam),
                "output_bigwig": output_name,
                "sample_bam_path": sample_bam,
                "control_bam_path": control_bam,
                "output_bigwig_path": output_path
            }
            conversion_items.append(conversion_item)
            
        yaml_data["conversion_plan"] = conversion_items
        
        # 表形式データをDataFrameとして作成し、TSVも出力する
        table_rows = []
        for item in conversion_items:
            table_row = {
                "sample_name": item["sample_name"],
                "control_name": item["control_name"],
                "sample_bam": item["sample_bam"],
                "control_bam": item["control_bam"],
                "output_bigwig": item["output_bigwig"],
                "sample_bam_path": item["sample_bam_path"],
                "control_bam_path": item["control_bam_path"],
                "output_bigwig_path": item["output_bigwig_path"]
            }
            table_rows.append(table_row)
        
        table_df = pd.DataFrame(table_rows)
        
        # YAMLファイルの出力ディレクトリを確保
        yaml_dir = os.path.dirname(output_yaml_path)
        if yaml_dir and not os.path.exists(yaml_dir):
            os.makedirs(yaml_dir, exist_ok=True)
            logger.debug(f"Created directory for YAML file: {yaml_dir}")
        
        # YAMLファイルとして保存
        with open(output_yaml_path, 'w', encoding='utf-8') as f:
            yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
            
        # 同じ場所にTSVも保存
        tsv_path = os.path.splitext(output_yaml_path)[0] + ".tsv"
        table_df.to_csv(tsv_path, sep='\t', index=False)
        
        logger.info(f"Conversion plan saved as YAML: {output_yaml_path}")
        logger.info(f"Conversion plan table saved as TSV: {tsv_path}")
        
        return output_yaml_path
        
    except Exception as e:
        logger.error(f"Error creating YAML conversion plan: {e}", exc_info=True)
        return None
    

def load_yaml_conversion_plan(
    yaml_path: str, 
    logger: logging.Logger
) -> Optional[Dict[str, Any]]:
    """
    Loads the conversion plan from a YAML file.
    
    Args:
        yaml_path (str): Path to the YAML file
        logger (logging.Logger): Logger object for logging
        
    Returns:
        Optional[Dict[str, Any]]: Dictionary containing YAML data, or None if an error occurred
    """
    logger.info(f"Loading YAML conversion plan from: {yaml_path}")
    
    try:
        if not os.path.exists(yaml_path):
            logger.error(f"YAML file not found: {yaml_path}")
            return None
            
        with open(yaml_path, 'r', encoding='utf-8') as f:
            yaml_data = yaml.safe_load(f)
            
        # 必要なキーが存在するか確認
        required_keys = ["parameters", "conversion_plan"]
        missing_keys = [key for key in required_keys if key not in yaml_data]
        if missing_keys:
            logger.error(f"Required keys missing in YAML file: {missing_keys}")
            return None
            
        logger.info(f"Successfully loaded conversion plan with {len(yaml_data.get('conversion_plan', []))} items")
        return yaml_data
        
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML file {yaml_path}: {e}")
        return None
    except Exception as e:
        logger.error(f"Unexpected error loading YAML file {yaml_path}: {e}", exc_info=True)
        return None


def execute_conversion_plan(
    yaml_path: str,
    logger: logging.Logger,
    output_dir: Optional[str] = None,
    bamcompare_path: str = 'bamCompare',
    force_overwrite: bool = False,
    max_concurrent: int = 1,
    additional_options: Optional[List[str]] = None,
    write_results_yaml: bool = False,
    results_yaml_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Loads the conversion plan from a YAML file and executes bamCompare operations.
    
    Args:
        yaml_path (str): Path to the conversion plan YAML file
        logger (logging.Logger): Logger object for logging
        output_dir (Optional[str]): Override output directory (if specified, overrides the output destinations in YAML)
        bamcompare_path (str): Path to the bamCompare executable
        force_overwrite (bool): Whether to overwrite existing files
        max_concurrent (int): Maximum number of concurrent jobs (currently only max_concurrent=1 is supported as parallel processing is not implemented)
        additional_options (Optional[List[str]]): Additional options to pass to bamCompare
        write_results_yaml (bool): Whether to write the results as a YAML file
        results_yaml_path (Optional[str]): Path for the results YAML file (if None, appends .results to the original YAML)
        
    Returns:
        Dict[str, Any]: Dictionary containing the extended YAML data with results
            - Original configuration parameters and conversion plan
            - results: Detailed processing results
                - summary: Processing summary
                    - total_jobs: Total number of jobs
                    - completed_jobs: Number of completed jobs
                    - skipped_jobs: Number of skipped jobs
                    - failed_jobs: Number of failed jobs
                - successful_files: List of successful bigWig files
                - failed_pairs: List of failed pair IDs
                - timestamp: Timestamp of processing completion
    """
    
    # YAMLファイルを読み込む
    yaml_data = load_yaml_conversion_plan(yaml_path, logger)
    if yaml_data is None:
        logger.error("Failed to load YAML conversion plan. Cannot proceed.")
        empty_result = {
            "parameters": {},
            "conversion_plan": [],
            "results": {
                "summary": {
                    "total_jobs": 0,
                    "completed_jobs": 0,
                    "skipped_jobs": 0,
                    "failed_jobs": 0
                },
                "successful_files": [],
                "failed_pairs": [],
                "timestamp": datetime.datetime.now().isoformat()
            }
        }
        
        if write_results_yaml:
            _write_results_yaml(empty_result, yaml_path, results_yaml_path, logger)
            
        return empty_result
    
    # パラメータの取得
    parameters = yaml_data.get("parameters", {})
    conversion_plan = yaml_data.get("conversion_plan", [])
    
    if not conversion_plan:
        logger.warning("No conversion jobs found in the YAML file.")
        empty_result = {
            "parameters": parameters,
            "conversion_plan": [],
            "results": {
                "summary": {
                    "total_jobs": 0,
                    "completed_jobs": 0,
                    "skipped_jobs": 0,
                    "failed_jobs": 0
                },
                "successful_files": [],
                "failed_pairs": [],
                "timestamp": datetime.datetime.now().isoformat()
            }
        }
        
        if write_results_yaml:
            _write_results_yaml(empty_result, yaml_path, results_yaml_path, logger)
            
        return empty_result
    
    # bamCompareのデフォルトパラメータ設定
    operation = parameters.get("operation", "log2")
    bin_size = parameters.get("bin_size", 50)
    pseudocount = parameters.get("pseudocount", 1.0)
    effective_genome_size = parameters.get("effective_genome_size", 2913022398)
    threads = parameters.get("threads", 1)
    
    # 結果追跡用の変数
    total_jobs = len(conversion_plan)
    completed_jobs = 0
    skipped_jobs = 0
    failed_jobs = 0
    successful_files = []
    failed_pairs = []
    job_results = []
    
    logger.info(f"Starting execution of {total_jobs} bamCompare jobs from conversion plan")
    
    # 各変換ジョブを処理
    for i, job in enumerate(conversion_plan, 1):
        pair_id = job.get("pair_id", f"job_{i}")
        sample_bam = job.get("sample_bam_path", "")
        control_bam = job.get("control_bam_path", "")
        output_bigwig = job.get("output_bigwig", "")
        
        # 出力パスの調整（もしoutput_dirが指定されている場合）
        if output_dir:
            output_bigwig_filename = os.path.basename(output_bigwig)
            output_bigwig_path = os.path.join(output_dir, output_bigwig_filename)
        else:
            output_bigwig_path = job.get("output_bigwig_path", "")
            # パスが指定されていない場合はファイル名だけの可能性
            if not os.path.dirname(output_bigwig_path):
                parent_dir = os.path.dirname(yaml_path)  # YAMLファイルと同じディレクトリを使用
                output_bigwig_path = os.path.join(parent_dir, output_bigwig_path)
        
        logger.info(f"Processing job {i}/{total_jobs}: {pair_id}")
        
        # ジョブごとの結果を追跡
        job_result = {
            "pair_id": pair_id,
            "sample_bam": sample_bam,
            "control_bam": control_bam,
            "output_bigwig": output_bigwig_path,
            "status": "pending"
        }
        
        # 入力ファイルの存在確認
        if not os.path.exists(sample_bam):
            logger.error(f"Sample BAM file not found: {sample_bam}")
            failed_jobs += 1
            failed_pairs.append(pair_id)
            job_result["status"] = "failed"
            job_result["error"] = "Sample BAM file not found"
            job_results.append(job_result)
            continue
            
        if not os.path.exists(control_bam):
            logger.error(f"Control BAM file not found: {control_bam}")
            failed_jobs += 1
            failed_pairs.append(pair_id)
            job_result["status"] = "failed"
            job_result["error"] = "Control BAM file not found"
            job_results.append(job_result)
            continue
            
        # 出力ファイルが既に存在するかチェック
        if os.path.exists(output_bigwig_path) and not force_overwrite:
            logger.info(f"Output file already exists, skipping: {output_bigwig_path}")
            skipped_jobs += 1
            # 既存ファイルも成功扱いとする
            successful_files.append(output_bigwig_path)
            job_result["status"] = "skipped"
            job_result["note"] = "Output file already exists"
            job_results.append(job_result)
            continue
            
        # 出力ディレクトリの作成
        output_dir_path = os.path.dirname(output_bigwig_path)
        try:
            if output_dir_path:
                os.makedirs(output_dir_path, exist_ok=True)
        except OSError as e:
            logger.error(f"Failed to create output directory {output_dir_path}: {e}")
            failed_jobs += 1
            failed_pairs.append(pair_id)
            job_result["status"] = "failed"
            job_result["error"] = f"Failed to create output directory: {str(e)}"
            job_results.append(job_result)
            continue
            
        # ログファイルパスの設定
        log_file = os.path.splitext(output_bigwig_path)[0] + ".bamCompare.log"
        
        # 処理開始時刻の記録
        start_time = datetime.datetime.now()
        job_result["start_time"] = start_time.isoformat()
        
        # run_bamcompareの実行
        result = run_bamcompare(
            bam1_file=sample_bam,
            bam2_file=control_bam,
            output_dir=output_dir_path,
            logger=logger,
            output_filename=os.path.basename(output_bigwig_path),
            operation=operation,
            pseudocount=pseudocount,
            bin_size=bin_size,
            effective_genome_size=effective_genome_size,
            bamcompare_path=bamcompare_path,
            threads=threads,
            log_file=log_file,
            additional_options=additional_options
        )
        
        # 処理終了時刻の記録
        end_time = datetime.datetime.now()
        job_result["end_time"] = end_time.isoformat()
        job_result["duration_seconds"] = (end_time - start_time).total_seconds()
        
        if result:
            logger.info(f"Successfully completed job {i}/{total_jobs}: {pair_id}")
            completed_jobs += 1
            successful_files.append(result)
            job_result["status"] = "completed"
            job_result["output_file"] = result
        else:
            logger.error(f"Failed to complete job {i}/{total_jobs}: {pair_id}")
            failed_jobs += 1
            failed_pairs.append(pair_id)
            job_result["status"] = "failed"
            job_result["error"] = "bamCompare execution failed"
        
        job_results.append(job_result)
    
    # 処理結果サマリーのログ出力
    logger.info(f"Conversion plan execution completed:")
    logger.info(f"  Total jobs: {total_jobs}")
    logger.info(f"  Completed: {completed_jobs}")
    logger.info(f"  Skipped (existing files): {skipped_jobs}")
    logger.info(f"  Failed: {failed_jobs}")
    
    # 結果YAMLの作成
    results_yaml = {
        # meta_infoセクションを追加
        "meta_info": {
            "generated_date": datetime.datetime.now().isoformat(),
            "pipeline": "bamCompare processing pipeline",
            "host": os.uname().nodename,
            "working_directory": os.getcwd()
        },
        "parameters": parameters,
        "conversion_plan": conversion_plan,
        "results": {
            "summary": {
                "total_jobs": total_jobs,
                "completed_jobs": completed_jobs,
                "skipped_jobs": skipped_jobs,
                "failed_jobs": failed_jobs
            },
            "job_results": job_results,
            "successful_files": successful_files,
            "failed_pairs": failed_pairs,
            "timestamp": datetime.datetime.now().isoformat()
        }
    }
        
    # 結果YAMLの保存（オプション）
    if write_results_yaml:
        _write_results_yaml(results_yaml, yaml_path, results_yaml_path, logger)
    
    return results_yaml

def _write_results_yaml(results_yaml, original_yaml_path, results_yaml_path=None, logger=None):
    """
    Helper function for writing results YAML to a file
    
    Args:
        results_yaml (dict): Results data to write to YAML
        original_yaml_path (str): Path to the original YAML file
        results_yaml_path (Optional[str]): Path for saving the results YAML (auto-generated if None)
        logger (Optional[logging.Logger]): Logger
    
    Returns:
        str: Path to the written YAML file
    """

    
    # 結果YAML保存先が指定されていなければ自動生成
    if results_yaml_path is None:
        base_name, ext = os.path.splitext(original_yaml_path)
        results_yaml_path = f"{base_name}.results{ext}"
    
    try:
        with open(results_yaml_path, 'w') as yaml_file:
            yaml.safe_dump(results_yaml, yaml_file, default_flow_style=False)
        
        if logger:
            logger.info(f"Results saved to YAML file: {results_yaml_path}")
        
        return results_yaml_path
    except Exception as e:
        if logger:
            logger.error(f"Failed to write results YAML: {e}")
        return None

def save_bamdiff_results(
    results_yaml: Dict[str, Any],
    output_dir: str,
    logger: logging.Logger,
    output_basename: str = "bamdiff_to_bigwig_details"
) -> Tuple[Optional[str], Optional[str]]:
    """
    Saves the output results of the enhanced execute_conversion_plan as TSV and YAML files.
    These files are formatted to be compatible with bigwig_peakprep_summary.py as input,
    maintaining compatibility with the output format of bigwig_peakprep_bamconvert.py.
    
    Args:
        results_yaml (Dict[str, Any]): Output result YAML from the enhanced execute_conversion_plan
        output_dir (str): Output directory
        logger (logging.Logger): Logger object for logging messages
        output_basename (str): Base name for output files
        
    Returns:
        Tuple[Optional[str], Optional[str]]: 
            - Path to the TSV file, or None if an error occurred
            - Path to the YAML file, or None if an error occurred
    """

    logger.info(f"Saving bamCompare analysis results to: {output_dir}")
    
    # 実行時間の計算用
    current_time = time.time()
    
    try:
        # 出力ディレクトリの確認
        os.makedirs(output_dir, exist_ok=True)
        
        # 改良版execute_conversion_planの出力構造から必要な情報を抽出
        parameters = results_yaml.get("parameters", {})
        conversion_plan = results_yaml.get("conversion_plan", [])
        
        # 結果セクションへのアクセス
        results = results_yaml.get("results", {})
        summary = results.get("summary", {})
        
        # 成功したファイルリストの取得
        successful_files = results.get("successful_files", [])
        # 失敗したペアIDのリスト
        failed_pairs = results.get("failed_pairs", [])
        # 詳細なジョブ結果情報
        job_results = results.get("job_results", [])
        
        if not successful_files:
            logger.warning("No successful bigWig files to report.")
            return None, None
        
        # 成功したジョブを特定（改良版では job_results から直接取得できる）
        successful_pairs = []
        
        # job_resultsがある場合はそれを使用（改良版からの出力）
        if job_results:
            for job_result in job_results:
                if job_result.get("status") in ["completed", "skipped"]:
                    # job_resultsには元のjobの情報が含まれていないため、conversion_planから補完
                    pair_id = job_result.get("pair_id", "")
                    for job in conversion_plan:
                        if job.get("pair_id") == pair_id:
                            # job_resultの情報で一部更新
                            merged_job = job.copy()
                            merged_job["output_bigwig_path"] = job_result.get("output_file", job.get("output_bigwig_path", ""))
                            successful_pairs.append(merged_job)
                            break
        else:
            # 従来の方法（job_resultsがない場合）
            for job in conversion_plan:
                pair_id = job.get("pair_id", "")
                if pair_id in failed_pairs:
                    continue  # 失敗したペアはスキップ
                    
                output_path = job.get("output_bigwig_path", "")
                if output_path in successful_files or any(f.endswith(job.get("output_bigwig", "")) for f in successful_files):
                    successful_pairs.append(job)
        
        # 整理されたデータを準備
        rows = []
        for job in successful_pairs:
            sample_name = job.get("sample_name", "")
            bigwig_filename = job.get("output_bigwig", "")
            bigwig_path = job.get("output_bigwig_path", "")
            
            # サンプル名がない場合はペアIDを使用
            if not sample_name:
                sample_name = job.get("pair_id", "unknown")
                
            row = {
                "Sample_name": sample_name,
                "BigWig_filename": bigwig_filename,
                "BigWig_fullpath": bigwig_path
            }
            rows.append(row)
        
        # DataFrame作成
        df = pd.DataFrame(rows)
        
        # --- TSVファイルとして保存 ---
        tsv_path = os.path.join(output_dir, f"{output_basename}.tsv")
        
        with open(tsv_path, 'w', encoding='utf-8') as f:
            # コメント行（メタデータ）
            f.write(f"# bigwig_peakprep_bamdiff.py - BigWig files generated from BAM differential analysis\n")
            f.write(f"# Generation date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Host: {os.uname().nodename}\n")
            f.write(f"# Working directory: {os.getcwd()}\n")
            f.write(f"# Output directory: {os.path.abspath(output_dir)}\n")
            f.write(f"# Processing results: {len(successful_files)} successful bigWig files\n")
            f.write(f"# already_log2_transformed: true\n")
            f.write(f"# bamCompare operation: {parameters.get('operation', 'log2')}\n")
            f.write(f"# Pseudocount: {parameters.get('pseudocount', 1.0)}\n")
            f.write(f"# bin_size: {parameters.get('bin_size', 50)}\n")
            f.write(f"# effective_genome_size: {parameters.get('effective_genome_size', 2913022398)}\n")
            f.write(f"# This file can be used as input for bigwig_peakprep_summary.py\n")
            f.write(f"# Format: TSV (Sample_name, BigWig_filename, BigWig_fullpath)\n")
            f.write("#\n")
            
        # TSV本体の保存
        df.to_csv(tsv_path, sep='\t', index=False, mode='a')
        logger.info(f"Saved bigWig details to TSV file: {tsv_path}")
        
        # --- YAMLファイルとして詳細情報を保存 ---
        yaml_path = os.path.join(output_dir, f"{output_basename}.yaml")
        
        # YAMLデータの構築
        yaml_output = {
            "meta_info": {
                "generated_date": time.strftime("%Y-%m-%d %H:%M:%S"),
                "pipeline": "bigwig_peakprep_bamdiff pipeline",
                "host": os.uname().nodename,
                "working_directory": os.getcwd(),
                "output_directory": os.path.abspath(output_dir),
                "already_log2_transformed": True
            },
            "parameters": parameters,
            "execution_results": {
                "total_jobs": summary.get("total_jobs", 0),
                "completed_jobs": summary.get("completed_jobs", 0),
                "skipped_jobs": summary.get("skipped_jobs", 0),
                "failed_jobs": summary.get("failed_jobs", 0)
            },
            "successful_files": [
                {
                    "sample_name": job.get("sample_name", ""),
                    "control_name": job.get("control_name", ""),
                    "pair_id": job.get("pair_id", ""),
                    "bigwig_filename": job.get("output_bigwig", ""),
                    "bigwig_path": job.get("output_bigwig_path", ""),
                    "sample_bam": job.get("sample_bam", ""),
                    "control_bam": job.get("control_bam", ""),
                    "sample_bam_path": job.get("sample_bam_path", ""),
                    "control_bam_path": job.get("control_bam_path", "")
                }
                for job in successful_pairs
            ]
        }
        
        # 改良版の詳細なジョブ情報がある場合は追加
        if job_results:
            yaml_output["job_details"] = job_results
            
        # 元のYAMLにあったタイムスタンプ情報を追加
        if "timestamp" in results:
            yaml_output["meta_info"]["execution_timestamp"] = results["timestamp"]
        
        # YAMLファイル保存
        with open(yaml_path, 'w', encoding='utf-8') as f:
            yaml.dump(yaml_output, f, default_flow_style=False, sort_keys=False)
        
        logger.info(f"Saved detailed results to YAML file: {yaml_path}")
        
        return tsv_path, yaml_path
        
    except Exception as e:
        logger.error(f"Error saving bamCompare results: {e}", exc_info=True)
        return None, None



# # --- より高度なサンプル-コントロールマッピング関数 --- 
# """
# def map_sample_to_control(
#     sample_to_bam: Dict[str, str],
#     control_to_bam: Dict[str, str],
#     logger: logging.Logger,
#     mapping_file: Optional[str] = None,
#     mapping_pattern: Optional[str] = None,
#     use_same_name: bool = False
# ) -> Dict[str, Tuple[str, str, str, str]]:
#     """
#     Map sample BAM files to control BAM files with advanced matching options
    
#     Args:
#         sample_to_bam (Dict[str, str]): Dictionary mapping sample names to BAM paths
#         control_to_bam (Dict[str, str]): Dictionary mapping control names to BAM paths
#         logger (logging.Logger): Logger object
#         mapping_file (Optional[str]): Path to tab-delimited mapping file
#         mapping_pattern (Optional[str]): Regex pattern to extract control name from sample name
#         use_same_name (bool): Use same base name for sample and control
        
#     Returns:
#         Dict[str, Tuple[str, str, str, str]]: Dictionary mapping pair_id to (sample_bam, control_bam, sample_name, control_name)
#     """
#     logger.info("Creating sample to control BAM mapping")
    
#     pairs = {}
    
#     # マッピングファイルを優先的に使用
#     if mapping_file and os.path.exists(mapping_file):
#         logger.info(f"Using mapping file: {mapping_file}")
#         try:
#             with open(mapping_file, 'r') as f:
#                 for line in f:
#                     if line.startswith('#') or not line.strip():
#                         continue
#                     parts = line.strip().split('\t')
#                     if len(parts) >= 2:
#                         sample_name, control_name = parts[0], parts[1]
#                         if sample_name in sample_to_bam and control_name in control_to_bam:
#                             pair_id = f"{sample_name}_vs_{control_name}"
#                             pairs[pair_id] = (
#                                 sample_to_bam[sample_name],
#                                 control_to_bam[control_name],
#                                 sample_name,
#                                 control_name
#                             )
#                             logger.debug(f"Mapped: {sample_name} -> {control_name} (Pair ID: {pair_id})")
#                         else:
#                             if sample_name not in sample_to_bam:
#                                 logger.warning(f"Sample '{sample_name}' from mapping file not found in BAM files")
#                             if control_name not in control_to_bam:
#                                 logger.warning(f"Control '{control_name}' from mapping file not found in BAM files")
            
#             logger.info(f"Created {len(pairs)} sample-control pairs from mapping file")
#             return pairs
#         except Exception as e:
#             logger.error(f"Error reading mapping file: {e}")
#             logger.warning("Falling back to automatic mapping")
    
#     # 正規表現パターンを使用
#     if mapping_pattern:
#         pattern = re.compile(mapping_pattern)
#         for sample_name, sample_bam in sample_to_bam.items():
#             match = pattern.match(sample_name)
#             if match and match.groupdict().get('control'):
#                 control_name = match.groupdict().get('control')
#                 if control_name in control_to_bam:
#                     pair_id = f"{sample_name}_vs_{control_name}"
#                     pairs[pair_id] = (
#                         sample_bam,
#                         control_to_bam[control_name],
#                         sample_name,
#                         control_name
#                     )
#                     logger.debug(f"Regex mapped: {sample_name} -> {control_name} (Pair ID: {pair_id})")
#                 else:
#                     logger.warning(f"Control '{control_name}' extracted from sample name '{sample_name}' not found")
        
#         if pairs:
#             logger.info(f"Created {len(pairs)} sample-control pairs using regex pattern")
#             return pairs
#         else:
#             logger.warning("No pairs created with regex pattern, trying other methods")
    
#     # 同じ名前を使用（例：sample1.bam -> input1.bam）
#     if use_same_name:
#         # サンプル名から共通のベース名を抽出
#         for sample_name, sample_bam in sample_to_bam.items():
#             base_name = None
#             # 数字が含まれる場合は数字部分を抽出
#             match = re.search(r'(\d+)', sample_name)
#             if match:
#                 base_name = match.group(1)
            
#             # 対応するコントロールを探す
#             matching_controls = []
#             for control_name, control_bam in control_to_bam.items():
#                 if base_name and base_name in control_name:
#                     matching_controls.append((control_name, control_bam))
            
#             # 一致するコントロールがあればペアを作成
#             if matching_controls:
#                 # 複数ある場合は "Input" または "IgG" を含むものを優先
#                 best_match = None
#                 for control_name, control_bam in matching_controls:
#                     if "input" in control_name.lower() or "igg" in control_name.lower():
#                         best_match = (control_name, control_bam)
#                         break
                
#                 # 優先順位が高いものがなければ最初のものを使用
#                 if not best_match:
#                     best_match = matching_controls[0]
                
#                 control_name, control_bam = best_match
#                 pair_id = f"{sample_name}_vs_{control_name}"
#                 pairs[pair_id] = (
#                     sample_bam,
#                     control_bam,
#                     sample_name,
#                     control_name
#                 )
#                 logger.debug(f"Name matched: {sample_name} -> {control_name} (Pair ID: {pair_id})")
        
#         if pairs:
#             logger.info(f"Created {len(pairs)} sample-control pairs using name matching")
#             return pairs
    
#     # 何も指定がない場合、単一のコントロールに全てのサンプルをマッピング
#     if len(control_to_bam) == 1:
#         control_name, control_bam = next(iter(control_to_bam.items()))
#         logger.info(f"Using single control '{control_name}' for all samples")
        
#         for sample_name, sample_bam in sample_to_bam.items():
#             pair_id = f"{sample_name}_vs_{control_name}"
#             pairs[pair_id] = (
#                 sample_bam,
#                 control_bam,
#                 sample_name,
#                 control_name
#             )
#             logger.debug(f"Single control: {sample_name} -> {control_name} (Pair ID: {pair_id})")
        
#         logger.info(f"Created {len(pairs)} sample-control pairs with single control")
#         return pairs
    
#     # それでもマッピングができない場合
#     if not pairs:
#         logger.warning("Could not automatically map samples to controls. Please provide a mapping method.")
#         logger.info(f"Available samples: {', '.join(sample_to_bam.keys())}")
#         logger.info(f"Available controls: {', '.join(control_to_bam.keys())}")
    
#     return pairs
# """

# """
# def format_sample_control_pairs(
#     sample_control_pairs: Dict[str, Tuple[str, str, str, str]],
#     operation: str = "log2"
# ) -> Dict[str, Tuple[str, str, str]]:
#     """
#     サンプル-コントロールのペア情報を、prepare_yaml関数で使用可能な形式に変換します
    
#     Args:
#         sample_control_pairs (Dict[str, Tuple[str, str, str, str]]): 
#             map_sample_to_control関数の返り値
#             pair_id -> (sample_bam, control_bam, sample_name, control_name)
#         operation (str): 操作タイプ（log2, ratio, etc.）
        
#     Returns:
#         Dict[str, Tuple[str, str, str]]: prepare_yaml関数用の形式
#             pair_id -> (sample_bam, control_bam, output_filename)
#     """
#     formatted_pairs = {}
#     for pair_id, (sample_bam, control_bam, sample_name, control_name) in sample_control_pairs.items():
#         output_filename = f"{pair_id}.{operation}ratio.bigwig"
#         formatted_pairs[pair_id] = (sample_bam, control_bam, output_filename)
#     return formatted_pairs
# """