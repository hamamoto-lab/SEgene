# === cpm_peakprep.py の修正案 ===
# 以下のコードは、--calcnormオプション使用時も両方の出力を保存できるようにするための修正です

import argparse
import logging
import os
import sys
import pandas as pd
import time
import json
import tempfile

# --- cpm_peakprep_utils.py から関数をインポート ---
try:
    from cpm_peakprep_utils import (
        get_sample_data,
        run_samtools_flagstat,
        run_featurecounts,
        calculate_logcpm,
        convert_bed_to_saf,
        convert_mergese_to_saf,
        invert_dictionary
    )
except ImportError as e:
    print(f"Error importing functions from cpm_peakprep_utils.py: {e}", file=sys.stderr)
    print("Please ensure 'cpm_peakprep_utils.py' exists and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- cpm_peakprep_calcnorm_utils.py から関数をインポート(高度な正規化計算用) ---
try:
    from cpm_peakprep_calcnorm_utils import (
        apply_edger_normalization,
        validate_counts_file,
        clean_sample_names,
        clean_count_file_sample_names,
        convert_to_tsv_format,
        load_pattern_rules,
        FEATURE_COUNTS_METADATA_COLS,
        setup_logging,
        cleanup_temp_files
    )
    CALCNORM_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import functions from cpm_peakprep_calcnorm_utils.py: {e}", file=sys.stderr)
    print("Advanced normalization calculation will not be available.", file=sys.stderr)
    CALCNORM_AVAILABLE = False

# --- もしライブラリ側を修正できない場合のバックアップ実装 ---
if CALCNORM_AVAILABLE:
    # 既存関数の引数を確認
    import inspect
    convert_to_tsv_params = inspect.signature(convert_to_tsv_format).parameters
    if len(convert_to_tsv_params) < 3 or 'output_id_name' not in convert_to_tsv_params:
        # 関数シグネチャに output_id_name がない場合、ローカルに同名関数を再定義
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

# --- 引数パーサーの設定 ---
def parse_arguments():
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Simplified ChIP-seq Quantification Pipeline: Runs flagstat, featureCounts, calculates logCPM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数 ---
    parser.add_argument("-b", "--bam_dir", required=True,
                        help="Directory containing BAM files.")
    parser.add_argument("-a", "--annotation_file", required=True,
                        help="Path to the annotation file (SAF, BED, or merge_SE format). Specify format with --is_mergese_format if needed.")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files.")

    # --- オプション引数: ファイルフォーマット ---
    parser.add_argument("--is_mergese_format", action='store_true',
                        help="Specify if the annotation file is in 'merge_SE.tsv' format (chr_start_end in the first column).")

    # --- オプション引数: ファイル/サンプル選択 ---
    parser.add_argument("--filename_pattern", default="*.bam",
                        help="Wildcard pattern for BAM filenames.")
    parser.add_argument("--sample_delimiter", default=None,
                        help="Delimiter string to extract sample names from BAM filenames.")

    # --- オプション引数: ツールパスとスレッド数 ---
    parser.add_argument("--samtools_path", default="samtools", help="Path to samtools.")
    parser.add_argument("--featurecounts_path", default="featureCounts", help="Path to featureCounts.")
    parser.add_argument("-T", "--threads", type=int, default=4, help="Threads for samtools/featureCounts.")
    parser.add_argument("--single_end", action='store_true', help="Run in single-end mode (default: paired-end).")

    # --- オプション引数: featureCounts関連 ---
    parser.add_argument("--fc_basename", default=None,
                        help="Base name for featureCounts output/log files (default: ANNOTATION_fc).")
    parser.add_argument("--fc_options", nargs='+', default=None,
                        help="Additional options for featureCounts (e.g., --fc_options --minOverlap 10).")

    # --- オプション引数: CPM計算関連（通常方式） ---
    cpm_group = parser.add_argument_group('Regular CPM calculation options (used when --calcnorm is not specified)')
    cpm_group.add_argument("--add_peakid", action='store_true', help="Replace original IDs with sequential PeakIDs.")
    cpm_group.add_argument("--id_column", default="Geneid", help="ID column name in SAF/featureCounts input.")
    cpm_group.add_argument("--output_id_name", default="PeakID", help="ID column name in the final logCPM table.")
    cpm_group.add_argument("--log_base", type=int, default=2, help="Log base (2 or 10). <=0 for no log.")
    cpm_group.add_argument("--pseudocount", type=float, default=1.0, help="Pseudocount for log transform.")

    # --- 高度な正規化計算オプション（新規追加）---
    calcnorm_group = parser.add_argument_group('Advanced normalization options')
    calcnorm_group.add_argument("--calcnorm", action='store_true', 
                        help="Use advanced statistical normalization instead of the default CPM calculation method.")
    calcnorm_group.add_argument("--min_cpm", type=float, default=1.0,
                        help="Minimum CPM threshold for filtering (calcnorm mode only).")
    calcnorm_group.add_argument("--min_samples", type=int, default=0,
                        help="Minimum number of samples with CPM > threshold. Default 0 means no filtering (calcnorm mode only).")
    calcnorm_group.add_argument("--calcnorm_method", default="upperquartile",
                        choices=["upperquartile", "TMM", "RLE", "none"],
                        help="Normalization method for calcnorm (calcnorm mode only).")
    calcnorm_group.add_argument("--remove_extensions", action="store_true",
                        help="Remove common BAM file extensions (.bam, .sorted.bam, etc) when using calcnorm.")
    calcnorm_group.add_argument("--pattern_rules", default=None,
                        help="JSON file with pattern rules for additional sample name cleaning (calcnorm mode only).")
    calcnorm_group.add_argument("--full_metadata_output", nargs="?", const="", default=None, metavar="PATH",
                        help="Path to save complete output with all featureCounts metadata columns. "
                             "If specified without a value, defaults to '{output_dir}/{logcpm_base}_full_metadata.txt'.")
    # 新規追加：calcnorm用出力ファイル名オプション
    calcnorm_group.add_argument("--calcnorm_output_name", default="calcnorm.tsv",
                        help="Filename for the calcnorm-normalized CPM table when using --calcnorm. "
                             "Set to '' to skip saving. Defaults to 'calcnorm.tsv'.")

    # --- オプション引数: 出力ファイル名 ---
    parser.add_argument("--logcpm_output_name", default="logCPM.tsv",
                        help="Filename for the logCPM table.")

    # --- オプション引数: ログ設定 ---
    parser.add_argument("--log_level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], help="Logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log", help="Filename for script execution log. Set to '' to disable.")

    # --- 引数の解析と検証 ---
    args = parser.parse_args()
    
    # 高度な正規化計算関連オプションのチェック
    if args.calcnorm and not CALCNORM_AVAILABLE:
        print("Error: --calcnorm is specified but cpm_peakprep_calcnorm_utils.py is not available.", file=sys.stderr)
        print("Please ensure cpm_peakprep_calcnorm_utils.py exists and necessary libraries (rpy2, edgeR) are installed.", file=sys.stderr)
        sys.exit(1)
    
    return args

# --- 高度な正規化計算処理関数 ---
def process_with_calcnorm(
    counts_file_path: str,
    output_dir: str,
    logger: logging.Logger,
    args: argparse.Namespace,
    bam_to_sample_map: dict = None
) -> tuple:
    """
    高度な正規化計算を実行（edgeRを利用）
    
    Args:
        counts_file_path (str): featureCounts出力ファイルのパス
        output_dir (str): 出力ディレクトリ
        logger (logging.Logger): ロガーオブジェクト
        args (argparse.Namespace): コマンドライン引数
        bam_to_sample_map (dict): BAMパスからサンプル名へのマッピング
        
    Returns:
        tuple: (normalized_df, output_path, full_metadata_path)
    """
    logger.info("--- Using advanced normalization calculation ---")
    
    # 入力ファイルの検証を追加
    logger.info(f"Validating counts file: {counts_file_path}")
    if not validate_counts_file(counts_file_path, logger):
        logger.critical("Counts file validation failed")
        return None, None, None
    
    # デリミターの準備
    sample_delimiter = args.sample_delimiter
    if sample_delimiter:
        logger.info(f"Using sample_delimiter: {sample_delimiter}")
    
    # パターンルールの読み込み
    pattern_rules = load_pattern_rules(args.pattern_rules, logger, args.remove_extensions)
    
    # 一時ファイルのパスを一意に生成
    try:
        cleaned_output_fd, cleaned_output_path = tempfile.mkstemp(
            dir=output_dir, prefix="temp_", suffix="_cleaned_counts.txt")
        os.close(cleaned_output_fd)
        logger.debug(f"Created temporary file: {cleaned_output_path}")
    except Exception as e:
        logger.error(f"Failed to create temporary file: {e}")
        cleaned_output_path = os.path.join(output_dir, f"temp_cleaned_counts_{int(time.time())}.txt")
        logger.warning(f"Using fallback temporary file path: {cleaned_output_path}")
    
    # サンプル名クリーニングの実行
    logger.info(f"Cleaning sample names in {counts_file_path}...")
    if not clean_count_file_sample_names(counts_file_path, cleaned_output_path, logger, sample_delimiter, pattern_rules):
        logger.critical("Sample name cleaning failed")
        return None, None, None
    
    # クリーニングされたファイルの読み込み
    try:
        counts_df = pd.read_csv(cleaned_output_path, sep='\t', comment='#')
        logger.info(f"Loaded counts data with {len(counts_df)} features and {len(counts_df.columns) - 6} samples")
    except Exception as e:
        logger.critical(f"Failed to read cleaned counts file: {e}")
        return None, None, None
    
    # フィルタリング情報をログに記録
    if args.min_samples > 0:
        logger.info(f"Feature filtering is ENABLED (min_cpm={args.min_cpm}, min_samples={args.min_samples})")
    else:
        logger.info("Feature filtering is DISABLED (min_samples=0)")
    
    # 高度な正規化計算の適用（edgeR利用）
    logger.info(f"Applying advanced normalization with method: {args.calcnorm_method}...")
    normalized_df = apply_edger_normalization(
        counts_df, 
        args.min_cpm, 
        args.min_samples, 
        logger,
        method=args.calcnorm_method
    )
    
    if normalized_df is None:
        logger.critical("Advanced normalization calculation failed")
        return None, None, None
    
    # --- 出力ファイル名の設定 ---
    # 標準TSV出力 - ここを修正: args.calcnorm_output_nameを使用
    if args.calcnorm_output_name.lower().endswith('.tsv'):
        tsv_output_name = args.calcnorm_output_name
    else:
        base_name, ext = os.path.splitext(args.calcnorm_output_name)
        tsv_output_name = f"{base_name}.tsv"
        logger.info(f"Changed output file extension to .tsv: {tsv_output_name}")
    
    tsv_output_path = os.path.join(output_dir, tsv_output_name)
    
    # フルメタデータ出力（オプション）
    full_metadata_path = None
    if args.full_metadata_output is not None:
        if args.full_metadata_output == "":
            # 値なしで指定された場合はデフォルトのパスを生成
            base_name = os.path.splitext(tsv_output_name)[0]
            full_metadata_name = f"{base_name}_full_metadata.txt"
            full_metadata_path = os.path.join(output_dir, full_metadata_name)
        else:
            # パスが指定された場合はそれを使用
            full_metadata_path = args.full_metadata_output
            if not os.path.isabs(full_metadata_path):
                full_metadata_path = os.path.join(output_dir, full_metadata_path)
                
        logger.info(f"Will save complete data with full metadata to: {full_metadata_path}")
    
    # --- 結果の保存 ---
    # フルメタデータ形式の保存（オプション）
    if full_metadata_path:
        try:
            normalized_df.to_csv(full_metadata_path, sep='\t', index=False, float_format='%.6f')
            logger.info(f"Complete data with full metadata saved to: {full_metadata_path}")
        except Exception as e:
            logger.error(f"Failed to save full metadata output: {e}")
            logger.warning("Continuing with standard TSV output only")
    
    # calcnorm出力が空文字列の場合はスキップ
    if not args.calcnorm_output_name:
        logger.info("Skipping calcnorm output file (--calcnorm_output_name is empty)")
        return normalized_df, None, full_metadata_path
    
    # 標準TSV形式に変換 (output_id_nameパラメータを追加)
    tsv_df = convert_to_tsv_format(normalized_df, logger, output_id_name=args.output_id_name)
    
    # 標準形式の保存
    try:
        tsv_df.to_csv(tsv_output_path, sep='\t', index=False, float_format='%.6f')
        logger.info(f"Standard TSV format saved to: {tsv_output_path}")
    except Exception as e:
        logger.critical(f"Failed to save TSV output: {e}")
        return normalized_df, None, full_metadata_path
    
    # 一時ファイルの削除
    try:
        os.remove(cleaned_output_path)
        logger.debug(f"Removed temporary cleaned file: {cleaned_output_path}")
    except Exception as e:
        logger.warning(f"Failed to remove temporary file {cleaned_output_path}: {e}")
    
    return normalized_df, tsv_output_path, full_metadata_path

# --- メイン実行関数 ---
def main():
    """pipeline main function"""
    start_time = time.time()
    args = parse_arguments()

    # --- ロガー設定 (StreamHandlerを明示的に追加) ---
    log_level_str = args.log_level.upper()
    log_level = getattr(logging, log_level_str, logging.INFO)
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    # ルートロガーを取得し、既存ハンドラをクリア
    root_logger = logging.getLogger()
    # 現在のハンドラをリストのコピーに対してループ処理して安全に削除
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close() # ハンドラを閉じる

    # ルートロガー自体のレベルを設定
    root_logger.setLevel(log_level)

    # 新しいコンソールハンドラ(StreamHandler)を追加し、レベルを設定
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level) # ハンドラにもレベルを設定
    formatter = logging.Formatter(log_format, datefmt=date_format)
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)

    # このスクリプト用の（子）ロガー取得
    logger = logging.getLogger("PipelineMain")

    # 出力ディレクトリ作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Output directory: {args.output_dir}")
    except OSError as e:
        logger.critical(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1)

    # ファイルハンドラ設定
    if args.script_log_file:
        log_file_path = os.path.join(args.output_dir, args.script_log_file)
        try:
            fh = logging.FileHandler(log_file_path, mode='w')
            fh.setLevel(log_level)
            fh.setFormatter(formatter)
            root_logger.addHandler(fh)
            logger.info(f"Logging script execution to file: {log_file_path}")
        except Exception as e:
            # ファイル書き込みエラーが出ても処理は止めない
            logger.error(f"Failed to configure file logging to {log_file_path}: {e}")

    logger.debug(f"Parsed arguments: {vars(args)}")

    # --- 変数初期化 ---
    logcpm_df = None
    counts_file_path = None
    temp_saf_path = None # 一時SAFファイルのパス
    saf_file_to_use = None # featureCountsに実際に渡すSAFパス
    
    # 高度な正規化計算使用フラグの記録
    # この時点で args.calcnorm が True かつ CALCNORM_AVAILABLE が False の場合は
    # parse_arguments() 内で既に sys.exit(1) されているはずなので、
    # ここでは単純に args.calcnorm の値を使用
    using_calcnorm = args.calcnorm

    try: # メイン処理全体を try...finally で囲み、一時ファイルを確実に削除

        # --- アノテーションファイルの準備 ---
        logger.info("--- Preparing annotation file ---")
        annotation_file = args.annotation_file
        if not os.path.exists(annotation_file):
             logger.critical(f"Annotation file not found: {annotation_file}")
             sys.exit(1)

        # ファイル形式を判定し、必要ならSAFに変換
        if args.is_mergese_format:
            logger.info("Input annotation specified as merge_SE format. Converting to temporary SAF...")
            temp_saf_path = convert_mergese_to_saf(annotation_file, logger)
            if temp_saf_path is None: sys.exit(1)
            saf_file_to_use = temp_saf_path
        else:
            annot_lower = annotation_file.lower()
            if annot_lower.endswith(('.bed', '.bed.gz')):
                logger.info("Input annotation appears to be a BED file. Converting to temporary SAF...")
                temp_saf_path = convert_bed_to_saf(annotation_file, logger)
                if temp_saf_path is None: sys.exit(1)
                saf_file_to_use = temp_saf_path
            elif annot_lower.endswith('.saf'):
                 logger.info("Input annotation appears to be a SAF file. Using it directly.")
                 saf_file_to_use = annotation_file
            else:
                 logger.warning(f"Unknown extension for annotation file: {annotation_file}. Assuming it is in SAF format.")
                 saf_file_to_use = annotation_file

        # featureCounts出力のベース名を決定
        annot_base = os.path.splitext(os.path.basename(annotation_file))[0]
        logger.info(f"Using annotation for featureCounts: {saf_file_to_use}")
        
        # =============================
        # === BAMファイル処理 ===
        # =============================
        logger.info("--- Processing BAM files ---")
        try:
            # 1. Get sample data
            logger.info("Step 1: Getting sample info for BAM files...")
            sample_dict, bam_list = get_sample_data(
                args.bam_dir, logger, args.sample_delimiter, args.filename_pattern
            )
            if not bam_list: raise ValueError("No BAM files found.")
            bam_to_sample_map = invert_dictionary(sample_dict)
            logger.info(f"Found {len(sample_dict)} samples.")
            
            # 検出されたファイルリストをログ出力
            if bam_list:
                logger.debug(f"BAM files included ({len(bam_list)}):")
                for f_path in bam_list: logger.debug(f"  - {os.path.basename(f_path)}")
            
            # 検出されたファイルリストをファイル保存
            bam_list_path = os.path.join(args.output_dir, "bam_list.txt")
            try:
                with open(bam_list_path, 'w', encoding='utf-8') as f_out:
                    for bam_path in bam_list: f_out.write(bam_path + '\n')
                logger.info(f"Saved list of processed BAM files to: {bam_list_path}")
            except IOError as e:
                logger.warning(f"Could not write BAM list file {bam_list_path}: {e}")

            # 2. Run samtools flagstat
            logger.info("Step 2: Getting total mapped reads for BAM files...")
            flagstat_dir = os.path.join(args.output_dir, "flagstat")
            total_reads_map = run_samtools_flagstat(
                bam_list, flagstat_dir, logger, args.samtools_path, args.threads
            )
            logger.info(f"Extracted reads for {len(total_reads_map)} files.")
            if len(total_reads_map) != len(bam_list): 
                logger.warning("Could not get reads for all BAMs.")

            # 3. Run featureCounts
            logger.info("Step 3: Running featureCounts for BAM files...")
            fc_base = args.fc_basename if args.fc_basename else f"{annot_base}_fc"
            fc_output = os.path.join(args.output_dir, f"{fc_base}_featureCounts.txt")
            fc_log = os.path.join(args.output_dir, f"{fc_base}_featureCounts.log")
            counts_file_path = run_featurecounts(
                bam_list, saf_file_to_use, logger, fc_output, fc_log,
                args.threads, (not args.single_end), args.featurecounts_path, args.fc_options
            )
            if not counts_file_path: raise RuntimeError("featureCounts failed.")
            logger.info(f"featureCounts output: {counts_file_path}")
            
            # featureCountsの結果を読み込む
            try:
                counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=0)
                if counts_df.columns[0].startswith('#'): 
                    counts_df = pd.read_csv(counts_file_path, sep='\t', comment='#', header=1)
            except Exception as e: 
                raise RuntimeError(f"Failed to read counts file {counts_file_path}") from e

            # === 分岐: 高度な正規化計算を使用するかで処理を切り替え ===
            if using_calcnorm:
                # === 高度な正規化計算処理の実行 ===
                logger.info("Step 4: Calculating CPM values using advanced normalization...")
                
                # calcnorm結果を計算・保存
                calcnorm_df, calcnorm_output_path, full_metadata_path = process_with_calcnorm(
                    counts_file_path, args.output_dir, logger, args, bam_to_sample_map
                )
                
                if calcnorm_df is None:
                    logger.critical("Advanced normalization calculation failed. Try using the standard CPM calculation instead.")
                    return
                
                # 標準のlogCPMも計算・保存（両方出力するため）
                if args.logcpm_output_name:
                    logger.info("Step 5: Also calculating standard logCPM for comparison...")
                    log_base_val = args.log_base if args.log_base in [2, 10] else None
                    logcpm_df = calculate_logcpm(
                        counts_df=counts_df, 
                        total_reads_dict=total_reads_map, 
                        logger=logger,
                        bam_path_to_sample_name=bam_to_sample_map, 
                        add_peakid_column=args.add_peakid,
                        id_column_name=args.id_column, 
                        output_id_column_name=args.output_id_name,
                        log_transform_base=log_base_val, 
                        pseudocount=args.pseudocount
                    )
                    
                    if logcpm_df is None:
                        logger.warning("Standard logCPM calculation failed, but advanced normalization succeeded.")
                    else:
                        # 標準logCPMを保存
                        std_output_path = os.path.join(args.output_dir, args.logcpm_output_name)
                        logcpm_df.to_csv(std_output_path, sep='\t', index=False, float_format='%.4f')
                        logger.info(f"Saved standard logCPM table to: {std_output_path}")
                else:
                    logger.info("Skipping standard logCPM output (--logcpm_output_name is empty)")
                
                # 結果ファイルの記録
                logger.info("--- Results Summary ---")
                logger.info(f"featureCounts output: {counts_file_path}")
                
                if args.logcpm_output_name:
                    logger.info(f"Standard logCPM values: {os.path.join(args.output_dir, args.logcpm_output_name)}")
                
                if args.calcnorm_output_name:
                    logger.info(f"Advanced normalized CPM values: {calcnorm_output_path}")
                
                if full_metadata_path:
                    logger.info(f"Full metadata output: {full_metadata_path}")
            else:
                # === オリジナルのCPM計算処理 ===
                # 4. Calculate logCPM (original method)
                logger.info("Step 4: Calculating logCPM (standard method)...")

                log_base_val = args.log_base if args.log_base in [2, 10] else None
                logcpm_df = calculate_logcpm(
                    counts_df=counts_df, 
                    total_reads_dict=total_reads_map, 
                    logger=logger,
                    bam_path_to_sample_name=bam_to_sample_map, 
                    add_peakid_column=args.add_peakid,
                    id_column_name=args.id_column, 
                    output_id_column_name=args.output_id_name,
                    log_transform_base=log_base_val, 
                    pseudocount=args.pseudocount
                )
                if logcpm_df is None: raise RuntimeError("logCPM calculation failed.")
                logger.info("logCPM calculation completed.")

                # 5. Save logCPM results
                output_path = os.path.join(args.output_dir, args.logcpm_output_name)
                logcpm_df.to_csv(output_path, sep='\t', index=False, float_format='%.4f')
                logger.info(f"Saved logCPM table to: {output_path}")
            # === 分岐ここまで ===

        except Exception as e:
            # 処理中のエラーは致命的として終了
            logger.critical(f"Critical error during data processing: {e}", exc_info=True)
            sys.exit(1)

        # --- パイプライン終了 ---
        end_time = time.time()
        logger.info(f"Pipeline finished successfully in {end_time - start_time:.2f} seconds.")

    finally:
        # 一時ファイルを確実に削除 - ログレベル変更：INFO→DEBUG
        if temp_saf_path and os.path.exists(temp_saf_path):
            try:
                os.remove(temp_saf_path)
                logger.debug(f"Removed temporary annotation file: {temp_saf_path}")
            except Exception as e:
                logger.warning(f"Failed to remove temporary annotation file {temp_saf_path}: {e}")

# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    main()