#!/usr/bin/env python3
# === bigwig_peakprep_summary.py ===
# BigWig Summary and Processing Pipeline: 
# Runs multiBigwigSummary, processes the counts, applies log2 transformation,
# and adds PeakID information from the BED file.

import argparse
import logging
import os
import sys
import time
import pandas as pd
from typing import List, Dict, Optional, Tuple

# --- bigwig_peakprep_utils.py から関数をインポート ---
try:
    from bigwig_peakprep_utils import (
        run_multibigwigsummary,
        read_multibigwig_counts,
        dataframe_to_pyranges,
        sort_pyranges,
        pyranges_to_dataframe,
        apply_log2_transform_to_counts,
        read_bed_to_dataframe,
        add_peak_ids_from_bed,
        # 新規追加の変換関数
        convert_saf_to_bed,
        convert_mergese_to_bed
    )
except ImportError as e:
    print(f"Error importing functions from utility modules: {e}", file=sys.stderr)
    print("Please ensure 'bigwig_peakprep_utils.py' exists and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- 引数パーサーの設定 ---
def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ChIP-seq BigWig Summary Pipeline: Creates summaries from bigWig files and prepares normalized data tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数 ---
    parser.add_argument("-i", "--bigwig_details_file", required=True,
                        help="Path to the bam_to_bigwig_details.tsv file from bamconvert process.")
    parser.add_argument("-a", "--annotation_file", required=True,
                        help="Path to the annotation file (BED, SAF, or merge_SE format).")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files.")

    # --- オプション引数: ファイルフォーマット ---
    parser.add_argument("--annotation_format", default="auto",
                        choices=["auto", "bed", "saf", "mergese"],
                        help="Format of the annotation file. 'auto' detects based on extension.")
    parser.add_argument("--is_mergese_format", action='store_true',
                        help="Specify if the annotation file is in 'merge_SE.tsv' format (chr_start_end in the first column).")

    # --- オプション引数: ツールパス ---
    parser.add_argument("--tools_dir", default=None,
                        help="Directory containing deeptools executables. If provided, will look for multiBigwigSummary in this directory.")
    parser.add_argument("--multibigwigsummary_path", default="multiBigwigSummary",
                        help="Path to multiBigwigSummary executable. Ignored if --tools_dir is specified.")

    # --- オプション引数: multiBigwigSummary関連 ---
    parser.add_argument("--summary_dir", default=None,
                        help="Directory to save multiBigwigSummary files. Defaults to output_dir/summary.")
    parser.add_argument("--summary_basename", default="bigwig_summary",
                        help="Base name for multiBigwigSummary output files.")
    parser.add_argument("--multibigwigsummary_options", nargs='+', default=None,
                        help="Additional options for multiBigwigSummary as a list of strings.")

    # --- オプション引数: データテーブル関連 ---
    parser.add_argument("--force_chr_start_end_ids", action='store_true',
                        help="Force using chr_start_end format as PeakID even if BED file has name column. Default is to use name column if available.")
    parser.add_argument("--pseudocount", type=float, default=1.0,
                        help="Pseudocount for log2 transform.")
    parser.add_argument("--raw_counts_name", default="multibigwig_counts.tsv",
                        help="Filename for the raw counts table.")
    parser.add_argument("--log2_counts_name", default="multibigwig_counts_log2.tsv",
                        help="Filename for the log2-transformed counts table.")

    # --- オプション引数: スレッド数 ---
    parser.add_argument("-T", "--threads", type=int, default=4,
                        help="Number of threads to use for deeptools commands.")

    # --- オプション引数: ログ設定 ---
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for script execution log. Set to '' to disable.")


    # log変換済み対応
    parser.add_argument("--skip_log2_transform", action='store_true',
                        help="Skip log2 transformation (for pre-transformed data like bamCompare output)")

    parser.add_argument("--skip_log2_file", action='store_true',
                  help="Skip generation of log2-transformed counts file")

    return parser.parse_args()

# --- bigwig_details_file から bigWig ファイルとサンプル情報を読み込む関数 ---
# def load_bigwig_details(
#     bigwig_details_file: str,
#     logger: logging.Logger
# ) -> Tuple[List[str], List[str]]:
#     """
#     Load bigWig file paths and sample names from the
#     bam_to_bigwig_details.tsv file.
    
#     Args:
#         bigwig_details_file (str): Path to the bigWig details file
#         logger (logging.Logger): Logger object
        
#     Returns:
#         Tuple[List[str], List[str]]: 
#             - List of bigWig file paths
#             - List of sample names (in the same order as bigWig files)
#     """
#     logger.info(f"Reading bigWig details file: {bigwig_details_file}")
    
#     bigwig_files = []
#     sample_names = []
    
#     try:
#         # ファイルの存在確認
#         if not os.path.exists(bigwig_details_file):
#             logger.error(f"BigWig details file not found: {bigwig_details_file}")
#             return [], []
        
#         # TSVファイルを読み込み（#で始まる行はコメント扱い）
#         df = pd.read_csv(bigwig_details_file, sep='\t', comment='#')
        
#         # 必要なカラムが存在するか確認
#         required_columns = ['Sample_name', 'BigWig_fullpath']
#         missing_columns = [col for col in required_columns if col not in df.columns]
#         if missing_columns:
#             logger.error(f"Required columns missing in BigWig details file: {missing_columns}")
#             return [], []
        
#         # bigWigファイルパスとサンプル名を抽出
#         bigwig_files = df['BigWig_fullpath'].tolist()
#         sample_names = df['Sample_name'].tolist()
        
#         if not bigwig_files:
#             logger.warning("No bigWig files found in details file.")
#             return [], []
        
#         logger.info(f"Loaded {len(bigwig_files)} bigWig files and their sample names")
#         logger.debug(f"Sample names: {', '.join(sample_names[:5])}{'...' if len(sample_names) > 5 else ''}")
#         return bigwig_files, sample_names
        
#     except Exception as e:
#         logger.error(f"Error reading BigWig details file: {e}", exc_info=True)
#         return [], []

#メタデータ対応の為に差し替え
def load_bigwig_details(
    bigwig_details_file: str,
    logger: logging.Logger
) -> Tuple[List[str], List[str], bool]:  # 戻り値にbool型を追加
    """
    Load bigWig file paths and sample names from the
    bam_to_bigwig_details.tsv file.
    
    Args:
        bigwig_details_file (str): Path to the bigWig details file
        logger (logging.Logger): Logger object
        
    Returns:
        Tuple[List[str], List[str], bool]: 
            - List of bigWig file paths
            - List of sample names (in the same order as bigWig files)
            - Boolean indicating if data is already log2 transformed
    """
    logger.info(f"Reading bigWig details file: {bigwig_details_file}")
    
    bigwig_files = []
    sample_names = []
    already_log2_transformed = False  # デフォルトはfalse
    
    try:
        # ファイルの存在確認
        if not os.path.exists(bigwig_details_file):
            logger.error(f"BigWig details file not found: {bigwig_details_file}")
            return [], [], False
        
        # メタデータの確認（ファイルの先頭の#コメント行）
        with open(bigwig_details_file, 'r', encoding='utf-8') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                if 'already_log2_transformed: true' in line.lower():
                    already_log2_transformed = True
                    logger.info("Detected already log2 transformed data from file metadata")
                elif 'bamcompare operation: log2' in line.lower():
                    # bamCompareの操作がlog2の場合も変換済みと判断
                    already_log2_transformed = True
                    logger.info("Detected log2 operation from bamCompare metadata")
        
        # TSVファイルを読み込み（#で始まる行はコメント扱い）
        df = pd.read_csv(bigwig_details_file, sep='\t', comment='#')
        
        # 必要なカラムが存在するか確認
        required_columns = ['Sample_name', 'BigWig_fullpath']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            logger.error(f"Required columns missing in BigWig details file: {missing_columns}")
            return [], [], False
        
        # bigWigファイルパスとサンプル名を抽出
        bigwig_files = df['BigWig_fullpath'].tolist()
        sample_names = df['Sample_name'].tolist()
        
        if not bigwig_files:
            logger.warning("No bigWig files found in details file.")
            return [], [], False
        
        logger.info(f"Loaded {len(bigwig_files)} bigWig files and their sample names")
        logger.debug(f"Sample names: {', '.join(sample_names[:5])}{'...' if len(sample_names) > 5 else ''}")
        if already_log2_transformed:
            logger.info("BigWig files are already log2 transformed, will skip transformation during processing")
        
        return bigwig_files, sample_names, already_log2_transformed
        
    except Exception as e:
        logger.error(f"Error reading BigWig details file: {e}", exc_info=True)
        return [], [], False

def main():
    """
    Main function to generate bigWig summaries and create peak count tables.
    Returns the paths of the generated result files.
    """
    start_time = time.time()
    args = parse_arguments()

    # --- ロガー設定 ---
    log_level_str = args.log_level.upper()
    log_level = getattr(logging, log_level_str, logging.INFO)
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'

    # ルートロガーを取得し、既存ハンドラをクリア
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    # ルートロガー自体のレベルを設定
    root_logger.setLevel(log_level)

    # 新しいコンソールハンドラを追加
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    formatter = logging.Formatter(log_format, datefmt=date_format)
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)

    # このスクリプト用のロガー取得
    logger = logging.getLogger("BigWigSummary")

    # 出力ディレクトリ作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.debug(f"Created or verified output directory: {args.output_dir}")
    except OSError as e:
        logger.critical(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1)

    # ファイルログハンドラ設定
    if args.script_log_file:
        log_file_path = os.path.join(args.output_dir, args.script_log_file)
        
        # 既存のログファイルがあれば、タイムスタンプ付きでリネーム
        if os.path.exists(log_file_path):
            timestamp = time.strftime("%Y%m%d-%H%M%S")
            backup_log_path = f"{log_file_path}.{timestamp}.bak"
            try:
                os.rename(log_file_path, backup_log_path)
                logger.debug(f"Existing log file renamed to: {os.path.basename(backup_log_path)}")
            except Exception as e:
                logger.warning(f"Failed to rename existing log file {log_file_path}: {e}")
        
        try:
            fh = logging.FileHandler(log_file_path, mode='w')
            fh.setLevel(log_level)
            fh.setFormatter(formatter)
            root_logger.addHandler(fh)
            logger.debug(f"Logging to file: {log_file_path}")
        except Exception as e:
            logger.error(f"Failed to configure file logging to {log_file_path}: {e}")

    logger.debug(f"Parsed arguments: {vars(args)}")

    # --- ツールパスの決定 ---
    if args.tools_dir:
        multibigwigsummary_path = os.path.join(args.tools_dir, "multiBigwigSummary")
        logger.debug(f"Using multiBigwigSummary from directory: {args.tools_dir}")
    else:
        multibigwigsummary_path = args.multibigwigsummary_path
        logger.debug(f"Using specified multiBigwigSummary path: {multibigwigsummary_path}")

    # --- ディレクトリパスの設定 ---
    summary_dir = args.summary_dir if args.summary_dir else os.path.join(args.output_dir, "summary")
    
    try:
        os.makedirs(summary_dir, exist_ok=True)
        logger.debug(f"Created or verified summary directory: {summary_dir}")
    except OSError as e:
        logger.critical(f"Failed to create summary directory {summary_dir}: {e}")
        sys.exit(1)

    # --- アノテーションファイルの準備 ---
    logger.info("--- Preparing annotation file ---")
    annotation_file = args.annotation_file
    if not os.path.exists(annotation_file):
        logger.critical(f"Annotation file not found: {annotation_file}")
        sys.exit(1)

    # 変換されたBEDファイルのパス（クリーンアップ不要）
    bed_file_to_use = None
    is_temporary_file = False  # 変換ファイルが一時的かどうかの判定用

    # ファイル形式の決定
    annotation_format = args.annotation_format
    if args.is_mergese_format:
        annotation_format = "mergese" 

    if annotation_format == "auto":
        # 拡張子で自動判定
        annot_lower = annotation_file.lower()
        if annot_lower.endswith(('.bed', '.bed.gz')):
            annotation_format = "bed"
        elif annot_lower.endswith('.saf'):
            annotation_format = "saf"
        else:
            logger.warning(f"Could not determine format from extension of {annotation_file}. Assuming BED format.")
            annotation_format = "bed"

    # 形式に応じた処理
    if annotation_format == "bed":
        logger.info("Using BED format annotation file directly.")
        bed_file_to_use = annotation_file
    elif annotation_format == "saf":
        logger.info("Converting SAF format to BED file.")
        # 出力ディレクトリを指定して変換（ログとして保存）
        bed_file_to_use = convert_saf_to_bed(annotation_file, logger, args.output_dir)
        if bed_file_to_use is None:
            logger.critical("Failed to convert SAF file to BED. Cannot proceed.")
            sys.exit(1)
    elif annotation_format == "mergese":
        logger.info("Converting merge_SE format to BED file.")
        # 出力ディレクトリを指定して変換（ログとして保存）
        bed_file_to_use = convert_mergese_to_bed(annotation_file, logger, args.output_dir)
        if bed_file_to_use is None:
            logger.critical("Failed to convert merge_SE file to BED. Cannot proceed.")
            sys.exit(1)
    else:
        logger.critical(f"Unsupported annotation format: {annotation_format}")
        sys.exit(1)

    logger.info(f"Using annotation file for multiBigwigSummary: {bed_file_to_use}")

    # --- BEDファイルの読み込み ---
    logger.info(f"Reading annotation BED file: {bed_file_to_use}")
    bed_df = read_bed_to_dataframe(bed_file_to_use, logger)
    if bed_df is None:
        logger.critical(f"Failed to read BED file: {bed_file_to_use}")
        sys.exit(1)
    logger.info(f"Loaded BED file with {len(bed_df)} regions")

    # --- bigwig_details_file から情報を取得 ---
    bigwig_files, sample_names, already_log2_transformed = load_bigwig_details(args.bigwig_details_file, logger)
    if not bigwig_files:
        logger.critical("No bigWig files found in details file. Cannot proceed.")
        sys.exit(1)

    # --- multiBigwigSummary の実行 ---
    logger.info("Running multiBigwigSummary with BED regions")
    logger.debug(f"Summary output directory: {summary_dir}")
    logger.debug(f"Summary basename: {args.summary_basename}")
    
    summary_basename = args.summary_basename
    summary_npz = os.path.join(summary_dir, f"{summary_basename}.npz")
    summary_tab = os.path.join(summary_dir, f"{summary_basename}.tab")
    summary_log = os.path.join(summary_dir, f"{summary_basename}.log")
    
    # multiBigwigSummary の実行 - ラベル指定
    summary_result = run_multibigwigsummary(
        bigwig_files=bigwig_files,
        bed_file=bed_file_to_use,  # 変換後のBEDファイルを使用
        output_dir=summary_dir,
        logger=logger,
        output_filename=summary_basename,
        labels=sample_names,  # サンプル名をラベルとして使用
        multibigwigsummary_path=multibigwigsummary_path,
        threads=args.threads,
        log_file=summary_log,
        additional_options=args.multibigwigsummary_options
    )
    
    if not summary_result:
        logger.critical("Failed to generate multiBigwigSummary. Cannot proceed.")
        sys.exit(1)
    
    logger.info(f"MultiBigwigSummary completed successfully")
    logger.debug(f"Summary files: {summary_npz}, {summary_tab}")
    
    # --- multiBigwigSummary の読み込みと処理 ---
    logger.info("Processing multiBigwigSummary data")
    
    # tabファイルの読み込み
    logger.debug(f"Reading multiBigwigSummary tab file: {summary_tab}")
    multibigwig_df = read_multibigwig_counts(summary_tab, logger)
    if multibigwig_df is None:
        logger.critical("Failed to read multiBigwigSummary tab file")
        sys.exit(1)
    
    # PyRangesを使ったソート
    logger.debug("Converting to PyRanges and sorting")
    pyranges_obj = dataframe_to_pyranges(multibigwig_df, logger)
    if pyranges_obj is None:
        logger.critical("Failed to convert DataFrame to PyRanges")
        sys.exit(1)
    
    sorted_pyranges = sort_pyranges(pyranges_obj, logger)
    if sorted_pyranges is None:
        logger.critical("Failed to sort PyRanges object")
        sys.exit(1)
    
    sorted_df = pyranges_to_dataframe(sorted_pyranges, logger)
    if sorted_df is None:
        logger.critical("Failed to convert sorted PyRanges back to DataFrame")
        sys.exit(1)
    
    logger.debug(f"Successfully sorted multiBigwigSummary data with {len(sorted_df)} regions")
    
    # コマンドラインオプションでのスキップが優先
    skip_log2 = args.skip_log2_transform or already_log2_transformed
    
    if skip_log2:
        logger.info("Skipping log2 transformation as data is already transformed or skip flag is set")
        log2_df = sorted_df.copy()  # 変換せずにコピーするだけ
    else:
        # Log2変換の適用
        logger.info("Applying log2 transformation to counts")
        log2_df = apply_log2_transform_to_counts(
            sorted_df, logger, pseudocount=args.pseudocount
        )
        if log2_df is None:
            logger.critical("Failed to apply log2 transformation")
            sys.exit(1)
    
    logger.debug(f"Data preparation completed" + (" (log2 transformation skipped)" if skip_log2 else f" with pseudocount: {args.pseudocount}"))
    
    # --- PeakID情報の追加 ---
    logger.info("Adding PeakID information from BED file")
    
    # BEDファイルのname列が存在するか確認
    has_name_column = 'name' in bed_df.columns and not bed_df['name'].isnull().all()
    if has_name_column:
        logger.debug("BED file has 'name' column (4th column), will use as PeakID by default")
        if args.force_chr_start_end_ids:
            logger.debug("--force_chr_start_end_ids specified, will override BED name column and use chr_start_end format")
    else:
        logger.debug("BED file does not have valid 'name' column, will use chr_start_end format for PeakIDs")
    
    # 通常のカウントデータにPeakID追加
    counts_with_peaks = add_peak_ids_from_bed(
        sorted_df, bed_df, logger, 
        force_chr_start_end_ids=args.force_chr_start_end_ids
    )
    if counts_with_peaks is None:
        logger.critical("Failed to add PeakID information to raw counts")
        sys.exit(1)
    
    # Log2変換したデータにPeakID追加
    log2_counts_with_peaks = add_peak_ids_from_bed(
        log2_df, bed_df, logger,
        force_chr_start_end_ids=args.force_chr_start_end_ids
    )
    if log2_counts_with_peaks is None:
        logger.critical("Failed to add PeakID information to log2 counts")
        sys.exit(1)
    
    logger.debug("Successfully added PeakID information to count tables")
    
    # --- 結果の保存 ---
    logger.info("Saving result tables")
    
    # 通常のカウントデータを保存
    raw_counts_path = os.path.join(args.output_dir, args.raw_counts_name)
    counts_with_peaks.to_csv(raw_counts_path, sep='\t', index=False, float_format='%.6f')
    logger.info(f"Raw counts saved to: {raw_counts_path}")
    
    # # Log2変換したデータを保存
    # log2_counts_path = os.path.join(args.output_dir, args.log2_counts_name)
    # log2_counts_with_peaks.to_csv(log2_counts_path, sep='\t', index=False, float_format='%.6f')
    # # Log2変換済みデータかどうかをログに表示
    # if skip_log2:
    #     logger.info(f"Pre-transformed data saved to: {log2_counts_path} (log2 transformation was skipped)")
    # else:
    #     logger.info(f"Log2 counts saved to: {log2_counts_path}")


    # Log2変換したデータを保存（skip_log2_fileが設定されていない場合のみ）
    if not args.skip_log2_file:
        log2_counts_path = os.path.join(args.output_dir, args.log2_counts_name)
        log2_counts_with_peaks.to_csv(log2_counts_path, sep='\t', index=False, float_format='%.6f')
        # Log2変換済みデータかどうかをログに表示
        if skip_log2:
            logger.info(f"Pre-transformed data saved to: {log2_counts_path} (log2 transformation was skipped)")
        else:
            logger.info(f"Log2 counts saved to: {log2_counts_path}")
    else:
        logger.info("Skipping generation of log2-transformed counts file as requested")
        log2_counts_path = None  # 生成していないのでNoneを設定

    # --- パイプラインの終了 ---
    end_time = time.time()
    logger.info(f"BigWig summary processing finished in {end_time - start_time:.2f} seconds")
    
    # 変換されたファイルはログとして保存されるため、クリーンアップは不要
    # 変換ファイルを含む conversion_logs ディレクトリの場所を表示
    if annotation_format != "bed" and bed_file_to_use:
        logger.info(f"Converted annotation file was saved to: {bed_file_to_use}")
        logger.info(f"The directory contains all conversion logs for reference.")
    
    # 生成されたファイルパスをタプルで返す
    return raw_counts_path, log2_counts_path


# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    try:
        raw_counts_file, log2_counts_file = main()
        print(f"BigWig summary processing completed successfully.")
        print(f"Raw counts file: {raw_counts_file}")
        print(f"Log2 counts file: {log2_counts_file}")
    except Exception as e:
        print(f"Error during execution: {e}", file=sys.stderr)
        sys.exit(1)