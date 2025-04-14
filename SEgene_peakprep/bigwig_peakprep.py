#!/usr/bin/env python3
# === bigwig_peakprep.py ===
# Wrapper script for the ChIP-seq BigWig Processing Pipeline
# This script calls the separate step scripts:
# 1. bigwig_peakprep_bamconvert.py - BAM to bigWig conversion
# 2. bigwig_peakprep_summary.py - multiBigwigSummary and downstream analysis

import argparse
import logging
import os
import sys
import time
import subprocess
from typing import List, Optional, Tuple

# --- bigwig_peakprep_utils.py から関数をインポート ---
try:
    from bigwig_peakprep_utils import (
        run_bamcoverage,
        run_multibigwigsummary,
        read_multibigwig_counts,
        dataframe_to_pyranges,
        sort_pyranges,
        pyranges_to_dataframe,
        apply_log2_transform_to_counts,
        read_bed_to_dataframe,
        add_peak_ids_from_bed,
        _get_bamcoverage_version,
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
        description="ChIP-seq BigWig Processing Pipeline: Generates bigWig files, creates summaries, and prepares normalized data tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数 ---
    parser.add_argument("-b", "--bam_dir", required=False,
                        help="Directory containing BAM files. Required for BAM to bigWig conversion.")
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
                        help="Directory containing deeptools executables. If provided, will look for bamCoverage and multiBigwigSummary in this directory.")
    parser.add_argument("--bamcoverage_path", default="bamCoverage",
                        help="Path to bamCoverage executable. Ignored if --tools_dir is specified.")
    parser.add_argument("--multibigwigsummary_path", default="multiBigwigSummary",
                        help="Path to multiBigwigSummary executable. Ignored if --tools_dir is specified.")

    # --- オプション引数: ファイル/サンプル選択 ---
    parser.add_argument("--filename_pattern", default="*.bam",
                        help="Wildcard pattern for BAM filenames.")
    parser.add_argument("--sample_delimiter", default=None, metavar="DELIMITER",
                        help="Delimiter string to extract sample names from BAM filenames.")

    # --- オプション引数: bamCoverage関連 ---
    parser.add_argument("--bigwig_dir", default=None,
                        help="Directory to save bigWig files. Defaults to output_dir/bigwig.")
    parser.add_argument("--bin_size", type=int, default=50,
                        help="Size of the bins in bases for bamCoverage.")
    parser.add_argument("--normalize_using", default="RPGC",
                        help="Normalization method for bamCoverage (e.g., RPGC, CPM).")
    parser.add_argument("--effective_genome_size", type=int, default=2913022398,
                        help="Effective genome size for RPGC normalization. Default is human genome.")
    parser.add_argument("--bamcoverage_options", nargs='+', default=None,
                        help="Additional options for bamCoverage as a list of strings.")

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
    parser.add_argument("--force_overwrite", action='store_true',
                        help="Force overwriting existing bigWig files. By default, existing files are skipped.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for script execution log. Set to '' to disable.")

    # --- オプション引数: ステップ設定 ---
    parser.add_argument("--run_bamconvert_only", action='store_true',
                        help="Run only BAM to bigWig conversion.")
    parser.add_argument("--run_summary_only", action='store_true',
                        help="Run only bigWig summary and downstream analysis.")
    parser.add_argument("--bigwig_details_file", default=None,
                        help="Path to bam_to_bigwig_details.tsv file. Required if --run_summary_only is specified.")
    
    # --- オプション引数: スクリプトパス ---
    parser.add_argument("--script_dir", default=None,
                        help="Directory containing the step scripts. If not specified, uses the same directory as this script.")

    return parser.parse_args()

# --- 引数を分離してステップスクリプト用のコマンドラインを構築する関数 ---
def build_bamconvert_command(args, script_path: str) -> List[str]:
    """
    Build command line for BAM to bigWig conversion
    
    Args:
        args: Parsed arguments
        script_path (str): Path to the bamconvert script
        
    Returns:
        List[str]: Command line argument list
    """
    # ベースコマンド
    cmd = [sys.executable, script_path]
    
    # 必須引数 - イコール渡し形式に変更
    cmd.append(f"--bam_dir={args.bam_dir}")
    cmd.append(f"--output_dir={args.output_dir}")
    
    # オプション引数: ツールパス - イコール渡し形式に変更
    if args.tools_dir:
        cmd.append(f"--tools_dir={args.tools_dir}")
    else:
        cmd.append(f"--bamcoverage_path={args.bamcoverage_path}")
    
    # オプション引数: ファイル/サンプル選択 - イコール渡し形式に変更
    if args.filename_pattern != "*.bam":
        cmd.append(f"--filename_pattern={args.filename_pattern}")
    if args.sample_delimiter:
        cmd.append(f"--sample_delimiter={args.sample_delimiter}")
    
    # オプション引数: bamCoverage関連 - イコール渡し形式に変更
    if args.bigwig_dir:
        cmd.append(f"--bigwig_dir={args.bigwig_dir}")
    if args.bin_size != 50:
        cmd.append(f"--bin_size={args.bin_size}")
    if args.normalize_using != "RPGC":
        cmd.append(f"--normalize_using={args.normalize_using}")
    if args.effective_genome_size != 2913022398:
        cmd.append(f"--effective_genome_size={args.effective_genome_size}")
    
    # bamcoverage_optionsは複数値なのでそのまま追加
    if args.bamcoverage_options:
        cmd.append("--bamcoverage_options")
        cmd.extend(args.bamcoverage_options)
    
    # オプション引数: スレッド数 - イコール渡し形式に変更
    if args.threads != 4:
        cmd.append(f"--threads={args.threads}")
    
    # オプション引数: ログ設定 - イコール渡し形式に変更
    if args.log_level != "INFO":
        cmd.append(f"--log_level={args.log_level}")
    if args.force_overwrite:
        cmd.append("--force_overwrite")
    if args.script_log_file != "pipeline_run.log":
        cmd.append(f"--script_log_file={args.script_log_file}")
    
    return cmd

def build_summary_command(args, script_path: str, bigwig_details_file: str) -> List[str]:
    """
    Build command line for bigWig summary processing
    
    Args:
        args: Parsed arguments
        script_path (str): Path to the summary script
        bigwig_details_file (str): Path to the bigWig details file
        
    Returns:
        List[str]: Command line argument list
    """
    # ベースコマンド
    cmd = [sys.executable, script_path]
    
    # 必須引数 - イコール渡し形式に変更
    cmd.append(f"--bigwig_details_file={bigwig_details_file}")
    cmd.append(f"--annotation_file={args.annotation_file}")
    cmd.append(f"--output_dir={args.output_dir}")
    
    # オプション引数: ツールパス - イコール渡し形式に変更
    if args.tools_dir:
        cmd.append(f"--tools_dir={args.tools_dir}")
    else:
        cmd.append(f"--multibigwigsummary_path={args.multibigwigsummary_path}")
    
    # オプション引数: multiBigwigSummary関連 - イコール渡し形式に変更
    if args.summary_dir:
        cmd.append(f"--summary_dir={args.summary_dir}")
    if args.summary_basename != "bigwig_summary":
        cmd.append(f"--summary_basename={args.summary_basename}")
    
    # multibigwigsummary_optionsは複数値なのでそのまま追加
    if args.multibigwigsummary_options:
        cmd.append("--multibigwigsummary_options")
        cmd.extend(args.multibigwigsummary_options)
    
    # オプション引数: データテーブル関連 - イコール渡し形式に変更
    if args.force_chr_start_end_ids:
        cmd.append("--force_chr_start_end_ids")
    if args.pseudocount != 1.0:
        cmd.append(f"--pseudocount={args.pseudocount}")
    if args.raw_counts_name != "multibigwig_counts.tsv":
        cmd.append(f"--raw_counts_name={args.raw_counts_name}")
    if args.log2_counts_name != "multibigwig_counts_log2.tsv":
        cmd.append(f"--log2_counts_name={args.log2_counts_name}")
    
    # オプション引数: スレッド数 - イコール渡し形式に変更
    if args.threads != 4:
        cmd.append(f"--threads={args.threads}")
    
    # オプション引数: ログ設定 - イコール渡し形式に変更
    if args.log_level != "INFO":
        cmd.append(f"--log_level={args.log_level}")
    if args.script_log_file != "pipeline_run.log":
        cmd.append(f"--script_log_file={args.script_log_file}")
    
    return cmd

# --- サブプロセスを実行し、結果を処理する関数 ---
def run_subprocess(cmd: List[str], step_name: str, logger: logging.Logger) -> Tuple[bool, Optional[str]]:
    """
    Execute a command as a subprocess and process the results
    
    Args:
        cmd (List[str]): List of command and arguments to execute
        step_name (str): Step name (for logging)
        logger (logging.Logger): Logger object
        
    Returns:
        Tuple[bool, Optional[str]]:
            - Whether execution was successful
            - Last line of standard output (typically the output file path)
    """
    logger.debug(f"Running {step_name} with command: {' '.join(cmd)}")
    
    try:
        # サブプロセスを実行し、標準出力と標準エラーをキャプチャ
        process = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8'
        )
        
        # 標準出力と標準エラーをログに出力
        important_info = []
        for line in process.stdout.strip().split('\n'):
            # 重要な情報を抽出してINFOで表示
            if any(marker in line for marker in ["completed successfully", "Output file:", "saved to:"]):
                important_info.append(line)
                logger.info(f"{step_name}: {line}")
            else:
                logger.debug(f"{step_name} stdout: {line}")
        
        for line in process.stderr.strip().split('\n'):
            if line:  # 空行はスキップ
                if "warning" in line.lower():
                    logger.warning(f"{step_name}: {line}")
                else:
                    logger.debug(f"{step_name} stderr: {line}")
        
        # 標準出力の最後の行を取得（通常は出力ファイルのパス）
        last_line = process.stdout.strip().split('\n')[-1] if process.stdout.strip() else None
        
        logger.info(f"{step_name} completed successfully")
        return True, last_line
        
    except subprocess.CalledProcessError as e:
        logger.error(f"{step_name} failed with return code {e.returncode}")
        if e.stdout:
            for line in e.stdout.strip().split('\n'):
                if any(marker in line for marker in ["error", "Error", "failed", "Failed"]):
                    logger.error(f"{step_name}: {line}")
                else:
                    logger.debug(f"{step_name} stdout: {line}")
        if e.stderr:
            for line in e.stderr.strip().split('\n'):
                logger.error(f"{step_name}: {line}")
        return False, None
    
    except Exception as e:
        logger.error(f"Error running {step_name}: {e}", exc_info=True)
        return False, None

# --- メイン実行関数 ---
def main():
    """Wrapper pipeline main function"""
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
    logger = logging.getLogger("BigWigPipeline")

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
    
    # --- スクリプトパスの設定 ---
    # スクリプトディレクトリが指定されていなければ、このスクリプトと同じディレクトリを使用
    if args.script_dir:
        script_dir = args.script_dir
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    
    bamconvert_script = os.path.join(script_dir, "bigwig_peakprep_bamconvert.py")
    summary_script = os.path.join(script_dir, "bigwig_peakprep_summary.py")
    
    # スクリプトの存在確認
    if not os.path.exists(bamconvert_script):
        logger.critical(f"BAM convert script not found: {bamconvert_script}")
        sys.exit(1)
    
    if not os.path.exists(summary_script):
        logger.critical(f"Summary script not found: {summary_script}")
        sys.exit(1)
    
    logger.debug(f"Using BAM convert script: {bamconvert_script}")
    logger.debug(f"Using summary script: {summary_script}")
    
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
        annotation_format = "mergese"  # --is_mergese_format フラグが指定されていれば優先

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

    logger.info(f"Using annotation file for processing: {bed_file_to_use}")
    
    # --- 実行するステップの決定 ---
    run_bamconvert = not args.run_summary_only
    run_summary = not args.run_bamconvert_only
    
    if args.run_bamconvert_only and args.run_summary_only:
        logger.critical("Cannot specify both --run_bamconvert_only and --run_summary_only at the same time.")
        sys.exit(1)
    
    if args.run_summary_only and not args.bigwig_details_file:
        logger.critical("--bigwig_details_file must be specified when --run_summary_only is used.")
        sys.exit(1)
    
    if run_bamconvert and not args.bam_dir:
        logger.critical("--bam_dir is required for BAM to bigWig conversion.")
        sys.exit(1)
    
    # --- BAMからbigWigへの変換 ---
    bigwig_details_file = None
    
    if run_bamconvert:
        logger.info("Starting BAM to bigWig conversion")
        
        # bamconvertのコマンドを構築
        bamconvert_cmd = build_bamconvert_command(args, bamconvert_script)
        
        # bamconvertを実行
        bamconvert_success, bamconvert_last_line = run_subprocess(bamconvert_cmd, "BAM convert", logger)
        
        if not bamconvert_success:
            logger.critical("BAM to bigWig conversion failed. Cannot proceed to summary processing.")
            if run_summary:
                logger.info("Exiting without running summary processing due to BAM convert failure.")
                sys.exit(1)
            else:
                logger.info("Only BAM conversion was requested, but it failed.")
                sys.exit(1)
        
        # 出力ファイルのパスを取得
        # 最後の行が出力ファイルのパスのはず
        if bamconvert_last_line and "Output file:" in bamconvert_last_line:
            bigwig_details_file = bamconvert_last_line.split("Output file:")[-1].strip()
        else:
            # 出力ファイルが見つからない場合はデフォルトのパスを使用
            bigwig_details_file = os.path.join(args.output_dir, "bam_to_bigwig_details.tsv")
            if os.path.exists(bigwig_details_file):
                logger.warning(f"Could not find bigWig details file path from stdout. Using default path: {bigwig_details_file}")
            else:
                logger.critical(f"bigWig details file not found at default path: {bigwig_details_file}")
                sys.exit(1)
        
        logger.info(f"BAM to bigWig conversion completed. Output file: {bigwig_details_file}")
    else:
        bigwig_details_file = args.bigwig_details_file
        logger.info(f"Skipping BAM to bigWig conversion as requested. Using provided details file: {bigwig_details_file}")
    
    # --- bigWigサマリーとピークカウントテーブルの作成 ---
    if run_summary:
        logger.info("Starting bigWig summary processing")
        
        # summaryのコマンドを構築 - bed_file_to_use を渡す
        summary_cmd = build_summary_command(args, summary_script, bigwig_details_file)
        
        # bed_file_to_use が annotation_file と異なる場合、コマンドラインの annotation_file を置き換える
        if bed_file_to_use != annotation_file:
            # "-a" または "--annotation_file=" で始まる引数を探して置き換え
            for i, arg in enumerate(summary_cmd):
                if arg.startswith("-a=") or arg.startswith("--annotation_file="):
                    summary_cmd[i] = f"{arg.split('=')[0]}={bed_file_to_use}"
                    break
                elif (arg == "-a" or arg == "--annotation_file") and i + 1 < len(summary_cmd):
                    summary_cmd[i + 1] = bed_file_to_use
                    break
        
        # 変換元の形式情報を渡す
        if annotation_format != "bed":
            summary_cmd.append(f"--annotation_format={annotation_format}")
        
        # summaryを実行
        summary_success, summary_last_line = run_subprocess(summary_cmd, "Summary", logger)
        
        if not summary_success:
            logger.critical("bigWig summary processing failed.")
            sys.exit(1)
        
        # 処理成功のログ出力
        logger.info("bigWig summary processing completed successfully")
        
    else:
        logger.info("Skipping bigWig summary processing as requested.")
    
    # --- パイプラインの終了 ---
    end_time = time.time()
    logger.info(f"Pipeline finished successfully in {end_time - start_time:.2f} seconds")

    # 変換されたファイルはログとして保存されるため、finally ブロックでの削除は不要
    # 変換ファイルを含む conversion_logs ディレクトリの場所を表示
    if annotation_format != "bed" and bed_file_to_use:
        logger.info(f"Converted annotation file was saved to: {bed_file_to_use}")
        logger.info(f"The directory contains all conversion logs for reference.")

# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    main()