#!/usr/bin/env python3
# === bigwig_peakprep_bamdiff.py ===
# ChIP-seq BAM Differential Analysis and Processing Pipeline:
# Generates log2 ratio bigWig files from sample/control BAM files and processes summary counts.

import argparse
import logging
import os
import sys
import time
import subprocess
import yaml
from typing import List, Dict, Optional, Tuple, Any

# --- 必要な関数をbigwig_peakprep_bamdiff_utils.pyからインポート ---
try:
    from bigwig_peakprep_bamdiff_utils import (
        run_bamcompare,
        scan_bam_files,
        prepare_yaml,
        execute_conversion_plan,
        save_bamdiff_results
    )
    # cpm_peakprep_utils からサンプル処理関数をインポート
    from cpm_peakprep_utils import get_sample_data, invert_dictionary
    from bigwig_peakprep_bamdiff_generatediff import create_sample_control_pairs
except ImportError as e:
    print(f"Error importing functions from utility modules: {e}", file=sys.stderr)
    print("Please ensure 'bigwig_peakprep_bamdiff_utils.py' and 'cpm_peakprep_utils.py' exist and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ChIP-seq BAM Differential Analysis Pipeline: Complete workflow from BAM to summary tables.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数: ファイル検索条件 ---
    parser.add_argument("--sample_bam_dir", required=True,
                        help="Directory containing sample BAM files")
    parser.add_argument("--control_bam_dir", required=True,
                        help="Directory containing control BAM files (Input/IgG)")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files (bigWig files and logs)")
    
    # --- オプション引数: ファイル検索条件 ---
    parser.add_argument("--sample_pattern", default="*.bam",
                        help="Glob pattern for sample BAM files")
    parser.add_argument("--control_pattern", default="*Input*.bam",
                        help="Glob pattern for control BAM files")
    parser.add_argument("--sample_delimiter", default="_",
                        help="Delimiter for sample name extraction")
                        
    # --- オプション引数: ツールパス ---
    parser.add_argument("--bamcompare_path", default="bamCompare",
                        help="Path to bamCompare executable.")
    parser.add_argument("--tools_dir", default=None,
                        help="Directory containing deeptools executables. If provided, bamCompare will be looked for in this directory.")
    
    # --- オプション引数: bamCompare関連 ---
    parser.add_argument("--bin_size", type=int, default=50,
                        help="Size of the bins in bases for bamCompare.")
    parser.add_argument("--operation", default="log2",
                        choices=["log2", "ratio", "subtract", "add", "mean", "reciprocal_ratio", "first", "second"],
                        help="Operation to perform on the two BAM files.")
    parser.add_argument("--pseudocount", type=float, default=1.0,
                        help="Pseudocount to add to avoid division by zero. Used in log2 and ratio operations.")
    parser.add_argument("--effective_genome_size", type=int, default=2913022398,
                        help="Effective genome size. Default is human genome.")
    parser.add_argument("--bamcompare_options", nargs='+', default=None,
                        help="Additional options for bamCompare as a list of strings.")
                        
    # --- オプション引数: スレッド数 ---
    parser.add_argument("-T", "--threads", type=int, default=4,
                        help="Number of threads to use for deeptools commands.")
                        
    # --- オプション引数: 実行制御 ---
    parser.add_argument("--force_overwrite", action='store_true',
                        help="Force overwriting existing bigWig files. By default, existing files are skipped.")

    # --- オプション引数: ログ設定 ---
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for script execution log. Set to '' to disable.")

    # --- サマリー関連の引数 ---
    parser.add_argument("-a", "--annotation_file", 
                        help="Path to the annotation file (BED, SAF, or merge_SE format) for summary processing.")
    parser.add_argument("--run_diff_only", action='store_true',
                        help="Run only BAM to differential bigWig conversion, skip summary processing.")
    parser.add_argument("--bigwig_peakprep_path", default=None,
                        help="Path to bigwig_peakprep.py script. If not specified, looks in the same directory.")
    parser.add_argument("--summary_script_path", default=None,
                        help="Path to bigwig_peakprep_summary.py script. If not specified, looks in the same directory.")

    # --- summary.py のオプション ---
    parser.add_argument("--force_chr_start_end_ids", action='store_true',
                        help="Force using chr_start_end format as PeakID even if BED file has name column.")
    parser.add_argument("--annotation_format", default="auto",
                        choices=["auto", "bed", "saf", "mergese"],
                        help="Format of the annotation file. 'auto' detects based on extension.")
    parser.add_argument("--raw_counts_name", default="multibigwig_counts.tsv",
                        help="Filename for the raw counts table.")
    parser.add_argument("--log2_counts_name", default="multibigwig_counts_log2.tsv",
                        help="Filename for the log2-transformed counts table.")
    
    return parser.parse_args()

def setup_logger(log_level_str: str, output_dir: str, script_log_file: Optional[str]) -> logging.Logger:
    """
    Set up the logger for the script
    
    Args:
        log_level_str (str): String log level (DEBUG, INFO, etc.)
        output_dir (str): Output directory for log file
        script_log_file (Optional[str]): Log file name, None to disable file logging
        
    Returns:
        logging.Logger: Configured logger
    """
    log_level = getattr(logging, log_level_str.upper(), logging.INFO)
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
    logger = logging.getLogger("BAMDiffAnalyzer")
    
    # ファイルログハンドラ設定（必要な場合）
    if script_log_file:
        log_file_path = os.path.join(output_dir, script_log_file)
        
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
            
    return logger

def main():
    """Main function for the BAM differential analysis pipeline."""
    start_time = time.time()
    args = parse_arguments()
    
    # 出力ディレクトリの作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory: {e}", file=sys.stderr)
        sys.exit(1)
    
    # ロガー設定
    logger = setup_logger(args.log_level, args.output_dir, args.script_log_file)
    
    # ツールパスの設定
    bamcompare_path = args.bamcompare_path
    if args.tools_dir:
        bamcompare_path = os.path.join(args.tools_dir, "bamCompare")
        logger.debug(f"Using bamCompare from tools directory: {bamcompare_path}")
    
    # --- BAM差分解析ステップの実行 ---
    logger.info("Starting BAM differential analysis pipeline")
    
    # bigWigディレクトリの作成
    bigwig_dir = os.path.join(args.output_dir, "bigwig")
    try:
        os.makedirs(bigwig_dir, exist_ok=True)
        logger.debug(f"Created or verified bigWig directory: {bigwig_dir}")
    except OSError as e:
        logger.critical(f"Failed to create bigWig directory {bigwig_dir}: {e}")
        sys.exit(1)
    
    # サンプルとコントロールのBAMファイルをスキャン
    sample_to_bam, control_to_bam = scan_bam_files(
        sample_dir=args.sample_bam_dir,
        control_dir=args.control_bam_dir,
        logger=logger,
        sample_pattern=args.sample_pattern,
        control_pattern=args.control_pattern,
        sample_delimiter=args.sample_delimiter
    )
    
    if not sample_to_bam:
        logger.critical("No sample BAM files found. Cannot proceed.")
        sys.exit(1)
        
    if not control_to_bam:
        logger.critical("No control BAM files found. Cannot proceed.")
        sys.exit(1)
    
    # サンプルとコントロールのペア化・マッピング

    
    sample_control_pairs = create_sample_control_pairs(
        sample_to_bam=sample_to_bam,
        control_to_bam=control_to_bam,
        logger=logger
    )
    
    if not sample_control_pairs:
        logger.critical("Failed to create sample-control pairs. Cannot proceed.")
        sys.exit(1)
    
    # bamCompareパラメータ設定
    bamcompare_params = {
        "operation": args.operation,
        "bin_size": args.bin_size,
        "pseudocount": args.pseudocount,
        "effective_genome_size": args.effective_genome_size,
        "threads": args.threads
    }
    
    # 変換プランの作成
    yaml_path = os.path.join(args.output_dir, "conversion_plan.yaml")
    
    yaml_file = prepare_yaml(
        sample_bams=sample_to_bam,
        control_bams=control_to_bam,
        sample_control_pairs=sample_control_pairs,
        output_dir=bigwig_dir,
        bamcompare_params=bamcompare_params,
        output_yaml_path=yaml_path,
        logger=logger
    )
    
    if not yaml_file:
        logger.critical("Failed to create conversion plan. Cannot proceed.")
        sys.exit(1)
    
    # 変換プランの実行
    results = execute_conversion_plan(
        yaml_path=yaml_file,
        logger=logger,
        output_dir=bigwig_dir,
        bamcompare_path=bamcompare_path,
        force_overwrite=args.force_overwrite,
        additional_options=args.bamcompare_options,
        write_results_yaml=True
    )
    
    # 結果の保存
    tsv_path, yaml_path = save_bamdiff_results(
        results_yaml=results,
        output_dir=args.output_dir,
        logger=logger,
        output_basename="bamdiff_to_bigwig_details"
    )
    
    # サマリー処理をスキップするかチェック
    if args.run_diff_only:
        logger.info("Running only differential analysis as requested (--run_diff_only)")
        if tsv_path:
            logger.info(f"Results saved to TSV: {tsv_path}")
            logger.info(f"Results saved to YAML: {yaml_path}")
            logger.info(f"BAM differential analysis pipeline completed in {time.time() - start_time:.2f} seconds")
            return 0
        else:
            logger.error("Failed to save results")
            return 1
    
    # --- サマリー処理を実行 ---
    if not args.annotation_file:
        logger.error("Annotation file (--annotation_file) is required for summary processing.")
        logger.error("Either provide annotation file or use --run_diff_only to skip summary processing.")
        return 1
    
    # サマリースクリプトのパスを決定
    if args.summary_script_path:
        summary_script = args.summary_script_path
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        summary_script = os.path.join(script_dir, "bigwig_peakprep_summary.py")
    
    # サマリースクリプトの存在確認
    if not os.path.exists(summary_script):
        logger.error(f"Summary script not found at: {summary_script}")
        logger.error("Please specify correct path with --summary_script_path")
        return 1
    
    logger.info(f"Running summary processing using: {summary_script}")
    
    # サマリースクリプトを呼び出すコマンドを構築
    cmd = [
        sys.executable,
        summary_script,
#        "--run_summary_only",
        "--skip_log2_transform",  # 重要: log2変換をスキップ
        "--skip_log2_file",  # log2ファイルの生成をスキップ
        f"--bigwig_details_file={tsv_path}",
        f"--annotation_file={args.annotation_file}",
        f"--output_dir={args.output_dir}",
        f"--threads={args.threads}",
        f"--log_level={args.log_level}"
    ]

    # tools_dir引数を渡す
    if args.tools_dir:
        cmd.append(f"--tools_dir={args.tools_dir}")
    
    # オプションフラグを条件付きで追加
    if args.force_chr_start_end_ids:
        cmd.append("--force_chr_start_end_ids")
    
    if args.annotation_format != "auto":
        cmd.append(f"--annotation_format={args.annotation_format}")
    
    if args.raw_counts_name != "multibigwig_counts.tsv":
        cmd.append(f"--raw_counts_name={args.raw_counts_name}")
    
    if args.log2_counts_name != "multibigwig_counts_log2.tsv":
        cmd.append(f"--log2_counts_name={args.log2_counts_name}")
    
    # コマンド実行をログに記録
    logger.info(f"Running summary processing with command: {' '.join(cmd)}")
    
    try:
        # サブプロセスとしてbigwig_peakprep_summary.pyを実行
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8'
        )
        
        # 標準出力と標準エラーをログに記録
        for line in result.stdout.strip().split('\n'):
            if line:  # 空行はスキップ
                logger.info(f"Summary: {line}")
            
        for line in result.stderr.strip().split('\n'):
            if line:  # 空行はスキップ
                logger.warning(f"Summary stderr: {line}")
        
        logger.info("Summary processing completed successfully")
        
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running summary processing: exit code {e.returncode}")
        if e.stdout:
            for line in e.stdout.strip().split('\n'):
                if line:
                    logger.info(f"Summary stdout: {line}")
        if e.stderr:
            for line in e.stderr.strip().split('\n'):
                if line:
                    logger.error(f"Summary stderr: {line}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error running summary processing: {e}")
        return 1
    
    # --- パイプライン完了メッセージ ---
    end_time = time.time()
    logger.info(f"Complete BAM differential analysis pipeline finished in {end_time - start_time:.2f} seconds")
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)