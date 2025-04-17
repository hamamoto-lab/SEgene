#!/usr/bin/env python3
# === bigwig_peakprep_bamdiff_generatediff.py ===
# ChIP-seq BAM Differential Analysis Pipeline: 
# Generates log2 ratio bigWig files from sample and control BAM files.

import argparse
import logging
import os
import sys
import time
import yaml
import glob
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
    from cpm_peakprep_utils import get_sample_data
except ImportError as e:
    print(f"Error importing functions from utility modules: {e}", file=sys.stderr)
    print("Please ensure 'bigwig_peakprep_bamdiff_utils.py' and 'cpm_peakprep_utils.py' exist and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- 引数parser の設定 ---
def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ChIP-seq BAM Differential Analysis Pipeline: Generates log2 ratio bigWig files from sample and control BAM files.",
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
                        
    # --- オプション引数: tool path ---
    parser.add_argument("--bamcompare_path", default="bamCompare",
                        help="Path to bamCompare executable.")
    
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
                        
    # --- オプション引数: thread 数 ---
    parser.add_argument("-T", "--threads", type=int, default=4,
                        help="Number of threads to use for deeptools commands.")
                        
    # --- オプション引数: 実行制御 ---
    parser.add_argument("--force_overwrite", action='store_true',
                        help="Force overwriting existing bigWig files. By default, existing files are skipped.")

    # --- オプション引数: log 設定 ---
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for script execution log. Set to '' to disable.")

    return parser.parse_args()

# --- Logger設定関数 ---
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

    # Root loggerを取得し、既存handlerをクリア
    root_logger = logging.getLogger()
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()

    # Root logger自体のレベルを設定
    root_logger.setLevel(log_level)

    # 新しいコンソールhandlerを追加
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    formatter = logging.Formatter(log_format, datefmt=date_format)
    ch.setFormatter(formatter)
    root_logger.addHandler(ch)
    
    # このスクリプト用のlogger取得
    logger = logging.getLogger("BAMDiffAnalyzer")
    
    # File log handler設定（必要な場合）
    if script_log_file:
        log_file_path = os.path.join(output_dir, script_log_file)
        
        # 既存のlog fileがあれば、タイムスタンプ付きでリネーム
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

# --- サンプルとコントロールをペア化する関数（ファイル名ベース） ---
def create_sample_control_pairs(
    sample_to_bam: Dict[str, str],
    control_to_bam: Dict[str, str],
    logger: logging.Logger
) -> Dict[str, Tuple[str, str, str]]:
    """
    Pairs samples with controls and creates conversion pairs for bamCompare
    Uses the original filename (without extension) for the output filename
    
    Args:
        sample_to_bam: Mapping from sample names to BAM paths
        control_to_bam: Mapping from control names to BAM paths
        logger: Logger
        
    Returns:
        Mapping from pair IDs to (sample BAM path, control BAM path, output filename)
    """
    pairs = {}
    
    # サンプルとコントロールの元ファイル名を取得
    sample_path_to_filename = {path: os.path.splitext(os.path.basename(path))[0] for path in sample_to_bam.values()}
    control_path_to_filename = {path: os.path.splitext(os.path.basename(path))[0] for path in control_to_bam.values()}
    
    # サンプル名->パス->ファイル名のマッピングも作成
    sample_name_to_filename = {name: sample_path_to_filename[path] for name, path in sample_to_bam.items()}
    control_name_to_filename = {name: control_path_to_filename[path] for name, path in control_to_bam.items()}
    
    # 単一コントロールの場合
    if len(control_to_bam) == 1:
        control_name, control_bam = next(iter(control_to_bam.items()))
        control_filename = control_path_to_filename[control_bam]
        logger.info(f"Using single control '{control_name}' for all samples")
        
        for sample_name, sample_bam in sample_to_bam.items():
            sample_filename = sample_path_to_filename[sample_bam]
            
            # デリミタ未適用のファイル名でペアIDを作成
            pair_id = f"{sample_name}_vs_{control_name}"  # ログ表示用のペアID
            output_name = f"{sample_filename}_vs_{control_filename}.log2ratio.bigwig"
            
            pairs[pair_id] = (sample_bam, control_bam, output_name)
            logger.debug(f"Created pair: {sample_name} -> {control_name}")
    
    # 同じ接頭辞のファイル名でマッチング（例: T1-H3K27ac.bam と T1-Input.bam）
    else:
        logger.info("Matching samples and controls by common prefix")
        for sample_name, sample_bam in sample_to_bam.items():
            sample_filename = sample_path_to_filename[sample_bam]
            
            # サンプル名から接頭辞を抽出（例: "T1-H3K27ac" -> "T1"）
            sample_prefix = sample_name.split('-')[0] if '-' in sample_name else sample_name
            
            # マッチするコントロールを探す
            for control_name, control_bam in control_to_bam.items():
                control_filename = control_path_to_filename[control_bam]
                control_prefix = control_name.split('-')[0] if '-' in control_name else control_name
                
                # 接頭辞が一致すればペア作成
                if sample_prefix and sample_prefix == control_prefix:
                    pair_id = f"{sample_name}_vs_{control_name}"  # ログ表示用のペアID
                    output_name = f"{sample_filename}_vs_{control_filename}.log2ratio.bigwig"
                    
                    pairs[pair_id] = (sample_bam, control_bam, output_name)
                    logger.debug(f"Matched by prefix: {sample_name} -> {control_name}")
                    break
    
    # マッチングできなかったサンプルの処理
    unmatched_samples = []
    for sample_name, sample_bam in sample_to_bam.items():
        if not any(sample_bam == sbam for _, (sbam, _, _) in pairs.items()):
            unmatched_samples.append((sample_name, sample_bam))
    
    if unmatched_samples:
        # コントロールが1つもマッチしなかったサンプルがある場合
        if control_to_bam:
            # 最初のコントロールを使用
            fallback_control_name, fallback_control_bam = next(iter(control_to_bam.items()))
            fallback_control_filename = control_path_to_filename[fallback_control_bam]
            logger.warning(f"Using '{fallback_control_name}' as fallback control for {len(unmatched_samples)} unmatched samples")
            
            for sample_name, sample_bam in unmatched_samples:
                sample_filename = sample_path_to_filename[sample_bam]
                
                pair_id = f"{sample_name}_vs_{fallback_control_name}"  # ログ表示用のペアID
                output_name = f"{sample_filename}_vs_{fallback_control_filename}.log2ratio.bigwig"
                
                pairs[pair_id] = (sample_bam, fallback_control_bam, output_name)
                logger.debug(f"Using fallback control for: {sample_name}")
    
    if not pairs:
        logger.error("Could not create any sample-control pairs")
    else:
        logger.info(f"Created {len(pairs)} sample-control pairs")
    
    return pairs

# --- メイン実行関数 ---
def main():
    """
    Main function for the BAM differential analysis pipeline.
    """
    start_time = time.time()
    args = parse_arguments()
    
    # 出力ディレクトリの作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory: {e}", file=sys.stderr)
        sys.exit(1)
    
    # --- Logger設定 ---
    logger = setup_logger(args.log_level, args.output_dir, args.script_log_file)
    
    # デバッグ用に引数を出力
    logger.debug(f"Command line arguments: {vars(args)}")
    
    # --- bigwigディレクトリの作成 ---
    bigwig_dir = os.path.join(args.output_dir, "bigwig")
    try:
        os.makedirs(bigwig_dir, exist_ok=True)
        logger.debug(f"Created or verified bigWig directory: {bigwig_dir}")
    except OSError as e:
        logger.critical(f"Failed to create bigWig directory {bigwig_dir}: {e}")
        sys.exit(1)
    
    # --- BAMファイルのスキャン ---
    logger.info("Starting BAM differential analysis pipeline")
    
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
    
    # --- サンプルとコントロールのペア化 ---
    sample_control_pairs = create_sample_control_pairs(
        sample_to_bam=sample_to_bam,
        control_to_bam=control_to_bam,
        logger=logger
    )
    
    if not sample_control_pairs:
        logger.critical("Failed to create sample-control pairs. Cannot proceed.")
        sys.exit(1)
    
    # --- bamcompareパラメータ設定 ---
    bamcompare_params = {
        "operation": args.operation,
        "bin_size": args.bin_size,
        "pseudocount": args.pseudocount,
        "effective_genome_size": args.effective_genome_size,
        "threads": args.threads
    }
    
    # --- 変換プランの作成 ---
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
    
    # --- 変換プランの実行 ---
    results = execute_conversion_plan(
        yaml_path=yaml_file,
        logger=logger,
        output_dir=bigwig_dir,
        bamcompare_path=args.bamcompare_path,
        force_overwrite=args.force_overwrite,
        additional_options=args.bamcompare_options,
        write_results_yaml=True
    )
    
    # --- 結果の保存 ---
    tsv_path, yaml_path = save_bamdiff_results(
        results_yaml=results,
        output_dir=args.output_dir,
        logger=logger,
        output_basename="bamdiff_to_bigwig_details"
    )
    
    # --- パイプライン完了メッセージ ---
    end_time = time.time()
    logger.info(f"BAM differential analysis pipeline completed in {end_time - start_time:.2f} seconds")
    
    if tsv_path:
        logger.info(f"Results saved to TSV: {tsv_path}")
        logger.info(f"Results saved to YAML: {yaml_path}")
        return 0
    else:
        logger.error("Failed to save results")
        return 1

# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)