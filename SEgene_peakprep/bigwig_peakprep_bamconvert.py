#!/usr/bin/env python3
# === bigwig_peakprep_bamconvert.py ===
# BAM to BigWig Conversion Pipeline: 
# Generates bigWig files from BAM files and creates detailed output for downstream analysis.

import argparse
import logging
import os
import sys
import time
import glob
import pandas as pd
from typing import List, Dict, Optional, Tuple

# --- bigwig_peakprep_utils.py と cpm_peakprep_utils.py から関数をインポート ---
try:
    # bigwig_peakprep_utils からの関数インポート
    from bigwig_peakprep_utils import (
        run_bamcoverage,
        _get_bamcoverage_version
    )
    # cpm_peakprep_utils からサンプル処理関数をインポート
    from cpm_peakprep_utils import get_sample_data, invert_dictionary
except ImportError as e:
    print(f"Error importing functions from utility modules: {e}", file=sys.stderr)
    print("Please ensure 'bigwig_peakprep_utils.py' and 'cpm_peakprep_utils.py' exist and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- 引数parser の設定 ---
def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="ChIP-seq BAM to BigWig Conversion Pipeline: Generates bigWig files from BAM files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # --- 必須引数 ---
    parser.add_argument("-b", "--bam_dir", required=True,
                        help="Directory containing BAM files.")
    # アノテーションファイルは当面使用しないのでコメントアウト
    # parser.add_argument("-a", "--annotation_file", required=True,
    #                     help="Path to the BED format annotation file (required for metadata).")
    parser.add_argument("-o", "--output_dir", required=True,
                        help="Directory to save all output files.")

    # --- オプション引数: tool path ---
    parser.add_argument("--tools_dir", default=None,
                        help="Directory containing deeptools executables. If provided, will look for bamCoverage in this directory.")
    parser.add_argument("--bamcoverage_path", default="bamCoverage",
                        help="Path to bamCoverage executable. Ignored if --tools_dir is specified.")

    # --- オプション引数: file/sample 選択 ---
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

    # --- オプション引数: thread 数 ---
    parser.add_argument("-T", "--threads", type=int, default=4,
                        help="Number of threads to use for deeptools commands.")

    # --- オプション引数: log 設定 ---
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level.")
    parser.add_argument("--force_overwrite", action='store_true',
                        help="Force overwriting existing bigWig files. By default, existing files are skipped.")
    parser.add_argument("--script_log_file", default="pipeline_run.log",
                        help="Filename for script execution log. Set to '' to disable.")

    return parser.parse_args()

# --- BAM file をスキャンして情報を抽出する関数 ---
def scan_bam_files(
    bam_dir: str,
    logger: logging.Logger,
    filename_pattern: str = "*.bam",
    sample_delimiter: Optional[str] = None
) -> Tuple[List[str], Dict[str, str], Dict[str, str]]:
    """
    Scan BAM files from directory using get_sample_data function,
    and retrieve the correspondence between paths and sample names.
    
    Args:
        bam_dir (str): Directory containing BAM files
        logger (logging.Logger): Logger object
        filename_pattern (str): Filename pattern (e.g., "*.bam")
        sample_delimiter (Optional[str]): Delimiter for extracting sample names
        
    Returns:
        Tuple[List[str], Dict[str, str], Dict[str, str]]: 
            - List of BAM file paths
            - Dictionary mapping paths to sample names
            - Dictionary mapping sample names to paths
    """
    logger.info(f"Scanning BAM files in directory: {bam_dir}")
    logger.debug(f"Using pattern '{filename_pattern}' and delimiter: {sample_delimiter}")
    
    try:
        # get_sample_data関数を使用してsample dataを取得
        sample_to_bam, bam_files = get_sample_data(
            bam_folder=bam_dir,
            logger=logger,
            sample_name_delimiter=sample_delimiter,
            filename_pattern=filename_pattern
        )
        
        if not bam_files:
            logger.warning(f"No BAM files found matching pattern '{filename_pattern}' in {bam_dir}")
            return [], {}, {}
        
        logger.info(f"Found {len(bam_files)} BAM files with {len(sample_to_bam)} distinct samples")
        
        # 逆dictionary（BAM path→sample name）も作成
        bam_to_sample = invert_dictionary(sample_to_bam)
        
        # デバッグ情報
        for i, (sample, bam_path) in enumerate(sample_to_bam.items()):
            if i < 5:  # 最初の5つのsampleだけlog出力
                logger.debug(f"Sample: {sample} -> BAM: {os.path.basename(bam_path)}")
        
        return bam_files, bam_to_sample, sample_to_bam
        
    except Exception as e:
        logger.error(f"Error scanning BAM files: {e}", exc_info=True)
        return [], {}, {}

# --- BAM fileからbigWig fileを生成する関数 ---
def generate_bigwig_files(
    bam_files: List[str],
    bam_to_sample: Dict[str, str],
    bigwig_dir: str,
    logger: logging.Logger,
    bamcoverage_path: str,
    bin_size: int,
    normalize_using: str,
    effective_genome_size: int,
    threads: int,
    force_overwrite: bool,
    bamcoverage_options: Optional[List[str]] = None
) -> Tuple[List[str], Dict[str, str], List[str], List[str]]:
    """
    Generate bigWig files from BAM files
    
    Args:
        bam_files (List[str]): List of BAM file paths
        bam_to_sample (Dict[str, str]): Dictionary mapping BAM paths to sample names
        bigwig_dir (str): Directory to save bigWig files
        logger (logging.Logger): Logger object
        bamcoverage_path (str): Path to bamCoverage executable
        bin_size (int): Bin size
        normalize_using (str): Normalization method
        effective_genome_size (int): Effective genome size
        threads (int): Number of threads to use
        force_overwrite (bool): Whether to overwrite existing files
        bamcoverage_options (Optional[List[str]]): Additional options for bamCoverage
    
    Returns:
        Tuple[List[str], Dict[str, str], List[str], List[str]]:
            - List of generated bigWig files
            - Dictionary mapping BAM paths to bigWig paths
            - List of skipped files
            - List of failed BAM files
    """
    logger.info(f"Starting bigWig generation for {len(bam_files)} BAM files")
    logger.debug(f"BigWig output directory: {bigwig_dir}")
    logger.debug(f"Using bamCoverage path: {bamcoverage_path}")
    logger.debug(f"Normalization method: {normalize_using}, bin size: {bin_size}, genome size: {effective_genome_size}")
    
    # 出力directory作成
    try:
        os.makedirs(bigwig_dir, exist_ok=True)
        logger.debug(f"Created or verified bigWig directory: {bigwig_dir}")
    except OSError as e:
        logger.critical(f"Failed to create bigWig directory {bigwig_dir}: {e}")
        return [], {}, [], []
    
    bigwig_files = []
    failed_bams = []
    # BAM pathからBigWig pathへのmapping dictionary
    bam_to_bigwig_map = {}
    
    # Mapping statusの追跡情報
    skipped_files = []
    generated_files = []
    
    # 進捗counter
    total_files = len(bam_files)
    processed_count = 0
    
    for bam_path in bam_files:
        processed_count += 1
        bam_filename = os.path.basename(bam_path)
        sample_name = bam_to_sample.get(bam_path, os.path.splitext(bam_filename)[0])
        output_filename = f"{sample_name}.{normalize_using}.bigwig"
        output_path = os.path.join(bigwig_dir, output_filename)
        log_file = os.path.join(bigwig_dir, f"{sample_name}.bamCoverage.log")
        
        # 進捗状況の表示（10ファイルごと）
        if processed_count % 10 == 0 or processed_count == total_files:
            logger.info(f"Processing BAM files: {processed_count}/{total_files}")
        
        # 既存fileの確認
        if os.path.exists(output_path) and not force_overwrite:
            logger.debug(f"Skipping existing bigWig file for {sample_name}: {output_filename}")
            bigwig_files.append(output_path)
            bam_to_bigwig_map[bam_path] = output_path
            skipped_files.append(output_path)
            continue
        
        logger.debug(f"Processing BAM file: {bam_filename} (Sample: {sample_name})")
        
        bigwig_path = run_bamcoverage(
            bam_file=bam_path,
            output_dir=bigwig_dir,
            logger=logger,
            output_filename=output_filename,
            bin_size=bin_size,
            normalize_using=normalize_using,
            effective_genome_size=effective_genome_size,
            bamcoverage_path=bamcoverage_path,
            threads=threads,
            log_file=log_file,
            additional_options=bamcoverage_options
        )
        
        if bigwig_path:
            bigwig_files.append(bigwig_path)
            # BAM pathとBigWig pathの対応を保存
            bam_to_bigwig_map[bam_path] = bigwig_path
            generated_files.append(bigwig_path)
            logger.debug(f"Successfully generated bigWig for {sample_name}: {os.path.basename(bigwig_path)}")
        else:
            failed_bams.append(bam_path)
            logger.error(f"Failed to generate bigWig for {sample_name}")
    
    # 生成結果の確認
    if not bigwig_files:
        logger.critical("Failed to generate any bigWig files.")
        return [], {}, [], []
    
    # 処理結果のsummaryを表示
    logger.info(f"BigWig generation completed: {len(generated_files)} created, {len(skipped_files)} skipped, {len(failed_bams)} failed")
    
    return bigwig_files, bam_to_bigwig_map, skipped_files, failed_bams

# --- Mapping関係の詳細情報をTSVとして保存する関数 ---
def save_bam_to_bigwig_mapping(
    bam_to_bigwig_map: Dict[str, str],
    bam_to_sample: Dict[str, str],
    output_dir: str,
    skipped_files: List[str],
    failed_bams: List[str],
    logger: logging.Logger
) -> bool:
    """
    Save the mapping relationship between BAM and BigWig files
    
    Args:
        bam_to_bigwig_map (Dict[str, str]): Dictionary mapping BAM paths to BigWig paths
        bam_to_sample (Dict[str, str]): Dictionary mapping BAM paths to sample names
        output_dir (str): Output directory
        skipped_files (List[str]): List of skipped files
        failed_bams (List[str]): List of failed BAM files
        logger (logging.Logger): Logger object
    
    Returns:
        bool: Whether the save was successful
    """
    mapping_file_path = os.path.join(output_dir, "bam_to_bigwig_mapping.txt")
    try:
        with open(mapping_file_path, 'w', encoding='utf-8') as f:
            f.write("# BAM file to BigWig file mapping\n")
            f.write("# Format: Sample_name\tBAM_file\tBigWig_file\tStatus\n")
            for bam_path, bw_path in bam_to_bigwig_map.items():
                sample_name = bam_to_sample.get(bam_path, "Unknown")
                bam_filename = os.path.basename(bam_path)
                bw_filename = os.path.basename(bw_path)
                
                # fileのstatusを追加（新規作成/スキップ）
                status = "SKIPPED" if bw_path in skipped_files else "GENERATED"
                if bam_path in failed_bams:
                    status = "FAILED"
                
                f.write(f"{sample_name}\t{bam_filename}\t{bw_filename}\t{status}\n")
        logger.debug(f"Saved BAM to BigWig mapping to: {mapping_file_path}")
        return True
    except Exception as e:
        logger.error(f"Failed to save BAM to BigWig mapping file: {e}")
        return False

# --- 詳細な変換情報をTSVとして保存する関数 ---
def save_bam_to_bigwig_details(
    bam_to_bigwig_map: Dict[str, str],
    bam_to_sample: Dict[str, str],
    output_dir: str,
    bam_dir: str,
    bigwig_files: List[str],
    skipped_files: List[str],
    failed_bams: List[str],
    normalize_using: str,
    bin_size: int,
    effective_genome_size: int,
    threads: int,
    bamcoverage_path: str,
    start_time: float,
    logger: logging.Logger
) -> str:
    """
    Save detailed output information as a TSV file
    
    Args:
        bam_to_bigwig_map (Dict[str, str]): Dictionary mapping BAM paths to BigWig paths
        bam_to_sample (Dict[str, str]): Dictionary mapping BAM paths to sample names
        output_dir (str): Output directory
        bam_dir (str): Directory containing BAM files
        bigwig_files (List[str]): List of generated bigWig files
        skipped_files (List[str]): List of skipped files
        failed_bams (List[str]): List of failed BAM files
        normalize_using (str): Normalization method used
        bin_size (int): Bin size
        effective_genome_size (int): Effective genome size
        threads (int): Number of threads used
        bamcoverage_path (str): Path to bamCoverage executable
        start_time (float): Processing start time
        logger (logging.Logger): Logger object
    
    Returns:
        str: Path to the saved TSV file, or empty string if failed
    """
    # bamCoverageのversionを取得
    bamcoverage_version = _get_bamcoverage_version(bamcoverage_path, logger) or "Unknown"
    
    # 実行時間の計算
    execution_time = time.time() - start_time
    
    # 詳細情報をTSVとして保存
    details_file_path = os.path.join(output_dir, "bam_to_bigwig_details.tsv")
    try:
        with open(details_file_path, 'w', encoding='utf-8') as f:
            # Comment行（metadata）
            f.write(f"# bigwig_peakprep_bamconvert.py - BAM to BigWig Conversion Details\n")
            f.write(f"# Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"# Host: {os.uname().nodename}\n")
            f.write(f"# Working directory: {os.getcwd()}\n")
            f.write(f"# Output directory: {os.path.abspath(output_dir)}\n")
            f.write(f"# BAM directory: {os.path.abspath(bam_dir)}\n")
            f.write(f"# Processed: {len(bigwig_files)} bigWig files ({len(skipped_files)} skipped, {len(failed_bams)} failed)\n")
            f.write(f"# Normalization method: {normalize_using}\n")
            f.write(f"# bin_size: {bin_size}\n")
            f.write(f"# effective_genome_size: {effective_genome_size}\n")
            f.write(f"# threads: {threads}\n")
            f.write(f"# bamCoverage version: {bamcoverage_version}\n")
            f.write(f"# Execution time: {execution_time:.2f} seconds\n")
            f.write(f"# Command: bigwig_peakprep_bamconvert.py --bam_dir=\"{bam_dir}\" --output_dir=\"{output_dir}\" ...\n")
            f.write(f"# This file can be used as input for downstream analysis\n")
            f.write(f"# Format: TSV with 5 columns (Sample_name, BAM_filename, BigWig_filename, BAM_fullpath, BigWig_fullpath)\n")
            f.write(f"# Note: Lines starting with '#' are comments and should be skipped when parsing\n")
            f.write("#\n")
            
            # Header行
            f.write("Sample_name\tBAM_filename\tBigWig_filename\tBAM_fullpath\tBigWig_fullpath\n")
            
            # 各行のdata
            for bam_path, bw_path in bam_to_bigwig_map.items():
                if bam_path in failed_bams:
                    continue  # 失敗したものは含めない
                    
                sample_name = bam_to_sample.get(bam_path, "Unknown")
                bam_filename = os.path.basename(bam_path)
                bw_filename = os.path.basename(bw_path)
                
                # フルpathとfilenameをタブ区切りで出力
                f.write(f"{sample_name}\t{bam_filename}\t{bw_filename}\t{bam_path}\t{bw_path}\n")
        
        logger.info(f"Results saved to: {details_file_path}")
        return os.path.abspath(details_file_path)
    except Exception as e:
        logger.error(f"Failed to save details file: {e}")
        return ""

# --- メイン実行関数 ---
def main():
    """
    Main function for the BAM to bigWig conversion pipeline.
    Returns the path to the bam_to_bigwig_details.tsv file.
    """
    start_time = time.time()
    args = parse_arguments()
    
    # --- Logger設定 ---
    log_level_str = args.log_level.upper()
    log_level = getattr(logging, log_level_str, logging.INFO)
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
    logger = logging.getLogger("BigWigConverter")

    # 出力directory作成
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.debug(f"Created or verified output directory: {args.output_dir}")
    except OSError as e:
        logger.critical(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1)

    # File log handler設定
    if args.script_log_file:
        log_file_path = os.path.join(args.output_dir, args.script_log_file)
        
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

    logger.debug(f"Parsed arguments: {vars(args)}")

    # --- Tool pathの決定 ---
    if args.tools_dir:
        bamcoverage_path = os.path.join(args.tools_dir, "bamCoverage")
        logger.debug(f"Using bamCoverage from directory: {args.tools_dir}")
    else:
        bamcoverage_path = args.bamcoverage_path
        logger.debug(f"Using specified bamCoverage path: {bamcoverage_path}")

    # --- bigWig directoryの設定 ---
    bigwig_dir = args.bigwig_dir if args.bigwig_dir else os.path.join(args.output_dir, "bigwig")
    logger.info(f"BAM to bigWig conversion started")
    logger.debug(f"BAM directory: {args.bam_dir}")
    logger.debug(f"BigWig output directory: {bigwig_dir}")
    
    # --- BAM fileのスキャン ---
    bam_files, bam_to_sample, sample_to_bam = scan_bam_files(
        args.bam_dir, logger, args.filename_pattern, args.sample_delimiter
    )
    if not bam_files:
        logger.critical("No BAM files found. Cannot proceed.")
        sys.exit(1)

    # --- BAM fileからbigWig生成 ---
    bigwig_files, bam_to_bigwig_map, skipped_files, failed_bams = generate_bigwig_files(
        bam_files=bam_files,
        bam_to_sample=bam_to_sample,
        output_dir=args.output_dir,
        bigwig_dir=bigwig_dir,
        logger=logger,
        bamcoverage_path=bamcoverage_path,
        bin_size=args.bin_size,
        normalize_using=args.normalize_using,
        effective_genome_size=args.effective_genome_size,
        threads=args.threads,
        force_overwrite=args.force_overwrite,
        bamcoverage_options=args.bamcoverage_options
    )

    # --- BAMとBigWigのmapping関係をfileに保存 ---
    save_bam_to_bigwig_mapping(
        bam_to_bigwig_map=bam_to_bigwig_map,
        bam_to_sample=bam_to_sample,
        output_dir=args.output_dir,
        skipped_files=skipped_files,
        failed_bams=failed_bams,
        logger=logger
    )
    
    # --- 詳細な変換情報をTSVとして保存 ---
    details_file_path = save_bam_to_bigwig_details(
        bam_to_bigwig_map=bam_to_bigwig_map,
        bam_to_sample=bam_to_sample,
        output_dir=args.output_dir,
        bam_dir=args.bam_dir,
        bigwig_files=bigwig_files,
        skipped_files=skipped_files,
        failed_bams=failed_bams,
        normalize_using=args.normalize_using,
        bin_size=args.bin_size,
        effective_genome_size=args.effective_genome_size,
        threads=args.threads,
        bamcoverage_path=bamcoverage_path,
        start_time=start_time,
        logger=logger
    )
    
    # --- パイプラインの終了 ---
    end_time = time.time()
    logger.info(f"BAM to bigWig conversion finished in {end_time - start_time:.2f} seconds")
    
    # 詳細出力のTSV file pathを返す
    return details_file_path

# --- スクリプト実行のエントリポイント ---
if __name__ == "__main__":
    output_file = main()
    if output_file:
        print(f"BAM to bigWig conversion completed successfully. Output file: {output_file}")
    else:
        print("BAM to bigWig conversion completed with errors. Check the log for details.")
        sys.exit(1)