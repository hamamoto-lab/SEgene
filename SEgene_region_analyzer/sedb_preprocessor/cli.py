import argparse
import logging
import sys
import os
from pathlib import Path

from sedb_common.logging import LoggerManager

def setup_logger(log_level=logging.INFO, log_file=None):
    """Setup logger"""
    manager = LoggerManager(
        logger_name="sedb_preprocessor",
        enable_console_logging=True,
        console_verbose=(log_level == logging.DEBUG),
        enable_file_logging=(log_file is not None),
        log_file_name=log_file,
        log_level=log_level
    )
    return manager.get_logger()


def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="SEdb Data Preprocessing Tool - Processes BED files and metadata to generate a standardized package"
    )
    
    # 必須引数
    parser.add_argument('--bed-file', '-b', required=True, help='Path to the SEDB BED file')
    parser.add_argument('--metadata', '-m', required=True, help='Path to the SEDB sample metadata file')
    parser.add_argument('--output', '-o', required=True, help='Path to the output directory')
    
    # オプション引数
    parser.add_argument('--partition-by-chr', '-p', action='store_true', 
                       help='Partition BED data by chromosome (for large-scale datasets)')
    parser.add_argument('--filter-chromosomes', '-f', action='store_true',
                       help='Retain only standard human chromosomes (chr1-22, X, Y)')
    parser.add_argument('--log-file', '-l', help='Path to the log file (console only if omitted)')
    parser.add_argument('--debug', '-d', action='store_true', help='Enable debug mode')
    
    return parser.parse_args()

def validate_args(args, logger):
    """Validate arguments"""
    # 入力ファイルの存在確認
    if not os.path.exists(args.bed_file):
        logger.error(f"BED file not found: {args.bed_file}")
        return False
    
    if not os.path.exists(args.metadata):
        logger.error(f"Metadata file not found: {args.metadata}")
        return False
    
    # 出力ディレクトリの確認と作成
    output_path = Path(args.output)
    if output_path.exists() and not output_path.is_dir():
        logger.error(f"Output path is not a directory: {args.output}")
        return False
    
    if output_path.exists() and any(output_path.iterdir()):
        logger.warning(f"Output directory is not empty: {args.output}")
    
    return True
