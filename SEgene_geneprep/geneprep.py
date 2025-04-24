#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gene Position and Expression Merger

This script extracts gene position information from a GTF file and merges it with
gene expression data from a TPM file to create a combined dataset.
"""

import argparse
import logging
import sys
import os
from pathlib import Path

# ユーティリティ関数をインポート (Import utility functions)
from geneprep_utils import (
    setup_logging,
    extract_gene_info_from_gtf,
    read_tpm_file,
    save_unmatched_lines,
    save_id_name_mismatches,
    save_undefined_strand_lines,
    merge_gene_data,
    print_statistics,
    validate_files,
    ensure_dir_exists
)

def parse_args():
    """
    Parse command line arguments.
    
    Returns
    -------
    argparse.Namespace
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Merge gene position information from GTF with expression data from TPM file"
    )
    parser.add_argument("--gtf", required=True, help="Path to GTF file")
    parser.add_argument("--tpm", required=True, help="Path to TPM file")
    parser.add_argument(
        "--output", default="gene_tpm_with_position.csv", 
        help="Output file name (default: gene_tpm_with_position.csv)"
    )
    parser.add_argument(
        "--output-dir", 
        help="Output directory for all files (overrides individual file paths)"
    )
    parser.add_argument(
        "--unmatched", default="unmatched_gtf_lines.tsv", 
        help="Output file for unmatched GTF lines (default: unmatched_gtf_lines.tsv)"
    )
    parser.add_argument(
        "--mismatch", default="id_name_mismatches.tsv", 
        help="Output file for gene_id/gene_name mismatches (default: id_name_mismatches.tsv)"
    )
    parser.add_argument(
        "--undefined-strand", default="undefined_strand_lines.tsv",
        help="Output file for lines with undefined strand (default: undefined_strand_lines.tsv)"
    )
    parser.add_argument(
        "--join", choices=["inner", "left", "outer"], default="inner",
        help=("Join method: inner (only genes in both files), "
              "left (all genes in GTF), outer (all genes in either file) (default: inner)")
    )
    parser.add_argument(
        "--log", help="Path to log file (if specified, logs will be saved to this file)"
    )
    parser.add_argument(
        "--no-create-dirs", action="store_true",
        help="Do not create output directories if they don't exist (by default, directories are created automatically)"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", 
        help="Enable verbose output with detailed information"
    )
    return parser.parse_args()

def main():
    """
    Main function that coordinates the merging process.
    
    Returns
    -------
    int
        Exit code (0 for success, 1 for failure)
    """
    # コマンドライン引数の解析 (Parse command line arguments)
    args = parse_args()
    
    # 出力ディレクトリの一括指定がある場合、個別のファイルパスを調整
    if args.output_dir:
        # メイン出力ファイルのファイル名のみを取得
        output_filename = os.path.basename(args.output)
        args.output = os.path.join(args.output_dir, output_filename)
        
        # 他のファイルパスもデフォルト値の場合は調整
        if args.unmatched == "unmatched_gtf_lines.tsv":
            args.unmatched = os.path.join(args.output_dir, "unmatched_gtf_lines.tsv")
        
        if args.mismatch == "id_name_mismatches.tsv":
            args.mismatch = os.path.join(args.output_dir, "id_name_mismatches.tsv")
        
        if args.undefined_strand == "undefined_strand_lines.tsv":
            args.undefined_strand = os.path.join(args.output_dir, "undefined_strand_lines.tsv")
    else:
        # 出力ディレクトリの取得 (Get output directory)
        output_dir = os.path.dirname(args.output) if os.path.dirname(args.output) else '.'
        
        # オプショナルファイルの出力パスを調整 (Adjust optional file paths)
        # 個別に指定されていない場合のみ、メイン出力と同じディレクトリを使用
        if args.unmatched == "unmatched_gtf_lines.tsv":
            args.unmatched = os.path.join(output_dir, "unmatched_gtf_lines.tsv")
        
        if args.mismatch == "id_name_mismatches.tsv":
            args.mismatch = os.path.join(output_dir, "id_name_mismatches.tsv")
        
        if args.undefined_strand == "undefined_strand_lines.tsv":
            args.undefined_strand = os.path.join(output_dir, "undefined_strand_lines.tsv")
    
    # ディレクトリ作成の設定 (デフォルトで自動作成、--no-create-dirsで無効化)
    create_dirs = not args.no_create_dirs
    
    # ログファイルのディレクトリが存在するか確認
    if args.log and not ensure_dir_exists(args.log, create_dirs):
        return 1
    
    # ロギングのセットアップ (Setup logging)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logging(log_level=log_level, log_file=args.log)
    
    # もし詳細モードなら、コマンドライン引数を表示 (If verbose, display command line arguments)
    if args.verbose:
        logging.debug("Running with arguments:")
        for arg, value in vars(args).items():
            logging.debug(f"  {arg}: {value}")
    
    # 出力ファイルのディレクトリが存在するか確認
    output_files = [args.output, args.unmatched, args.mismatch, args.undefined_strand]
    for output_file in output_files:
        if not ensure_dir_exists(output_file, create_dirs):
            return 1
    
    # ファイルの存在確認 (Validate file existence)
    if not validate_files(args.gtf, args.tpm):
        return 1
    
    # GTFファイルから遺伝子情報を抽出 (Extract gene information from GTF file)
    logging.info(f"Extracting gene information from GTF file {args.gtf}...")
    gene_df, unmatched_lines, id_name_mismatches, undefined_strand_lines = extract_gene_info_from_gtf(args.gtf, verbose=args.verbose)
    logging.info(f"Extracted {len(gene_df)} gene records from GTF")
    logging.info(f"Unmatched lines: {len(unmatched_lines)}")
    logging.info(f"Gene_id/gene_name mismatches: {len(id_name_mismatches)}")
    logging.info(f"Lines with undefined strand: {len(undefined_strand_lines)}")
    
    # マッチしなかった行を出力 (Save unmatched lines)
    save_unmatched_lines(unmatched_lines, args.unmatched)
    
    # gene_idとgene_nameが一致しない遺伝子を出力 (Save gene_id/gene_name mismatches)
    save_id_name_mismatches(id_name_mismatches, args.mismatch)
    
    # ストランドが未定義の行を出力 (Save lines with undefined strand)
    save_undefined_strand_lines(undefined_strand_lines, args.undefined_strand)
    
    # TPMファイルを読み込む (Read TPM file)
    logging.info(f"Reading TPM file {args.tpm}...")
    tpm_df = read_tpm_file(args.tpm, verbose=args.verbose)
    if tpm_df is None:
        logging.error("Failed to read TPM file. Aborting.")
        return 1
    
    logging.info(f"Read {len(tpm_df)} gene expression records from TPM file")
    
    # データを結合 (Merge data)
    result_df = merge_gene_data(gene_df, tpm_df, args.join, verbose=args.verbose)
    
    # CSVファイルとして出力 (Save as CSV)
    try:
        result_df.to_csv(args.output, index=False)
        logging.info(f"Saved results to {args.output}")
    except Exception as e:
        logging.error(f"Error saving results: {e}")
        return 1
    
    # 統計情報を表示 (Display statistics)
    print_statistics(gene_df, tpm_df, result_df, args.join, verbose=args.verbose)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())