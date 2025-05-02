#!/usr/bin/env python3
# === cpm_peakprep_calcnorm.py ===

# --- 標準ライブラリ ---
import argparse
import os
import sys
import time
from typing import Optional

# --- サードパーティライブラリ ---
import pandas as pd

# --- ユーティリティ関数のインポート ---
try:
    from cpm_peakprep_calcnorm_utils import (
        apply_edger_normalization,
        validate_counts_file,
        clean_count_file_sample_names,
        setup_logging,
        create_output_directories,
        get_temp_cleaned_file_path,
        load_pattern_rules,
        cleanup_temp_files,
        convert_to_tsv_format,
        ensure_tsv_extension
    )
except ImportError as e:
    print(f"Error importing functions from cpm_peakprep_calcnorm_utils.py: {e}", file=sys.stderr)
    print("Please ensure 'cpm_peakprep_calcnorm_utils.py' exists and necessary libraries are installed.", file=sys.stderr)
    sys.exit(1)

# --- 引数パーサー ---
def parse_arguments() -> argparse.Namespace:
    """
    Parses command line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="Apply edgeR-style normalization to featureCounts output and generate log-CPM values in TSV format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 必須引数
    parser.add_argument("-c", "--counts_file", required=True,
                        help="Path to featureCounts output file (raw counts)")
    parser.add_argument("-o", "--output_file", required=True,
                        help="Output file for normalized log-CPM values in TSV format (PeakID, Chr, Start, End, ...)")
    
    # サンプル名クリーニングオプション
    clean_group = parser.add_argument_group('Sample name cleaning options')
    clean_group.add_argument("--cleaned_output", default=None,
                        help="Path to save the intermediate counts file with cleaned sample names. "
                             "If not specified, a temporary file will be used.")
    clean_group.add_argument("--pattern_rules", default=None,
                        help="JSON file with pattern rules for additional sample name cleaning. "
                             "Format: [{\"pattern\": \"REGEX_PATTERN\", \"replace\": \"REPLACEMENT\"}]")
    clean_group.add_argument("--remove_extensions", action="store_true",
                        help="Remove common BAM file extensions (.bam, .sorted.bam, etc)")
    clean_group.add_argument("--sample_delimiter", default=None,
                        help="Delimiter string to extract sample names from BAM filenames. "
                             "If specified, the part of the filename before this delimiter will be used as the sample name.")
    
    # フィルタリングオプション
    filtering_group = parser.add_argument_group('Filtering options')
    filtering_group.add_argument("--min_cpm", type=float, default=1.0,
                        help="Minimum CPM threshold for filtering")
    filtering_group.add_argument("--min_samples", type=int, default=0,
                        help="Minimum number of samples with CPM > threshold. Default 0 means no filtering.")
    
    # 正規化手法
    parser.add_argument("--method", default="upperquartile",
                        choices=["upperquartile", "TMM", "RLE", "none"],
                        help="Normalization method (upperquartile, TMM, RLE, none)")
    
    # 出力形式オプション
    parser.add_argument("--full_metadata_output", default=None,
                        help="Path to save output in complete format with all featureCounts metadata columns (optional)")
    
    # ロギングオプション
    parser.add_argument("--log_level", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging level")
    parser.add_argument("--script_log_file", default="calcNormFactors_run.log",
                        help="Filename for script execution log. Set to '' to disable.")
    
    return parser.parse_args()

# --- メイン関数 ---
def main() -> None:
    """Main function to calculate CPM values from featureCounts file."""
    start_time = time.time()
    args = parse_arguments()
    
    # ログディレクトリの決定
    log_dir = os.path.dirname(args.output_file) or '.'
    
    # ロギングのセットアップ
    logger = setup_logging(args.log_level, args.script_log_file, log_dir)
    logger.debug(f"Parsed arguments: {vars(args)}")
    
    try:
        # 出力ファイルパスのTSV拡張子を保証
        output_file = ensure_tsv_extension(args.output_file, logger)
        
        # 出力ディレクトリの作成
        if not create_output_directories(output_file, args.cleaned_output, args.script_log_file, log_dir, args.full_metadata_output, logger):
            sys.exit(1)
        
        # 入力ファイルの検証
        logger.info("Validating counts file...")
        if not validate_counts_file(args.counts_file, logger):
            logger.critical("Counts file validation failed")
            sys.exit(1)
        
        # パターンルールの読み込み
        pattern_rules = load_pattern_rules(args.pattern_rules, logger, args.remove_extensions)
        
        # クリーニングファイルパスの決定
        cleaned_file_path = args.cleaned_output
        if not cleaned_file_path:
            cleaned_file_path = get_temp_cleaned_file_path(output_file, logger)
        
        # サンプル名クリーニングの実行
        logger.info(f"Cleaning sample names in {args.counts_file}...")
        if args.sample_delimiter:
            logger.info(f"Using delimiter '{args.sample_delimiter}' to extract sample names")
        if pattern_rules:
            logger.info(f"Using {len(pattern_rules)} pattern rules for sample name cleaning")
        else:
            logger.info("Using only basic sample name cleaning (basename extraction)")
            
        if not clean_count_file_sample_names(args.counts_file, cleaned_file_path, logger, args.sample_delimiter, pattern_rules):
            logger.critical("Sample name cleaning failed")
            sys.exit(1)
        
        logger.info(f"Sample names cleaned and saved to: {cleaned_file_path}")
        
        # クリーニングされたカウントファイルを読み込む
        logger.info(f"Reading cleaned counts file: {cleaned_file_path}")
        try:
            counts_df = pd.read_csv(cleaned_file_path, sep='\t', comment='#')
            logger.info(f"Loaded counts data with {len(counts_df)} features and {len(counts_df.columns) - 6} samples")
        except Exception as e:
            logger.critical(f"Failed to read cleaned counts file: {e}")
            sys.exit(1)
        
        # フィルタリング設定をログに記録
        if args.min_samples > 0:
            logger.info(f"Feature filtering is ENABLED (min_cpm={args.min_cpm}, min_samples={args.min_samples})")
        else:
            logger.info("Feature filtering is DISABLED (min_samples=0)")
        
        # 正規化を適用
        logger.info(f"Applying edgeR-style normalization with method: {args.method}...")
        normalized_df = apply_edger_normalization(
            counts_df, 
            args.min_cpm, 
            args.min_samples, 
            logger,
            method=args.method
        )
        
        if normalized_df is None:
            logger.critical("Normalization failed")
            sys.exit(1)
        
        # フルメタデータ形式の出力（オプション）
        if args.full_metadata_output:
            logger.info(f"Saving complete data with full featureCounts metadata to: {args.full_metadata_output}")
            try:
                normalized_df.to_csv(args.full_metadata_output, sep='\t', index=False, float_format='%.6f')
                logger.info(f"Complete data with full metadata saved successfully: {args.full_metadata_output}")
            except Exception as e:
                logger.error(f"Failed to save full metadata output: {e}")
                logger.warning("Continuing with standard TSV output only")
        
        # 標準TSV形式に変換
        tsv_df = convert_to_tsv_format(normalized_df, logger)
        
        # 標準形式で保存
        logger.info(f"Saving standard TSV format to: {output_file}")
        try:
            tsv_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
            logger.info(f"Standard TSV format saved successfully: {output_file}")
        except Exception as e:
            logger.critical(f"Failed to save TSV output: {e}")
            sys.exit(1)
        
        # 一時ファイルのクリーンアップ
        is_temp_cleaned_file = not args.cleaned_output
        cleanup_temp_files(cleaned_file_path, is_temp_cleaned_file, logger)
        
        # 完了をログに記録
        end_time = time.time()
        logger.info(f"Processing completed in {end_time - start_time:.2f} seconds")
        logger.info(f"Output file: {output_file}")
        if args.full_metadata_output:
            logger.info(f"Full metadata output: {args.full_metadata_output}")
        
    except Exception as e:
        logger.critical(f"Error during processing: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()