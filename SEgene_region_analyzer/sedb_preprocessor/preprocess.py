#!/usr/bin/env python3
"""SEdb Data Preprocessing Program"""

import sys
import logging
import time
from pathlib import Path

from .loader import SEdbDataLoader
from .processor import SEdbDataProcessor
from .packager import SEdbDataPackager
from .cli import setup_logger, parse_arguments, validate_args

def main():
    # コマンドライン引数の解析
    args = parse_arguments()
    
    # ロガーのセットアップ
    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = setup_logger(log_level=log_level, log_file=args.log_file)
    
    logger.info("Starting SEdb data preprocessing")
    start_time = time.time()
    
    # 引数の検証
    if not validate_args(args, logger):
        return 1
    
    try:
        # 1. データローダーの初期化と実行
        loader = SEdbDataLoader(logger=logger)
        
        logger.info(f"Loading BED file: {args.bed_file}")
        bed_df, original_columns, bed_md5 = loader.load_bed_file(args.bed_file)  # MD5も受け取る
        logger.info(f"BED file loaded: {len(bed_df)} rows")
        logger.info(f"MD5 hash of BED file: {bed_md5}")
        
        logger.info(f"Loading sample metadata: {args.metadata}")
        sample_info, metadata_md5 = loader.load_sample_info(args.metadata)  # MD5も受け取る
        logger.info(f"Metadata loaded: {len(sample_info)} samples")
        logger.info(f"MD5 hash of metadata file: {metadata_md5}")
        
        # 2. データプロセッサーの初期化と実行
        processor = SEdbDataProcessor(logger=logger)
        
        # サンプルIDのクリーニング
        sample_info = processor.remove_na_sample_ids(sample_info)
        
        # cell_idの追加
        sample_info = processor.add_cell_id_to_sample_info(sample_info)
        
        # 染色体フィルタリング（オプション）
        if args.filter_chromosomes:
            bed_df, original_bed_df, removed_bed_df = processor.filter_human_chromosomes(bed_df)
            logger.info(f"BED data after chromosome filtering: {len(bed_df)} rows")
        
        # cell_id出現回数のカウント
        cell_id_counts = processor.count_cell_id_occurrences(bed_df, sample_info)
        
        # 基本統計の計算
        statistics = processor.calculate_basic_statistics(bed_df, sample_info)
        
        # ソースファイル情報の追加
        source_files = {
            "bed_file": {
                "path": str(args.bed_file),
                "md5": bed_md5
            },
            "metadata_file": {
                "path": str(args.metadata),
                "md5": metadata_md5
            }
        }
        
        # 3. データパッケージャーの初期化と実行
        packager = SEdbDataPackager(logger=logger)
        
        # 処理済みデータの保存
        package_info = packager.save_processed_package(
            args.output,
            bed_df=bed_df,
            sample_info=sample_info,
            cell_id_counts=cell_id_counts,
            statistics=statistics,
            original_columns=original_columns,
            partition_by_chr=args.partition_by_chr,
            source_files=source_files  # ソースファイル情報を追加
        )
        
        # 処理時間の計算と表示
        elapsed_time = time.time() - start_time
        logger.info(f"Data preprocessing completed in: {elapsed_time:.2f} seconds")
        logger.info(f"Processed data package saved to: {package_info['output_dir']}")
        
        return 0
        
    except Exception as e:
        logger.exception(f"An error occurred during data preprocessing: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
