"""
Main processing module for the SEdb region extraction tool
"""

import logging
import sys
import time
from pathlib import Path

from sedb_preprocessor.loader import SEdbDataLoader
from sedb_preprocessor.processor import SEdbDataProcessor
from .cli import parse_args, validate_args, parse_region
from .extractor import RegionOverlapAnalyzer
from .packager import RegionPackager
from sedb_common.logging import get_module_logger


# print("DEBUG: extractor/main.py execution started.", file=sys.stderr) 


def main():
    """Execute main processing"""
    # コマンドライン引数の解析
    args = parse_args()


    # print(f"DEBUG: Logger 'sedb_extractor' handlers BEFORE get_module_logger: {logging.getLogger('sedb_extractor').handlers}", file=sys.stderr)

    
    # # ロガー設定
    # log_level = logging.DEBUG if args.debug else logging.INFO
    # # logger = get_module_logger("sedb_extractor", file_output=(args.log_file is not None), log_level=log_level)
    

    # logger = get_module_logger(
    #     module_name="sedb_extractor",
    #     # file_output は LoggerManager 側で log_file_name の有無で判断するので不要かも
    #     # file_output=(args.log_file is not None), # この行は削除またはコメントアウトしても良い
    #     log_level=log_level,
    #     log_file=args.log_file # args.log_file を log_file 引数として渡す
    # )


    # ロガー設定
    log_level = logging.DEBUG if args.debug else logging.INFO
    logger = get_module_logger(
        module_name="sedb_extractor",
        file_output=(args.log_file is not None),
        log_level=log_level,
        log_file=args.log_file
    )


    # print(f"DEBUG: Logger 'sedb_extractor' handlers AFTER get_module_logger: {logger.handlers}", file=sys.stderr)
    # print(f"DEBUG: Logger level: {logger.level}, Effective level: {logger.getEffectiveLevel()}", file=sys.stderr)
    # print(f"DEBUG: Attempting to log first message...", file=sys.stderr)


    logger.info("Starting SEdb region extraction")
    start_time = time.time()
    
    # 引数の検証
    if not validate_args(args, logger):
        logger.error("Argument validation failed")
        return 1
    
    # 出力ディレクトリが存在することを確認し、なければ作成する
    try:
        output_path = Path(args.output)
        output_path.mkdir(parents=True, exist_ok=True) # exist_ok=True で既に存在してもエラーにしない
        logger.info(f"Output directory confirmed/created: {args.output}")
    except Exception as e:
        logger.error(f"Failed to confirm/create output directory: {args.output} - {e}")
        return 1 # ディレクトリが確保できない場合は処理を中断





    try:
        # 領域のパース
        region_chr, region_start, region_end = parse_region(args.region)
        region_name = args.region_name or f"{region_chr}:{region_start}-{region_end}"
        
        logger.info(f"Target region: {region_name} ({region_chr}:{region_start}-{region_end})")
        
        # データの読み込み
        loader = SEdbDataLoader(logger=logger)
        
        logger.info(f"Loading BED file: {args.bed_file}")
        bed_df, original_columns, bed_md5 = loader.load_bed_file(args.bed_file)
        logger.info(f"Loaded BED file: {len(bed_df)} rows")
        
        logger.info(f"Loading sample metadata: {args.metadata}")
        sample_info, metadata_md5 = loader.load_sample_info(args.metadata)
        logger.info(f"Loaded metadata: {len(sample_info)} samples")

        processor = SEdbDataProcessor(logger=logger)
        # 必要であれば NaN の Sample ID を除去 (preprocessor でもやっているが念のため)
        # sample_info = processor.remove_na_sample_ids(sample_info)
        sample_info = processor.add_cell_id_to_sample_info(sample_info) # 'cell_id' カラムを追加
        if 'cell_id' in sample_info.columns:
             logger.info("Added 'cell_id' column to metadata.")
        else:
             # もしカラム追加に失敗したらエラーを出すなど（通常は成功するはず）
             logger.error("Failed to add 'cell_id' column to metadata!")
             return 1 # エラー終了させるなど


        # ソースファイル情報
        source_files = {
            "bed_file": {"path": str(args.bed_file), "md5": bed_md5},
            "metadata_file": {"path": str(args.metadata), "md5": metadata_md5}
        }
        
        # 領域抽出の実行
        analyzer = RegionOverlapAnalyzer(logger=logger)
        # extraction_results = analyzer.extract_region(
        #     region_chr, region_start, region_end, region_name,
        #     bed_df, sample_info, 
        #     extraction_method=args.extraction_method,
        #     analysis_level=args.analysis_level
        # )

        extraction_results = analyzer.extract_region(
            region_chr, region_start, region_end, region_name,
            bed_df, sample_info,
            extraction_method=args.extraction_method,
            analysis_level=args.analysis_level,
            output_dir=args.output # <-- この行を追加して出力ディレクトリを渡す
        )


        # 結果の保存
        packager = RegionPackager(logger=logger)
        package_info = packager.save_region_package(
            extraction_results,
            args.output,
            include_distributions=args.include_distributions,
            source_files=source_files,
            export_bed=args.export_bed,
            bed_format=args.bed_format
        )
        
        # 処理時間と結果の表示
        elapsed_time = time.time() - start_time
        logger.info(f"Region extraction completed: {elapsed_time:.2f} seconds")
        logger.info(f"Region data package saved: {package_info['output_dir']}")
        
        # 統計情報の表示
        stats = extraction_results['stats']
        logger.info(f"Number of overlapping SE regions: {stats['total_regions']}")
        logger.info(f"Number of unique samples: {stats.get('unique_samples', 'N/A')}")
        
        if 'tissue_types' in stats and stats['tissue_types'] > 0:
            logger.info(f"Number of tissue types: {stats['tissue_types']}")
        
        if 'biosample_types' in stats and stats['biosample_types'] > 0:
            logger.info(f"Number of biosample types: {stats['biosample_types']}")
        
        # BEDファイル出力結果の表示
        if args.export_bed:
            bedfile_path = Path(package_info['output_dir']) / "overlapping_regions.bed"
            if bedfile_path.exists():
                logger.info(f"BED file output: {bedfile_path}")
        
        return 0
        
    except Exception as e:
        logger.exception(f"An error occurred during region extraction: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
