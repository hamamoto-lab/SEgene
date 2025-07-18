#!/usr/bin/env python3
"""
SEdb Region Analyzer CLI with caching support

Command-line interface for analyzing super-enhancer regions using SEdb data.
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from datetime import datetime

# SEdbRegionAnalyzerとキャッシュマネージャーをインポート
from .sedb_tools_analyzer import SEdbRegionAnalyzer
from .cache_manager import SEdbCacheManager

def setup_logging(verbose=False):
    """ロギングの設定"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)


def validate_file_exists(filepath, file_type="File"):
    """ファイルの存在確認"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{file_type} not found: {filepath}")
    return filepath


def create_parser():
    """引数パーサーの作成"""
    parser = argparse.ArgumentParser(
        prog='sedb-analyzer',
        description='SEdb Region Analyzer - Analyze super-enhancer regions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Processing Flow:
  1. Load and filter SEdb database
  2. Calculate database-wide statistics (tissue distribution, etc.)
  3. Use these statistics as background for region-specific analysis

Examples:
  # Prepare statistics cache (recommended for large-scale analysis)
  %(prog)s prepare-stats -b SE.bed -s sample.txt -o cache/stats.pkl

  # Batch analysis with cached statistics
  %(prog)s batch -b SE.bed -s sample.txt -r regions.tsv --use-cached-stats cache/stats.pkl

  # Batch analysis without cache (calculates statistics on-the-fly)
  %(prog)s batch -b SE.bed -s sample.txt -r regions.tsv -o results/

  # Analyze single region
  %(prog)s single -b SE.bed -s sample.txt --chr chr7 --start 1000000 --end 2000000

  # Generate SEdb metadata report
  %(prog)s report -b SE.bed -s sample.txt -o report/
        """
    )
    
    # 共通の親パーサー
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-b', '--bed-file',
        required=True,
        help='Path to SE BED file (e.g., SE_package_hg38.bed)'
    )
    parent_parser.add_argument(
        '-s', '--sample-info',
        required=True,
        help='Path to sample information file'
    )
    parent_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parent_parser.add_argument(
        '--enable-peak-analysis',
        action='store_true',
        help='Enable peak analysis (experimental feature)'
    )
    
    # サブコマンド
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # 1. 統計準備コマンド（新規）
    prepare_parser = subparsers.add_parser(
        'prepare-stats',
        parents=[parent_parser],
        help='Pre-calculate and cache database statistics'
    )
    prepare_parser.add_argument(
        '-o', '--output',
        default='sedb_stats_cache.pkl',
        help='Output cache file (default: sedb_stats_cache.pkl)'
    )
    prepare_parser.add_argument(
        '--no-filter-chromosomes',
        action='store_true',
        help='Do not filter for standard human chromosomes'
    )
    
    # 2. バッチ解析コマンド（更新）
    batch_parser = subparsers.add_parser(
        'batch',
        parents=[parent_parser],
        help='Batch analyze multiple regions from file'
    )
    batch_parser.add_argument(
        '-r', '--regions-file',
        required=True,
        help='Path to regions file (TSV or BED format)'
    )
    batch_parser.add_argument(
        '-o', '--output-dir',
        default='sedb_results',
        help='Output directory (default: sedb_results)'
    )
    batch_parser.add_argument(
        '--use-cached-stats',
        help='Path to cached statistics file (speeds up analysis)'
    )
    batch_parser.add_argument(
        '--save-stats-cache',
        help='Save calculated statistics to cache file for future use'
    )
    batch_parser.add_argument(
        '--max-rows',
        type=int,
        help='Maximum number of regions to process'
    )
    batch_parser.add_argument(
        '--no-filter-chromosomes',
        action='store_true',
        help='Do not filter for standard human chromosomes'
    )
    batch_parser.add_argument(
        '--figure-formats',
        nargs='+',
        default=['png', 'svg'],
        choices=['png', 'svg', 'pdf', 'eps'],
        help='Figure output formats (default: png svg)'
    )
    batch_parser.add_argument(
        '--peak-threshold',
        type=float,
        default=0.2,
        help='Peak detection threshold (default: 0.2)'
    )
    batch_parser.add_argument(
        '--pvalue-threshold',
        type=float,
        default=0.05,
        help='P-value threshold for enrichment (default: 0.05)'
    )
    batch_parser.add_argument(
        '--fdr-threshold',
        type=float,
        default=0.05,
        help='FDR threshold for enrichment (default: 0.05)'
    )
    batch_parser.add_argument(
        '--fc-threshold',
        type=float,
        default=1.5,
        help='Fold change threshold (default: 1.5)'
    )
    batch_parser.add_argument(
        '--no-tables',
        action='store_true',
        help='Do not save data tables as CSV'
    )
    
    # 3. 単一領域解析コマンド（更新）
    single_parser = subparsers.add_parser(
        'single',
        parents=[parent_parser],
        help='Analyze a single genomic region'
    )
    single_parser.add_argument(
        '--chr',
        required=True,
        help='Chromosome (e.g., chr7)'
    )
    single_parser.add_argument(
        '--start',
        type=int,
        required=True,
        help='Start position'
    )
    single_parser.add_argument(
        '--end',
        type=int,
        required=True,
        help='End position'
    )
    single_parser.add_argument(
        '--region-name',
        help='Optional region name'
    )
    single_parser.add_argument(
        '-o', '--output-dir',
        default='sedb_results',
        help='Output directory (default: sedb_results)'
    )
    single_parser.add_argument(
        '--use-cached-stats',
        help='Path to cached statistics file'
    )
    single_parser.add_argument(
        '--figure-formats',
        nargs='+',
        default=['png', 'svg'],
        choices=['png', 'svg', 'pdf', 'eps'],
        help='Figure output formats (default: png svg)'
    )
    single_parser.add_argument(
        '--peak-threshold',
        type=float,
        default=0.2,
        help='Peak detection threshold (default: 0.2)'
    )
    single_parser.add_argument(
        '--pvalue-threshold',
        type=float,
        default=0.05,
        help='P-value threshold for enrichment (default: 0.05)'
    )
    
    # 4. レポート生成コマンド
    report_parser = subparsers.add_parser(
        'report',
        parents=[parent_parser],
        help='Generate SEdb metadata report'
    )
    report_parser.add_argument(
        '-o', '--output-dir',
        default='sedb_report',
        help='Output directory (default: sedb_report)'
    )
    report_parser.add_argument(
        '--top-n',
        type=int,
        default=20,
        help='Number of top categories to show (default: 20)'
    )
    report_parser.add_argument(
        '--dpi',
        type=int,
        default=300,
        help='DPI for figures (default: 300)'
    )
    
    return parser


def prepare_database_statistics(analyzer, args, logger):
    """
    データベース統計を事前計算
    
    Parameters:
    -----------
    analyzer : SEdbRegionAnalyzer
        初期化済みのアナライザー
    args : argparse.Namespace
        コマンドライン引数
    logger : Logger
        ロガー
    """
    logger.info("Calculating database-wide statistics...")
    
    # 染色体フィルタリング
    if not args.no_filter_chromosomes:
        logger.info("Filtering for standard human chromosomes...")
        analyzer.filter_human_chromosomes()
    
    # cell_idカウント
    logger.info("Counting cell_id occurrences...")
    analyzer.count_cell_id_occurrences()
    
    # 組織分布の計算
    logger.info("Calculating tissue distribution...")
    analyzer.visualize_tissue_distribution(store_all=True)
    
    # 統計情報が計算されたか確認
    if not hasattr(analyzer, '_tissue_distribution') or analyzer._tissue_distribution is None:
        raise RuntimeError("Failed to calculate tissue distribution")
    
    logger.info(f"Calculated statistics for {len(analyzer._tissue_distribution)} tissue types")


def run_prepare_stats(args, logger):
    """統計情報の事前計算とキャッシュ"""
    logger.info("Preparing database statistics...")
    
    # ファイル検証
    validate_file_exists(args.bed_file, "BED file")
    validate_file_exists(args.sample_info, "Sample info file")
    
    # アナライザーの初期化
    analyzer = SEdbRegionAnalyzer(
        logger=logger,
        enable_peak_analysis=False  # 統計計算には不要
    )
    
    # ファイル設定とデータベースロード
    analyzer.set_se_bed_file(args.bed_file)
    analyzer.set_sample_info_file(args.sample_info)
    analyzer.load_databases()
    
    # 統計を計算
    prepare_database_statistics(analyzer, args, logger)
    
    # キャッシュに保存
    cache_manager = SEdbCacheManager(logger=logger)
    metadata = cache_manager.save_statistics(
        analyzer,
        args.output,
        args.bed_file,
        args.sample_info
    )
    
    logger.info("Statistics preparation completed:")
    logger.info(f"  - Total samples: {metadata['total_samples']}")
    logger.info(f"  - Total SE records: {metadata['total_se_records']}")
    logger.info(f"  - Cache saved to: {args.output}")
    
    return args.output


def run_batch_analysis(args, logger):
    """バッチ解析の実行（キャッシュ対応版）"""
    logger.info("Starting batch analysis...")
    logger.info(f"  BED file: {args.bed_file}")
    logger.info(f"  Sample info: {args.sample_info}")
    logger.info(f"  Regions file: {args.regions_file}")
    logger.info(f"  Output directory: {args.output_dir}")
    
    # ファイル検証
    validate_file_exists(args.bed_file, "BED file")
    validate_file_exists(args.sample_info, "Sample info file")
    validate_file_exists(args.regions_file, "Regions file")
    
    # ファイル形式の検出（TSV or BED）
    from segene_analyzer.batch.file_converters import FileConverter
    converter = FileConverter(logger=logger)
    file_format = converter.detect_file_format(args.regions_file)
    logger.info(f"  Regions file format: {file_format.upper()}")
    
    # アナライザーの初期化
    analyzer = SEdbRegionAnalyzer(
        results_dir=args.output_dir,
        logger=logger,
        enable_peak_analysis=args.enable_peak_analysis
    )
    
    # ファイル設定とデータベースロード
    analyzer.set_se_bed_file(args.bed_file)
    analyzer.set_sample_info_file(args.sample_info)
    analyzer.load_databases()
    
    # キャッシュの処理
    cache_manager = SEdbCacheManager(logger=logger)
    
    if args.use_cached_stats:
        # キャッシュから統計を読み込み
        logger.info(f"Loading cached statistics from: {args.use_cached_stats}")
        success = cache_manager.load_statistics(
            analyzer,
            args.use_cached_stats,
            verify_files=True,
            bed_file=args.bed_file,
            sample_info_file=args.sample_info
        )
        
        if not success:
            logger.warning("Failed to load cached statistics. Calculating from scratch...")
            prepare_database_statistics(analyzer, args, logger)
    else:
        # 統計を新規計算
        prepare_database_statistics(analyzer, args, logger)
    
    # 統計をキャッシュに保存（指定された場合）
    if args.save_stats_cache and not args.use_cached_stats:
        logger.info(f"Saving statistics to cache: {args.save_stats_cache}")
        cache_manager.save_statistics(
            analyzer,
            args.save_stats_cache,
            args.bed_file,
            args.sample_info
        )
    
    # 出力ディレクトリの設定
    data_output_dir = os.path.join(args.output_dir, 'data')
    report_output_dir = os.path.join(args.output_dir, 'reports')
    
    # バッチ解析の実行
    logger.info("Running batch analysis...")
    result_file = analyzer.batch_analyze_regions_from_file(
        bed_file=args.bed_file,
        metadata_file=args.sample_info,
        data_output_dir=data_output_dir,
        report_output_dir=report_output_dir,
        regions_file=args.regions_file,  # TSV/BED両対応
        max_rows=args.max_rows,
        figure_formats=args.figure_formats,
        save_tables=not args.no_tables,
        peak_threshold=args.peak_threshold,
        pvalue_threshold=args.pvalue_threshold,
        fdr_threshold=args.fdr_threshold,
        fc_threshold=args.fc_threshold,
        filter_chromosomes=not args.no_filter_chromosomes,
        count_cell_id=True
    )
    
    logger.info(f"Batch analysis completed.")
    logger.info(f"Results saved to: {result_file}")
    return result_file


def run_single_analysis(args, logger):
    """単一領域解析の実行（キャッシュ対応版）"""
    logger.info(f"Analyzing single region: {args.chr}:{args.start}-{args.end}")
    
    # ファイル検証
    validate_file_exists(args.bed_file, "BED file")
    validate_file_exists(args.sample_info, "Sample info file")
    
    # アナライザーの初期化
    analyzer = SEdbRegionAnalyzer(
        results_dir=args.output_dir,
        logger=logger,
        enable_peak_analysis=args.enable_peak_analysis
    )
    
    # ファイル設定とデータベースロード
    analyzer.set_se_bed_file(args.bed_file)
    analyzer.set_sample_info_file(args.sample_info)
    analyzer.load_databases()
    
    # キャッシュの処理
    if args.use_cached_stats:
        cache_manager = SEdbCacheManager(logger=logger)
        logger.info(f"Loading cached statistics from: {args.use_cached_stats}")
        success = cache_manager.load_statistics(
            analyzer,
            args.use_cached_stats,
            verify_files=True,
            bed_file=args.bed_file,
            sample_info_file=args.sample_info
        )
        
        if not success:
            logger.warning("Failed to load cached statistics. Calculating from scratch...")
            analyzer.filter_human_chromosomes()
            analyzer.count_cell_id_occurrences()
            analyzer.visualize_tissue_distribution(store_all=True)
    else:
        # 統計を新規計算
        analyzer.filter_human_chromosomes()
        analyzer.count_cell_id_occurrences()
        logger.info("Calculating tissue distribution...")
        analyzer.visualize_tissue_distribution(store_all=True)
    
    # 単一領域の解析
    result = analyzer.analyze_and_report_region(
        chrom=args.chr,
        start=args.start,
        end=args.end,
        region_name=args.region_name,
        output_dir=args.output_dir,
        figure_formats=args.figure_formats,
        peak_threshold=args.peak_threshold,
        pvalue_threshold=args.pvalue_threshold
    )
    
    # 結果のサマリー表示
    logger.info("Analysis completed:")
    logger.info(f"  - Report: {result['report_path']}")
    logger.info(f"  - Overlapping samples: {result['overlap_sample_count']}")
    logger.info(f"  - Significant tissues: {len(result['significant_tissues'])}")
    
    # 有意な組織の表示（上位5つ）
    if result['significant_tissues']:
        logger.info("  - Top significant tissues:")
        for i, tissue in enumerate(result['significant_tissues'][:5]):
            logger.info(f"    {i+1}. {tissue['tissue']}: FC={tissue['fold_change']:.2f}, FDR={tissue['fdr']:.2e}")
    
    return result


def run_report_generation(args, logger):
    """SEdbメタデータレポートの生成"""
    logger.info("Generating SEdb metadata report...")
    
    # ファイル検証
    validate_file_exists(args.bed_file, "BED file")
    validate_file_exists(args.sample_info, "Sample info file")
    
    # アナライザーの初期化
    analyzer = SEdbRegionAnalyzer(logger=logger)
    
    # ファイル設定とデータベースロード
    analyzer.set_se_bed_file(args.bed_file)
    analyzer.set_sample_info_file(args.sample_info)
    analyzer.load_databases()
    
    # 標準的な前処理
    logger.info("Preprocessing data...")
    analyzer.filter_human_chromosomes()
    analyzer.count_cell_id_occurrences()
    
    # レポート生成
    report_path = analyzer.generate_sedb_report(
        output_dir=args.output_dir,
        top_n=args.top_n,
        dpi=args.dpi
    )
    
    logger.info(f"Report generated: {report_path}")
    return report_path


def main():
    """メイン関数"""
    parser = create_parser()
    args = parser.parse_args()
    
    # コマンドが指定されていない場合
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # ロギング設定
    logger = setup_logging(args.verbose)
    
    try:
        # コマンドに応じた処理
        if args.command == 'prepare-stats':
            run_prepare_stats(args, logger)
            
        elif args.command == 'batch':
            run_batch_analysis(args, logger)
            
        elif args.command == 'single':
            run_single_analysis(args, logger)
            
        elif args.command == 'report':
            run_report_generation(args, logger)
            
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)
    
    logger.info("Process completed successfully")


if __name__ == '__main__':
    main()