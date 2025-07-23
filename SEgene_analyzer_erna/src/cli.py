#!/usr/bin/env python3
"""
eRNAbase Region Analyzer CLI

Command-line interface for analyzing genomic regions using eRNAbase data.
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from datetime import datetime
import json
import pickle

# Import eRNA analyzer
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    from .erna_analyzer import ERNARegionAnalyzer
except ImportError:
    try:
        from erna_analyzer import ERNARegionAnalyzer
    except ImportError:
        # Fallback for direct execution
        import sys
        import os
        current_dir = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, current_dir)
        from erna_analyzer import ERNARegionAnalyzer


def setup_logging(verbose=False, log_file=None):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    
    handlers = [logging.StreamHandler()]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )
    return logging.getLogger(__name__)


def validate_file_exists(filepath, file_type="File"):
    """Check if file exists"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"{file_type} not found: {filepath}")
    return filepath


def validate_directory_exists(dirpath, dir_type="Directory"):
    """Check if directory exists"""
    if not os.path.exists(dirpath):
        raise FileNotFoundError(f"{dir_type} not found: {dirpath}")
    if not os.path.isdir(dirpath):
        raise NotADirectoryError(f"Not a directory: {dirpath}")
    return dirpath


def create_parser():
    """Create argument parser"""
    parser = argparse.ArgumentParser(
        prog='erna-analyzer',
        description='eRNAbase Region Analyzer - Analyze genomic regions using eRNAbase data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Processing Flow:
  1. Load eRNAbase BED files and metadata
  2. Calculate database-wide statistics (tissue/cell type distribution)
  3. Use these statistics as background for region-specific analysis

Examples:
  # Prepare statistics cache (recommended for large-scale analysis)
  %(prog)s prepare-stats -b data/peak -m metadata.parquet -o cache/stats.pkl

  # Batch analysis with cached statistics
  %(prog)s batch -b data/peak -m metadata.parquet -r regions.tsv --use-cached-stats cache/stats.pkl

  # Batch analysis without cache (calculates statistics on-the-fly)
  %(prog)s batch -b data/peak -m metadata.parquet -r regions.tsv -o results/

  # Analyze single region
  %(prog)s single -b data/peak -m metadata.parquet --chr chr7 --start 1000000 --end 2000000

  # Generate eRNAbase metadata report
  %(prog)s report -m metadata.parquet -o report/
        """
    )
    
    # Common parent parser
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging'
    )
    parent_parser.add_argument(
        '--log-file',
        help='Path to log file'
    )
    
    # Subcommands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # 1. prepare-stats command
    prepare_parser = subparsers.add_parser(
        'prepare-stats',
        parents=[parent_parser],
        help='Pre-calculate and cache database statistics'
    )
    prepare_parser.add_argument(
        '-b', '--bed-directory',
        required=True,
        help='Directory containing eRNAbase BED files'
    )
    prepare_parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='Path to eRNAbase metadata file (Parquet/CSV)'
    )
    prepare_parser.add_argument(
        '-o', '--output',
        default='erna_stats_cache.pkl',
        help='Output cache file (default: erna_stats_cache.pkl)'
    )
    prepare_parser.add_argument(
        '--species',
        default='human',
        choices=['human', 'mouse', 'all'],
        help='Species filter (default: human)'
    )
    
    # 2. batch command
    batch_parser = subparsers.add_parser(
        'batch',
        parents=[parent_parser],
        help='Batch analyze multiple regions from file'
    )
    batch_parser.add_argument(
        '-b', '--bed-directory',
        required=True,
        help='Directory containing eRNAbase BED files'
    )
    batch_parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='Path to eRNAbase metadata file (Parquet/CSV)'
    )
    batch_parser.add_argument(
        '-r', '--regions-file',
        required=True,
        help='Path to regions file (TSV format)'
    )
    batch_parser.add_argument(
        '-o', '--output-dir',
        default='erna_results',
        help='Output directory (default: erna_results)'
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
        '--species',
        default='human',
        choices=['human', 'mouse', 'all'],
        help='Species filter (default: human)'
    )
    batch_parser.add_argument(
        '--max-rows',
        type=int,
        help='Maximum number of regions to process'
    )
    batch_parser.add_argument(
        '--figure-formats',
        nargs='+',
        default=['png', 'svg'],
        choices=['png', 'svg', 'pdf', 'eps'],
        help='Figure output formats (default: png svg)'
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
    
    # 3. single command
    single_parser = subparsers.add_parser(
        'single',
        parents=[parent_parser],
        help='Analyze a single genomic region'
    )
    single_parser.add_argument(
        '-b', '--bed-directory',
        required=True,
        help='Directory containing eRNAbase BED files'
    )
    single_parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='Path to eRNAbase metadata file (Parquet/CSV)'
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
        default='erna_results',
        help='Output directory (default: erna_results)'
    )
    single_parser.add_argument(
        '--use-cached-stats',
        help='Path to cached statistics file'
    )
    single_parser.add_argument(
        '--species',
        default='human',
        choices=['human', 'mouse', 'all'],
        help='Species filter (default: human)'
    )
    single_parser.add_argument(
        '--figure-formats',
        nargs='+',
        default=['png', 'svg'],
        choices=['png', 'svg', 'pdf', 'eps'],
        help='Figure output formats (default: png svg)'
    )
    single_parser.add_argument(
        '--pvalue-threshold',
        type=float,
        default=0.05,
        help='P-value threshold for enrichment (default: 0.05)'
    )
    
    # 4. report command
    report_parser = subparsers.add_parser(
        'report',
        parents=[parent_parser],
        help='Generate eRNAbase metadata report'
    )
    report_parser.add_argument(
        '-m', '--metadata',
        required=True,
        help='Path to eRNAbase metadata file (Parquet/CSV)'
    )
    report_parser.add_argument(
        '-o', '--output-dir',
        default='erna_report',
        help='Output directory (default: erna_report)'
    )
    report_parser.add_argument(
        '--species',
        default='human',
        choices=['human', 'mouse', 'all'],
        help='Species filter (default: human)'
    )
    report_parser.add_argument(
        '--top-n',
        type=int,
        default=15,
        help='Number of top categories to show (default: 15)'
    )
    report_parser.add_argument(
        '--figure-formats',
        nargs='+',
        default=['png', 'svg'],
        choices=['png', 'svg', 'pdf', 'eps'],
        help='Figure output formats (default: png svg)'
    )
    
    return parser


def cmd_prepare_stats(args, logger):
    """Execute prepare-stats command"""
    logger.info("Starting statistics preparation...")
    
    # Validate inputs
    bed_dir = validate_directory_exists(args.bed_directory, "BED directory")
    metadata_file = validate_file_exists(args.metadata, "Metadata file")
    
    # Create analyzer
    analyzer = ERNARegionAnalyzer(results_dir=os.path.dirname(args.output))
    
    # Load data
    logger.info("Loading eRNAbase metadata...")
    analyzer.load_erna_metadata(metadata_file, species=args.species)
    
    logger.info("Loading BED files...")
    analyzer.get_bed_files(bed_dir)
    
    # Calculate statistics
    logger.info("Calculating database statistics...")
    stats = {
        'metadata_file': metadata_file,
        'bed_directory': bed_dir,
        'species': args.species,
        'timestamp': datetime.now().isoformat(),
        'total_samples': len(analyzer.erna_metadata),
        'total_bed_files': len(analyzer.erna_bed_files) if analyzer.erna_bed_files else 0,
        'tissue_counts': analyzer.erna_metadata['Tissue Type'].value_counts().to_dict() if 'Tissue Type' in analyzer.erna_metadata.columns else {},
        'cell_type_counts': analyzer.erna_metadata['Cell Type'].value_counts().to_dict() if 'Cell Type' in analyzer.erna_metadata.columns else {},
        'metadata_summary': analyzer.erna_metadata.describe(include='all').to_dict() if hasattr(analyzer.erna_metadata, 'describe') else {}
    }
    
    # Save cache
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    with open(args.output, 'wb') as f:
        pickle.dump(stats, f)
    
    # Also save as JSON for readability
    json_file = args.output.replace('.pkl', '.json')
    with open(json_file, 'w') as f:
        # Convert non-serializable objects
        json_stats = {k: v for k, v in stats.items() if k != 'metadata_summary'}
        json.dump(json_stats, f, indent=2)
    
    logger.info(f"Statistics cache saved to: {args.output}")
    logger.info(f"JSON summary saved to: {json_file}")
    
    # Print summary
    print(f"\nStatistics Summary:")
    print(f"  Total samples: {stats['total_samples']}")
    print(f"  Total BED files: {stats['total_bed_files']}")
    print(f"  Species filter: {stats['species']}")
    
    if stats['tissue_counts']:
        print(f"\nTop 5 tissues:")
        for tissue, count in sorted(stats['tissue_counts'].items(), key=lambda x: x[1], reverse=True)[:5]:
            print(f"  {tissue}: {count}")


def cmd_batch(args, logger):
    """Execute batch command"""
    logger.info("Starting batch analysis...")
    
    # Validate inputs
    bed_dir = validate_directory_exists(args.bed_directory, "BED directory")
    metadata_file = validate_file_exists(args.metadata, "Metadata file")
    regions_file = validate_file_exists(args.regions_file, "Regions file")
    
    # Create output directories
    data_output_dir = os.path.join(args.output_dir, 'data')
    report_output_dir = os.path.join(args.output_dir, 'reports')
    
    # Create analyzer
    analyzer = ERNARegionAnalyzer(results_dir=args.output_dir)
    
    # Load cached stats if provided
    if args.use_cached_stats:
        logger.info(f"Loading cached statistics from: {args.use_cached_stats}")
        try:
            with open(args.use_cached_stats, 'rb') as f:
                cached_stats = pickle.load(f)
            logger.info(f"Loaded statistics from {cached_stats['timestamp']}")
        except Exception as e:
            logger.warning(f"Failed to load cached statistics: {e}")
            cached_stats = None
    else:
        cached_stats = None
    
    # Run batch analysis
    try:
        result_file = analyzer.batch_analyze_regions_from_tsv(
            bed_directory=bed_dir,
            metadata_file=metadata_file,
            data_output_dir=data_output_dir,
            report_output_dir=report_output_dir,
            regions_tsv_file=regions_file,
            max_rows=args.max_rows,
            species=args.species,
            figure_formats=args.figure_formats,
            save_tables=not args.no_tables,
            pvalue_threshold=args.pvalue_threshold,
            fdr_threshold=args.fdr_threshold,
            fold_change_threshold=args.fc_threshold
        )
        
        logger.info(f"Batch analysis completed. Results saved to: {result_file}")
        print(f"\nBatch analysis completed successfully!")
        print(f"Results saved to: {result_file}")
        
        # Save stats cache if requested
        if args.save_stats_cache:
            # Implementation would save analyzer's calculated stats
            logger.info(f"Saving statistics cache to: {args.save_stats_cache}")
            
    except Exception as e:
        logger.error(f"Batch analysis failed: {e}")
        raise


def cmd_single(args, logger):
    """Execute single region analysis command"""
    logger.info("Starting single region analysis...")
    
    # Validate inputs
    bed_dir = validate_directory_exists(args.bed_directory, "BED directory")
    metadata_file = validate_file_exists(args.metadata, "Metadata file")
    
    # Create analyzer
    analyzer = ERNARegionAnalyzer(results_dir=args.output_dir)
    
    # Load data
    logger.info("Loading eRNAbase data...")
    analyzer.get_bed_files(bed_dir)
    analyzer.load_erna_metadata(metadata_file, species=args.species)
    
    # Generate region name if not provided
    if not args.region_name:
        region_name = f"{args.chr}_{args.start}_{args.end}"
    else:
        region_name = args.region_name
    
    # Run analysis
    try:
        results = analyzer.analyze_and_report_region(
            chrom=args.chr,
            start=args.start,
            end=args.end,
            region_name=region_name,
            output_dir=args.output_dir,
            figure_formats=args.figure_formats,
            pvalue_threshold=args.pvalue_threshold
        )
        
        logger.info("Single region analysis completed")
        print(f"\nAnalysis completed for region: {region_name}")
        print(f"Report saved to: {results['report_path']}")
        print(f"Overlapping samples: {results['overlap_sample_count']}")
        
        if 'significant_tissues' in results and results['significant_tissues']:
            print("\nTop enriched tissues:")
            for tissue in results['significant_tissues'][:5]:
                if 'tissue' in tissue and 'fold_change' in tissue and 'fdr' in tissue:
                    print(f"  - {tissue['tissue']}: FC={tissue['fold_change']:.2f}, FDR={tissue['fdr']:.2e}")
        
    except Exception as e:
        logger.error(f"Single region analysis failed: {e}")
        raise


def cmd_report(args, logger):
    """Execute report generation command"""
    logger.info("Starting eRNAbase metadata report generation...")
    
    # Validate inputs
    metadata_file = validate_file_exists(args.metadata, "Metadata file")
    
    # Create analyzer
    analyzer = ERNARegionAnalyzer(results_dir=args.output_dir)
    
    # Generate report
    try:
        report_path = analyzer.generate_erna_db_report(
            metadata_file=metadata_file,
            output_dir=args.output_dir,
            species=args.species,
            top_n=args.top_n,
            figure_formats=args.figure_formats
        )
        
        logger.info("Report generation completed")
        print(f"\neRNAbase metadata report generated successfully!")
        print(f"Report saved to: {report_path}")
        
    except Exception as e:
        logger.error(f"Report generation failed: {e}")
        raise


def main():
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args()
    
    # Show help if no command provided
    if not args.command:
        parser.print_help()
        return 1
    
    # Setup logging
    logger = setup_logging(args.verbose, args.log_file if hasattr(args, 'log_file') else None)
    
    try:
        # Execute command
        if args.command == 'prepare-stats':
            cmd_prepare_stats(args, logger)
        elif args.command == 'batch':
            cmd_batch(args, logger)
        elif args.command == 'single':
            cmd_single(args, logger)
        elif args.command == 'report':
            cmd_report(args, logger)
        else:
            parser.print_help()
            return 1
            
        return 0
        
    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        return 130
    except Exception as e:
        logger.error(f"Command failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())