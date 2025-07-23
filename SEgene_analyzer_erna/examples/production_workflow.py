#!/usr/bin/env python3
"""
Production workflow example for SEgene_analyzer_erna

This script demonstrates a complete analysis workflow from data preparation 
to final result generation, following the patterns observed in actual usage.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
from datetime import datetime

# Add src to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from erna_analyzer import ERNARegionAnalyzer

def setup_logging(log_dir='logs'):
    """Setup logging configuration"""
    os.makedirs(log_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(log_dir, f'erna_analysis_{timestamp}.log')
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)

def validate_input_files(bed_dir, metadata_file, regions_file):
    """Validate that all required input files exist"""
    errors = []
    
    if not os.path.exists(bed_dir):
        errors.append(f"BED directory not found: {bed_dir}")
    elif not os.listdir(bed_dir):
        errors.append(f"BED directory is empty: {bed_dir}")
    
    if not os.path.exists(metadata_file):
        errors.append(f"Metadata file not found: {metadata_file}")
    
    if not os.path.exists(regions_file):
        errors.append(f"Regions file not found: {regions_file}")
    
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        sys.exit(1)

def run_database_overview(analyzer, output_dir, species='human'):
    """Generate eRNAbase database overview report"""
    logger = logging.getLogger(__name__)
    
    logger.info("Starting database overview analysis...")
    
    db_report_dir = os.path.join(output_dir, 'database_overview')
    os.makedirs(db_report_dir, exist_ok=True)
    
    try:
        report_path = analyzer.generate_erna_db_report(
            output_dir=db_report_dir,
            species=species
        )
        logger.info(f"Database overview report generated: {report_path}")
        return report_path
    except Exception as e:
        logger.error(f"Failed to generate database overview: {e}")
        raise

def run_batch_analysis(analyzer, bed_dir, metadata_file, regions_file, 
                      output_dir, max_rows=None, species='human'):
    """Run batch analysis on specified regions"""
    logger = logging.getLogger(__name__)
    
    logger.info("Starting batch analysis...")
    
    data_output_dir = os.path.join(output_dir, 'region_analysis', 'data')
    report_output_dir = os.path.join(output_dir, 'region_analysis', 'reports')
    
    os.makedirs(data_output_dir, exist_ok=True)
    os.makedirs(report_output_dir, exist_ok=True)
    
    try:
        result_file = analyzer.batch_analyze_regions_from_tsv(
            bed_directory=bed_dir,
            metadata_file=metadata_file,
            data_output_dir=data_output_dir,
            report_output_dir=report_output_dir,
            regions_tsv_file=regions_file,
            max_rows=max_rows,
            species=species,
            figure_formats=['png', 'svg'],
            save_tables=True,
            # Statistical parameters
            # pvalue_threshold=0.05,
            # fdr_threshold=0.05
        )
        
        logger.info(f"Batch analysis completed. Results saved to: {result_file}")
        return result_file
        
    except Exception as e:
        logger.error(f"Batch analysis failed: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(
        description='eRNAbase region analysis production workflow',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--bed-dir', required=True,
                       help='Directory containing eRNAbase BED files')
    parser.add_argument('--metadata', required=True,
                       help='eRNAbase metadata file (CSV or Parquet)')
    parser.add_argument('--regions', required=True,
                       help='TSV file containing regions to analyze')
    parser.add_argument('--output-dir', default='results',
                       help='Output directory for results')
    parser.add_argument('--max-rows', type=int, default=None,
                       help='Maximum number of regions to process')
    parser.add_argument('--species', default='human',
                       help='Species to analyze')
    parser.add_argument('--log-dir', default='logs',
                       help='Directory for log files')
    parser.add_argument('--skip-overview', action='store_true',
                       help='Skip database overview generation')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.log_dir)
    logger.info("Starting eRNAbase analysis workflow")
    logger.info(f"Arguments: {vars(args)}")
    
    # Validate input files
    validate_input_files(args.bed_dir, args.metadata, args.regions)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Initialize analyzer
        analyzer = ERNARegionAnalyzer(results_dir=args.output_dir)
        
        # Load metadata
        logger.info(f"Loading metadata from: {args.metadata}")
        analyzer.load_erna_metadata(args.metadata, species=args.species)
        
        # Generate database overview (optional)
        if not args.skip_overview:
            db_report = run_database_overview(
                analyzer, args.output_dir, args.species
            )
        
        # Run batch analysis
        result_file = run_batch_analysis(
            analyzer=analyzer,
            bed_dir=args.bed_dir,
            metadata_file=args.metadata,
            regions_file=args.regions,
            output_dir=args.output_dir,
            max_rows=args.max_rows,
            species=args.species
        )
        
        logger.info("Analysis workflow completed successfully!")
        logger.info(f"Main results file: {result_file}")
        
        # Summary of outputs
        print("\n" + "="*60)
        print("ANALYSIS COMPLETED SUCCESSFULLY")
        print("="*60)
        print(f"Results directory: {args.output_dir}")
        print(f"Main results file: {result_file}")
        if not args.skip_overview:
            print(f"Database overview: {args.output_dir}/database_overview/")
        print(f"Region analysis: {args.output_dir}/region_analysis/")
        print(f"Log file: {args.log_dir}/erna_analysis_*.log")
        print("="*60)
        
    except Exception as e:
        logger.error(f"Analysis workflow failed: {e}")
        logger.exception("Full traceback:")
        sys.exit(1)

if __name__ == "__main__":
    main()