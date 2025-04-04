"""
sedb_analyzer.cli - Command-line interface module for sedb_analyzer

This module provides CLI functionality for the SEdb analysis tools.
"""

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, Any, Optional, Union, List, Tuple

# Import common logging utility
from sedb_common.logging import get_module_logger

def parse_args():
    """Parse command line arguments for the SEdb analyzer tool"""
    parser = argparse.ArgumentParser(
        description="SEdb Analyzer Tool - Performs statistical analysis and enrichment testing on extracted SE regions"
    )
    
    # Create argument groups for better organization
    input_group = parser.add_argument_group('Input/Output Arguments')
    enrich_group = parser.add_argument_group('Enrichment Analysis Options')
    output_group = parser.add_argument_group('Output Options')
    other_group = parser.add_argument_group('Other Options')
    
    # Required input/output arguments
    input_group.add_argument('--input', '-i', required=True, 
                      help='Directory path containing extracted region data package')
    input_group.add_argument('--background', '-b', required=True, 
                      help='Directory path containing preprocessed background data package')
    input_group.add_argument('--output', '-o', required=True, 
                      help='Directory path for output analysis results')
    
    # Enrichment analysis options
    enrich_group.add_argument('--enrichment-method', '--method', choices=['fisher', 'chi2'], default='fisher',
                        help='Statistical method for enrichment testing (default: fisher)')
    enrich_group.add_argument('--correction-method', '--correction', default='fdr_bh',
                        help='Multiple testing correction method (default: fdr_bh)')
    enrich_group.add_argument('--pvalue-threshold', '--pvalue', type=float, default=0.05,
                        help='P-value threshold for significance (default: 0.05)')
    enrich_group.add_argument('--fc-threshold', '--fc', type=float, default=1.5,
                        help='Fold-change threshold for significance (default: 1.5)')
    enrich_group.add_argument('--min-region-count', '--min-count', type=int, default=1,
                        help='Minimum region count threshold for analysis (default: 1)')
    
    # Output options
    output_group.add_argument('--figures', action='store_true', default=True,
                       help='Generate and save figures (default)')
    output_group.add_argument('--no-figures', action='store_false', dest='figures',
                       help='Do not generate figures')
    output_group.add_argument('--html-report', action='store_true', default=True,
                       help='Generate HTML report (default)')
    output_group.add_argument('--no-html-report', action='store_false', dest='html_report',
                       help='Do not generate HTML report')
    output_group.add_argument('--data-formats-machine', nargs='+', 
                       choices=['parquet', 'json', 'pickle'], default=['parquet', 'json'],
                       help='Output formats for machine-readable data (default: parquet json)')
    output_group.add_argument('--data-formats-human', nargs='+',
                       choices=['csv', 'txt', 'tsv'], default=['csv', 'txt'],
                       help='Output formats for human-readable data (default: csv txt)')
    output_group.add_argument('--image-formats', nargs='+',
                       choices=['png', 'svg', 'pdf'], default=['png', 'svg'],
                       help='Image formats for figures (default: png svg)')
    
    # Other options
    other_group.add_argument('--log-file', '-l', 
                      help='Path to log file (if not provided, logs to console only)')
    other_group.add_argument('--debug', '-d', action='store_true', 
                      help='Enable debug mode with detailed logging')
    
    return parser.parse_args()

def setup_analyzer_logger(log_level: int = logging.INFO, log_file: Optional[str] = None) -> logging.Logger:
    """
    Set up a logger for the SEdb analyzer
    
    Parameters
    ----------
    log_level : int
        Logging level (default: logging.INFO)
    log_file : str, optional
        Path to log file (default: None, logs to console only)
        
    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    # Use the common logging utility
    console_output = True
    file_output = log_file is not None
    
    logger = get_module_logger(
        module_name="sedb_analyzer",
        console=console_output,
        file_output=file_output,
        log_level=log_level
    )
    
    # If log_file is provided, we need to manually add a file handler since
    # get_module_logger doesn't allow specifying the file path directly
    if file_output:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def validate_args(args, logger: logging.Logger) -> bool:
    """
    Validate command line arguments
    
    Parameters
    ----------
    args : Namespace
        Parsed command line arguments
    logger : logging.Logger
        Logger instance
        
    Returns
    -------
    bool
        True if arguments are valid, False otherwise
    """
    # Check input directory
    if not os.path.exists(args.input):
        logger.error(f"Input directory does not exist: {args.input}")
        return False
    
    if not os.path.isdir(args.input):
        logger.error(f"Input path is not a directory: {args.input}")
        return False
    
    # Check background directory
    if not os.path.exists(args.background):
        logger.error(f"Background directory does not exist: {args.background}")
        return False
    
    if not os.path.isdir(args.background):
        logger.error(f"Background path is not a directory: {args.background}")
        return False
    
    # Prepare output directory
    output_path = Path(args.output)
    if output_path.exists() and not output_path.is_dir():
        logger.error(f"Output path exists but is not a directory: {args.output}")
        return False
    
    # Create output directory if it doesn't exist
    if not output_path.exists():
        try:
            output_path.mkdir(parents=True)
            logger.info(f"Created output directory: {args.output}")
        except Exception as e:
            logger.error(f"Failed to create output directory: {e}")
            return False
    
    # Validate correction method
    valid_correction_methods = [
        'bonferroni', 'sidak', 'holm-sidak', 'holm', 
        'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 
        'fdr_tsbh', 'fdr_tsbky'
    ]
    
    if args.correction_method not in valid_correction_methods:
        logger.warning(f"Warning: Correction method '{args.correction_method}' may not be supported. "
                      f"Valid methods: {', '.join(valid_correction_methods)}")
    
    # Validate p-value threshold
    if not (0 < args.pvalue_threshold < 1):
        logger.error(f"P-value threshold must be between 0 and 1: {args.pvalue_threshold}")
        return False
    
    # Validate fold-change threshold
    if args.fc_threshold <= 0:
        logger.error(f"Fold-change threshold must be positive: {args.fc_threshold}")
        return False
    
    # Validate min region count
    if args.min_region_count < 0:
        logger.error(f"Minimum region count must be non-negative: {args.min_region_count}")
        return False
    
    return True