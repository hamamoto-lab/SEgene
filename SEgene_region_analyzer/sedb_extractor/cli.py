"""
Command-line interface for the SEDb Region Extraction Tool
"""

import argparse
import os
import sys
from pathlib import Path

def parse_region(region_str):
    """Parse region string (e.g., "chr1:1000-2000")"""
    try:
        if ':' not in region_str:
            raise ValueError(f"Invalid region format: {region_str}")
        
        chrom, pos = region_str.split(':')
        if '-' not in pos:
            raise ValueError(f"Invalid position format: {pos}")
            
        start, end = map(int, pos.split('-'))
        if start >= end:
            raise ValueError(f"Start position must be less than end position: {start}-{end}")
            
        return chrom, start, end
    except Exception as e:
        raise ValueError(f"Invalid region format: {region_str}. Use format 'chr:start-end'") from e

def parse_args(args=None):
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(
        description="SEdb Region Extraction Tool - Extract and save super-enhancers from a specified genomic region"
    )
    
    # Required arguments
    parser.add_argument('--bed-file', '-b', required=True, help='Path to SEDb BED file')
    parser.add_argument('--metadata', '-m', required=True, help='Path to SEDb sample metadata file')
    parser.add_argument('--region', '-r', required=True, help='Target region (e.g., "chr1:1000-2000")')
    parser.add_argument('--output', '-o', required=True, help='Path to output directory')
    
    # Optional arguments
    parser.add_argument('--region-name', '-n', help='Region name (if not specified, coordinates will be used)')
    parser.add_argument('--extraction-method', '-e', choices=['pybedtools', 'pandas'], 
                       default='pybedtools', help='Extraction method (default: pybedtools)')
    parser.add_argument('--include-distributions', '-d', action='store_true',
                       help='Include distribution data (e.g., tissue type, biosample type)')
    parser.add_argument('--export-bed', '-x', action='store_true',
                       help='Also output extracted regions as a BED file')
    parser.add_argument('--bed-format', choices=['classic', 'extended'], default='classic',
                       help='Output BED file format (classic: basic format, extended: includes additional information)')
    parser.add_argument('--log-file', '-l', help='Path to log file')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')
    
    # Analysis level option
    parser.add_argument('--analysis-level', choices=['none', 'basic'], default='basic',
                       help='Analysis level (none: extraction only, basic: includes basic statistics)')
    
    return parser.parse_args(args)

def validate_args(args, logger):
    """Validate arguments"""
    # Check if input files exist
    if not os.path.exists(args.bed_file):
        logger.error(f"BED file not found: {args.bed_file}")
        return False
    
    if not os.path.exists(args.metadata):
        logger.error(f"Metadata file not found: {args.metadata}")
        return False
    
    # Validate region
    try:
        parse_region(args.region)
    except ValueError as e:
        logger.error(str(e))
        return False
    
    # Check output directory
    output_path = Path(args.output)

    
    return True
