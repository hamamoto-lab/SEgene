#!/usr/bin/env python3
"""
Basic usage example for SEgene_analyzer_erna

This example demonstrates how to use the ERNARegionAnalyzer to analyze
a genomic region against the eRNAbase database.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from erna_analyzer import ERNARegionAnalyzer

def main():
    """Basic usage example"""
    
    # Initialize the analyzer
    analyzer = ERNARegionAnalyzer(results_dir='example_output')
    
    # Load eRNAbase metadata (CSV or Parquet format)
    # Replace with actual path to your metadata file
    metadata_path = 'path/to/eRNAbase_metadata.csv'
    # analyzer.load_erna_metadata(metadata_path)
    
    # Load BED files directory
    # Replace with actual path to your BED files directory
    bed_files_dir = 'path/to/bed_files/'
    # analyzer.get_bed_files(bed_files_dir)
    
    # Analyze a specific genomic region
    region_results = analyzer.analyze_region(
        chrom='chr1',
        start=1000000,
        end=2000000,
        region_name='example_region'
    )
    
    # Generate comprehensive HTML report
    # analyzer.generate_region_analysis_report()
    
    print("Analysis completed! Check the 'example_output' directory for results.")

if __name__ == "__main__":
    main()