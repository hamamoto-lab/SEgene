#!/usr/bin/env python3
"""
Batch analysis example for SEgene_analyzer_erna

This example demonstrates how to analyze multiple genomic regions
using a TSV file input.
"""

import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from erna_analyzer import ERNARegionAnalyzer

def create_example_regions_file():
    """Create an example regions TSV file"""
    regions_content = """region_name\tchrom\tstart\tend
region_1\tchr1\t1000000\t2000000
region_2\tchr2\t5000000\t6000000
region_3\tchr3\t10000000\t11000000
"""
    
    with open('example_regions.tsv', 'w') as f:
        f.write(regions_content)
    
    print("Created example_regions.tsv")

def main():
    """Batch analysis example"""
    
    # Create example regions file
    create_example_regions_file()
    
    # Initialize the analyzer
    analyzer = ERNARegionAnalyzer(results_dir='batch_output')
    
    # Load eRNAbase metadata (uncomment when you have actual data)
    # analyzer.load_erna_metadata('path/to/eRNAbase_metadata.csv')
    # analyzer.get_bed_files('path/to/bed_files/')
    
    # Perform batch analysis
    # analyzer.batch_analyze_from_tsv('example_regions.tsv')
    
    print("Batch analysis setup completed!")
    print("To run with real data:")
    print("1. Uncomment the load_erna_metadata and get_bed_files lines")
    print("2. Update the paths to your actual data files")
    print("3. Uncomment the batch_analyze_from_tsv line")

if __name__ == "__main__":
    main()