#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script reads a pickle file containing region information, processes it to generate a
featureCounts-compatible GTF file, and writes the final GTF file.
It also saves the top-n (default 20) regions DataFrame (raw, sorted by se_count or sample_count)
for later inspection.

Command line arguments:
  --input_pickle : Path to the input pickle file.
                   (Default: "temp_full_df_info.pkl", assumed to be in the same directory as the script)
  --top_n        : Number of top regions to use (default: 20)
  --sort_by      : Sorting method: "se" (default) for se_count, "sample" for sample_count.
  --output_gtf   : Path to the final output GTF file (default: output/modified_for_featurecounts.gtf)
  --unsorted_gtf : Temporary file path for the intermediate GTF file (default: output/unsorted_temp.gtf)
  --output_top_df: Path to save the top regions DataFrame (default: output/top_regions_df.pkl)
"""

import argparse
import os
import pandas as pd
from gtf_processing import (
    get_top_regions_df,
    generate_gtf_df_from_top_df,
    sort_gtf_df,
    write_gtf_from_df,
    modify_gtf_file
)

def main():
    parser = argparse.ArgumentParser(
        description="Generate a featureCounts-compatible GTF file and save the top-n regions DataFrame."
    )
    parser.add_argument("--input_pickle", type=str, default="temp_full_df_info.pkl",
                        help='Path to input pickle file (default: "temp_full_df_info.pkl", assumed to be in the same directory as this script)')
    parser.add_argument("--top_n", type=int, default=20,
                        help="Number of top regions to use (default: 20)")
    parser.add_argument("--sort_by", type=str, choices=["se", "sample"], default="se",
                        help='Sorting method: "se" (default) for se_count, "sample" for sample_count')
    parser.add_argument("--output_gtf", type=str, default="output/modified_for_featurecounts.gtf",
                        help="Path to output GTF file (default: output/modified_for_featurecounts.gtf)")
    parser.add_argument("--unsorted_gtf", type=str, default="output/unsorted_temp.gtf",
                        help="Temporary file for intermediate GTF (default: output/unsorted_temp.gtf)")
    parser.add_argument("--output_top_df", type=str, default="output/top_regions_df.pkl",
                        help="Path to save the top regions DataFrame (default: output/top_regions_df.pkl)")
    
    args = parser.parse_args()
    
    # Resolve input file path relative to the script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    input_pickle_path = os.path.join(script_dir, args.input_pickle)
    
    # Read the input pickle file
    info_df = pd.read_pickle(input_pickle_path)
    
    # Get the top regions DataFrame (raw, sorted by the specified method)
    raw_top_df = get_top_regions_df(info_df, top_n=args.top_n, sort_by=args.sort_by)
    
    # Save the raw top regions DataFrame for inspection
    output_top_df_path = args.output_top_df
    os.makedirs(os.path.dirname(output_top_df_path), exist_ok=True)
    raw_top_df.to_pickle(output_top_df_path)
    print(f"Top regions DataFrame (n={args.top_n}) saved to {output_top_df_path}")
    
    # Convert the raw top regions DataFrame into a GTF-formatted DataFrame
    gtf_df = generate_gtf_df_from_top_df(raw_top_df)
    
    # Optionally sort the GTF DataFrame by chr, start, and end
    sorted_gtf_df = sort_gtf_df(gtf_df)
    
    # Ensure output directories exist
    unsorted_gtf_path = args.unsorted_gtf
    output_gtf_path = args.output_gtf
    os.makedirs(os.path.dirname(unsorted_gtf_path), exist_ok=True)
    os.makedirs(os.path.dirname(output_gtf_path), exist_ok=True)
    
    # Write the sorted (but not yet attribute-modified) GTF to a temporary file
    write_gtf_from_df(sorted_gtf_df, unsorted_gtf_path)
    print(f"Temporary GTF file written to {unsorted_gtf_path}")
    
    # Modify the temporary GTF file (adjust attribute fields) to produce final GTF output
    modify_gtf_file(unsorted_gtf_path, output_gtf_path)
    print(f"Final GTF file written to {output_gtf_path}")

if __name__ == "__main__":
    main()
