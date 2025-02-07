#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script reads an already generated merged CPM pickle file,
adjusts the 'Start' coordinate by subtracting 1 to convert from 1-indexed to 0-indexed,
and recreates the 'Geneid' column (formatted as "Chr_Start_End").
The adjusted data is then saved as a new pickle file.
Optionally, the data can also be exported as CSV and/or TSV.

Usage:
    python convert_featurecount_CPM_zero_indexed.py [--input_pickle PATH]
                                                    [--output_pickle PATH]
                                                    [--export_csv]
                                                    [--export_tsv]

Default paths:
    --input_pickle:  ./output/merge_featurecount_CPM.pkl
    --output_pickle: ./output/merge_featurecount_CPM_zero_indexed.pkl
"""

import pandas as pd
import os
import re
import argparse

def adjust_index_string(index_str):
    """
    Given an index string of the format "chrX_start_end",
    subtract 1 from the start value to convert from 1-indexed to 0-indexed.

    Parameters:
        index_str (str): The original index string.

    Returns:
        str: The adjusted index string.
    """
    m = re.match(r'^(chr\S+)_(\d+)_(\d+)$', index_str)
    if m:
        chr_val = m.group(1)
        start_val = int(m.group(2))
        end_val = m.group(3)
        # Subtract 1 from the start value for 0-indexing
        return f"{chr_val}_{start_val - 1}_{end_val}"
    else:
        # If pattern does not match, return the original string
        return index_str

def convert_to_zero_indexed(input_pickle, output_pickle):
    """
    Reads the merged CPM pickle file, converts the 'Start' column from 1-indexed to 0-indexed,
    adjusts the index strings, and reconstitutes the 'Geneid' column based on 'Chr', 'Start', and 'End'.
    The result is saved as a pickle file at output_pickle.
    
    Parameters:
        input_pickle (str): Path to the input pickle file.
        output_pickle (str): Path to the output pickle file.
    
    Returns:
        pd.DataFrame: The adjusted DataFrame.
    """
    # 1. Read the merged CPM pickle file
    df = pd.read_pickle(input_pickle)
    
    # 2. Adjust the 'Start' column from 1-indexed to 0-indexed if it exists
    if 'Start' in df.columns:
        df['Start'] = df['Start'] - 1
    else:
        print("Warning: 'Start' column not found.")
    
    # 3. If the index is of string type, adjust it using the adjust_index_string function
    if df.index.dtype == object:
        df.index = df.index.map(adjust_index_string)
    
    # 4. Reconstruct the 'Geneid' column using 'Chr', 'Start', 'End' if available
    if set(['Chr', 'Start', 'End']).issubset(df.columns):
        df['Geneid'] = df['Chr'] + '_' + df['Start'].astype(str) + '_' + df['End'].astype(str)
        df.set_index('Geneid', inplace=True)
    else:
        print("Warning: One or more of 'Chr', 'Start', 'End' columns are missing. Skipping index reconstruction.")
    
    # 5. Create the output directory if it does not exist
    out_dir = os.path.dirname(output_pickle)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    # 6. Save the adjusted dataframe as a pickle file
    df.to_pickle(output_pickle)
    print(f"Saved zero-indexed data to {output_pickle}")
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Convert merged featureCounts CPM pickle data from 1-indexed to 0-indexed coordinates."
    )
    parser.add_argument("--input_pickle", type=str, default="./output/merge_featurecount_CPM.pkl",
                        help="Path to input merged CPM pickle file (default: ./output/merge_featurecount_CPM.pkl)")
    parser.add_argument("--output_pickle", type=str, default="./output/merge_featurecount_CPM_zero_indexed.pkl",
                        help="Path to output zero-indexed CPM pickle file (default: ./output/merge_featurecount_CPM_zero_indexed.pkl)")
    parser.add_argument("--export_csv", action="store_true",
                        help="If set, also export the zero-indexed data as a CSV file (same basename as output_pickle)")
    parser.add_argument("--export_tsv", action="store_true",
                        help="If set, also export the zero-indexed data as a TSV file (same basename as output_pickle)")
    
    args = parser.parse_args()
    
    # Convert the data to 0-indexed and obtain the adjusted DataFrame
    df = convert_to_zero_indexed(args.input_pickle, args.output_pickle)
    
    # Optionally export as CSV
    if args.export_csv:
        output_csv = os.path.splitext(args.output_pickle)[0] + ".csv"
        df.to_csv(output_csv)
        print(f"Exported zero-indexed data to CSV: {output_csv}")
    
    # Optionally export as TSV
    if args.export_tsv:
        output_tsv = os.path.splitext(args.output_pickle)[0] + ".tsv"
        df.to_csv(output_tsv, sep='\t')
        print(f"Exported zero-indexed data to TSV: {output_tsv}")

if __name__ == "__main__":
    main()
