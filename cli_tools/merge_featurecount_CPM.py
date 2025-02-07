#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script merges multiple featureCounts output files and calculates CPM values based on the
BAM read counts provided in a CSV file. The final CPM data is saved as a pickle file.
Optionally, the data can also be exported as CSV and/or TSV.

Usage:
    python merge_featurecount_CPM.py [--featurecounts_dir PATH]
                                     [--bam_read_counts_csv PATH]
                                     [--output_pickle PATH]
                                     [--export_csv]
                                     [--export_tsv]

Default paths (if not provided):
    --featurecounts_dir:   ./output/featurecounts_output
    --bam_read_counts_csv: ./output/bam_read_counts.csv
    --output_pickle:       ./output/merge_featurecount_CPM.pkl
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np

def merge_featurecounts(featurecounts_dir):
    """
    Reads all featureCounts files (matching "featurecounts_*.txt") in the specified directory,
    extracts the required columns, renames the count column with the sample name, and merges
    all dataframes on the key columns (Geneid, Chr, Start, End, Strand, Length).

    Parameters:
        featurecounts_dir (str): Directory containing featureCounts output files.

    Returns:
        pd.DataFrame: Merged dataframe of featureCounts data, or None if no valid files are found.
    """
    file_pattern = os.path.join(featurecounts_dir, "featurecounts_*.txt")
    files = glob.glob(file_pattern)
    if not files:
        print("No featureCounts files found in", featurecounts_dir)
        return None

    dfs = []
    for file_path in files:
        # Extract sample name from the file name.
        # For example: "featurecounts_TM_input_R1.mLb.clN.sorted.txt" -> "TM_input"
        basename = os.path.basename(file_path)
        sample_name = basename.replace("featurecounts_", "").split("_R1")[0]
        
        # Read the file excluding comment lines (lines starting with '#')
        df = pd.read_csv(file_path, sep='\t', comment='#')
        
        # Identify the column that contains ".bam" in its name
        count_columns = [col for col in df.columns if ".bam" in col]
        if not count_columns:
            print(f"No .bam column found in {file_path}. Skipping.")
            continue
        count_column = count_columns[0]
        
        # Extract only the necessary columns
        df = df[['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', count_column]]
        # Rename the count column to the sample name
        df.rename(columns={count_column: sample_name}, inplace=True)
        
        dfs.append(df)
    
    if not dfs:
        print("No valid featureCounts dataframes were obtained.")
        return None

    # Merge all dataframes on the key columns
    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(merged_df, df,
                             on=['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length'],
                             how='outer')
    return merged_df

def calculate_CPM(merged_df, bam_read_counts_csv):
    """
    For each sample in merged_df, calculates the CPM value using the total read counts from
    bam_read_counts_csv. The calculation is: CPM = (count / total_reads) * 1e6.
    A new column (sample_name + '_CPM') is added for each sample.
    Finally, only meta information and CPM columns are extracted and returned.

    Parameters:
        merged_df (pd.DataFrame): Merged featureCounts dataframe.
        bam_read_counts_csv (str): Path to the CSV file containing BAM read counts.

    Returns:
        pd.DataFrame: DataFrame containing meta information and CPM columns,
                      with 'Geneid' set as the index.
    """
    # Read BAM read counts CSV
    read_counts_df = pd.read_csv(bam_read_counts_csv)
    # Remove trailing '_R1' and beyond from sample names, if present
    read_counts_df['Sample'] = read_counts_df['Sample'].str.replace(r'_R1.*$', '', regex=True)
    read_counts_dict = dict(zip(read_counts_df['Sample'], read_counts_df['ReadCount']))
    
    meta_columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']
    cpm_data = {}
    # For each sample column in merged_df, if the sample exists in bam_read_counts,
    # calculate CPM = (count / total_reads) * 1e6
    for col in merged_df.columns:
        if col in meta_columns:
            continue
        if col in read_counts_dict:
            total_reads = read_counts_dict[col]
            cpm_data[col + '_CPM'] = (merged_df[col] / total_reads) * 1e6
        else:
            print(f"Sample '{col}' not found in bam read counts. Skipping CPM calculation for it.")
    
    # Create a DataFrame from the CPM dictionary and concatenate with merged_df
    cpm_df = pd.DataFrame(cpm_data)
    merged_df = pd.concat([merged_df, cpm_df], axis=1)
    
    # Select meta columns and CPM columns only
    cpm_columns = [col for col in merged_df.columns if col.endswith('_CPM')]
    selected_columns = meta_columns + cpm_columns
    merged_df_cpm_only = merged_df[selected_columns]
    
    # Optionally simplify column names by removing '_CPM'
    new_column_names = {col: col.replace('_CPM', '') for col in cpm_columns}
    merged_df_cpm_only = merged_df_cpm_only.rename(columns=new_column_names)
    
    # Set 'Geneid' as the index
    merged_df_cpm_only.set_index('Geneid', inplace=True)
    
    return merged_df_cpm_only

def main():
    parser = argparse.ArgumentParser(
        description="Merge featureCounts files and calculate CPM values based on BAM read counts."
    )
    parser.add_argument("--featurecounts_dir", type=str, default="./output/featurecounts_output",
                        help="Directory containing featureCounts output files (default: ./output/featurecounts_output)")
    parser.add_argument("--bam_read_counts_csv", type=str, default="./output/bam_read_counts.csv",
                        help="CSV file with BAM read counts (default: ./output/bam_read_counts.csv)")
    parser.add_argument("--output_pickle", type=str, default="./output/merge_featurecount_CPM.pkl",
                        help="Output pickle file for merged CPM data (default: ./output/merge_featurecount_CPM.pkl)")
    parser.add_argument("--export_csv", action="store_true",
                        help="If set, also export the merged CPM data as a CSV file (same basename as output_pickle)")
    parser.add_argument("--export_tsv", action="store_true",
                        help="If set, also export the merged CPM data as a TSV file (same basename as output_pickle)")
    
    args = parser.parse_args()
    
    featurecounts_dir = args.featurecounts_dir
    bam_read_counts_csv = args.bam_read_counts_csv
    output_pickle = args.output_pickle
    
    # 1. Merge featureCounts files
    print("Merging featureCounts files …")
    merged_df = merge_featurecounts(featurecounts_dir)
    if merged_df is None:
        print("Merging featureCounts files failed.")
        return
    print("Merged featureCounts data (first 5 rows):")
    print(merged_df.head())
    
    # 2. Calculate CPM values using BAM read counts
    print("Calculating CPM values …")
    merged_df_cpm_only = calculate_CPM(merged_df, bam_read_counts_csv)
    print("CPM-only merged data (first 5 rows):")
    print(merged_df_cpm_only.head())
    
    # 3. Save the final CPM data as a pickle file
    os.makedirs(os.path.dirname(output_pickle), exist_ok=True)
    merged_df_cpm_only.to_pickle(output_pickle)
    print(f"Saved merged CPM data to {output_pickle}")

    # Optionally export as CSV
    if args.export_csv:
        output_csv = os.path.splitext(output_pickle)[0] + ".csv"
        merged_df_cpm_only.to_csv(output_csv)
        print(f"Exported merged CPM data to CSV: {output_csv}")

    # Optionally export as TSV
    if args.export_tsv:
        output_tsv = os.path.splitext(output_pickle)[0] + ".tsv"
        merged_df_cpm_only.to_csv(output_tsv, sep='\t')
        print(f"Exported merged CPM data to TSV: {output_tsv}")

if __name__ == "__main__":
    main()
