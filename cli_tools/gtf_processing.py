#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module contains functions for processing a pickle file containing region information
and converting it into a featureCounts-compatible GTF file.
It also provides a function to extract the top-n rows (sorted by se_count or sample_count)
from the input DataFrame.
"""

import pandas as pd
import csv

### 基本関数群 ###

def add_region_data_to_info(info_df):
    """
    Extracts region information (chr, start, end) from the DataFrame index
    (assumed format: "chr_start_end") and adds them as new columns.
    Returns a DataFrame with 'chr', 'start', and 'end' as the first columns.
    """
    df = info_df.copy()
    df['chr'] = df.index.str.split('_').str[0]
    df['start'] = df.index.str.split('_').str[1].astype(int)
    df['end'] = df.index.str.split('_').str[2].astype(int)
    columns = ['chr', 'start', 'end'] + [col for col in df.columns if col not in ['chr', 'start', 'end']]
    df = df[columns]
    return df

def sort_se_merge_SE_count_df(df):
    """
    Sorts the DataFrame by the 'se_count' column in descending order, if present.
    Otherwise returns the DataFrame as is.
    """
    if 'se_count' in df.columns:
        return df.sort_values(by='se_count', ascending=False)
    else:
        return df

def sort_se_merge_SAMPLE_count_df(df):
    """
    Sorts the DataFrame by the 'sample_count' column in descending order, if present.
    Otherwise returns the DataFrame as is.
    """
    if 'sample_count' in df.columns:
        return df.sort_values(by='sample_count', ascending=False)
    else:
        return df

def get_top_regions_df(info_df, top_n=20, sort_by="se"):
    """
    Processes the input DataFrame (from pickle) to extract the top_n rows.
    
    Steps:
      1. Apply add_region_data_to_info to extract region information.
      2. Sort the DataFrame by:
         - 'se_count' (if sort_by=="se", default) or
         - 'sample_count' (if sort_by=="sample")
      3. Take the top_n rows.
    
    Returns:
      The raw top regions DataFrame (with all original columns) without further GTF conversion.
    """
    df_regions = add_region_data_to_info(info_df)
    if sort_by == "sample":
        sorted_df = sort_se_merge_SAMPLE_count_df(df_regions)
    else:
        sorted_df = sort_se_merge_SE_count_df(df_regions)
    top_df = sorted_df.head(top_n).copy()
    return top_df

def generate_gtf_df_from_top_df(top_df):
    """
    Converts the top regions DataFrame into a GTF-formatted DataFrame.
    
    Note: In the raw top_df, the 'start' column is 0-indexed (or as in the input);
    for GTF format, we need 1-based coordinates. Therefore, add 1 to the 'start' values.
    
    The resulting DataFrame will contain the following columns:
      [chr, source, feature, start, end, score, strand, frame, attribute]
    
    - source is fixed to "SuperEnhancer"
    - feature is fixed to "exon"
    - score, strand, frame are set to '.'
    - attribute is constructed using the index (used as gene_id and transcript_id)
    """
    # Create a copy of top_df to avoid modifying the raw data
    top = top_df.copy()
    # Adjust start to 1-based for GTF
    top['start'] = top['start'].astype(int) + 1
    gtf_df = top[['chr', 'start', 'end']].copy()
    gtf_df['source'] = 'SuperEnhancer'
    gtf_df['feature'] = 'exon'
    gtf_df['score'] = '.'
    gtf_df['strand'] = '.'
    gtf_df['frame'] = '.'
    # Use the DataFrame index as gene_id and transcript_id
    gtf_df['attribute'] = 'gene_id "' + top.index.astype(str) + '"; transcript_id "' + top.index.astype(str) + '";'
    # Reorder columns to match GTF format
    gtf_df = gtf_df[['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']]
    return gtf_df

def sort_gtf_df(gtf_df):
    """
    Sorts the GTF DataFrame by chromosome, start, and end.
    """
    return gtf_df.sort_values(by=['chr', 'start', 'end'])

def write_gtf_from_df(gtf_df, output_path, quoting=csv.QUOTE_NONE):
    """
    Writes the GTF DataFrame to a file in tab-separated format without header or index.
    """
    gtf_df.to_csv(output_path, sep='\t', header=False, index=False, quoting=quoting)

def modify_gtf_line(line):
    """
    Modifies a single GTF line by forcing the feature to 'exon' and reconstructing the attribute field.
    """
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return line.strip()
    fields[2] = 'exon'
    chrom = fields[0]
    start = fields[3]
    end = fields[4]
    gene_id = f"{chrom}_{start}_{end}"
    gene_name = f"{chrom}_{start}"
    new_attrs = f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{gene_id}";'
    fields[8] = new_attrs
    return '\t'.join(fields)

def modify_gtf_file(input_path, output_path):
    """
    Reads a GTF file line by line, modifies each line using modify_gtf_line,
    and writes the result to the output file.
    """
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            modified_line = modify_gtf_line(line)
            outfile.write(modified_line + '\n')
