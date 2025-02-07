#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script provides:
1) create_gene_se_correlation_chart_log10 (Pearson correlation)
2) process_se_gene_correlation

It generates scatter plots of log10(TPM) vs. log10(CPM) for each SE–Gene pair
and saves them in SVG format by default. If --save_png is specified, PNG files will be
saved as well.

Usage example:
    python correlation_analysis.py \
        --top_regions_df <path_to_top_regions_df.pkl> \
        --se_count 10 \
        --output_dir <output_directory> \
        --merge_featurecount_pkl <path_to_merge_featurecount_CPM_zero_indexed.pkl> \
        --rna_csv <path_to_gene_tpm.csv> \
        --full_info_pkl <path_to_full_info.pkl> \
        --use_all_samples \
        --save_png
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def create_gene_se_correlation_chart_log10(
    input_se,
    input_gene,
    merge_featurecount_pkl,
    rna_csv,
    input_df_full_info,
    use_all_samples=False
):
    """
    Creates a log10 scatter plot (TPM vs. CPM) for a given SE–Gene pair,
    and computes the Pearson correlation coefficient.

    Parameters:
        input_se (str): SE name (peak identifier).
        input_gene (str): Gene name.
        merge_featurecount_pkl (str): Path to a pickle file containing CPM data for SE.
        rna_csv (str): Path to a CSV file containing TPM data for Genes.
        input_df_full_info (str): Path to a pickle with peak/sample info (SAMPLE_distinct, etc.).
        use_all_samples (bool): If True, uses all samples in common between TPM/CPM data;
                                otherwise only uses those listed in SAMPLE_distinct.

    Returns:
        (fig, ax, pearson_corr):
          fig (matplotlib Figure),
          ax (matplotlib Axes),
          pearson_corr (float).
        If no data or an error occurs, returns (None, None, None).
    """
    temp_full_df_info = pd.read_pickle(input_df_full_info)

    if input_se not in temp_full_df_info.index:
        print(f"[Warning] SE '{input_se}' not found in full info pickle. Skipping.")
        return None, None, None

    if "SAMPLE_distinct" not in temp_full_df_info.columns:
        print("[Warning] 'SAMPLE_distinct' column not found. Skipping.")
        return None, None, None

    sample_list_str = temp_full_df_info.loc[input_se, "SAMPLE_distinct"]
    if not isinstance(sample_list_str, str):
        print("[Warning] SAMPLE_distinct is not a string. Skipping.")
        return None, None, None
    sample_list = sample_list_str.split(",")

    df_rna_csv = pd.read_csv(rna_csv, index_col=0)
    df_SE_CPM = pd.read_pickle(merge_featurecount_pkl)

    if input_gene not in df_rna_csv.index:
        print(f"[Warning] Gene '{input_gene}' not found in RNA CSV. Skipping.")
        return None, None, None
    if input_se not in df_SE_CPM.index:
        print(f"[Warning] SE '{input_se}' not found in SE CPM PKL. Skipping.")
        return None, None, None

    print(f"Processing correlation: SE = {input_se}, Gene = {input_gene}")

    series_tmp_x_TPM = df_rna_csv.loc[input_gene][4:]
    series_tmp_y_CPM = df_SE_CPM.loc[input_se][5:]

    if use_all_samples:
        common_samples = list(set(series_tmp_x_TPM.index) & set(series_tmp_y_CPM.index))
    else:
        common_samples = list(
            set(sample_list) & set(series_tmp_x_TPM.index) & set(series_tmp_y_CPM.index)
        )

    if not common_samples:
        print(f"No common samples for SE '{input_se}' and gene '{input_gene}'. Skipping.")
        return None, None, None

    filtered_x_TPM = series_tmp_x_TPM[common_samples].astype(float)
    filtered_y_CPM = series_tmp_y_CPM[common_samples].astype(float)
    log_filtered_x_TPM = np.log10(filtered_x_TPM + 1)
    log_filtered_y_CPM = np.log10(filtered_y_CPM + 1)

    pearson_corr = log_filtered_x_TPM.corr(log_filtered_y_CPM)
    print(f"Pearson correlation: {pearson_corr:.3f}")

    try:
        slope, intercept = np.polyfit(log_filtered_x_TPM, log_filtered_y_CPM, 1)
        reg_line = slope * log_filtered_x_TPM + intercept
    except Exception as e:
        print(f"Error calculating regression line: {e}")
        return None, None, None

    fig, ax = plt.subplots(figsize=(8, 6))
    # Removed label arguments for scatter and plot
    ax.scatter(log_filtered_x_TPM, log_filtered_y_CPM)
    ax.plot(log_filtered_x_TPM, reg_line, color="red", linestyle="-")

    ax.set_title(f"Correlation: {input_gene} (TPM) vs. {input_se} (CPM)")
    ax.set_xlabel(f"Log10 TPM of {input_gene}")
    ax.set_ylabel(f"Log10 CPM of {input_se}")

    # Removed legend call here
    ax.text(
        0.05, 0.95, f"Pearson: {pearson_corr:.2f}",
        transform=ax.transAxes, fontsize=12, verticalalignment="top"
    )

    return fig, ax, pearson_corr


def process_se_gene_correlation(
    se_data_list,
    gene_list,
    merge_featurecount_pkl,
    rna_csv,
    path_to_full_info_pickle,
    output_dir,
    use_all_samples=False,
    save_png=False
):
    """
    For each SE in se_data_list and each gene in the corresponding list of genes,
    generates a scatter plot (Pearson correlation), saves the figure (SVG and optionally PNG),
    and compiles the results into a TSV file.

    Parameters:
        se_data_list (list[str]): SE names (peak identifiers).
        gene_list (list[list[str]]): Each element is a list of gene names to evaluate.
        merge_featurecount_pkl (str): Path to a PKL file containing CPM data for SE.
        rna_csv (str): Path to a CSV file with TPM data for genes.
        path_to_full_info_pickle (str): Path to a pickle containing sample/peak info.
        output_dir (str): Directory to save correlation plots and results.
        use_all_samples (bool): If True, uses all overlapping samples. Otherwise uses SAMPLE_distinct.
        save_png (bool): If True, also save PNG files.

    Returns:
        A pandas DataFrame with columns [SE, Gene, Pearson].
    """
    results = []
    os.makedirs(output_dir, exist_ok=True)

    counter = 1
    for se_data, genes in zip(se_data_list, gene_list):
        for gene in genes:
            result = create_gene_se_correlation_chart_log10(
                se_data,
                gene,
                merge_featurecount_pkl,
                rna_csv,
                path_to_full_info_pickle,
                use_all_samples=use_all_samples
            )

            if result[0] is None:
                print(f"Skipping SE '{se_data}' and gene '{gene}' (no data).")
                continue

            fig, ax, pearson_corr = result
            file_index = str(counter).zfill(3)
            svg_filename = f"{file_index}_correlation_{se_data}_{gene}.svg"
            svg_path = os.path.join(output_dir, svg_filename)

            # Save SVG
            fig.savefig(svg_path, format="svg", dpi=300)

            # Optionally save PNG
            if save_png:
                png_filename = f"{file_index}_correlation_{se_data}_{gene}.png"
                png_path = os.path.join(output_dir, png_filename)
                fig.savefig(png_path, format="png", dpi=300)

            plt.close(fig)
            results.append((se_data, gene, pearson_corr))
            counter += 1

    columns = ["SE", "Gene", "Pearson"]
    results_df = pd.DataFrame(results, columns=columns)
    results_csv_path = os.path.join(output_dir, "correlation_results.csv")
    results_df.to_csv(results_csv_path, index=False, sep="\t")
    print(f"Correlation results saved to {results_csv_path}")

    return results_df


def main():
    parser = argparse.ArgumentParser(
        description="Command-line script for SE-Gene correlation analysis (Pearson only)."
    )
    parser.add_argument("--top_regions_df", type=str, required=True,
                        help="Path to a top_regions_df.pkl file with SE index and a 'gene_list' column.")
    parser.add_argument("--se_count", type=int, default=10,
                        help="Number of SE entries to process from top_regions_df (default: 10).")
    parser.add_argument("--output_dir", type=str, default="correlation_output",
                        help="Directory to save correlation plots and results.")
    parser.add_argument("--merge_featurecount_pkl", type=str,
                        default="output/merge_featurecount_CPM_zero_indexed.pkl",
                        help="Path to a PKL file with SE CPM data.")
    parser.add_argument("--rna_csv", type=str, required=True,
                        help="Path to a CSV file with gene TPM data.")
    parser.add_argument("--full_info_pkl", type=str, required=True,
                        help="Path to a pickle file with sample/peak info.")
    parser.add_argument("--use_all_samples", action="store_true",
                        help="If set, uses all samples in common between TPM/CPM data.")
    parser.add_argument("--save_png", action="store_true",
                        help="If set, also save PNG files in addition to SVG.")

    args = parser.parse_args()

    if not os.path.isfile(args.top_regions_df):
        print(f"[Error] top_regions_df file not found: {args.top_regions_df}")
        sys.exit(1)

    top_df = pd.read_pickle(args.top_regions_df)
    print(f"Loaded top_regions_df: {args.top_regions_df}, shape: {top_df.shape}")

    if "gene_list" not in top_df.columns:
        print("[Error] 'gene_list' column not found in top_regions_df.")
        sys.exit(1)

    se_data_list = top_df.index[:args.se_count].tolist()
    gene_list_raw = top_df["gene_list"].iloc[:args.se_count]

    processed_gene_lists = []
    for item in gene_list_raw:
        if isinstance(item, str):
            import ast
            try:
                parsed = ast.literal_eval(item)
                if not isinstance(parsed, list):
                    parsed = [str(parsed)]
                processed_gene_lists.append(parsed)
            except Exception:
                processed_gene_lists.append([item])
        elif isinstance(item, list):
            processed_gene_lists.append(item)
        else:
            processed_gene_lists.append([item])

    print(f"Using the top {args.se_count} SE entries. Each with a corresponding list of genes.")

    results_df = process_se_gene_correlation(
        se_data_list=se_data_list,
        gene_list=processed_gene_lists,
        merge_featurecount_pkl=args.merge_featurecount_pkl,
        rna_csv=args.rna_csv,
        path_to_full_info_pickle=args.full_info_pkl,
        output_dir=args.output_dir,
        use_all_samples=args.use_all_samples,
        save_png=args.save_png
    )

    print("Correlation analysis completed.")


if __name__ == "__main__":
    main()
