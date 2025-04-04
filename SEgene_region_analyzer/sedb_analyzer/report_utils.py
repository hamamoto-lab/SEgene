"""
sedb_analyzer.report_utils - Report Generation Utility Module

This module provides functionality for generating text reports from SEdb analysis results.
"""

import pandas as pd
import numpy as np
import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

def dict_to_human_readable(data_dict, indent=0):
    """
    Convert dictionary data into a human-readable text format.
    
    Parameters
    ----------
    data_dict : dict
        The dictionary data to be converted.
    indent : int, optional
        The level of indentation (used for recursive calls).
        
    Returns
    -------
    str
        The formatted text string.
    """
    if not isinstance(data_dict, dict):
        return str(data_dict)
    
    lines = []
    indent_str = "  " * indent
    
    for key, value in data_dict.items():
        # 値の型によって処理を分ける
        if isinstance(value, dict):
            # 辞書の場合は再帰的に処理
            if not value:  # 空辞書の場合
                lines.append(f"{indent_str}{key}: {{}}")
            else:
                lines.append(f"{indent_str}{key}:")
                lines.append(dict_to_human_readable(value, indent + 1))
        elif isinstance(value, (list, tuple)):
            # リストやタプルの場合
            if not value:  # 空リストの場合
                lines.append(f"{indent_str}{key}: []")
            else:
                lines.append(f"{indent_str}{key}:")
                for item in value:
                    if isinstance(item, dict):
                        lines.append(f"{indent_str}  -")
                        lines.append(dict_to_human_readable(item, indent + 2))
                    else:
                        lines.append(f"{indent_str}  - {item}")
        elif isinstance(value, pd.Series):
            # Seriesの場合
            lines.append(f"{indent_str}{key}: Series with {len(value)} elements")
            # 先頭5要素を表示
            top_items = value.head(5)
            for idx, val in top_items.items():
                lines.append(f"{indent_str}  {idx}: {val}")
            if len(value) > 5:
                lines.append(f"{indent_str}  ... ({len(value) - 5} more items)")
        elif isinstance(value, pd.DataFrame):
            # DataFrameの場合
            lines.append(f"{indent_str}{key}: DataFrame with {value.shape[0]} rows and {value.shape[1]} columns")
            # 列名を表示
            lines.append(f"{indent_str}  Columns: {', '.join(value.columns)}")
            # 先頭3行をテキスト形式で表示
            if not value.empty:
                lines.append(f"{indent_str}  First {min(3, len(value))} rows:")
                for i, row in value.head(3).iterrows():
                    row_items = [f"{col}={val}" for col, val in row.items()]
                    lines.append(f"{indent_str}    Row {i}: {', '.join(row_items)}")
        else:
            # その他の値はそのまま表示
            lines.append(f"{indent_str}{key}: {value}")
    
    return "\n".join(lines)

def dataframe_to_human_readable(df, max_rows=20):
    """
    Convert a DataFrame into a human-readable text format.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to be converted.
    max_rows : int, optional
        Maximum number of rows to display.
        
    Returns
    -------
    str
        The formatted text string.
    """
    if df is None or df.empty:
        return "Empty DataFrame"
    
    # Basic information
    header = f"DataFrame with {df.shape[0]} rows and {df.shape[1]} columns\n"
    header += f"Columns: {', '.join(df.columns)}\n"
    
    # Convert DataFrame to text format
    with pd.option_context('display.max_rows', max_rows, 'display.max_columns', None):
        df_text = df.to_string()
    
    # If the DataFrame has many rows, display only the head
    if len(df) > max_rows:
        df_text = df.head(max_rows).to_string()
        df_text += f"\n... ({len(df) - max_rows} more rows)"
    
    return header + df_text

def save_figure(figure_object, output_base_path_without_ext, formats=['png', 'svg']):
    """
    Save a figure object in the specified formats.
    
    Parameters
    ----------
    figure_object : matplotlib.figure.Figure
        The figure object to be saved.
    output_base_path_without_ext : str or Path
        The base path for output files without extension.
    formats : list, optional
        List of formats to save the figure in (default is ['png', 'svg']).
        
    Returns
    -------
    list
        A list of file paths for the saved figures.
    """
    import matplotlib.pyplot as plt
    saved_paths = []
    
    for fmt in formats:
        try:
            output_path = f"{output_base_path_without_ext}.{fmt}"
            figure_object.savefig(output_path, bbox_inches='tight', dpi=300)
            saved_paths.append(output_path)
        except Exception as e:
            print(f"Failed to save figure ({output_base_path_without_ext}.{fmt}): {e}")
    
    # Release figure resources
    plt.close(figure_object)
    
    return saved_paths

def generate_analysis_summary(region_info, enrichment_results=None, 
                           sample_analysis=None, gene_analysis=None):
    """
    Generate a comprehensive summary text for the analysis results.
    
    Parameters
    ----------
    region_info : dict
        Region information.
    enrichment_results : pd.DataFrame, optional
        Enrichment analysis results.
    sample_analysis : dict, optional
        Sample distribution analysis results.
    gene_analysis : dict, optional
        Gene distribution analysis results.
        
    Returns
    -------
    str
        The generated summary text.
    """
    summary_lines = []
    
    # Header
    summary_lines.append("=" * 80)
    summary_lines.append("SEdb Analysis Results Summary".center(80))
    summary_lines.append("=" * 80)
    summary_lines.append(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append("-" * 80)
    
    # Region Information
    summary_lines.append("\n[Target Region for Analysis]")
    summary_lines.append("-" * 40)
    region_name = region_info.get('name', 'Unknown')
    chrom = region_info.get('chr', region_info.get('chrom', 'Unknown'))
    start = region_info.get('start', 0)
    end = region_info.get('end', 0)
    size = region_info.get('size', end - start)
    
    summary_lines.append(f"Region Name:     {region_name}")
    summary_lines.append(f"Chromosome:     {chrom}")
    summary_lines.append(f"Start Position:   {start:,}")
    summary_lines.append(f"End Position:   {end:,}")
    summary_lines.append(f"Size:     {size:,} bp")
    
    # Tissue Enrichment Analysis
    if enrichment_results is not None and not enrichment_results.empty:
        summary_lines.append("\n[Tissue Enrichment Analysis]")
        summary_lines.append("-" * 40)
        
        # Filter significant results
        if 'Significant' in enrichment_results.columns:
            significant = enrichment_results[enrichment_results['Significant'] == True]
            depleted = enrichment_results[enrichment_results['Depleted'] == True] if 'Depleted' in enrichment_results.columns else pd.DataFrame()
            
            summary_lines.append(f"Total tissue types analyzed:       {len(enrichment_results)}")
            summary_lines.append(f"Significantly enriched tissue types: {len(significant)}")
            if not depleted.empty:
                summary_lines.append(f"Significantly depleted tissue types:       {len(depleted)}")
            
            # Display top 5 entries
            if len(significant) > 0:
                summary_lines.append("\nTop Enriched Tissues:")
                summary_lines.append(f"{'Tissue Type':30} {'Fold Change':12} {'p-value':12} {'Adjusted p-value':16}")
                summary_lines.append("-" * 70)
                
                for i, (_, row) in enumerate(significant.sort_values('p-value').head(5).iterrows()):
                    tissue = row.get('Tissue Type', 'Unknown')
                    fc = row.get('Fold Change', 0)
                    pval = row.get('p-value', 1)
                    adj_pval = row.get('adjusted p-value', pval)
                    summary_lines.append(f"{tissue:30} {fc:12.2f} {pval:12.2e} {adj_pval:16.2e}")
        else:
            summary_lines.append(f"Number of tissue types analyzed: {len(enrichment_results)}")
    
    # Sample Distribution Analysis
    if sample_analysis:
        summary_lines.append("\n[Sample Distribution]")
        summary_lines.append("-" * 40)
        
        sample_count = sample_analysis.get('sample_count', 0)
        summary_lines.append(f"Total samples: {sample_count}")
        
        # Statistics of SE Counts
        se_stats = sample_analysis.get('se_count_stats', {})
        if se_stats:
            summary_lines.append("\nStatistics of SE Counts:")
            summary_lines.append(f"  Total:    {se_stats.get('total', 0):,}")
            summary_lines.append(f"  Mean:    {se_stats.get('mean', 0):.2f}")
            summary_lines.append(f"  Median:  {se_stats.get('median', 0):.2f}")
            summary_lines.append(f"  Standard Deviation: {se_stats.get('std', 0):.2f}")
            summary_lines.append(f"  Minimum:    {se_stats.get('min', 0)}")
            summary_lines.append(f"  Maximum:    {se_stats.get('max', 0)}")
        
        # Distribution of tissue types
        tissue_counts = sample_analysis.get('tissue_counts')
        if tissue_counts is not None and len(tissue_counts) > 0:
            summary_lines.append(f"\nTissue Type Distribution (Total: {len(tissue_counts)}):")
            summary_lines.append(f"{'Tissue Type':30} {'Sample Count':10} {'Percentage (%)':10}")
            summary_lines.append("-" * 50)
            
            # Display top 5 entries
            total_samples = tissue_counts.sum()
            for tissue, count in tissue_counts.nlargest(5).items():
                percent = (count / total_samples * 100) if total_samples > 0 else 0
                summary_lines.append(f"{tissue:30} {count:10} {percent:10.2f}")
            
            if len(tissue_counts) > 5:
                summary_lines.append(f"... and {len(tissue_counts) - 5} more types")
        
        # Distribution of biosample types
        biosample_counts = sample_analysis.get('biosample_counts')
        if biosample_counts is not None and len(biosample_counts) > 0:
            summary_lines.append(f"\nBiosample Type Distribution (Total: {len(biosample_counts)}):")
            summary_lines.append(f"{'Biosample Type':30} {'Sample Count':10} {'Percentage (%)':10}")
            summary_lines.append("-" * 50)
            
            # Display top 5 entries
            total_samples = biosample_counts.sum()
            for biosample, count in biosample_counts.nlargest(5).items():
                percent = (count / total_samples * 100) if total_samples > 0 else 0
                summary_lines.append(f"{biosample:30} {count:10} {percent:10.2f}")
            
            if len(biosample_counts) > 5:
                summary_lines.append(f"... and {len(biosample_counts) - 5} more types")
    
    # Gene Distribution Analysis
    if gene_analysis:
        summary_lines.append("\n[Gene Distribution]")
        summary_lines.append("-" * 40)
        
        gene_count = gene_analysis.get('gene_count', 0)
        summary_lines.append(f"Total genes detected: {gene_count}")
        
        # Statistics of Occurrence Counts
        occ_stats = gene_analysis.get('occurrence_stats', {})
        if occ_stats:
            summary_lines.append("\nStatistics of Occurrence Counts:")
            summary_lines.append(f"  Total:    {occ_stats.get('total', 0):,}")
            summary_lines.append(f"  Mean:    {occ_stats.get('mean', 0):.2f}")
            summary_lines.append(f"  Median:  {occ_stats.get('median', 0):.2f}")
            summary_lines.append(f"  Standard Deviation: {occ_stats.get('std', 0):.2f}")
            summary_lines.append(f"  Minimum:    {occ_stats.get('min', 0)}")
            summary_lines.append(f"  Maximum:    {occ_stats.get('max', 0)}")
        
        # Most Frequently Occurring Genes
        top_genes = gene_analysis.get('top_genes')
        if top_genes is not None and len(top_genes) > 0:
            summary_lines.append("\nMost Frequently Occurring Genes:")
            summary_lines.append(f"{'Gene Name':20} {'Occurrence Count':10}")
            summary_lines.append("-" * 30)
            
            for gene, count in top_genes.nlargest(10).items():
                summary_lines.append(f"{gene:20} {count:10}")
    
    # Footer
    summary_lines.append("\n" + "=" * 80)
    summary_lines.append("This analysis was generated automatically.")
    summary_lines.append("For detailed results, please refer to the individual files in the directory.")
    summary_lines.append("=" * 80)
    
    return "\n".join(summary_lines)
