#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions for processing GTF and TPM files for genomic data analysis.
"""

import pandas as pd
import re
import logging
import os
from pathlib import Path

def setup_logging(log_level=logging.INFO, log_file=None):
    """
    Configure logging settings.
    
    Parameters
    ----------
    log_level : int
        Logging level (default: logging.INFO)
    log_file : str, optional
        Path to log file. If provided, logs will be saved to this file
    """
    # ログのフォーマット設定 (Set log format)
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # ルートロガーの設定をリセット (Reset root logger config)
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    
    # 基本設定 (Basic configuration)
    logging.basicConfig(
        level=log_level,
        format=log_format,
        datefmt=date_format
    )
    
    # ファイルへのログ出力が指定されている場合 (If log file is specified)
    if log_file:
        # ファイルハンドラーを作成 (Create file handler)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_formatter = logging.Formatter(log_format, datefmt=date_format)
        file_handler.setFormatter(file_formatter)
        
        # ルートロガーにファイルハンドラーを追加 (Add file handler to root logger)
        logging.getLogger('').addHandler(file_handler)
        
        logging.info(f"Logging to file: {log_file}")

def validate_files(gtf_path, tpm_path):
    """
    Validate that the input files exist.
    
    Parameters
    ----------
    gtf_path : str or Path
        Path to the GTF file
    tpm_path : str or Path
        Path to the TPM file
        
    Returns
    -------
    bool
        True if both files exist, False otherwise
    """
    gtf_file = Path(gtf_path)
    tpm_file = Path(tpm_path)
    
    if not gtf_file.exists():
        logging.error(f"GTF file {gtf_path} not found.")
        return False
    
    if not tpm_file.exists():
        logging.error(f"TPM file {tpm_path} not found.")
        return False
    
    return True

def ensure_dir_exists(file_path, create_dirs=False):
    """
    Ensure the directory for a file path exists.
    
    Parameters
    ----------
    file_path : str or Path
        Path to the file
    create_dirs : bool
        If True, create the directory if it doesn't exist
        
    Returns
    -------
    bool
        True if directory exists or was created, False otherwise
    """
    directory = os.path.dirname(file_path)
    if not directory:  # カレントディレクトリの場合
        return True
    
    if not os.path.exists(directory):
        if create_dirs:
            try:
                os.makedirs(directory)
                logging.info(f"Created directory: {directory}")
                return True
            except Exception as e:
                logging.error(f"Failed to create directory {directory}: {e}")
                return False
        else:
            logging.error(f"Directory {directory} does not exist. Use --create-dirs to create it automatically.")
            return False
    return True

def parse_gtf_attribute(attribute_string, attribute_name):
    """
    Extract a specific attribute value from a GTF attribute string.
    
    Parameters
    ----------
    attribute_string : str
        The attribute section of a GTF line (9th column)
    attribute_name : str
        The name of the attribute to extract
        
    Returns
    -------
    str or None
        The attribute value if found, None otherwise
    """
    pattern = rf'{attribute_name} "([^"]+)"'
    match = re.search(pattern, attribute_string)
    return match.group(1) if match else None

def parse_gtf_line(line):
    """
    Parse a line from a GTF file and extract relevant information.
    
    Parameters
    ----------
    line : str
        A line from a GTF file
        
    Returns
    -------
    dict or None
        Dictionary containing gene information or None if parsing failed
    str or None
        Error message if parsing failed, None otherwise
    """
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return None, "Insufficient number of fields"
    
    # フィーチャータイプをチェック (Check feature type)
    feature = fields[2]
    if feature != "exon":
        return None, "Feature not exon"
    
    # 属性を抽出 (Extract attributes)
    attributes = fields[8]
    gene_id = parse_gtf_attribute(attributes, 'gene_id')
    
    if not gene_id:
        return None, "No gene_id found"
    
    gene_name = parse_gtf_attribute(attributes, 'gene_name') or gene_id
    
    # 位置情報を抽出 (Extract position information)
    chrom = fields[0]
    start = int(fields[3]) - 1  # 0-based開始位置に変換 (Convert to 0-based start position)
    end = int(fields[4])
    strand = fields[6]
    
    # Check for undefined strand and skip
    if strand == '.':
        return None, "Undefined strand: Gene skipped"
    
    gene_info = {
        'gene_id': gene_id,
        'gene_name': gene_name,
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand
    }
    
    return gene_info, None

def update_gene_coordinates(existing_gene, new_gene_info):
    """
    Update gene coordinates based on new information.
    
    Parameters
    ----------
    existing_gene : dict
        The existing gene information
    new_gene_info : dict
        New gene information to incorporate
        
    Returns
    -------
    dict
        Updated gene information
    """
    updated_gene = existing_gene.copy()
    updated_gene['start'] = min(existing_gene['start'], new_gene_info['start'])
    updated_gene['end'] = max(existing_gene['end'], new_gene_info['end'])
    
    # strandが異なる場合の処理 (Handle case where strands differ)
    if existing_gene['strand'] != new_gene_info['strand']:
        logging.warning(
            f"Warning: Gene {existing_gene['gene_id']} has inconsistent strands ({existing_gene['strand']} vs {new_gene_info['strand']})"
        )
    
    return updated_gene

def check_id_name_mismatch(gene_id, gene_name):
    """
    Check if gene_id and gene_name are considered mismatched.
    
    Parameters
    ----------
    gene_id : str
        Gene ID
    gene_name : str
        Gene name
        
    Returns
    -------
    bool
        True if mismatched, False otherwise
    """
    # gene_idがgene_nameと完全一致する場合はOK (No mismatch if gene_id equals gene_name)
    if gene_id == gene_name:
        return False
    
    # gene_idからハイフン以降を除いた部分がgene_nameに一致する場合はOK 
    # (No mismatch if base_id equals gene_name)
    base_id = gene_id.split('-')[0]
    if base_id == gene_name:
        return False
        
    # それ以外は不一致 (Otherwise, considered mismatch)
    return True

def read_gtf_file(gtf_file, report_progress=True, progress_interval=1000000, verbose=False):
    """
    Read a GTF file and process it line by line.
    
    Parameters
    ----------
    gtf_file : str
        Path to the GTF file
    report_progress : bool
        Whether to report progress while processing large files
    progress_interval : int
        How often to report progress (in number of lines)
    verbose : bool
        If True, provide more detailed log messages
        
    Returns
    -------
    dict
        Dictionary of gene information
    list
        List of unmatched lines
    list
        List of gene_id/gene_name mismatches
    list
        List of lines with undefined strand
    """
    try:
        genes = {}
        unmatched_lines = []
        id_name_mismatches = []
        undefined_strand_lines = []  # New list to store lines with undefined strand
        
        with open(gtf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                # 進捗報告 (Progress reporting)
                if report_progress and line_num % progress_interval == 0:
                    logging.info(f"Processed {line_num} lines...")
                
                # 行をパース (Parse line)
                gene_info, error = parse_gtf_line(line)
                
                if error:
                    # Check if this is an undefined strand error
                    if error == "Undefined strand: Gene skipped":
                        # Extract gene_id for better reporting
                        fields = line.strip().split('\t')
                        attributes = fields[8] if len(fields) > 8 else ""
                        gene_id = parse_gtf_attribute(attributes, 'gene_id') or "unknown"
                        undefined_strand_lines.append((line_num, gene_id, line.strip()))
                        if verbose:
                            logging.warning(f"Line {line_num}: Gene {gene_id} skipped due to undefined strand")
                    else:
                        unmatched_lines.append((line_num, line.strip(), error))
                        if verbose:
                            logging.debug(f"Line {line_num}: {error} - {line.strip()[:50]}...")
                    continue
                
                gene_id = gene_info['gene_id']
                gene_name = gene_info['gene_name']
                
                # gene_idとgene_nameの不一致チェック (Check for gene_id/gene_name mismatch)
                if check_id_name_mismatch(gene_id, gene_name):
                    gene_tuple = (gene_id, gene_name)
                    if gene_tuple not in [item[2:4] for item in id_name_mismatches]:
                        id_name_mismatches.append((line_num, line.strip(), gene_id, gene_name))
                        if verbose:
                            logging.debug(f"ID-name mismatch: {gene_id} vs {gene_name}")
                
                # 遺伝子情報を更新または追加 (Update or add gene information)
                if gene_id in genes:
                    old_start = genes[gene_id]['start']
                    old_end = genes[gene_id]['end']
                    genes[gene_id] = update_gene_coordinates(genes[gene_id], gene_info)
                    if verbose and (old_start != genes[gene_id]['start'] or old_end != genes[gene_id]['end']):
                        logging.debug(f"Updated coords for {gene_id}: {old_start}-{old_end} -> {genes[gene_id]['start']}-{genes[gene_id]['end']}")
                else:
                    genes[gene_id] = gene_info
                    if verbose:
                        logging.debug(f"Added new gene: {gene_id} ({gene_name}) at {gene_info['chrom']}:{gene_info['start']}-{gene_info['end']}")
        
        return genes, unmatched_lines, id_name_mismatches, undefined_strand_lines
    
    except FileNotFoundError:
        logging.error(f"Error: File {gtf_file} not found.")
        return {}, [], [], []
    except Exception as e:
        logging.error(f"Error: Problem processing GTF file: {e}")
        return {}, [], [], []

def convert_genes_to_dataframe(genes):
    """
    Convert a dictionary of gene information to a DataFrame.
    
    Parameters
    ----------
    genes : dict
        Dictionary containing gene information
        
    Returns
    -------
    pandas.DataFrame
        DataFrame of gene information
    """
    gene_info_list = list(genes.values())
    return pd.DataFrame(gene_info_list)

def extract_gene_info_from_gtf(gtf_file, verbose=False):
    """
    Extract gene position information from a GTF file.
    
    This is the main entry point for GTF processing.
    
    Parameters
    ----------
    gtf_file : str
        Path to the GTF file
    verbose : bool
        If True, provide more detailed log messages
        
    Returns
    -------
    pandas.DataFrame
        DataFrame of gene information
    list
        List of unmatched lines
    list
        List of gene_id/gene_name mismatches
    list
        List of lines with undefined strand
    """
    # GTFファイルを読み込み、行ごとに処理 (Read GTF file and process line by line)
    genes, unmatched_lines, id_name_mismatches, undefined_strand_lines = read_gtf_file(gtf_file, verbose=verbose)
    
    # 遺伝子情報をDataFrameに変換 (Convert gene information to DataFrame)
    gene_df = convert_genes_to_dataframe(genes)
    
    return gene_df, unmatched_lines, id_name_mismatches, undefined_strand_lines

def read_tpm_file(tpm_file, verbose=False):
    """
    Read a TPM file.
    
    Parameters
    ----------
    tpm_file : str
        Path to the TPM file
    verbose : bool
        If True, provide more detailed log messages
        
    Returns
    -------
    pandas.DataFrame or None
        DataFrame containing TPM data or None if reading failed
    """
    try:
        df = pd.read_csv(tpm_file, sep='\t')
        # 必須カラムの確認 (Check required columns)
        required_columns = ['gene_id']
        for col in required_columns:
            if col not in df.columns:
                raise ValueError(f"TPM file is missing required column '{col}'")
        
        if verbose:
            # サンプル数とカラム情報を表示 (Display sample count and column information)
            sample_columns = [col for col in df.columns if col not in ['gene_id', 'gene_name']]
            logging.debug(f"TPM file has {len(sample_columns)} samples: {', '.join(sample_columns)}")
            logging.debug(f"TPM value range: {df[sample_columns].min().min():.4f} - {df[sample_columns].max().max():.4f}")
        
        return df
    except FileNotFoundError:
        logging.error(f"Error: File {tpm_file} not found.")
        return None
    except Exception as e:
        logging.error(f"Error: Problem reading TPM file: {e}")
        return None

def save_unmatched_lines(unmatched_lines, output_file):
    """
    Save unmatched lines to a file.
    
    Parameters
    ----------
    unmatched_lines : list
        List of unmatched lines
    output_file : str
        Path to the output file
    """
    try:
        with open(output_file, 'w') as f:
            f.write("Line number\tContent\tReason\n")
            for line_num, content, reason in unmatched_lines:
                f.write(f"{line_num}\t{content}\t{reason}\n")
        logging.info(f"Saved unmatched lines to {output_file}")
    except Exception as e:
        logging.error(f"Error: Problem saving unmatched lines: {e}")

def save_id_name_mismatches(id_name_mismatches, output_file):
    """
    Save gene_id/gene_name mismatches to a file.
    
    Parameters
    ----------
    id_name_mismatches : list
        List of gene_id/gene_name mismatches
    output_file : str
        Path to the output file
    """
    try:
        with open(output_file, 'w') as f:
            f.write("Line number\tgene_id\tgene_name\tContent\n")
            for line_num, content, gene_id, gene_name in id_name_mismatches:
                f.write(f"{line_num}\t{gene_id}\t{gene_name}\t{content}\n")
        logging.info(f"Saved gene_id/gene_name mismatches to {output_file}")
    except Exception as e:
        logging.error(f"Error: Problem saving gene_id/gene_name mismatches: {e}")

def save_undefined_strand_lines(undefined_strand_lines, output_file):
    """
    Save lines with undefined strand to a file.
    
    Parameters
    ----------
    undefined_strand_lines : list
        List of lines with undefined strand
    output_file : str
        Path to the output file
    """
    try:
        with open(output_file, 'w') as f:
            f.write("Line number\tGene ID\tOriginal GTF line\n")
            for line_num, gene_id, content in undefined_strand_lines:
                f.write(f"{line_num}\t{gene_id}\t{content}\n")
        logging.info(f"Saved {len(undefined_strand_lines)} lines with undefined strand to {output_file}")
    except Exception as e:
        logging.error(f"Error: Problem saving undefined strand lines: {e}")

def merge_gene_data(gene_df, tpm_df, join_method='inner', verbose=False):
    """
    Merge gene position information with TPM data.
    
    Parameters
    ----------
    gene_df : pandas.DataFrame
        DataFrame containing gene position information
    tpm_df : pandas.DataFrame
        DataFrame containing TPM data
    join_method : str
        Method for joining the data ('inner', 'left', 'outer')
    verbose : bool
        If True, provide more detailed log messages
        
    Returns
    -------
    pandas.DataFrame
        Merged DataFrame
    """
    # gene_idをキーにして結合 (Join using gene_id as key)
    logging.info(f"Merging data (join method: {join_method})...")
    merged_df = pd.merge(gene_df, tpm_df, on='gene_id', how=join_method)
    logging.info(f"Merged data has {len(merged_df)} rows")
    
    # サンプルカラムを特定 (Identify sample columns)
    sample_columns = [col for col in tpm_df.columns if col not in ['gene_id', 'gene_name']]
    
    if verbose:
        # 結合結果の詳細情報 (Detailed information about merge result)
        logging.debug(f"Number of unique chromosomes in result: {merged_df['chrom'].nunique()}")
        logging.debug(f"Number of unique strands in result: {merged_df['strand'].nunique()}")
        logging.debug(f"Chromosome distribution: {merged_df['chrom'].value_counts().head(5).to_dict()}")
    
    # 結果を指定形式で出力 (Format result as specified)
    # gene_idをsymbolとして使用 (Use gene_id as symbol)
    result_df = merged_df.rename(columns={'gene_id': 'symbol'})
    ordered_columns = ['symbol', 'chrom', 'strand', 'start', 'end'] + sample_columns
    
    # 指定されたカラムが存在するかチェック (Check if specified columns exist)
    existing_columns = [col for col in ordered_columns if col in result_df.columns]
    result_df = result_df[existing_columns]
    
    return result_df

def print_statistics(gene_df, tpm_df, merged_df, join_method='inner', verbose=False):
    """
    Print statistics about the merging operation.
    
    Parameters
    ----------
    gene_df : pandas.DataFrame
        DataFrame containing gene position information
    tpm_df : pandas.DataFrame
        DataFrame containing TPM data
    merged_df : pandas.DataFrame
        Merged DataFrame
    join_method : str
        Method used for joining the data
    verbose : bool
        If True, provide more detailed statistics
    """
    total_genes_in_gtf = len(gene_df)
    total_genes_in_tpm = len(tpm_df)
    
    if join_method == 'inner':
        matched_genes = len(merged_df)
    else:
        matched_genes = sum(merged_df['symbol'].isin(gene_df['gene_id']) & 
                        merged_df['symbol'].isin(tpm_df['gene_id']))
    
    print("\nStatistics:")
    print(f"Number of genes in GTF: {total_genes_in_gtf}")
    print(f"Number of genes in TPM file: {total_genes_in_tpm}")
    print(f"Number of matched genes: {matched_genes}")
    print(f"Number of genes in GTF but not in TPM: {total_genes_in_gtf - matched_genes}")
    print(f"Number of genes in TPM but not in GTF: {total_genes_in_tpm - matched_genes}")
    
    if verbose:
        # より詳細な統計情報 (More detailed statistics)
        print("\nDetailed Statistics:")
        
        # クロモソーム分布 (Chromosome distribution)
        if 'chrom' in gene_df.columns:
            print("\nChromosome distribution in GTF:")
            chrom_counts = gene_df['chrom'].value_counts().head(10)
            for chrom, count in chrom_counts.items():
                print(f"  {chrom}: {count} genes")
        
        # TPMの統計 (TPM statistics)
        if 'gene_id' in tpm_df.columns:
            sample_columns = [col for col in tpm_df.columns if col not in ['gene_id', 'gene_name']]
            if sample_columns:
                print("\nTPM statistics:")
                for sample in sample_columns[:5]:  # 最初の5サンプルのみ表示 (Show only first 5 samples)
                    print(f"  {sample}: mean={tpm_df[sample].mean():.2f}, median={tpm_df[sample].median():.2f}, max={tpm_df[sample].max():.2f}")
                if len(sample_columns) > 5:
                    print(f"  ... and {len(sample_columns) - 5} more samples")
        
        # 結合結果の詳細 (Merge result details)
        if 'start' in merged_df.columns and 'end' in merged_df.columns:
            gene_lengths = merged_df['end'] - merged_df['start']
            print("\nGene length statistics:")
            print(f"  Mean: {gene_lengths.mean():.1f} bp")
            print(f"  Median: {gene_lengths.median():.1f} bp")
            print(f"  Min: {gene_lengths.min():.1f} bp")
            print(f"  Max: {gene_lengths.max():.1f} bp")