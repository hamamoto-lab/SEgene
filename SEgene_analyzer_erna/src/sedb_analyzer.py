import pandas as pd

# Set non-interactive backend before ANY matplotlib import  
import matplotlib
matplotlib.use('Agg', force=True)  # Force non-interactive backend for WSL/headless environments
import matplotlib.pyplot as plt

# Disable matplotlib warnings for non-interactive backend
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
warnings.filterwarnings('ignore', message='.*non-interactive.*cannot be shown.*')
import seaborn as sns
import numpy as np
from collections import Counter
import os
import pybedtools
from io import StringIO
from scipy import stats
from statsmodels.stats.multitest import multipletests
import logging
import warnings
import time
from typing import Dict, List, Tuple, Any, Self, Optional
import japanize_matplotlib
import networkx as nx
from datetime import datetime


class SEdbRegionAnalyzer:
    """
    Super-enhancer (SE) genomic region analyzer for SEdb data.
    
    This class loads SEdb database files and provides methods to analyze
    specific genomic regions for overlapping super enhancers.
    """
    
    def __init__(self, results_dir='results', logger=None, 
                 se_bed_file=None, sample_info_file=None):
        """
        Initialize the analyzer
        
        Parameters:
        results_dir (str): Directory to save results (default: 'results')
        logger (Logger): Optional logger instance to use. If None, creates a basic console logger.
        se_bed_file (str, optional): Path to the SE_package BED file
        sample_info_file (str, optional): Path to the sample information file
        """
        # Logger setup - accept external logger or create a simple one
        self.logger = logger or self._setup_default_logger()
        
        self.logger.info("Initializing SEdbRegionAnalyzer instance")
        
        # Store parameters
        self.se_bed_file = se_bed_file
        self.sample_info_file = sample_info_file
        self.results_dir = results_dir
        
        # Create results directory
        os.makedirs(self.results_dir, exist_ok=True)
        self.logger.info(f"Results directory: {self.results_dir}")
        
        # Initialize empty attributes
        self.bed_df = None
        self._original_bed_df = None  
        self._filtered_bed_df = None 
        self.sample_info = None
        self.region = None
        self.overlapping_se = None
        self.final_data = None
        self.unique_se_data = None
        self.extraction_method = None

        # Region-specific data containers
        self.region_tissue_counts = None
        self.region_biosample_counts = None
        self.region_gene_counts = None
        
        # Load files automatically if specified (for backward compatibility)
        if se_bed_file and sample_info_file:
            self.load_databases()
        
        self.logger.info("SEgeneRegionAnalyzer initialization completed successfully")

    def _remove_na_sample_ids(self) -> None:
        """
        Remove rows from sample_info where 'Sample ID' is NaN.
        This method is intended for internal use (private).
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Nothing to remove.")
            return

        original_count = len(self.sample_info)
        # Use .copy() to avoid modifying the original DataFrame
        self.sample_info = self.sample_info.dropna(subset=['Sample ID']).copy()
        filtered_count = len(self.sample_info)
        removed_count = original_count - filtered_count

        if removed_count > 0:
            self.logger.info(f"Removed rows with NaN 'Sample ID' from sample_info:")
            self.logger.info(f"  Original row count: {original_count}")
            self.logger.info(f"  Removed row count: {removed_count}")
            self.logger.info(f"  Filtered row count: {filtered_count}")
        else:
            self.logger.info("No rows with NaN 'Sample ID' found in sample_info.")

    def _add_cell_id_to_sample_info(self) -> None:
        """
        Add 'cell_id' column to sample_info, derived from 'Sample ID', and make it the first column.
        Internal method.
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot add cell_id.")
            return

        try:
            # Type check output
            print(f"Type of self.sample_info['Sample ID']: {type(self.sample_info['Sample ID'])}")
            
            # Convert using simple string replacement (safer method)
            self.sample_info['cell_id'] = self.sample_info['Sample ID'].str.replace('Sample_', 'SE_', regex=False)
            
            # Verify that 'cell_id' column was created correctly
            if 'cell_id' not in self.sample_info.columns:
                raise ValueError("Failed to create 'cell_id' column.")

        except (AttributeError, ValueError) as e:
            self.logger.error(f"Error creating 'cell_id' column: {e}.  Check 'Sample ID' format.")
            return

        try:
            # Rearrange columns to put 'cell_id' first
            cols = ['cell_id'] + [col for col in self.sample_info.columns if col != 'cell_id']
            self.sample_info = self.sample_info[cols]
        except KeyError as e:
            self.logger.error(f"Error reordering columns: {e}.  'cell_id' might not be present.")
            return
        except Exception as e:
            self.logger.error(f"An unexpected error occurred during column reordering")
            return

        self.logger.info("'cell_id' column added to sample_info and moved to the first position.")

    def _setup_default_logger(self):
        """Set up and return a default logger when none is provided"""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        
        # Avoid duplicate handlers
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        
        return logger
    
    def set_se_bed_file(self, file_path: str) -> Self:
        """
        Set the SE BED file path
        
        Parameters:
        file_path (str): Path to the SE BED file
        
        Returns:
        SEgeneRegionAnalyzer: self for method chaining
        """
        self.logger.info(f"Setting SE BED file: {file_path}")
        self.se_bed_file = file_path
        # Reset associated data
        self.bed_df = None
        return self
    
    def set_sample_info_file(self, file_path: str) -> Self:
        """
        Set the sample information file path
        
        Parameters:
        file_path (str): Path to the sample information file
        
        Returns:
        SEgeneRegionAnalyzer: self for method chaining
        """
        self.logger.info(f"Setting sample information file: {file_path}")
        self.sample_info_file = file_path
        # Reset associated data
        self.sample_info = None
        return self

    def load_databases(self) -> Self:
        """
        Load and prepare the SEdb databases, including filtering and sample ID processing.

        Returns:
            SEgeneRegionAnalyzer: self for method chaining
        """
        if not self.se_bed_file or not self.sample_info_file:
            self.logger.warning("Database files not set. Use set_se_bed_file() and set_sample_info_file() first.")
            return self

        self.logger.info("Loading SEdb databases")

        try:
            # Load BED file
            start_time = time.time()
            self.bed_df = self._read_bed_file(self.se_bed_file)
            self.logger.info(f"BED file loaded with {len(self.bed_df)} rows - elapsed time: {time.time() - start_time:.2f}s")

            # Load sample information
            start_time = time.time()
            self.sample_info = self._read_sample_info(self.sample_info_file)
            self.logger.info(f"Sample information loaded with {len(self.sample_info)} rows - elapsed time: {time.time() - start_time:.2f}s")

            # Remove rows with NaN 'Sample ID' values *after* loading.
            self._remove_na_sample_ids()

            # Add and reorder 'cell_id' column in sample_info
            self._add_cell_id_to_sample_info()

            self.logger.info("Database loading completed successfully")
        except Exception as e:
            self.logger.exception(f"Error loading databases: {e}")
            raise

        return self

    def _read_bed_file(self, file_path):
        """
        Read BED file with original column names - internal method

        Parameters:
        file_path (str): Path to the BED file

        Returns:
        pandas.DataFrame: Loaded BED data with original column names
        """
        self.logger.info(f"Reading BED file: {file_path}")
        try:
            # Read file with header, preserving original column names
            se_df = pd.read_csv(file_path, sep='\s+', header=0)

            # Store original column names for reference
            self._original_bed_columns = list(se_df.columns)
            self.logger.debug(f"BED file original columns: {', '.join(self._original_bed_columns)}")

            # Log basic information about the data
            self.logger.info(f"BED file loaded with {len(se_df)} rows and {len(se_df.columns)} columns")
            self.logger.debug(f"First few sample IDs: {', '.join(se_df['cell_id'].head().tolist())}")

            # Log column data types
            column_types = {col: str(se_df[col].dtype) for col in se_df.columns}
            self.logger.debug(f"Column data types: {column_types}")

            return se_df
        except Exception as e:
            self.logger.exception(f"Error reading BED file: {e}")
            raise
            
    def _read_sample_info(self, file_path):
        """
        Read sample information file - internal method (with UTF-16 encoding support)
        
        Parameters:
        file_path (str): Path to the sample information file
        
        Returns:
        pandas.DataFrame: Loaded sample information
        """
        self.logger.info(f"Reading sample information file: {file_path}")
        try:
            # Try UTF-16 first
            try:
                df = pd.read_csv(file_path, sep='\t', encoding='utf-16')
                self.logger.debug("Successfully read sample info with UTF-16 encoding")
                return df
            except Exception as e:
                self.logger.debug(f"UTF-16 encoding failed: {str(e)}")
                
                # Try UTF-16-LE
                try:
                    df = pd.read_csv(file_path, sep='\t', encoding='utf-16-le')
                    self.logger.debug("Successfully read sample info with UTF-16-LE encoding")
                    return df
                except Exception as e:
                    self.logger.debug(f"UTF-16-LE encoding failed: {str(e)}")
                    
                    # Try UTF-8
                    try:
                        df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
                        self.logger.debug("Successfully read sample info with UTF-8 encoding")
                        return df
                    except Exception as e:
                        self.logger.debug(f"UTF-8 encoding failed: {str(e)}")
                        
                        # Try UTF-8 with BOM
                        try:
                            df = pd.read_csv(file_path, sep='\t', encoding='utf-8-sig')
                            self.logger.debug("Successfully read sample info with UTF-8-sig encoding")
                            return df
                        except Exception as e:
                            self.logger.error("All encodings failed")
                            
                            # Check file header bytes
                            with open(file_path, 'rb') as f:
                                header = f.read(4)
                                self.logger.error(f"File header bytes: {[hex(b) for b in header]}")
                            raise e
        except Exception as e:
            self.logger.exception(f"Error reading sample information file: {e}")
            raise

    def clear_data(self) -> Self:
        """
        Clear all loaded data
        
        Returns:
        SEgeneRegionAnalyzer: self for method chaining
        """
        self.bed_df = None
        self.sample_info = None
        self.region = None
        self.overlapping_se = None
        self.final_data = None
        self.unique_se_data = None
        self.region_tissue_counts = None
        self.region_biosample_counts = None
        self.region_gene_counts = None
        self.logger.info("All data cleared")
        return self
    
    def display_data_summary(self, include_head=True, include_stats=True, include_missing=True, head_rows=5):
        """
        Display summary information about loaded data

        Parameters:
        include_head (bool): Whether to display the first few rows of data
        include_stats (bool): Whether to display basic statistics for numerical columns
        include_missing (bool): Whether to display information about missing values
        head_rows (int): Number of rows to display in head preview
        """
        self.logger.info("Displaying data summary")

        if self.bed_df is None:
            print("BED data not loaded. Call load_databases first.")
            return

        # Basic information
        print(f"=== BED Data Summary ===")
        print(f"Dimensions: {self.bed_df.shape[0]} rows × {self.bed_df.shape[1]} columns")

        # Column information
        print("\nColumn Information:")
        col_info = []
        for col in self.bed_df.columns:
            dtype = self.bed_df[col].dtype
            non_null = self.bed_df[col].count()
            missing = self.bed_df[col].isna().sum()
            missing_pct = (missing / len(self.bed_df)) * 100
            unique = self.bed_df[col].nunique()

            col_info.append({
                'Column': col,
                'Type': dtype,
                'Non-Null': non_null,
                'Missing': missing,
                'Missing %': f"{missing_pct:.2f}%",
                'Unique Values': unique
            })

        # Create a DataFrame for better display
        col_df = pd.DataFrame(col_info)
        print(col_df.to_string(index=False))

        # Data preview
        if include_head and len(self.bed_df) > 0:
            print(f"\nData Preview (first {head_rows} rows):")
            print(self.bed_df.head(head_rows).to_string())

        # Statistics for numerical columns
        if include_stats:
            print("\nNumerical Column Statistics:")
            numeric_cols = self.bed_df.select_dtypes(include=['number']).columns
            if len(numeric_cols) > 0:
                print(self.bed_df[numeric_cols].describe().to_string())
            else:
                print("No numerical columns found.")

        # Missing value information
        if include_missing:
            missing_count = self.bed_df.isna().sum().sum()
            missing_percent = (missing_count / (self.bed_df.shape[0] * self.bed_df.shape[1])) * 100
            print(f"\nMissing Values: {missing_count} ({missing_percent:.2f}% of all data)")

        # Sample information if available
        if self.sample_info is not None:
            print(f"\n=== Sample Information Summary ===")
            print(f"Dimensions: {self.sample_info.shape[0]} rows × {self.sample_info.shape[1]} columns")
            print(f"Columns: {', '.join(self.sample_info.columns)}")

            if include_head and len(self.sample_info) > 0:
                print(f"\nSample Data Preview (first {head_rows} rows):")
                print(self.sample_info.head(head_rows).to_string())
        else:
            print("\nSample information not loaded.")

        self.logger.info("Data summary displayed successfully")
        
    def display_gene_columns(self, max_entries=10):
        """
        Display information about gene-related columns in the BED data

        Parameters:
        max_entries (int): Maximum number of gene entries to display per column
        """
        self.logger.info("Displaying gene column information")

        if self.bed_df is None:
            print("BED data not loaded. Call load_databases first.")
            return

        # Find gene-related columns based on original column names
        gene_columns = [col for col in self.bed_df.columns if 'gene' in col.lower()]

        if not gene_columns:
            print("No gene-related columns found in the data")
            return

        print(f"=== Gene Column Information ({len(gene_columns)} columns) ===")

        for col in gene_columns:
            print(f"\n=== {col} ===")

            # Count number of entries that have data
            non_empty = self.bed_df[col].count()
            percent_filled = (non_empty / len(self.bed_df)) * 100
            print(f"Non-empty entries: {non_empty}/{len(self.bed_df)} ({percent_filled:.2f}%)")

            # Count unique values
            unique_values = self.bed_df[col].nunique()
            print(f"Unique values: {unique_values}")

            # Extract all genes and count them
            all_genes = []
            for entry in self.bed_df[col].dropna():
                if isinstance(entry, str):
                    genes = entry.split(',')
                    all_genes.extend([g.strip() for g in genes if g.strip()])

            # Count unique genes
            gene_counter = Counter(all_genes)
            unique_genes = len(gene_counter)
            print(f"Unique genes: {unique_genes}")

            # Display most common genes
            if gene_counter:
                print(f"\nTop {min(max_entries, len(gene_counter))} most common genes:")
                for gene, count in gene_counter.most_common(max_entries):
                    print(f"  {gene}: {count}")

        self.logger.info("Gene column information displayed successfully")

    def filter_human_chromosomes(self, chromosomes: List[str] | None = None) -> Self:
        """
        Filter the BED data to include only specified chromosomes (default: chr1-22, chrX, chrY).
        Keeps the original data in _original_bed_df, and filtered data in bed_df.

        Parameters:
            chromosomes (list of str, optional): List of chromosomes to keep.
                Defaults to standard human chromosomes (chr1-chr22, chrX, chrY).

        Returns:
            SEgeneRegionAnalyzer: self, for method chaining
        """
        if self.bed_df is None:
            self.logger.warning("BED data is not loaded. Nothing to filter.")
            return self

        # Default to standard human chromosomes if not specified
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

        self.logger.info(f"Filtering BED data for chromosomes: {', '.join(chromosomes)}")

        # Store the original data if not already stored
        if self._original_bed_df is None:
            self._original_bed_df = self.bed_df.copy()

        # Filter data, keep original in self._original_bed_df
        filtered_df = self.bed_df.loc[self.bed_df['se_chr'].isin(chromosomes)].copy()
        removed_df = self.bed_df.loc[~self.bed_df['se_chr'].isin(chromosomes)].copy()

        original_size = len(self._original_bed_df)  # Use original data
        filtered_size = len(filtered_df)
        removed_size = len(removed_df)

        self.logger.info(f"BED data filtering completed:")
        self.logger.info(f"  Original size: {original_size} rows")
        self.logger.info(f"  Removed size: {removed_size} rows")
        self.logger.info(f"  Filtered size: {filtered_size} rows")

        # Store filtered data in bed_df
        self.bed_df = filtered_df

        # Store removed data internally
        if not removed_df.empty:
            self._removed_bed_df = removed_df
            self.logger.info("Removed BED data saved as attribute: _removed_bed_df")

        return self

    def count_cell_id_occurrences(self) -> pd.DataFrame:
        """
        Count the occurrences of 'cell_id' in bed_df based on the 'cell_id' column in sample_info.
        For each cell_id, calculate how many times it appears in bed_df. Store the results as an internal 
        attribute self._cell_id_counts. If there are cell_ids that don't exist in bed_df, a warning message 
        is displayed.
        
        Returns:
            pandas.DataFrame: DataFrame with 'cell_id' and its occurrence 'count'
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot count cell_id occurrences.")
            return None
        if self.bed_df is None:
            self.logger.warning("BED data is not loaded. Cannot count cell_id occurrences.")
            return None

        # Get occurrence counts of 'cell_id' in bed_df
        bed_counts = self.bed_df['cell_id'].value_counts()

        result = []
        # Count based on cell_id in sample_info
        for cell_id in self.sample_info['cell_id']:
            count = bed_counts.get(cell_id, 0)
            result.append({'cell_id': cell_id, 'count': count})
            if count == 0:
                message = f"Warning: cell_id '{cell_id}' not found in BED data."
                self.logger.warning(message)
                warnings.warn(message)

        # Save results as internal attribute
        self._cell_id_counts = pd.DataFrame(result)
        self.logger.info("cell_id occurrence counting completed.")
        return self._cell_id_counts

    def display_cell_id_summary(self, top_n=10, bottom_n=10, show_histogram=True, bins=30, figsize=(15, 10)):
        """
        Display summary information about cell_id occurrences.
        
        Parameters:
        -----------
        top_n : int, optional
            Number of most frequent cell_ids to display
        bottom_n : int, optional
            Number of least frequent cell_ids to display
        show_histogram : bool, optional
            Whether to display a histogram
        bins : int, optional
            Number of bins for the histogram
        figsize : tuple, optional
            Figure size
        """
        if self.cell_id_counts is None:
            self.logger.warning("cell_id counts not available. Running count_cell_id_occurrences()...")
            self.count_cell_id_occurrences()
            if self.cell_id_counts is None:
                return
        
        counts_df = self.cell_id_counts
        
        # Display basic statistics
        print("=== Basic Statistics of cell_id Occurrences ===")
        stats = counts_df['count'].describe()
        print(stats)
        
        # Number of cell_ids with zero occurrences
        zero_counts = (counts_df['count'] == 0).sum()
        if zero_counts > 0:
            print(f"\ncell_ids with zero occurrences: {zero_counts} ({zero_counts/len(counts_df)*100:.2f}%)")
        
        # Most frequent cell_ids
        print(f"\n=== Most Frequent cell_ids (Top {top_n}) ===")
        top_counts = counts_df.nlargest(top_n, 'count')
        print(top_counts)
        
        # Least frequent cell_ids (excluding zeros)
        non_zero_df = counts_df[counts_df['count'] > 0]
        if len(non_zero_df) > 0:
            print(f"\n=== Least Frequent cell_ids (Bottom {bottom_n}, excluding zeros) ===")
            bottom_counts = non_zero_df.nsmallest(bottom_n, 'count')
            print(bottom_counts)
        
        # Display histogram
        if show_histogram:
            plt.figure(figsize=figsize)
            
            # 1. Histogram
            plt.subplot(2, 2, 1)
            plt.hist(counts_df['count'], bins=bins, edgecolor='black')
            plt.title('Histogram of cell_id Occurrences')
            plt.xlabel('Occurrence Count')
            plt.ylabel('Number of cell_ids')
            
            # 2. Log-scale histogram
            plt.subplot(2, 2, 2)
            plt.hist(counts_df['count'], bins=bins, edgecolor='black', log=True)
            plt.title('Histogram of cell_id Occurrences (Log Scale)')
            plt.xlabel('Occurrence Count')
            plt.ylabel('Number of cell_ids (log)')
            
            # 3. Box plot
            plt.subplot(2, 2, 3)
            plt.boxplot(counts_df['count'])
            plt.title('Box Plot of cell_id Occurrences')
            plt.ylabel('Occurrence Count')
            
            # 4. Cumulative distribution
            plt.subplot(2, 2, 4)
            counts_sorted = np.sort(counts_df['count'])
            p = 1. * np.arange(len(counts_df)) / (len(counts_df) - 1)
            plt.plot(counts_sorted, p)
            plt.title('Cumulative Distribution of cell_id Occurrences')
            plt.xlabel('Occurrence Count')
            plt.ylabel('Cumulative Probability')
            
            plt.tight_layout()
            plt.show()

    def analyze_counts_by_metadata(self, group_by_columns=['Tissue'], figsize=(12, 8)):
        """
        Analyze cell_id occurrences based on sample metadata.
        
        Parameters:
        -----------
        group_by_columns : list of str, optional
            List of column names in sample_info to group by
        figsize : tuple, optional
            Figure size
        """
        if self.cell_id_counts is None:
            self.logger.warning("cell_id counts not available. Running count_cell_id_occurrences()...")
            self.count_cell_id_occurrences()
            if self.cell_id_counts is None:
                return
        
        # Merge cell_id_counts with sample_info
        merged_df = pd.merge(self.cell_id_counts, self.sample_info, on='cell_id')
        
        # Group and visualize for each group_by_column
        for column in group_by_columns:
            if column not in merged_df.columns:
                self.logger.warning(f"Column '{column}' not found in sample_info. Skipping.")
                continue
            
            # Calculate statistics by group
            group_stats = merged_df.groupby(column)['count'].agg(['mean', 'median', 'std', 'min', 'max', 'count'])
            group_stats = group_stats.sort_values('mean', ascending=False)
            
            print(f"\n=== cell_id Occurrences by '{column}' ===")
            print(group_stats)
            
            # Display box plots
            plt.figure(figsize=figsize)
            
            # 1. Box plot
            plt.subplot(1, 2, 1)
            sns.boxplot(x=column, y='count', data=merged_df)
            plt.title(f'Distribution of cell_id Occurrences by {column}')
            plt.xticks(rotation=90)
            
            # 2. Violin plot
            plt.subplot(1, 2, 2)
            sns.violinplot(x=column, y='count', data=merged_df)
            plt.title(f'Distribution of cell_id Occurrences by {column} (Violin Plot)')
            plt.xticks(rotation=90)
            
            plt.tight_layout()
            plt.show()

    def analyze_sample_data(self, figsize=(15, 10)):
        """
        Display basic statistics and distribution of sample information.
        
        Parameters:
        -----------
        figsize : tuple, optional
            Figure size
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Please run load_databases().")
            return
        
        sample_df = self.sample_info
        
        # Display basic information
        print(f"=== Sample Information Summary ===")
        print(f"Total samples: {len(sample_df)}")
        
        # Check missing values in each column
        print("\n=== Missing Values by Column ===")
        missing_counts = sample_df.isnull().sum()
        missing_percent = (missing_counts / len(sample_df) * 100).round(2)
        missing_df = pd.DataFrame({
            'Column': missing_counts.index,
            'Missing Count': missing_counts.values,
            'Missing Percentage (%)': missing_percent.values
        })
        print(missing_df)
        
        # Display distribution of categorical data
        categorical_cols = ['Species', 'Data source', 'Biosample type', 'Tissue type']
        for col in categorical_cols:
            if col in sample_df.columns:
                print(f"\n=== {col} Distribution ===")
                value_counts = sample_df[col].value_counts()
                percent = (value_counts / len(sample_df) * 100).round(2)
                dist_df = pd.DataFrame({
                    f'{col}': value_counts.index,
                    'Sample Count': value_counts.values,
                    'Percentage (%)': percent.values
                })
                print(dist_df)
        
        # Visualization
        plt.figure(figsize=figsize)
        
        # Biosample type distribution
        if 'Biosample type' in sample_df.columns:
            plt.subplot(2, 2, 1)
            biosample_counts = sample_df['Biosample type'].value_counts()
            biosample_counts.plot(kind='bar')
            plt.title('Distribution of Biosample Types')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
        
        # Tissue type distribution
        if 'Tissue type' in sample_df.columns:
            plt.subplot(2, 2, 2)
            tissue_counts = sample_df['Tissue type'].value_counts()
            tissue_counts.plot(kind='bar')
            plt.title('Distribution of Tissue Types')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
        
        # Species distribution (usually just Human, but included for completeness)
        if 'Species' in sample_df.columns:
            plt.subplot(2, 2, 3)
            species_counts = sample_df['Species'].value_counts()
            species_counts.plot(kind='pie', autopct='%1.1f%%')
            plt.title('Species Distribution')
            plt.axis('equal')
        
        # Data source distribution
        if 'Data source' in sample_df.columns:
            plt.subplot(2, 2, 4)
            source_counts = sample_df['Data source'].value_counts()
            source_counts.plot(kind='pie', autopct='%1.1f%%')
            plt.title('Data Source Distribution')
            plt.axis('equal')
        
        plt.tight_layout()
        plt.show()
        
        return sample_df

    def visualize_tissue_distribution(self, top_n=20, horizontal=True, store_all=True, figsize=(10, 6),
                                save_dir=None, save_filename=None, 
                                save_png=False, save_svg=False, save_pdf=False, save_eps=False,
                                save_dpi=300, save_transparent=False):
        """
        Calculate and visualize tissue type distribution.
        
        Parameters:
        -----------
        top_n : int, optional
            表示する上位組織タイプの数（デフォルト: 20）
        horizontal : bool, optional
            水平バーチャートを使用するか（デフォルト: True）
        store_all : bool, optional
            全ての分布データを保存するか（デフォルト: True）
        figsize : tuple, optional
            図のサイズ（デフォルト: (10, 6)）
        save_dir : str, optional
            図を保存するディレクトリ（指定しない場合は保存しない）
        save_filename : str, optional
            保存するファイル名（拡張子なし）
        save_png, save_svg, save_pdf, save_eps : bool, optional
            各フォーマットで保存するかどうか（デフォルト: False）
        save_dpi : int, optional
            保存時の解像度（デフォルト: 300）
        save_transparent : bool, optional
            透明な背景で保存するかどうか（デフォルト: False）
            
        Returns:
        --------
        pandas.DataFrame
            組織タイプの分布データ
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot calculate tissue distribution.")
            return None
        
        print("組織型分布を計算...")
        
        # 組織タイプのカラムを確認
        tissue_col = 'Tissue type' if 'Tissue type' in self.sample_info.columns else None
        
        if tissue_col is None:
            self.logger.warning("Tissue type column not found in sample information.")
            return None
        
        # 組織タイプの分布を計算
        tissue_counts = self.sample_info[tissue_col].value_counts()
        total_samples = len(self.sample_info)
        
        # 計算結果をDataFrameに変換
        tissue_distribution = pd.DataFrame({
            'Tissue Type': tissue_counts.index,
            'Sample Count': tissue_counts.values,
            'Percentage (%)': (tissue_counts / total_samples * 100).round(2)
        })
        
        # 全データを保存（オプション）
        if store_all:
            self._tissue_distribution = tissue_distribution.copy()
        
        # トップNを表示用に選択
        if top_n > 0 and len(tissue_distribution) > top_n:
            display_data = tissue_distribution.head(top_n)
        else:
            display_data = tissue_distribution
        
        # 可視化用のフィギュアを作成
        plt.figure(figsize=figsize)
        
        if horizontal:
            # 水平バーチャート（長い名前のリストに適している）
            plt.barh(display_data['Tissue Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='skyblue')
            plt.xlabel('Sample Count')
            plt.ylabel('Tissue Type')
        else:
            # 垂直バーチャート
            plt.bar(display_data['Tissue Type'], display_data['Sample Count'], color='skyblue')
            plt.xlabel('Tissue Type')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
        
        plt.title(f'Tissue Type Distribution (Top {len(display_data)} of {len(tissue_distribution)})')
        plt.tight_layout()
        
        # 図を保存（オプション）
        if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
            import os
            from datetime import datetime
            
            # ディレクトリが存在しなければ作成
            os.makedirs(save_dir, exist_ok=True)
            
            # タイムスタンプを生成（ファイル名がない場合）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # ファイル名がなければデフォルト名を生成
            if save_filename is None:
                save_filename = f'tissue_distribution_{timestamp}'
                
            # 指定された各フォーマットで保存
            saved_files = []
            
            if save_png:
                png_path = os.path.join(save_dir, f'{save_filename}.png')
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_svg:
                svg_path = os.path.join(save_dir, f'{save_filename}.svg')
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_pdf:
                pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = os.path.join(save_dir, f'{save_filename}.eps')
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # 保存されたファイルをログに記録
            if saved_files:
                self.logger.info(f"Tissue distribution visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return tissue_distribution

    def visualize_tissue_distribution(self, top_n=20, horizontal=True, store_all=True, figsize=(10, 6),
                                save_dir=None, save_filename=None, 
                                save_png=False, save_svg=False, save_pdf=False, save_eps=False,
                                save_dpi=300, save_transparent=False):
        """
        Calculate and visualize tissue type distribution.
        
        Parameters:
        -----------
        top_n : int, optional
            表示する上位組織タイプの数（デフォルト: 20）
        horizontal : bool, optional
            水平バーチャートを使用するか（デフォルト: True）
        store_all : bool, optional
            全ての分布データを保存するか（デフォルト: True）
        figsize : tuple, optional
            図のサイズ（デフォルト: (10, 6)）
        save_dir : str, optional
            図を保存するディレクトリ（指定しない場合は保存しない）
        save_filename : str, optional
            保存するファイル名（拡張子なし）
        save_png, save_svg, save_pdf, save_eps : bool, optional
            各フォーマットで保存するかどうか（デフォルト: False）
        save_dpi : int, optional
            保存時の解像度（デフォルト: 300）
        save_transparent : bool, optional
            透明な背景で保存するかどうか（デフォルト: False）
            
        Returns:
        --------
        pandas.DataFrame
            組織タイプの分布データ
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot calculate tissue distribution.")
            return None
        
        print("組織型分布を計算...")
        
        # 組織タイプのカラムを確認
        tissue_col = 'Tissue type' if 'Tissue type' in self.sample_info.columns else None
        
        if tissue_col is None:
            self.logger.warning("Tissue type column not found in sample information.")
            return None
        
        # 組織タイプの分布を計算
        tissue_counts = self.sample_info[tissue_col].value_counts()
        total_samples = len(self.sample_info)
        
        # 計算結果をDataFrameに変換
        tissue_distribution = pd.DataFrame({
            'Tissue Type': tissue_counts.index,
            'Sample Count': tissue_counts.values,
            'Percentage (%)': (tissue_counts / total_samples * 100).round(2)
        })
        
        # 全データを保存（オプション）
        if store_all:
            self._tissue_distribution = tissue_distribution.copy()
        
        # トップNを表示用に選択
        if top_n > 0 and len(tissue_distribution) > top_n:
            display_data = tissue_distribution.head(top_n)
        else:
            display_data = tissue_distribution
        
        # 可視化用のフィギュアを作成
        plt.figure(figsize=figsize)
        
        if horizontal:
            # 水平バーチャート（長い名前のリストに適している）
            plt.barh(display_data['Tissue Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='skyblue')
            plt.xlabel('Sample Count')
            plt.ylabel('Tissue Type')
        else:
            # 垂直バーチャート
            plt.bar(display_data['Tissue Type'], display_data['Sample Count'], color='skyblue')
            plt.xlabel('Tissue Type')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
        
        plt.title(f'Tissue Type Distribution (Top {len(display_data)} of {len(tissue_distribution)})')
        plt.tight_layout()
        
        # 図を保存（オプション）
        if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
            import os
            from datetime import datetime
            
            # ディレクトリが存在しなければ作成
            os.makedirs(save_dir, exist_ok=True)
            
            # タイムスタンプを生成（ファイル名がない場合）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # ファイル名がなければデフォルト名を生成
            if save_filename is None:
                save_filename = f'tissue_distribution_{timestamp}'
                
            # 指定された各フォーマットで保存
            saved_files = []
            
            if save_png:
                png_path = os.path.join(save_dir, f'{save_filename}.png')
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_svg:
                svg_path = os.path.join(save_dir, f'{save_filename}.svg')
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_pdf:
                pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = os.path.join(save_dir, f'{save_filename}.eps')
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # 保存されたファイルをログに記録
            if saved_files:
                self.logger.info(f"Tissue distribution visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return tissue_distribution




    def analyze_tissue_biosample_relationship(self, figsize=(15, 12), min_count=5):
        """
        Analyze the relationship between tissue types and biosample types.
        
        Parameters:
        -----------
        figsize : tuple, optional
            Figure size
        min_count : int, optional
            Minimum count to include in the heatmap
        
        Returns:
        --------
        pandas.DataFrame
            Cross-tabulation of tissue types and biosample types
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Please run load_databases().")
            return None
        
        if 'Tissue type' not in self.sample_info.columns or 'Biosample type' not in self.sample_info.columns:
            self.logger.warning("'Tissue type' or 'Biosample type' column not found in sample information.")
            return None
        
        # Create cross-tabulation
        cross_tab = pd.crosstab(self.sample_info['Tissue type'], self.sample_info['Biosample type'])
        
        # Add row and column totals
        cross_tab['Total'] = cross_tab.sum(axis=1)
        cross_tab.loc['Total'] = cross_tab.sum()
        
        print(f"=== Tissue Type and Biosample Type Cross-tabulation ===")
        print(cross_tab)
        
        # Prepare data for heatmap (exclude total row and column)
        heatmap_data = cross_tab.iloc[:-1, :-1].copy()
        
        # Sort rows and columns by total count (descending)
        row_order = heatmap_data.sum(axis=1).sort_values(ascending=False).index
        col_order = heatmap_data.sum(axis=0).sort_values(ascending=False).index
        
        heatmap_data = heatmap_data.loc[row_order, col_order]
        
        # Filter out tissue types and biosample types with few samples
        if min_count > 0:
            row_mask = heatmap_data.sum(axis=1) >= min_count
            col_mask = heatmap_data.sum(axis=0) >= min_count
            filtered_data = heatmap_data.loc[row_mask, col_mask]
            
            if filtered_data.empty:
                print(f"No tissue type and biosample type combinations meet the min_count={min_count} criteria.")
                filtered_data = heatmap_data
        else:
            filtered_data = heatmap_data
        
        # Visualization
        plt.figure(figsize=figsize)
        
        # Color map (0 values displayed as gray)
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        mask = filtered_data == 0
        
        # Heatmap
        sns.heatmap(filtered_data, annot=True, fmt='d', cmap=cmap, mask=mask,
                linewidths=0.5, cbar_kws={'label': 'Sample Count'})
        
        plt.title('Relationship Between Tissue Type and Biosample Type')
        plt.xlabel('Biosample Type')
        plt.ylabel('Tissue Type')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.show()
        
        # Network graph data preparation
        plt.figure(figsize=figsize)
        
        # Create edge list
        edges = []
        for tissue in filtered_data.index:
            for biosample in filtered_data.columns:
                count = filtered_data.loc[tissue, biosample]
                if count > 0:
                    edges.append((tissue, biosample, count))
        
        # Create graph
        G = nx.Graph()
        
        # Add nodes
        tissue_nodes = list(filtered_data.index)
        biosample_nodes = list(filtered_data.columns)
        
        for tissue in tissue_nodes:
            G.add_node(tissue, bipartite=0)  # Tissue type nodes
        
        for biosample in biosample_nodes:
            G.add_node(biosample, bipartite=1)  # Biosample type nodes
        
        # Add edges
        for source, target, weight in edges:
            G.add_edge(source, target, weight=weight)
        
        # Calculate node positions (bipartite layout)
        pos = nx.spring_layout(G, k=0.3, iterations=50)
        
        # Edge width proportional to weight
        edge_width = [G[u][v]['weight'] / 5 for u, v in G.edges()]
        
        # Draw edges
        nx.draw_networkx_edges(G, pos, width=edge_width, alpha=0.7, edge_color='gray')
        
        # Draw tissue type nodes
        nx.draw_networkx_nodes(G, pos, nodelist=tissue_nodes, node_color='skyblue', 
                            node_size=300, alpha=0.8, label='Tissue Type')
        
        # Draw biosample type nodes
        nx.draw_networkx_nodes(G, pos, nodelist=biosample_nodes, node_color='lightgreen', 
                            node_size=300, alpha=0.8, label='Biosample Type')
        
        # Draw labels
        nx.draw_networkx_labels(G, pos, font_size=8)
        
        plt.title('Tissue Type and Biosample Type Network Diagram')
        plt.legend()
        plt.axis('off')
        plt.tight_layout()
        plt.show()
        
        return cross_tab

    def analyze_data_sources(self, figsize=(12, 8)):
        """
        Analyze data by source.
        
        Parameters:
        -----------
        figsize : tuple, optional
            Figure size
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Please run load_databases().")
            return
        
        if 'Data source' not in self.sample_info.columns:
            self.logger.warning("'Data source' column not found in sample information.")
            return
        
        # Calculate sample count by data source
        source_counts = self.sample_info['Data source'].value_counts()
        
        # Cross-tabulation of tissue type and data source
        tissue_source = None
        if 'Tissue type' in self.sample_info.columns:
            tissue_source = pd.crosstab(self.sample_info['Tissue type'], self.sample_info['Data source'])
            tissue_source['Total'] = tissue_source.sum(axis=1)
            tissue_source = tissue_source.sort_values('Total', ascending=False)
        
        # Cross-tabulation of biosample type and data source
        biosample_source = None
        if 'Biosample type' in self.sample_info.columns:
            biosample_source = pd.crosstab(self.sample_info['Biosample type'], self.sample_info['Data source'])
            biosample_source['Total'] = biosample_source.sum(axis=1)
            biosample_source = biosample_source.sort_values('Total', ascending=False)
        
        # Visualization
        plt.figure(figsize=figsize)
        
        # Sample count by data source
        plt.subplot(2, 2, 1)
        bars = plt.bar(source_counts.index, source_counts.values, 
                    color=plt.cm.Set3(np.arange(len(source_counts))))
        plt.title('Sample Count by Data Source')
        plt.ylabel('Sample Count')
        plt.xticks(rotation=45, ha='right')
        
        # Display values above bars
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height}', ha='center', va='bottom')
        
        # Percentage by data source (pie chart)
        plt.subplot(2, 2, 2)
        plt.pie(source_counts.values, labels=source_counts.index, autopct='%1.1f%%', 
            colors=plt.cm.Set3(np.arange(len(source_counts))))
        plt.title('Percentage by Data Source')
        plt.axis('equal')
        
        # Stacked bar chart (data source breakdown by tissue type)
        if tissue_source is not None:
            plt.subplot(2, 2, 3)
            # Show only top 10 tissue types
            top_tissues = tissue_source.iloc[:10, :-1]  # Exclude 'Total' column
            top_tissues.plot(kind='bar', stacked=True, ax=plt.gca())
            plt.title('Data Source Breakdown by Major Tissue Types')
            plt.ylabel('Sample Count')
            plt.xlabel('Tissue Type')
            plt.xticks(rotation=45, ha='right')
            plt.legend(loc='upper right')
        
        # Stacked bar chart (data source breakdown by biosample type)
        if biosample_source is not None:
            plt.subplot(2, 2, 4)
            # Show only top 10 biosample types
            top_biosamples = biosample_source.iloc[:10, :-1]  # Exclude 'Total' column
            top_biosamples.plot(kind='bar', stacked=True, ax=plt.gca())
            plt.title('Data Source Breakdown by Major Biosample Types')
            plt.ylabel('Sample Count')
            plt.xlabel('Biosample Type')
            plt.xticks(rotation=45, ha='right')
            plt.legend(loc='upper right')
        
        plt.tight_layout()
        plt.show()
        
        print("\n=== Sample Count by Data Source ===")
        source_df = pd.DataFrame({
            'Data Source': source_counts.index,
            'Sample Count': source_counts.values,
            'Percentage (%)': (source_counts / source_counts.sum() * 100).round(2).values
        })
        print(source_df)
        
        if tissue_source is not None:
            print("\n=== Cross-tabulation of Tissue Type and Data Source ===")
            print(tissue_source)
        
        if biosample_source is not None:
            print("\n=== Cross-tabulation of Biosample Type and Data Source ===")
            print(biosample_source)

    def test_tissue_enrichment(self, region_tissue_counts, method='fisher', correction='fdr_bh'):
        """
        Test for tissue type enrichment in a specific region.
        
        Parameters:
        -----------
        region_tissue_counts : dict or pandas.Series
            Tissue type counts in the specific region
        method : str, optional
            Statistical test method ('fisher' or 'chi2')
        correction : str, optional
            Multiple testing correction method. Same as the method parameter in statsmodels.stats.multitest.multipletests.
            
        Returns:
        --------
        pandas.DataFrame
            Test results for each tissue type (p-values, corrected p-values, odds ratios, etc.)
        """
        if not hasattr(self, '_tissue_distribution'):
            self.logger.warning("Overall tissue type distribution data not available. Please run visualize_tissue_distribution() first.")
            return None
        
        # Get overall sample distribution
        background = self._tissue_distribution.set_index('Tissue Type')['Sample Count']
        total_background = background.sum()
        
        # Convert region counts to Series (if it's a dictionary)
        if isinstance(region_tissue_counts, dict):
            region_counts = pd.Series(region_tissue_counts)
        else:
            region_counts = region_tissue_counts
        
        total_region = region_counts.sum()
        
        # List to store results
        results = []
        
        # Perform test for each tissue type
        for tissue in background.index:
            if tissue in region_counts.index:
                a = region_counts[tissue]  # Count of tissue in region
            else:
                a = 0
            
            b = total_region - a  # Count of other tissues in region
            c = background[tissue]  # Count of tissue in overall data
            d = total_background - c  # Count of other tissues in overall data
            
            # 2x2 contingency table
            contingency_table = np.array([[a, b], [c, d]])
            
            if method == 'fisher':
                # Fisher's exact test
                odds_ratio, p_value = stats.fisher_exact(contingency_table)
            elif method == 'chi2':
                # Chi-squared test
                chi2, p_value, dof, expected = stats.chi2_contingency(contingency_table)
                odds_ratio = (a * d) / (b * c) if (b * c != 0) else np.nan
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            # Observed and expected values
            observed = a
            expected = (c / total_background) * total_region
            
            # Fold change
            fold_change = (a / total_region) / (c / total_background) if (c != 0) else np.nan
            
            results.append({
                'Tissue Type': tissue,
                'Region Count': a,
                'Region Percentage': (a / total_region * 100) if total_region > 0 else 0,
                'Total Count': c,
                'Total Percentage': (c / total_background * 100),
                'Expected Value': expected,
                'Odds Ratio': odds_ratio,
                'Fold Change': fold_change,
                'p-value': p_value
            })
        
        # Convert results to DataFrame
        result_df = pd.DataFrame(results)
        
        # Sort by p-value
        result_df = result_df.sort_values('p-value')
        
        # Multiple testing correction
        if len(result_df) > 1:
            reject, p_corr, _, _ = multipletests(result_df['p-value'].values, method=correction)
            result_df['adjusted p-value'] = p_corr
            result_df['Significant'] = reject
        
        return result_df



    # def visualize_tissue_enrichment(self, enrichment_results=None, title=None, region_name=None,
    #                             plot_type='bar_fdr_horizontal', top_n=15,
    #                             pvalue_threshold=0.05, fc_threshold=1.5,
    #                             save_dir=None, save_filename=None, 
    #                             save_png=False, save_svg=False, save_pdf=False, save_eps=False,
    #                             save_dpi=300, save_transparent=False):
    #     """
    #     Visualize tissue enrichment results.
        
    #     Parameters:
    #     -----------
    #     enrichment_results : pandas.DataFrame, optional
    #         エンリッチメント分析の結果
    #     title : str, optional
    #         図のタイトル
    #     region_name : str, optional
    #         領域名（タイトルに含める）
    #     plot_type : str, optional
    #         プロットタイプ
    #     top_n : int, optional
    #         表示する上位組織タイプの数
    #     pvalue_threshold, fc_threshold : float, optional
    #         有意性とフォールドチェンジの閾値
    #     save_dir : str, optional
    #         図を保存するディレクトリ（指定しない場合は保存しない）
    #     save_filename : str, optional
    #         保存するファイル名（拡張子なし）
    #     save_png, save_svg, save_pdf, save_eps : bool, optional
    #         各フォーマットで保存するかどうか（デフォルト: False）
    #     save_dpi : int, optional
    #         保存時の解像度（デフォルト: 300）
    #     save_transparent : bool, optional
    #         透明な背景で保存するかどうか（デフォルト: False）
            
    #     Returns:
    #     --------
    #     matplotlib.figure.Figure
    #         生成された図
    #     """
    #     # 既存の実装 (省略) ...
    #     # 既存のコードに保存機能を追加する
        
    #     # 図を保存する部分の実装例
    #     if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
    #         import os
    #         from datetime import datetime
            
    #         # ディレクトリが存在しなければ作成
    #         os.makedirs(save_dir, exist_ok=True)
            
    #         # タイムスタンプを生成（ファイル名がない場合）
    #         timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
    #         # ファイル名がなければデフォルト名を生成
    #         if save_filename is None:
    #             plot_type_cleaned = plot_type.replace('_', '-')
    #             save_filename = f'tissue_enrichment_{plot_type_cleaned}_{timestamp}'
                
    #         # 指定された各フォーマットで保存
    #         saved_files = []
            
    #         if save_png:
    #             png_path = os.path.join(save_dir, f'{save_filename}.png')
    #             plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
    #             saved_files.append(f"PNG: {png_path}")
            
    #         if save_svg:
    #             svg_path = os.path.join(save_dir, f'{save_filename}.svg')
    #             plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
    #             saved_files.append(f"SVG: {svg_path}")
            
    #         if save_pdf:
    #             pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
    #             plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
    #             saved_files.append(f"PDF: {pdf_path}")
            
    #         if save_eps:
    #             eps_path = os.path.join(save_dir, f'{save_filename}.eps')
    #             plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
    #             saved_files.append(f"EPS: {eps_path}")
            
    #         # 保存されたファイルをログに記録
    #         if saved_files:
    #             self.logger.info(f"Tissue enrichment visualization saved in the following formats:")
    #             for file_info in saved_files:
    #                 self.logger.info(f"  - {file_info}")
        
    #     # 図のオブジェクトを返す（既存の実装と同様）
    #     return fig


#f3ced4bより復元
    def visualize_tissue_enrichment(self, enrichment_results=None, plot_type='volcano', 
                                title=None, figsize=(14, 8), top_n=15, 
                                pvalue_threshold=0.05, fc_threshold=1.5,
                                region_name=None, interactive=False,
                                enriched_color='#e66101', depleted_color='#5e3c99', 
                                nonsig_color='#b2b2b2',
                                save_dir=None, save_filename=None, 
                                save_png=False, save_svg=False, save_pdf=False, save_eps=False,
                                save_dpi=300, save_transparent=False):
        """
        Visualize tissue enrichment results based on unique samples.
        
        Parameters:
        -----------
        enrichment_results : pandas.DataFrame, optional
            Results from analyze_region_tissue_enrichment(). If None, uses stored results.
        plot_type : str, optional
            Plot type to use for visualization:
            
            FDR-based plots:
            - 'bar_fdr_horizontal' (was 'bar_fdr2'): Horizontal bar chart with -log10(FDR) as width
            - 'bar_fdr_vertical' (was 'bar_fdr'): Vertical bar chart with -log10(FDR) as height
            
            Fold change-based plots:
            - 'bar_fc_horizontal' (was 'bar'): Horizontal bar chart with fold change as width
            - 'bar_fc_vertical' (was 'bar2'): Vertical bar chart with fold change as height
            
            Other visualization types:
            - 'volcano': Scatter plot with fold change vs -log10(p-value)
            - 'bubble': Bubble plot where size represents sample count
            - 'heatmap': Heatmap visualization of multiple metrics
            - 'summary': Multi-panel figure with several plot types
        
        title : str, optional
            Plot title. Auto-generated if None.
        figsize : tuple, optional
            Figure size
        top_n : int, optional
            Number of top tissues to display
        pvalue_threshold : float, optional
            P-value threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        region_name : str, optional
            Region name for the plot title
        interactive : bool, optional
            Whether to use interactive plots (not currently supported)
        enriched_color : str, optional
            Color for enriched tissues (default: orange)
        depleted_color : str, optional
            Color for depleted tissues (default: purple)
        nonsig_color : str, optional
            Color for non-significant tissues (default: light gray)
        save_dir : str, optional
            Directory to save output files. If None, figure is not saved.
        save_filename : str, optional
            Base filename for saved files (without extension). If None, auto-generated.
        save_png, save_svg, save_pdf, save_eps : bool, optional
            Whether to save in respective formats (default: False)
        save_dpi : int, optional
            DPI for raster formats like PNG (default: 300)
        save_transparent : bool, optional
            Whether to save with transparent background (default: False)
        
        Returns:
        --------
        matplotlib.figure.Figure: Plot figure
        """
        # Handle old plot type names for backward compatibility
        plot_type_map = {
            'bar': 'bar_fc_horizontal',
            'bar2': 'bar_fc_vertical',
            'bar_fdr': 'bar_fdr_vertical',
            'bar_fdr2': 'bar_fdr_horizontal'
        }
        if plot_type in plot_type_map:
            original_type = plot_type
            plot_type = plot_type_map[plot_type]
            self.logger.info(f"Plot type '{original_type}' renamed to '{plot_type}' for clarity")
        
        # Check for results
        if enrichment_results is None:
            if hasattr(self, 'region_tissue_enrichment'):
                enrichment_results = self.region_tissue_enrichment
            else:
                self.logger.warning("No enrichment results provided. Run analyze_region_tissue_enrichment() first.")
                return None
        
        # Make a copy and handle missing values
        results_df = enrichment_results.copy()
        
        # Check for English column names and use them
        p_col = 'adjusted p-value' if 'adjusted p-value' in results_df.columns else 'p-value' if 'p-value' in results_df.columns else 'p値'
        fc_col = 'Fold Change' if 'Fold Change' in results_df.columns else 'フォールドチェンジ'
        tissue_col = 'Tissue Type' if 'Tissue Type' in results_df.columns else '組織型'
        count_col = 'Region Count' if 'Region Count' in results_df.columns else 'リージョンカウント'
        
        # Set plot title
        if title is None:
            title_prefix = "Tissue Enrichment Analysis" 
            if region_name:
                title_prefix += f" - {region_name}"
            title = f"{title_prefix} ({len(results_df)} tissues, based on unique samples)"
        
        # Calculate depleted threshold in a readable format (2 decimal places)
        depleted_threshold = round(1/fc_threshold, 2)
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        if plot_type == 'bar_fdr_horizontal':
            # Horizontal bar chart with -log10(FDR) as width
            # Get FDR column name
            if 'adjusted p-value' in results_df.columns:
                fdr_col = 'adjusted p-value'
            elif 'FDR' in results_df.columns:
                fdr_col = 'FDR'
            elif '補正p値' in results_df.columns:
                fdr_col = '補正p値'
            else:
                fdr_col = p_col
                self.logger.warning(f"No FDR column found, using p-value column '{p_col}' instead")
            
            # Calculate -log10(FDR) for visualization
            results_df['-log10(FDR)'] = -np.log10(results_df[fdr_col])
            
            # Sort by -log10(FDR) and get top tissues
            sorted_df = results_df.sort_values('-log10(FDR)', ascending=False)
            if top_n and len(sorted_df) > top_n:
                plot_data = sorted_df.head(top_n)
            else:
                plot_data = sorted_df
                
            # Reverse order for horizontal bar chart to show highest significance at top
            plot_data = plot_data.iloc[::-1].reset_index(drop=True)  
            
            # Create horizontal bar plot, colored by fold change
            bars = plt.barh(
                y=plot_data[tissue_col],
                width=plot_data['-log10(FDR)'],
                color=[
                    enriched_color if row[fc_col] >= fc_threshold else
                    depleted_color if row[fc_col] <= 1/fc_threshold else
                    nonsig_color
                    for _, row in plot_data.iterrows()
                ],
                alpha=0.7
            )
            
            # Add threshold line for significant FDR (e.g., 0.05)
            significance_line = -np.log10(pvalue_threshold)
            plt.axvline(x=significance_line, linestyle='--', color='black', alpha=0.7)
            
            # Make sure we have reasonable margins
            plt.subplots_adjust(right=0.85, left=0.15)
            
            # Add FDR, sample count, and fold change labels next to bars
            for i, bar in enumerate(bars):
                row = plot_data.iloc[i]
                fc = row[fc_col]
                count = row[count_col]
                fdr_value = row[fdr_col]
                bar_width = bar.get_width()
                y_pos = bar.get_y() + bar.get_height() * 0.5
                
                # Determine text position - just after the bar, 
                # but if the bar crosses the significance line, move text beyond the line
                if bar_width >= significance_line:
                    # Bar crosses the significance line, place text after the line
                    x_pos = max(bar_width + 0.1, significance_line + 0.3)
                else:
                    # Bar is shorter than significance line, place text right after the bar
                    x_pos = bar_width + 0.1
                
                # Create a compact string with all info (to reduce clutter)
                label_text = f"FDR={fdr_value:.1e}, n={int(count)}, FC={fc:.2f}"
                
                plt.text(
                    x=x_pos,
                    y=y_pos,
                    s=label_text,
                    va='center',
                    ha='left',
                    fontsize=9,
                    color='dimgray'
                )
            
            # Format plot
            plt.xlabel('-log10(FDR)')
            plt.ylabel('Tissue Type')
            plt.title(f'{title} - Statistical Significance (Horizontal Bar Plot)')
            
            # Add legend for fold change interpretation with significance threshold
            # Position in the lower right corner within the plot area
            from matplotlib.patches import Patch
            from matplotlib.lines import Line2D
            
            # Use simple ASCII symbols for >= and <= to avoid font issues
            legend_elements = [
                Patch(facecolor=enriched_color, alpha=0.7, label=f'Enriched (FC >= {fc_threshold})'),
                Patch(facecolor=depleted_color, alpha=0.7, label=f'Depleted (FC <= {depleted_threshold})'),
                Patch(facecolor=nonsig_color, alpha=0.7, label='Not enriched/depleted'),
                Line2D([0], [0], color='black', linestyle='--', 
                    label=f'FDR = {pvalue_threshold} threshold')
            ]
            plt.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.98, 0.02))
            
            # Adjust x-axis to ensure all bars and labels are visible
            max_width = max(plot_data['-log10(FDR)'])
            plt.xlim(0, max_width * 1.3)  # Add 30% extra space for labels
        
        elif plot_type == 'bar_fdr_vertical':
            # Bar plot with -log10(FDR) as height - vertical bar chart
            # Get FDR column name
            if 'adjusted p-value' in results_df.columns:
                fdr_col = 'adjusted p-value'
            elif 'FDR' in results_df.columns:
                fdr_col = 'FDR'
            elif '補正p値' in results_df.columns:
                fdr_col = '補正p値'
            else:
                fdr_col = p_col
                self.logger.warning(f"No FDR column found, using p-value column '{p_col}' instead")
            
            # Calculate -log10(FDR) for visualization
            results_df['-log10(FDR)'] = -np.log10(results_df[fdr_col])
            
            # Sort by -log10(FDR) and get top tissues
            sorted_df = results_df.sort_values('-log10(FDR)', ascending=False)
            if top_n and len(sorted_df) > top_n:
                plot_data = sorted_df.head(top_n)
            else:
                plot_data = sorted_df
                
            # Create bar plot, colored by fold change
            bars = plt.bar(
                plot_data[tissue_col],
                plot_data['-log10(FDR)'],
                color=[
                    enriched_color if row[fc_col] >= fc_threshold else
                    depleted_color if row[fc_col] <= 1/fc_threshold else
                    nonsig_color
                    for _, row in plot_data.iterrows()
                ],
                alpha=0.7
            )
            
            # Add threshold line for significant FDR (e.g., 0.05)
            significance_line = -np.log10(pvalue_threshold)
            plt.axhline(y=significance_line, linestyle='--', color='black', alpha=0.7, 
                    label=f'FDR={pvalue_threshold}')
            
            # Add fold change and count labels to bars
            for i, bar in enumerate(bars):
                row = plot_data.iloc[i]
                fc = row[fc_col]
                count = row[count_col]
                
                # Display fold change value inside or above bars
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.3,
                    f"FC={fc:.2f}",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )
                
                # Display count above the FC
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.7,
                    f"n={int(count)}",
                    ha='center',
                    va='bottom',
                    fontsize=8
                )
            
            # Format plot
            plt.xticks(rotation=45, ha='right')
            plt.xlabel('Tissue Type')
            plt.ylabel('-log10(FDR)')
            plt.title(f'{title} - Statistical Significance Bar Plot')
            
            # Add legend for fold change interpretation
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor=enriched_color, alpha=0.7, label=f'Enriched (FC >= {fc_threshold})'),
                Patch(facecolor=depleted_color, alpha=0.7, label=f'Depleted (FC <= {depleted_threshold})'),
                Patch(facecolor=nonsig_color, alpha=0.7, label='Not enriched/depleted')
            ]
            plt.legend(handles=legend_elements, loc='best')
            
            # Adjust y-axis to ensure all bars and labels are visible
            plt.ylim(0, max(plot_data['-log10(FDR)']) * 1.2)
        
        elif plot_type == 'bar_fc_horizontal':
            # Bar chart - ranking by fold change - horizontal layout
            # Extract top N tissue types
            if 'adjusted p-value' in results_df.columns:
                p_column = 'adjusted p-value'
            elif '補正p値' in results_df.columns:
                p_column = '補正p値'
            else:
                p_column = p_col
            
            # Sort by significance and fold change
            results_df['Significance'] = -np.log10(results_df[p_column])
            sorted_df = results_df.sort_values(by=['Significance', fc_col], ascending=[False, False])
            
            # Extract top N
            plot_df = sorted_df.head(top_n)
            
            # Create horizontal bar plot
            bars = plt.barh(
                y=plot_df[tissue_col],
                width=plot_df[fc_col],
                color=[
                    enriched_color if (row[fc_col] >= fc_threshold and row[p_column] <= pvalue_threshold) else
                    depleted_color if (row[fc_col] <= 1/fc_threshold and row[p_column] <= pvalue_threshold) else
                    nonsig_color 
                    for _, row in plot_df.iterrows()
                ],
                alpha=0.7
            )
            
            # Add reference line at fold change = 1
            plt.axvline(x=1, linestyle='--', color='gray')
            
            # Display p-values inside bars
            for i, bar in enumerate(bars):
                plt.text(
                    x=bar.get_width() * 0.5,
                    y=bar.get_y() + bar.get_height() * 0.5,
                    s=f"p={plot_df.iloc[i][p_column]:.1e}",
                    va='center',
                    ha='center',
                    color='white' if bar.get_width() > 2 else 'black',
                    fontweight='bold',
                    fontsize=8
                )
            
            # Format plot
            plt.xlabel('Fold Change')
            plt.ylabel('Tissue Type')
            plt.title(f'{title} - Tissue Enrichment (Top {len(plot_df)} tissues)')
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor=enriched_color, alpha=0.7, label=f'Enriched (p <= {pvalue_threshold}, FC >= {fc_threshold})'),
                Patch(facecolor=depleted_color, alpha=0.7, label=f'Depleted (p <= {pvalue_threshold}, FC <= {depleted_threshold})'),
                Patch(facecolor=nonsig_color, alpha=0.7, label='Not significant')
            ]
            plt.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.98, 0.02))
        
        elif plot_type == 'bar_fc_vertical':
            # Bar plot for tissue enrichment - vertical format based on fold change
            # Filter significant tissues if too many
            if 'adjusted p-value' in results_df.columns:
                p_column = 'adjusted p-value'
            elif '補正p値' in results_df.columns:
                p_column = '補正p値'
            else:
                p_column = p_col
                
            if top_n is not None and len(results_df) > top_n:
                # Get significant tissues first
                sig_tissues = results_df[(results_df[p_column] <= pvalue_threshold) & 
                                        (results_df[fc_col] >= fc_threshold)]
                
                # If we have enough significant tissues, use them
                if len(sig_tissues) >= top_n:
                    # Sort by fold change and get top ones
                    plot_data = sig_tissues.sort_values(fc_col, ascending=False).head(top_n)
                else:
                    # Otherwise, include some non-significant but high fold change tissues
                    remaining = top_n - len(sig_tissues)
                    non_sig = results_df[(results_df[p_column] > pvalue_threshold) | 
                                        (results_df[fc_col] < fc_threshold)]
                    non_sig = non_sig.sort_values(fc_col, ascending=False).head(remaining)
                    plot_data = pd.concat([sig_tissues, non_sig]).sort_values(fc_col, ascending=False)
            else:
                # Use all data, sorted by fold change
                plot_data = results_df.sort_values(fc_col, ascending=False)
            
            # Create bar plot, colored by significance
            bars = plt.bar(
                plot_data[tissue_col],
                plot_data[fc_col],
                color=[
                    enriched_color if (row[fc_col] >= fc_threshold and row[p_column] <= pvalue_threshold) else
                    nonsig_color 
                    for _, row in plot_data.iterrows()
                ],
                alpha=0.7
            )
            
            # Add reference line at fold change = 1
            plt.axhline(y=1, linestyle='--', color='black', alpha=0.5)
            
            # Add threshold line
            plt.axhline(y=fc_threshold, linestyle=':', color=enriched_color, alpha=0.5)
            
            # Add count and p-value annotations to bars
            for i, bar in enumerate(bars):
                row = plot_data.iloc[i]
                count = row[count_col]
                p_val = row[p_column]
                
                # Display counts inside or above bars
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    max(0.5, min(row[fc_col] - 0.5, row[fc_col] * 0.8)) if row[fc_col] > 2 else 0.2,
                    str(int(count)),
                    ha='center',
                    va='bottom',
                    color='white' if row[fc_col] > 2 else 'black',
                    fontweight='bold',
                    fontsize=9
                )
                
                # Display p-values above bars
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    row[fc_col] + 0.2,
                    f"p={p_val:.1e}",
                    ha='center',
                    va='bottom',
                    color='red' if p_val <= pvalue_threshold else 'gray',
                    fontsize=8,
                    rotation=90
                )
            
            # Format plot
            plt.xticks(rotation=45, ha='right')
            plt.xlabel('Tissue Type')
            plt.ylabel('Fold Change')
            plt.title(f'{title} - Fold Change Bar Plot')
            
            # Add legend
            from matplotlib.patches import Patch
            legend_elements = [
                Patch(facecolor=enriched_color, alpha=0.7, label=f'Significant (p <= {pvalue_threshold}, FC >= {fc_threshold})'),
                Patch(facecolor=nonsig_color, alpha=0.7, label='Not significant')
            ]
            plt.legend(handles=legend_elements, loc='best')
            
            # Adjust y-axis to include annotations
            y_max = max(results_df[fc_col]) * 1.2
            plt.ylim(0, y_max)
        
        elif plot_type == 'volcano':
            # Volcano plot
            if 'adjusted p-value' in results_df.columns:
                p_column = 'adjusted p-value'
            elif '補正p値' in results_df.columns:
                p_column = '補正p値'
            else:
                p_column = p_col
                
            plt.scatter(
                results_df[fc_col], 
                -np.log10(results_df[p_column]),
                s=results_df[count_col] * 10,  # Size proportional to count
                alpha=0.6,
                c=[
                    enriched_color if (row[fc_col] >= fc_threshold and row[p_column] <= pvalue_threshold) else
                    depleted_color if (row[fc_col] <= 1/fc_threshold and row[p_column] <= pvalue_threshold) else
                    nonsig_color 
                    for _, row in results_df.iterrows()
                ]
            )
            
            # Add labels for significant tissues
            for _, row in results_df.iterrows():
                if (row[p_column] <= pvalue_threshold and 
                    (row[fc_col] >= fc_threshold or row[fc_col] <= 1/fc_threshold)):
                    plt.annotate(
                        row[tissue_col],
                        xy=(row[fc_col], -np.log10(row[p_column])),
                        xytext=(5, 0),
                        textcoords='offset points',
                        ha='left',
                        va='center',
                        fontsize=8,
                        bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3)
                    )
            
            # Add threshold lines
            plt.axhline(-np.log10(pvalue_threshold), linestyle='--', color='gray', alpha=0.7)
            plt.axvline(fc_threshold, linestyle='--', color=enriched_color, alpha=0.7)
            plt.axvline(1/fc_threshold, linestyle='--', color=depleted_color, alpha=0.7)
            
            # Format plot
            plt.xlabel('Fold Change')
            plt.ylabel('-log10(p-value)')
            plt.title(f'{title} - Volcano Plot')
            plt.grid(True, linestyle='--', alpha=0.5)
            
            # Adjust x-axis if fold changes are extreme
            fc_max = results_df[fc_col].max()
            if fc_max > 10:
                plt.xlim(0, min(fc_max * 1.1, 20))
                
        elif plot_type == 'bubble':
            # Bubble plot implementation (omitted for brevity)
            pass
            
        elif plot_type == 'summary':
            # Summary plot implementation (omitted for brevity)
            pass
            
        elif plot_type == 'heatmap':
            # Heatmap implementation (omitted for brevity)
            pass
            
        else:
            self.logger.warning(f"Unknown plot type: {plot_type}")
            plt.close(fig)
            return None
        
        # plt.tight_layout()





        # plt.show()
        
        # # Store the last results
        # self._last_enrichment_results = enrichment_results
        
        # return fig





        plt.tight_layout()
        
        # Save figure if requested
        if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
            import os
            from datetime import datetime
            
            # Create directory if needed
            os.makedirs(save_dir, exist_ok=True)
            
            # Generate timestamp (for filename if not provided)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Filename should be descriptive if not provided
            if save_filename is None:
                plot_type_cleaned = plot_type.replace('_', '-')
                save_filename = f'tissue_enrichment_{plot_type_cleaned}_{timestamp}'
                
            # Save in each requested format
            saved_files = []
            
            if save_png:
                png_path = os.path.join(save_dir, f'{save_filename}.png')
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_svg:
                svg_path = os.path.join(save_dir, f'{save_filename}.svg')
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_pdf:
                pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = os.path.join(save_dir, f'{save_filename}.eps')
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Tissue enrichment visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        plt.show()
        
        # Store the last results
        self._last_enrichment_results = enrichment_results
        
        return fig



    def compare_multiple_regions(self, region_results_dict, plot_type='heatmap', 
                                figsize=(14, 10), p_threshold=0.05):
        """
        Compare tissue enrichment results across multiple regions.
        
        Parameters:
        -----------
        region_results_dict : dict
            Dictionary with region names as keys and results from test_tissue_enrichment() as values
        plot_type : str, optional
            Plot type ('heatmap', 'barplot', 'network')
        figsize : tuple, optional
            Figure size
        p_threshold : float, optional
            P-value threshold for significance
        
        Returns:
        --------
        matplotlib.figure.Figure
            Created figure object
        """
        if not region_results_dict or len(region_results_dict) == 0:
            self.logger.warning("Region results dictionary is empty.")
            return None
        
        # List of region names
        region_names = list(region_results_dict.keys())
        
        # Find all unique tissue types
        all_tissues = set()
        for results in region_results_dict.values():
            all_tissues.update(results['Tissue Type'].values)
        
        # Collect significant tissue types for each region
        significant_tissues = {}
        for region, results in region_results_dict.items():
            p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
            sig_tissues = results[results[p_col] <= p_threshold]['Tissue Type'].values
            significant_tissues[region] = set(sig_tissues)
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        if plot_type == 'heatmap':
            # Create dataframe for heatmap
            heatmap_data = []
            
            for tissue in sorted(all_tissues):
                row = {'Tissue Type': tissue}
                
                for region in region_names:
                    results = region_results_dict[region]
                    tissue_row = results[results['Tissue Type'] == tissue]
                    
                    if len(tissue_row) > 0:
                        fc = tissue_row['Fold Change'].values[0]
                        p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                        p_val = tissue_row[p_col].values[0]
                        
                        # Set fold change to 1 (no change) for non-significant results
                        if p_val > p_threshold:
                            row[region] = 1.0
                        else:
                            row[region] = fc
                    else:
                        row[region] = 1.0  # No data means fold change = 1 (no change)
                
                heatmap_data.append(row)
            
            # Convert to DataFrame
            heatmap_df = pd.DataFrame(heatmap_data).set_index('Tissue Type')
            
            # Filter by significance (significant in at least one region)
            sig_in_any = set()
            for sig_set in significant_tissues.values():
                sig_in_any.update(sig_set)
            
            if len(sig_in_any) > 0:
                filtered_df = heatmap_df.loc[list(sig_in_any)]
            else:
                filtered_df = heatmap_df
            
            # Sort by maximum fold change
            max_fc = filtered_df.max(axis=1)
            filtered_df = filtered_df.iloc[max_fc.argsort(ascending=False)]
            
            # Create heatmap
            cmap = sns.diverging_palette(220, 10, as_cmap=True)
            
            # Create mask for values = 1 (to display in gray)
            mask = filtered_df == 1.0
            
            # Draw heatmap
            sns.heatmap(filtered_df, cmap=cmap, center=1, annot=True, fmt='.2f', linewidths=0.5,
                    mask=mask, cbar_kws={'label': 'Fold Change'})
            
            plt.title(f'Tissue Enrichment Comparison Across Multiple Regions (p < {p_threshold})')
            
        elif plot_type == 'barplot':
            # Collect all significant tissue types
            all_sig_tissues = set()
            for sig_set in significant_tissues.values():
                all_sig_tissues.update(sig_set)
            
            # Limit to top 15 tissues
            if len(all_sig_tissues) > 15:
                # Get maximum fold change for each tissue across all regions
                max_fc_by_tissue = {}
                for tissue in all_sig_tissues:
                    max_fc = 0
                    for region, results in region_results_dict.items():
                        tissue_row = results[results['Tissue Type'] == tissue]
                        if len(tissue_row) > 0:
                            fc = abs(tissue_row['Fold Change'].values[0] - 1.0)  # Absolute change from 1.0
                            max_fc = max(max_fc, fc)
                    max_fc_by_tissue[tissue] = max_fc
                
                # Sort by maximum fold change and get top 15 tissues
                sorted_tissues = sorted(max_fc_by_tissue.items(), key=lambda x: x[1], reverse=True)
                all_sig_tissues = {t for t, _ in sorted_tissues[:15]}
            
            # Plot fold changes for significant tissue types in each region
            n_regions = len(region_names)
            region_colors = plt.cm.tab10(np.arange(n_regions) % 10)
            
            for i, (region, color) in enumerate(zip(region_names, region_colors)):
                results = region_results_dict[region]
                p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                
                # Filter to get significant tissue types only
                sig_results = results[
                    (results[p_col] <= p_threshold) & 
                    (results['Tissue Type'].isin(all_sig_tissues))
                ]
                
                if len(sig_results) > 0:
                    plt.subplot(n_regions, 1, i+1)
                    bars = plt.barh(
                        y=sig_results['Tissue Type'],
                        width=sig_results['Fold Change'],
                        color=color,
                        alpha=0.7
                    )
                    
                    plt.axvline(x=1, linestyle='--', color='gray')
                    plt.xlim(0, max(5, sig_results['Fold Change'].max() * 1.1))
                    plt.title(f'{region} - Significant Tissue Enrichment')
                    
                    if i == n_regions - 1:
                        plt.xlabel('Fold Change')
                    
                    # Display p-values
                    for bar, p_val in zip(bars, sig_results[p_col]):
                        plt.text(
                            x=max(bar.get_width() + 0.1, 1.1),
                            y=bar.get_y() + bar.get_height() * 0.5,
                            s=f"p={p_val:.1e}",
                            va='center',
                            ha='left',
                            fontsize=8
                        )
            
        elif plot_type == 'network':
            # Prepare nodes and edges for network plot
            import networkx as nx
            
            G = nx.Graph()
            
            # Add region nodes
            for region in region_names:
                G.add_node(region, type='region')
            
            # Add significant tissue type nodes
            all_sig_tissues = set()
            for sig_set in significant_tissues.values():
                all_sig_tissues.update(sig_set)
            
            # Add tissue nodes
            for tissue in all_sig_tissues:
                G.add_node(tissue, type='tissue')
            
            # Add edges (between regions and their significant tissue types)
            for region, results in region_results_dict.items():
                p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                
                for _, row in results[(results[p_col] <= p_threshold) & 
                                    (results['Tissue Type'].isin(all_sig_tissues))].iterrows():
                    tissue = row['Tissue Type']
                    fc = row['Fold Change']
                    G.add_edge(region, tissue, weight=fc)
            
            # Calculate node positions
            pos = nx.spring_layout(G, k=0.3, seed=42)
            
            # Separate region nodes and tissue nodes
            region_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'region']
            tissue_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'tissue']
            
            # Edge colors based on fold change
            edge_colors = []
            for u, v in G.edges():
                fc = G[u][v]['weight']
                if fc > 1:
                    edge_colors.append('red')
                else:
                    edge_colors.append('blue')
            
            # Edge width proportional to fold change
            edge_widths = [abs(G[u][v]['weight'] - 1) * 2 + 0.5 for u, v in G.edges()]
            
            # Draw network
            plt.figure(figsize=figsize)
            
            # Draw region nodes
            nx.draw_networkx_nodes(G, pos, nodelist=region_nodes, node_size=800,
                                node_color='lightgreen', alpha=0.8, 
                                node_shape='s', label='Region')
            
            # Draw tissue nodes
            nx.draw_networkx_nodes(G, pos, nodelist=tissue_nodes, node_size=600,
                                node_color='lightblue', alpha=0.8,
                                node_shape='o', label='Tissue Type')
            
            # Draw edges
            nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.7,
                                edge_color=edge_colors)
            
            # Draw labels
            nx.draw_networkx_labels(G, pos, font_size=8)
            
            plt.title('Network of Regions and Significant Tissue Types')
            plt.legend()
            plt.axis('off')
        
        else:
            self.logger.warning(f"Unknown plot type: {plot_type}")
            plt.close(fig)
            return None
        
        plt.tight_layout()
        plt.show()
        
        return fig

    def extract_overlapping_se(self, region_chr, region_start, region_end, region_name=None, 
                            extraction_method='pybedtools'):
        """
        Extract super enhancers that overlap with the specified genomic region
        
        Parameters:
        region_chr (str): Chromosome name (e.g., 'chr1')
        region_start (int): Start position of the region
        region_end (int): End position of the region
        region_name (str, optional): Name for this region (used in plots and reports)
        extraction_method (str): Method to use for extraction ('pybedtools' or 'pandas')
        
        Returns:
        SEgeneRegionAnalyzer: self for method chaining
        """
        if self.bed_df is None:
            self.logger.warning("BED data not loaded. Call load_databases first.")
            return self
        
        self.logger.info(f"Extracting SEs overlapping with region: {region_chr}:{region_start}-{region_end}")
        
        # Store region information
        self.region = {
            'chr': region_chr,
            'start': region_start,
            'end': region_end,
            'name': region_name or f"{region_chr}:{region_start}-{region_end}"
        }
        
        # Store extraction method
        self.extraction_method = extraction_method
        
        try:
            if extraction_method == 'pybedtools':
                self._extract_with_pybedtools()
            elif extraction_method == 'pandas':
                self._extract_with_pandas()
            else:
                self.logger.warning(f"Unknown extraction method: {extraction_method}. Using pybedtools.")
                self._extract_with_pybedtools()
            
            # If overlapping_se is successfully created, analyze it
            if hasattr(self, 'overlapping_se') and self.overlapping_se is not None:
                self._analyze_overlapping_se()
            
            return self
            
        except Exception as e:
            self.logger.exception(f"Error extracting overlapping SEs: {e}")
            return self

    def _extract_with_pybedtools(self):
        """
        Extract overlapping SEs using pybedtools (internal method)
        Uses pybedtools to find SEs that overlap with the specified region.
        """
        # Get region information
        region_chr = self.region['chr']
        region_start = self.region['start']
        region_end = self.region['end']
        
        # Prepare BED data for pybedtools
        required_columns = ['se_chr', 'se_start', 'se_end']
        
        # Check if all needed columns exist
        missing_cols = [col for col in required_columns if col not in self.bed_df.columns]
        if missing_cols:
            self.logger.warning(f"Missing columns in BED data: {missing_cols}")
            self.logger.warning("Trying to find alternative column names...")
            
            # Try to find alternative column names
            alt_cols = {'se_chr': ['chr'], 'se_start': ['start'], 'se_end': ['end']}
            col_mapping = {}
            
            for req_col in missing_cols:
                for alt_col in alt_cols.get(req_col, []):
                    if alt_col in self.bed_df.columns:
                        col_mapping[alt_col] = req_col
                        self.logger.info(f"Using '{alt_col}' for '{req_col}'")
                        break
            
            # Check if we found alternatives for all missing columns
            still_missing = [col for col in missing_cols if col not in col_mapping.values()]
            if still_missing:
                self.logger.error(f"Could not find alternatives for columns: {still_missing}. Cannot proceed.")
                return
            
            # Create a temporary copy with renamed columns
            temp_df = self.bed_df.copy()
            temp_df.rename(columns=col_mapping, inplace=True)
        else:
            temp_df = self.bed_df
        
        # Convert to BedTool object
        try:
            # Extract just the necessary columns for BED format
            bed_cols = ['se_chr', 'se_start', 'se_end', 'cell_id']
            se_df_bed = temp_df[bed_cols].copy()
            # Rename to standard BED format
            se_df_bed.columns = ['chrom', 'start', 'end', 'name']
            
            # Convert to BedTool
            se_bedtool = pybedtools.BedTool.from_dataframe(se_df_bed)
            
            # Define the region of interest
            interval_str = f"{region_chr}\t{region_start}\t{region_end}"
            interval_bedtool = pybedtools.BedTool(interval_str, from_string=True)
            
            # Find overlapping SEs
            intersect_result = se_bedtool.intersect(interval_bedtool, wa=True)
            
            # Convert result back to DataFrame
            result_df = pd.read_table(
                intersect_result.fn, 
                names=['chrom', 'start', 'end', 'name']
            )
            
            # Extract the overlapping cell_ids
            overlapping_cell_ids = set(result_df['name'])
            
            # Get the original rows that correspond to the overlapping SEs
            self.overlapping_se = temp_df[temp_df['cell_id'].isin(overlapping_cell_ids)].copy()
            
            self.logger.info(f"Found {len(self.overlapping_se)} SEs overlapping with the region")
            
        except Exception as e:
            self.logger.exception(f"Error using pybedtools: {e}")
            return

    def _extract_with_pandas(self):
        """
        Extract overlapping SEs using pandas (internal method)
        A pure pandas approach that doesn't require pybedtools
        """
        # Get region information
        region_chr = self.region['chr']
        region_start = self.region['start']
        region_end = self.region['end']
        
        # Check if we have the necessary columns
        required_columns = ['se_chr', 'se_start', 'se_end']
        alt_columns = {'se_chr': 'chr', 'se_start': 'start', 'se_end': 'end'}
        
        col_mapping = {}
        for req_col in required_columns:
            if req_col in self.bed_df.columns:
                col_mapping[req_col] = req_col
            elif alt_columns[req_col] in self.bed_df.columns:
                col_mapping[alt_columns[req_col]] = req_col
            else:
                self.logger.error(f"Could not find column for {req_col}. Cannot proceed.")
                return
        
        # Create a filtered dataframe with renamed columns
        chr_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_chr')
        start_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_start')
        end_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_end')
        
        # Filter SEs that overlap with the region
        self.overlapping_se = self.bed_df[
            (self.bed_df[chr_col] == region_chr) & 
            (self.bed_df[start_col] <= region_end) & 
            (self.bed_df[end_col] >= region_start)
        ].copy()
        
        self.logger.info(f"Found {len(self.overlapping_se)} SEs overlapping with the region")

    def _analyze_overlapping_se(self):
        """
        Analyze the overlapping SEs to calculate tissue, biosample, and gene counts.
        Counts are based on unique samples (each sample is counted only once).
        This is an internal method called by extract_overlapping_se().
        """
        if self.overlapping_se is None or self.sample_info is None:
            self.logger.warning("No overlapping SEs or sample info available.")
            return
        
        self.logger.info("Analyzing overlapping SEs based on unique samples")
        
        # Get unique cell_ids from overlapping SEs
        overlapping_cell_ids = set(self.overlapping_se['cell_id'])
        self.unique_sample_count = len(overlapping_cell_ids)
        self.logger.info(f"Found {self.unique_sample_count} unique samples with overlapping SEs")
        
        # Filter sample_info to get only rows corresponding to unique overlapping cell_ids
        # Each sample is represented exactly once in this filtered dataframe
        overlapping_samples = self.sample_info[self.sample_info['cell_id'].isin(overlapping_cell_ids)]
        
        # Count tissues - each sample contributes exactly one count
        if 'Tissue type' in overlapping_samples.columns:
            self.region_tissue_counts = overlapping_samples['Tissue type'].value_counts()
            self.logger.info(f"Found {len(self.region_tissue_counts)} tissue types across {self.unique_sample_count} unique samples")
        else:
            self.logger.warning("'Tissue type' column not found in sample information")
        
        # Count biosamples - each sample contributes exactly one count
        if 'Biosample type' in overlapping_samples.columns:
            self.region_biosample_counts = overlapping_samples['Biosample type'].value_counts()
            self.logger.info(f"Found {len(self.region_biosample_counts)} biosample types across {self.unique_sample_count} unique samples")
        else:
            self.logger.warning("'Biosample type' column not found in sample information")
        
        # Count genes (if gene column exists)
        # For gene counting, we're still using all overlapping SEs since multiple SEs from the same sample
        # can be associated with different genes
        gene_column = next((col for col in self.overlapping_se.columns if 'gene' in col.lower()), None)
        if gene_column:
            all_genes = []
            for entry in self.overlapping_se[gene_column].dropna():
                if isinstance(entry, str):
                    genes = entry.split(',')
                    all_genes.extend([g.strip() for g in genes if g.strip()])
            
            self.region_gene_counts = pd.Series(Counter(all_genes))
            self.logger.info(f"Found {len(self.region_gene_counts)} unique genes in overlapping SEs")
        else:
            self.logger.warning("No gene column found in BED data")
            
        # Create a consolidated dataframe with all overlapping SEs and their metadata
        try:
            # Merge overlapping_se with sample_info to get tissue types etc.
            self.final_data = pd.merge(
                self.overlapping_se, 
                self.sample_info,
                on='cell_id',
                how='left'
            )
            
            # Add region information
            for key, value in self.region.items():
                self.final_data[f'region_{key}'] = value
                
            self.logger.info(f"Created consolidated data with {len(self.final_data)} rows from {self.unique_sample_count} unique samples")
            
            # Also create a sample-based summary dataframe where each sample appears only once
            self.sample_summary = overlapping_samples.copy()
            
            # Add SE count per sample
            se_counts = self.overlapping_se['cell_id'].value_counts().reset_index()
            se_counts.columns = ['cell_id', 'se_count']
            self.sample_summary = pd.merge(
                self.sample_summary,
                se_counts,
                on='cell_id',
                how='left'
            )
            
            self.logger.info(f"Created sample summary data with {len(self.sample_summary)} rows (one per unique sample)")
            
        except Exception as e:
            self.logger.warning(f"Could not create consolidated data: {e}")
            
        # Store unique SE data
        self.get_unique_se_data()

    def analyze_samples_distribution(self, figsize=(12, 8)):
        """
        Analyze and visualize the distribution of samples and their SEs in the region.
        
        Parameters:
        figsize (tuple): Figure size for the plots
        
        Returns:
        None: Displays plots and prints summary statistics
        """
        if not hasattr(self, 'sample_summary') or self.sample_summary is None:
            self.logger.warning("Sample summary data not available. Run extract_overlapping_se first.")
            return
        
        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        
        print(f"=== Sample Distribution Analysis for Region: {region_name} ===")
        print(f"Total unique samples with overlapping SEs: {len(self.sample_summary)}")
        
        # SE count per sample statistics
        se_count_stats = self.sample_summary['se_count'].describe()
        print("\nSE Count per Sample Statistics:")
        print(se_count_stats)
        
        # Plot SE count distribution
        plt.figure(figsize=figsize)
        
        plt.subplot(2, 2, 1)
        sns.histplot(data=self.sample_summary, x='se_count', bins=20)
        plt.title('Distribution of SE Counts per Sample')
        plt.xlabel('Number of SEs per Sample')
        plt.ylabel('Number of Samples')
        
        # Pie chart of samples by tissue type
        if 'Tissue type' in self.sample_summary.columns:
            plt.subplot(2, 2, 2)
            tissue_counts = self.sample_summary['Tissue type'].value_counts()
            top_tissues = tissue_counts.head(8)
            
            # If there are more than 8 tissues, group the rest as "Other"
            if len(tissue_counts) > 8:
                other_count = tissue_counts[8:].sum()
                top_tissues = pd.concat([top_tissues, pd.Series({'Other': other_count})])
            
            top_tissues.plot.pie(autopct='%1.1f%%', startangle=90)
            plt.title('Distribution of Samples by Tissue Type')
            plt.ylabel('')
        
        # Top samples by SE count
        plt.subplot(2, 2, 3)
        top_samples = self.sample_summary.sort_values('se_count', ascending=False).head(10)
        sns.barplot(data=top_samples, y='cell_id', x='se_count', palette='viridis')
        plt.title('Top 10 Samples by SE Count')
        plt.xlabel('Number of SEs')
        plt.ylabel('Sample ID')
        
        # SE count vs Tissue Type
        if 'Tissue type' in self.sample_summary.columns:
            plt.subplot(2, 2, 4)
            sns.boxplot(data=self.sample_summary, y='Tissue type', x='se_count')
            plt.title('SE Count by Tissue Type')
            plt.xlabel('Number of SEs')
            plt.ylabel('Tissue Type')
        
        plt.tight_layout()
        plt.show()
        
        # Print additional statistics
        if 'Tissue type' in self.sample_summary.columns:
            print("\nSample Count by Tissue Type (each sample counted once):")
            tissue_counts = self.sample_summary['Tissue type'].value_counts()
            tissue_percent = (tissue_counts / len(self.sample_summary) * 100).round(2)
            tissue_df = pd.DataFrame({
                'Sample Count': tissue_counts,
                'Percentage': tissue_percent
            })
            print(tissue_df)
        
        if 'Biosample type' in self.sample_summary.columns:
            print("\nSample Count by Biosample Type (each sample counted once):")
            biosample_counts = self.sample_summary['Biosample type'].value_counts()
            biosample_percent = (biosample_counts / len(self.sample_summary) * 100).round(2)
            biosample_df = pd.DataFrame({
                'Sample Count': biosample_counts,
                'Percentage': biosample_percent
            })
            print(biosample_df)

    # def analyze_region_tissue_enrichment(self, method='fisher', correction='fdr_bh', plot=True, **kwargs):
    #     """
    #     Analyze tissue enrichment in the overlapping SEs compared to the background distribution.
    #     Counts are based on unique samples (each sample is counted only once).
        
    #     Parameters:
    #     method (str): Statistical test method ('fisher' or 'chi2')
    #     correction (str): Multiple testing correction method
    #     plot (bool): Whether to plot the results
    #     **kwargs: Additional parameters passed to visualize_tissue_enrichment
        
    #     Returns:
    #     pandas.DataFrame: Results of the enrichment analysis
    #     """
    #     if self.region_tissue_counts is None:
    #         self.logger.warning("No region tissue counts available. Call extract_overlapping_se first.")
    #         return None
        
    #     region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
    #     self.logger.info(f"Analyzing tissue enrichment for region {region_name} based on unique samples")
        
    #     # Make sure we have the overall tissue distribution
    #     if not hasattr(self, '_tissue_distribution'):
    #         self.logger.warning("Overall tissue distribution not available.")
    #         # Try to calculate it automatically
    #         self.visualize_tissue_distribution(store_all=True)
            
    #         if not hasattr(self, '_tissue_distribution'):
    #             self.logger.error("Could not calculate tissue distribution. Aborting.")
    #             return None
        
    #     # Perform enrichment test
    #     enrichment_results = self.test_tissue_enrichment(
    #         self.region_tissue_counts,
    #         method=method,
    #         correction=correction
    #     )
        
    #     # Rename columns to use English consistently
    #     english_column_mapping = {
    #         '組織型': 'Tissue Type',
    #         'リージョンカウント': 'Region Count',
    #         'リージョン内割合': 'Region Percentage',
    #         '全体カウント': 'Total Count',
    #         '全体内割合': 'Total Percentage',
    #         '期待値': 'Expected Value',
    #         'オッズ比': 'Odds Ratio',
    #         'フォールドチェンジ': 'Fold Change',
    #         'p値': 'p-value',
    #         '補正p値': 'adjusted p-value',
    #         '有意': 'Significant'
    #     }
        
    #     # Rename columns if they exist
    #     for old_col, new_col in english_column_mapping.items():
    #         if old_col in enrichment_results.columns:
    #             enrichment_results = enrichment_results.rename(columns={old_col: new_col})
        
    #     # Store results
    #     self.region_tissue_enrichment = enrichment_results
        
    #     # Visualize if requested
    #     if plot:
    #         self.visualize_tissue_enrichment(
    #             enrichment_results=enrichment_results,
    #             region_name=region_name,
    #             **kwargs
    #         )
        
    #     self.logger.info("Tissue enrichment analysis completed based on unique samples")
        
    #     return enrichment_results


    def analyze_region_tissue_enrichment(self, method='fisher', correction='fdr_bh', plot=True, 
                                        save_dir=None, save_filename=None, save_formats=None, 
                                        save_dpi=600, save_transparent=False, **kwargs):
        """
        Analyze tissue enrichment in the overlapping SEs compared to the background distribution.
        Counts are based on unique samples (each sample is counted only once).

        Parameters:
            method (str): Statistical test method ('fisher' or 'chi2')
            correction (str): Multiple testing correction method
            plot (bool): Whether to plot the results
            save_dir (str): Directory path to save the figure. If provided, the figure will be saved.
            save_filename (str): Base filename (without extension) for the saved figure.
                                If not provided, a default name will be generated.
            save_formats (list): List of formats to save the figure, e.g. ['png', 'svg']. 
                                If None, the figure is not saved.
            save_dpi (int): DPI for saving the figure (default: 600)
            save_transparent (bool): Whether to save the figure with transparent background.
            **kwargs: Additional parameters passed to visualize_tissue_enrichment

        Returns:
            pandas.DataFrame: Results of the enrichment analysis
        """
        if self.region_tissue_counts is None:
            self.logger.warning("No region tissue counts available. Call extract_overlapping_se first.")
            return None

        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        self.logger.info(f"Analyzing tissue enrichment for region {region_name} based on unique samples")

        # Ensure overall tissue distribution is available
        if not hasattr(self, '_tissue_distribution'):
            self.logger.warning("Overall tissue distribution not available.")
            self.visualize_tissue_distribution(store_all=True)
            if not hasattr(self, '_tissue_distribution'):
                self.logger.error("Could not calculate tissue distribution. Aborting.")
                return None

        # Perform enrichment test
        enrichment_results = self.test_tissue_enrichment(
            self.region_tissue_counts,
            method=method,
            correction=correction
        )

        # Rename columns to use English consistently
        english_column_mapping = {
            '組織型': 'Tissue Type',
            'リージョンカウント': 'Region Count',
            'リージョン内割合': 'Region Percentage',
            '全体カウント': 'Total Count',
            '全体内割合': 'Total Percentage',
            '期待値': 'Expected Value',
            'オッズ比': 'Odds Ratio',
            'フォールドチェンジ': 'Fold Change',
            'p値': 'p-value',
            '補正p値': 'adjusted p-value',
            '有意': 'Significant'
        }
        for old_col, new_col in english_column_mapping.items():
            if old_col in enrichment_results.columns:
                enrichment_results = enrichment_results.rename(columns={old_col: new_col})

        # Store results
        self.region_tissue_enrichment = enrichment_results

        # Visualize if requested
        if plot:
            fig = self.visualize_tissue_enrichment(
                enrichment_results=enrichment_results,
                region_name=region_name,
                **kwargs
            )
            # Save figure if saving options are provided
            if save_formats:
                import os
                if save_dir is None:
                    save_dir = os.getcwd()
                os.makedirs(save_dir, exist_ok=True)
                if save_filename is None:
                    from datetime import datetime
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    save_filename = f"tissue_enrichment_{region_name.replace(':', '_').replace('-', '_')}_{timestamp}"
                for fmt in save_formats:
                    file_path = os.path.join(save_dir, f"{save_filename}.{fmt}")
                    fig.savefig(file_path, format=fmt, dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                    self.logger.info(f"Figure saved to {file_path}")
            plt.show()

        self.logger.info("Tissue enrichment analysis completed based on unique samples")
        return enrichment_results





    def get_unique_se_data(self):
        """
        Get unique SE data for the region (one row per unique SE)
        
        Returns:
        pandas.DataFrame: Unique SE data
        """
        if self.overlapping_se is None:
            self.logger.warning("No overlapping SEs. Call extract_overlapping_se first.")
            return None
        
        # Group by SE coordinates and get unique SEs
        try:
            if 'se_id' in self.overlapping_se.columns:
                # If we have SE IDs, use those
                unique_se = self.overlapping_se.drop_duplicates('se_id')
            else:
                # Otherwise, use coordinates
                coords_cols = [col for col in ['se_chr', 'se_start', 'se_end'] 
                            if col in self.overlapping_se.columns]
                if len(coords_cols) < 3:
                    self.logger.warning("Missing coordinate columns, using fallback.")
                    coords_cols = [col for col in ['chr', 'start', 'end'] 
                                if col in self.overlapping_se.columns]
                
                if len(coords_cols) < 3:
                    self.logger.error("Cannot identify coordinate columns. Aborting.")
                    return None
                    
                unique_se = self.overlapping_se.drop_duplicates(coords_cols)
            
            self.unique_se_data = unique_se
            self.logger.info(f"Found {len(unique_se)} unique SEs in the region")
            
            return unique_se
            
        except Exception as e:
            self.logger.exception(f"Error getting unique SE data: {e}")
            return None

    def display_region_summary(self):
        """
        Display a summary of the region analysis based on unique samples.
        Each sample is counted only once in the statistics.
        """
        if self.region is None or self.overlapping_se is None:
            self.logger.warning("No region analysis results. Call extract_overlapping_se first.")
            return
        
        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        region_size = self.region['end'] - self.region['start']
        
        print(f"=== Region Analysis Summary: {region_name} ===")
        print(f"Coordinates: {self.region['chr']}:{self.region['start']}-{self.region['end']}")
        print(f"Region size: {region_size:,} bp")
        print(f"Extraction method: {self.extraction_method}")
        
        # SE statistics
        total_ses = len(self.overlapping_se)
        unique_ses = len(self.unique_se_data) if hasattr(self, 'unique_se_data') and self.unique_se_data is not None else 0
        unique_samples = getattr(self, 'unique_sample_count', 0)
        
        print(f"\nTotal overlapping SEs: {total_ses:,}")
        print(f"Unique SEs: {unique_ses:,}")
        print(f"Unique samples: {unique_samples:,}")
        
        if unique_samples > 0:
            avg_ses_per_sample = total_ses / unique_samples
            print(f"Average SEs per sample: {avg_ses_per_sample:.2f}")
        
        if self.region_tissue_counts is not None:
            print(f"\nTissue types (based on unique samples): {len(self.region_tissue_counts)}")
            top_tissues = self.region_tissue_counts.nlargest(5)
            print("\nTop Tissues (sample-based counts):")
            for tissue, count in top_tissues.items():
                print(f"  {tissue}: {count} samples")
        
        if self.region_biosample_counts is not None:
            print(f"\nBiosample types (based on unique samples): {len(self.region_biosample_counts)}")
            top_biosamples = self.region_biosample_counts.nlargest(5)
            print("\nTop Biosamples (sample-based counts):")
            for biosample, count in top_biosamples.items():
                print(f"  {biosample}: {count} samples")
        
        if self.region_gene_counts is not None:
            print(f"\nGenes: {len(self.region_gene_counts):,}")
            top_genes = self.region_gene_counts.nlargest(5)
            print("\nTop Genes:")
            for gene, count in top_genes.items():
                print(f"  {gene}: {count}")
            print("\nNote: Gene counts represent the frequency of genes associated with SEs.")
            print("      The same gene may be annotated to multiple SEs, even from the same sample.")
        
        if hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None:
            p_col = 'adjusted p-value' if 'adjusted p-value' in self.region_tissue_enrichment.columns else 'p-value'
            enriched_tissues = self.region_tissue_enrichment[
                (self.region_tissue_enrichment['Fold Change'] > 1.5) & 
                (self.region_tissue_enrichment[p_col] < 0.05)
            ].sort_values('Fold Change', ascending=False)
            
            if len(enriched_tissues) > 0:
                print("\nEnriched Tissues (p < 0.05, FC > 1.5):")
                for _, row in enriched_tissues.head(5).iterrows():
                    print(f"  {row['Tissue Type']}: {row['Region Count']} samples, "
                        f"{row['Fold Change']:.2f}-fold, p={row[p_col]:.1e}")




    def visualize_region_tissues(self, top_n=20, figsize=(12, 8), plot_type='bar',
                            save_dir=None, save_filename=None, 
                            save_svg=False, save_png=False, save_pdf=False, save_eps=False,
                            save_dpi=300, save_transparent=False):
        """
        Visualize the tissue distribution in the region based on unique samples.
        Each tissue is counted only once per sample.
        
        Parameters:
        top_n (int): Number of top tissues to show
        figsize (tuple): Figure size
        plot_type (str): Plot type ('bar', 'pie', 'horizontal')
        save_dir (str, optional): Directory to save the figure
        save_filename (str, optional): Base filename for saved figure (without extension)
        save_svg, save_png, save_pdf, save_eps (bool, optional): Whether to save in respective formats
        save_dpi (int, optional): DPI for saved images
        save_transparent (bool, optional): Whether to save with transparent background
        
        Returns:
        matplotlib.figure.Figure: The plot figure
        """
        if self.region_tissue_counts is None:
            self.logger.warning("No region tissue counts. Call extract_overlapping_se first.")
            return None
        
        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        
        # Get top N tissues
        top_tissues = self.region_tissue_counts.nlargest(top_n)
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        if plot_type == 'bar':
            bars = ax.bar(top_tissues.index, top_tissues.values, 
                        color=plt.cm.tab10(np.arange(len(top_tissues)) % 10))
            ax.set_ylabel('Count (unique samples)')
            ax.set_title(f'Top {len(top_tissues)} Tissues in Region: {region_name} (based on unique samples)')
            plt.xticks(rotation=45, ha='right')
            
            # Add counts above bars
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{height}', ha='center', va='bottom')
        
        elif plot_type == 'pie':
            ax.pie(top_tissues.values, labels=top_tissues.index, autopct='%1.1f%%',
                colors=plt.cm.tab10(np.arange(len(top_tissues)) % 10))
            ax.set_title(f'Tissue Distribution in Region: {region_name} (based on unique samples)')
            ax.axis('equal')
        
        elif plot_type == 'horizontal':
            bars = ax.barh(top_tissues.index, top_tissues.values, 
                        color=plt.cm.tab10(np.arange(len(top_tissues)) % 10))
            ax.set_xlabel('Count (unique samples)')
            ax.set_title(f'Top {len(top_tissues)} Tissues in Region: {region_name} (based on unique samples)')
            
            # Add counts next to bars
            for bar in bars:
                width = bar.get_width()
                ax.text(width + 0.1, bar.get_y() + bar.get_height()/2.,
                    f'{width}', ha='left', va='center')
        
        plt.tight_layout()

        if save_dir is not None and (save_svg or save_png or save_pdf or save_eps):
            import os
            os.makedirs(save_dir, exist_ok=True)
            
            if save_filename is None:
                save_filename = f"tissue_distribution_{region_name.replace(':', '_').replace('-', '_')}"
            
            base_path = os.path.join(save_dir, save_filename)
            
            saved_files = []
            
            if save_png:
                plt.savefig(f"{base_path}.png", format='png', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"PNG: {base_path}.png")
                
            if save_svg:
                plt.savefig(f"{base_path}.svg", format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {base_path}.svg")
                
            if save_pdf:
                plt.savefig(f"{base_path}.pdf", format='pdf', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"PDF: {base_path}.pdf")
                
            if save_eps:
                plt.savefig(f"{base_path}.eps", format='eps', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"EPS: {base_path}.eps")
                
            if saved_files:
                self.logger.info(f"Tissue distribution visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return fig




    def refine_tissue_enrichment_results(self, enrichment_results=None, fdr_threshold=0.05, fc_threshold=1.5, 
                                    top_n=None, rename_columns=True):
        """
        Refine tissue enrichment analysis results and add significance flags.
        
        Parameters:
        -----------
        enrichment_results : pandas.DataFrame, optional
            Analysis results DataFrame. If None, use the last analysis results
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        top_n : int, optional
            Number of top tissues to display. If None, display all
        rename_columns : bool, optional
            Whether to convert column names to English
            
        Returns:
        --------
        pandas.DataFrame
            Refined analysis results
        """
        # Get analysis results
        if enrichment_results is None:
            if hasattr(self, 'region_tissue_enrichment'):
                enrichment_results = self.region_tissue_enrichment
            else:
                self.logger.warning("No enrichment analysis results available. Please run analyze_region_tissue_enrichment() first.")
                return None
        
        # Copy the results for manipulation
        results = enrichment_results.copy()
        
        # Check and rename FDR column
        if '補正p値' in results.columns:
            fdr_col = '補正p値'
        elif 'adjusted p-value' in results.columns:
            fdr_col = 'adjusted p-value'
        else:
            # Use regular p-value if FDR column is not available
            fdr_col = 'p値' if 'p値' in results.columns else 'p-value'
            self.logger.warning(f"FDR-corrected p-value not found. Using {fdr_col} instead.")
        
        # Check p-value column
        if 'p値' in results.columns:
            p_col = 'p値'
        else:
            p_col = 'p-value'
        
        # Check fold change column
        if 'フォールドチェンジ' in results.columns:
            fc_col = 'フォールドチェンジ'
        else:
            fc_col = 'Fold Change'
        
        # Check tissue type column
        if '組織型' in results.columns:
            tissue_col = '組織型'
        else:
            tissue_col = 'Tissue Type'
        
        # Check count columns
        if 'リージョンカウント' in results.columns:
            region_count_col = 'リージョンカウント'
        else:
            region_count_col = 'Region Count'
            
        if '全体カウント' in results.columns:
            total_count_col = '全体カウント'
        else:
            total_count_col = 'Total Count'
        
        # Add significance column
        results['Significant'] = (results[fdr_col] < fdr_threshold) & (results[fc_col] > fc_threshold)
        
        # Check for negative enrichment (depletion)
        results['Depleted'] = (results[fdr_col] < fdr_threshold) & (results[fc_col] < (1/fc_threshold))
        
        # Sort for display (by FDR)
        results = results.sort_values(by=fdr_col)
        
        # Limit to top N entries
        if top_n is not None and top_n > 0:
            results = results.head(top_n)
        
        # Rename columns to English
        if rename_columns:
            column_mapping = {
                '組織型': 'Tissue Type',
                'リージョンカウント': 'Region Count',
                '全体カウント': 'Total Count',
                'フォールドチェンジ': 'Fold Change',
                '補正p値': 'FDR',
                'p値': 'p-value',
                '有意': 'Significant',
                'デプリート': 'Depleted'
            }
            
            # Only rename existing columns
            rename_cols = {col: column_mapping[col] for col in results.columns if col in column_mapping}
            results = results.rename(columns=rename_cols)
        
        # Select and order main columns for display
        main_columns = [
            tissue_col if not rename_columns else 'Tissue Type',
            region_count_col if not rename_columns else 'Region Count',
            total_count_col if not rename_columns else 'Total Count',
            fc_col if not rename_columns else 'Fold Change',
            fdr_col if not rename_columns else 'FDR',
            'Significant' if rename_columns else '有意',
            'Depleted' if rename_columns else 'デプリート'
        ]
        
        # Select only columns that exist
        available_columns = [col for col in main_columns if col in results.columns]
        results = results[available_columns]
        
        return results

    def save_enrichment_results_to_csv(self, results=None, filename=None, **kwargs):
        """
        Save tissue enrichment analysis results to a CSV file.
        
        Parameters:
        -----------
        results : pandas.DataFrame, optional
            Analysis results to save. If None, run refine_tissue_enrichment_results
        filename : str, optional
            Target file name. If None, generate using region name and timestamp
        **kwargs : dict
            Additional parameters to pass to refine_tissue_enrichment_results
            
        Returns:
        --------
        str
            Path to the saved CSV file
        """
        # Get or refine analysis results
        if results is None:
            results = self.refine_tissue_enrichment_results(**kwargs)
            
        if results is None or len(results) == 0:
            self.logger.warning("No analysis results to save.")
            return None
        
        # Set target file name
        if filename is None:
            # Get region name
            if hasattr(self, 'region') and self.region is not None:
                region_name = self.region.get('name', f"{self.region['chr']}_{self.region['start']}-{self.region['end']}")
            else:
                region_name = "unknown_region"
            
            # Generate timestamp
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Generate file name
            filename = os.path.join(self.results_dir, f"tissue_enrichment_{region_name}_{timestamp}.csv")
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(os.path.abspath(filename)), exist_ok=True)
        
        # Save as CSV
        try:
            results.to_csv(filename, index=False, encoding='utf-8-sig')
            self.logger.info(f"Enrichment analysis results saved: {filename}")
            return filename
        except Exception as e:
            self.logger.error(f"File save error: {e}")
            return None

    def visualize_enrichment_results(self, results=None, plot_type='barplot', figsize=(12, 8), 
                                    top_n=15, fdr_threshold=0.05, fc_threshold=1.5,
                                    show_depleted=False, save_fig=False, fig_filename=None,
                                    **kwargs):
        """
        Visualize tissue enrichment analysis results.
        
        Parameters:
        -----------
        results : pandas.DataFrame, optional
            Analysis results to visualize. If None, run refine_tissue_enrichment_results
        plot_type : str, optional
            Plot type ('barplot', 'horizontal', 'volcano', 'dotplot')
        figsize : tuple, optional
            Figure size
        top_n : int, optional
            Number of top tissues to display
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        show_depleted : bool, optional
            Whether to display depleted tissues
        save_fig : bool, optional
            Whether to save the figure
        fig_filename : str, optional
            Filename for the saved figure
        **kwargs : dict
            Additional parameters to pass to refine_tissue_enrichment_results
            
        Returns:
        --------
        matplotlib.figure.Figure
            Created figure
        """
        # Get or refine analysis results
        if results is None:
            results = self.refine_tissue_enrichment_results(
                fdr_threshold=fdr_threshold, 
                fc_threshold=fc_threshold,
                rename_columns=True,
                **kwargs
            )
            
        if results is None or len(results) == 0:
            self.logger.warning("No analysis results to visualize.")
            return None
        
        # Get region name (for title)
        if hasattr(self, 'region') and self.region is not None:
            region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        else:
            region_name = "Unknown Region"
        
        # Prepare data
        if 'Significant' in results.columns:
            # Extract enriched tissues only
            if not show_depleted:
                enriched = results[results['Significant']].sort_values('FDR')
            else:
                # Extract enriched or depleted tissues
                enriched = results[(results['Significant']) | (results['Depleted'])].sort_values('FDR')
        else:
            # No Significant column, use thresholds directly
            if 'FDR' in results.columns:
                fdr_col = 'FDR'
            elif 'adjusted p-value' in results.columns:
                fdr_col = 'adjusted p-value'
            else:
                fdr_col = 'p-value'
                
            if 'Fold Change' in results.columns:
                fc_col = 'Fold Change'
            else:
                fc_col = 'フォールドチェンジ'
                
            if not show_depleted:
                enriched = results[(results[fdr_col] < fdr_threshold) & 
                                (results[fc_col] > fc_threshold)].sort_values(fdr_col)
            else:
                enriched = results[((results[fdr_col] < fdr_threshold) & 
                                (results[fc_col] > fc_threshold)) |
                                ((results[fdr_col] < fdr_threshold) & 
                                (results[fc_col] < 1/fc_threshold))].sort_values(fdr_col)
        
        # Limit to top N
        if top_n is not None and top_n > 0 and len(enriched) > top_n:
            # Show warning
            if len(enriched) > top_n:
                self.logger.info(f"Limiting display to top {top_n} tissues (out of {len(enriched)} total)")
            plot_data = enriched.head(top_n)
        else:
            plot_data = enriched
            
        # Use default plot type if not specified
        if plot_type not in ['barplot', 'horizontal', 'volcano', 'dotplot']:
            self.logger.warning(f"Unknown plot type: {plot_type}. Using 'barplot' instead.")
            plot_type = 'barplot'
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # Check column names
        tissue_col = 'Tissue Type' if 'Tissue Type' in plot_data.columns else '組織型'
        fc_col = 'Fold Change' if 'Fold Change' in plot_data.columns else 'フォールドチェンジ'
        if 'FDR' in plot_data.columns:
            fdr_col = 'FDR'
        elif 'adjusted p-value' in plot_data.columns:
            fdr_col = 'adjusted p-value'
        else:
            fdr_col = 'p-value'
        region_count_col = 'Region Count' if 'Region Count' in plot_data.columns else 'リージョンカウント'
        
        if plot_type == 'barplot':
            # Vertical bar chart
            bars = plt.bar(
                plot_data[tissue_col], 
                plot_data[fc_col],
                color=[
                    'red' if row['Significant'] else 
                    'blue' if row.get('Depleted', False) else 
                    'gray' 
                    for _, row in plot_data.iterrows()
                ] if 'Significant' in plot_data.columns else plt.cm.tab10(range(len(plot_data)))
            )
            
            # Reference line (fold change = 1)
            plt.axhline(y=1, linestyle='--', color='gray', alpha=0.7)
            
            # Axis labels and title
            plt.xlabel('Tissue Type')
            plt.ylabel('Fold Change')
            plt.title(f'Tissue Enrichment in {region_name} (FDR < {fdr_threshold})')
            
            # Display FDR values on bars
            for bar, fdr in zip(bars, plot_data[fdr_col]):
                height = bar.get_height()
                plt.text(
                    bar.get_x() + bar.get_width()/2,
                    height,
                    f'FDR={fdr:.1e}',
                    ha='center',
                    va='bottom',
                    rotation=90,
                    fontsize=8
                )
            
            # Rotate x-axis labels
            plt.xticks(rotation=45, ha='right')
            
        elif plot_type == 'horizontal':
            # Horizontal bar chart (better for long tissue names)
            # Reverse order to display highest significance at top
            plot_data = plot_data.iloc[::-1].copy()
            
            bars = plt.barh(
                plot_data[tissue_col], 
                plot_data[fc_col],
                color=[
                    'red' if row['Significant'] else 
                    'blue' if row.get('Depleted', False) else 
                    'gray' 
                    for _, row in plot_data.iterrows()
                ] if 'Significant' in plot_data.columns else plt.cm.tab10(range(len(plot_data)))
            )
            
            # Reference line (fold change = 1)
            plt.axvline(x=1, linestyle='--', color='gray', alpha=0.7)
            
            # Axis labels and title
            plt.ylabel('Tissue Type')
            plt.xlabel('Fold Change')
            plt.title(f'Tissue Enrichment in {region_name} (FDR < {fdr_threshold})')
            
            # Display FDR values and counts on bars
            for i, (bar, fdr, count) in enumerate(zip(bars, plot_data[fdr_col], plot_data[region_count_col])):
                width = bar.get_width()
                plt.text(
                    width + 0.1,
                    i,
                    f'FDR={fdr:.1e} (n={count})',
                    va='center',
                    ha='left',
                    fontsize=8
                )
        
        elif plot_type == 'volcano':
            # Volcano plot (use all data)
            # Use original data
            all_data = results
            
            # Calculate -log10(FDR)
            log_fdr = -np.log10(all_data[fdr_col])
            
            # Point size proportional to count
            sizes = all_data[region_count_col] * 5
            
            # Color coding
            if 'Significant' in all_data.columns and 'Depleted' in all_data.columns:
                colors = [
                    'red' if row['Significant'] else 
                    'blue' if row['Depleted'] else 
                    'gray' 
                    for _, row in all_data.iterrows()
                ]
            else:
                colors = [
                    'red' if (row[fdr_col] < fdr_threshold and row[fc_col] > fc_threshold) else 
                    'blue' if (row[fdr_col] < fdr_threshold and row[fc_col] < 1/fc_threshold) else 
                    'gray' 
                    for _, row in all_data.iterrows()
                ]
            
            # Scatter plot
            plt.scatter(
                all_data[fc_col],
                log_fdr,
                s=sizes,
                c=colors,
                alpha=0.7
            )
            
            # Threshold lines
            plt.axhline(y=-np.log10(fdr_threshold), linestyle='--', color='gray', alpha=0.7)
            plt.axvline(x=fc_threshold, linestyle='--', color='red', alpha=0.7)
            plt.axvline(x=1/fc_threshold, linestyle='--', color='blue', alpha=0.7)
            plt.axvline(x=1, linestyle='-', color='black', alpha=0.3)
            
            # Axis labels and title
            plt.xlabel('Fold Change')
            plt.ylabel('-log10(FDR)')
            plt.title(f'Tissue Enrichment Volcano Plot for {region_name}')
            
            # Label significant points
            for _, row in plot_data.iterrows():
                tissue = row[tissue_col]
                x = row[fc_col]
                y = -np.log10(row[fdr_col])
                
                plt.annotate(
                    tissue,
                    xy=(x, y),
                    xytext=(5, 0),
                    textcoords='offset points',
                    ha='left',
                    va='center',
                    fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3', fc='yellow', alpha=0.3)
                )
            
            # Grid lines
            plt.grid(True, linestyle='--', alpha=0.3)
        
        elif plot_type == 'dotplot':
            # Dot plot (size and color represent enrichment and significance)
            # Reverse order to display highest significance at top
            plot_data = plot_data.iloc[::-1].copy()
            
            # Use simple index for X-axis
            plt.scatter(
                plot_data[fc_col],
                range(len(plot_data)),
                s=plot_data[region_count_col] * 10,  # Size proportional to count
                c=-np.log10(plot_data[fdr_col]),     # Color based on FDR
                cmap='plasma',
                alpha=0.8
            )
            
            # Reference line (fold change = 1)
            plt.axvline(x=1, linestyle='--', color='gray', alpha=0.7)
            
            # Axis labels and title
            plt.xlabel('Fold Change')
            plt.yticks(range(len(plot_data)), plot_data[tissue_col])
            plt.title(f'Tissue Enrichment in {region_name} (bubble size = count, color = -log10(FDR))')
            
            # Color bar
            cbar = plt.colorbar()
            cbar.set_label('-log10(FDR)')
            
            # Display FDR values
            for i, (fc, fdr, count) in enumerate(zip(plot_data[fc_col], plot_data[fdr_col], plot_data[region_count_col])):
                plt.text(
                    fc + 0.1,
                    i,
                    f'FDR={fdr:.1e} (n={count})',
                    va='center',
                    ha='left',
                    fontsize=8
                )
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure if requested
        if save_fig:
            if fig_filename is None:
                # Set target file name
                import datetime
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                fig_filename = os.path.join(self.results_dir, f"tissue_enrichment_{region_name}_{plot_type}_{timestamp}.png")
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(os.path.abspath(fig_filename)), exist_ok=True)
            
            try:
                plt.savefig(fig_filename, dpi=300, bbox_inches='tight')
                self.logger.info(f"Enrichment analysis figure saved: {fig_filename}")
            except Exception as e:
                self.logger.error(f"Figure save error: {e}")
        
        plt.show()
        return fig

    def visualize_tissue_significance(self, results=None, figsize=(12, 8), 
                                fdr_threshold=0.05, fc_threshold=1.5,
                                show_all_tissues=False, max_tissues=20,
                                sort_by='significance', # 'significance', 'fold_change', 'count'
                                horizontal=True, show_counts=True,
                                show_stars=True, save_fig=False, fig_filename=None):
        """
        Display a bar chart focusing on tissue significance.
        
        Parameters:
        -----------
        results : pandas.DataFrame, optional
            Analysis results to visualize. If None, run refine_tissue_enrichment_results
        figsize : tuple, optional
            Figure size
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        show_all_tissues : bool, optional
            Whether to show all tissues (False = show only significant ones)
        max_tissues : int, optional
            Maximum number of tissues to display
        sort_by : str, optional
            Sort method ('significance'=by significance, 'fold_change'=by fold change, 'count'=by count)
        horizontal : bool, optional
            Whether to use horizontal bar chart
        show_counts : bool, optional
            Whether to display counts
        show_stars : bool, optional
            Whether to display significance stars
        save_fig : bool, optional
            Whether to save the figure
        fig_filename : str, optional
            Filename for the saved figure
            
        Returns:
        --------
        matplotlib.figure.Figure
            Created figure
        """
        # Get or refine analysis results
        if results is None:
            results = self.refine_tissue_enrichment_results(
                fdr_threshold=fdr_threshold, 
                fc_threshold=fc_threshold,
                rename_columns=True
            )
            
        if results is None or len(results) == 0:
            self.logger.warning("No analysis results to visualize.")
            return None
        
        # Check column names
        tissue_col = 'Tissue Type' if 'Tissue Type' in results.columns else '組織型'
        fc_col = 'Fold Change' if 'Fold Change' in results.columns else 'フォールドチェンジ'
        fdr_col = 'FDR' if 'FDR' in results.columns else ('adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value')
        region_count_col = 'Region Count' if 'Region Count' in results.columns else 'リージョンカウント'
        total_count_col = 'Total Count' if 'Total Count' in results.columns else '全体カウント'
        sig_col = 'Significant' if 'Significant' in results.columns else '有意'
        
        # Add significance column if missing
        if sig_col not in results.columns:
            results[sig_col] = (results[fdr_col] < fdr_threshold) & (results[fc_col] > fc_threshold)
        
        # Prepare data - show only significant tissues or all
        if not show_all_tissues:
            plot_data = results[results[sig_col]].copy()
            if len(plot_data) == 0:
                self.logger.warning(f"No tissues meet the significance criteria: FDR < {fdr_threshold}, Fold Change > {fc_threshold}.")
                plot_data = results.head(min(max_tissues, len(results))).copy()
        else:
            plot_data = results.copy()
        
        # Sort
        if sort_by == 'significance':
            # Sort by significance + FDR (significant first, then by FDR)
            plot_data['_sort_key'] = ~plot_data[sig_col]  # Significant first
            plot_data = plot_data.sort_values(['_sort_key', fdr_col])
        elif sort_by == 'fold_change':
            # Sort by fold change (descending)
            plot_data = plot_data.sort_values(fc_col, ascending=False)
        elif sort_by == 'count':
            # Sort by count (descending)
            plot_data = plot_data.sort_values(region_count_col, ascending=False)
        else:
            self.logger.warning(f"Unknown sort method: {sort_by}. Sorting by significance.")
            plot_data['_sort_key'] = ~plot_data[sig_col]
            plot_data = plot_data.sort_values(['_sort_key', fdr_col])
        
        # Limit to maximum display count
        if max_tissues > 0 and len(plot_data) > max_tissues:
            plot_data = plot_data.head(max_tissues)
        
        # Get region name (for title)
        if hasattr(self, 'region') and self.region is not None:
            region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        else:
            region_name = "Unknown Region"
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Horizontal or vertical bar chart
        if horizontal:
            # Reverse index order (display from top to bottom)
            plot_data = plot_data.iloc[::-1].reset_index(drop=True)
            
            # Create bar chart
            bars = ax.barh(
                plot_data[tissue_col],
                plot_data[fc_col],
                color=[
                    'red' if row[sig_col] else 'gray'
                    for _, row in plot_data.iterrows()
                ],
                alpha=0.7
            )
            
            # Reference line (fold change = 1)
            ax.axvline(x=1, linestyle='--', color='black', alpha=0.5)
            
            # Display significance stars
            if show_stars:
                for i, (_, row) in enumerate(plot_data.iterrows()):
                    if row[sig_col]:
                        ax.text(
                            0.1,  # X position (slightly left of axis)
                            i,    # Y position (tissue index)
                            "★",  # Star symbol
                            ha='right',
                            va='center',
                            fontsize=14,
                            color='darkred',
                            weight='bold'
                        )
            
            # Display counts and FDR values
            if show_counts:
                for i, (_, row) in enumerate(plot_data.iterrows()):
                    fc = row[fc_col]
                    fdr = row[fdr_col]
                    count = row[region_count_col]
                    total = row[total_count_col]
                    
                    # Create text content
                    text = f"FDR={fdr:.1e}"
                    if show_counts:
                        text += f", {count}/{total} samples"
                    
                    ax.text(
                        max(fc + 0.1, 1.1),  # X position (after bar)
                        i,                   # Y position (tissue index)
                        text,
                        va='center',
                        ha='left',
                        fontsize=9
                    )
            
            # Axis labels and title
            ax.set_xlabel('Fold Change')
            ax.set_ylabel('Tissue Type')
            
        else:
            # Vertical bar chart
            bars = ax.bar(
                plot_data[tissue_col],
                plot_data[fc_col],
                color=[
                    'red' if row[sig_col] else 'gray'
                    for _, row in plot_data.iterrows()
                ],
                alpha=0.7
            )
            
            # Reference line (fold change = 1)
            ax.axhline(y=1, linestyle='--', color='black', alpha=0.5)
            
            # Display significance stars
            if show_stars:
                for i, (_, row) in enumerate(plot_data.iterrows()):
                    if row[sig_col]:
                        ax.text(
                            i,      # X position (bar center)
                            0.5,    # Y position (slightly above bottom)
                            "★",    # Star symbol
                            ha='center',
                            va='bottom',
                            fontsize=14,
                            color='darkred',
                            weight='bold'
                        )
            
            # Display counts and FDR values
            if show_counts:
                for i, (bar, (_, row)) in enumerate(zip(bars, plot_data.iterrows())):
                    fc = row[fc_col]
                    fdr = row[fdr_col]
                    count = row[region_count_col]
                    total = row[total_count_col]
                    
                    # Create text content
                    text = f"FDR={fdr:.1e}"
                    if show_counts:
                        text += f"\n{count}/{total}"
                    
                    ax.text(
                        bar.get_x() + bar.get_width()/2,  # X position (bar center)
                        fc + 0.1,                        # Y position (above bar)
                        text,
                        ha='center',
                        va='bottom',
                        fontsize=8,
                        rotation=90
                    )
            
            # Rotate x-axis labels
            plt.xticks(rotation=45, ha='right')
            
            # Axis labels and title
            ax.set_xlabel('Tissue Type')
            ax.set_ylabel('Fold Change')
        
        # Add significance threshold line
        if fc_threshold > 1:
            if horizontal:
                ax.axvline(x=fc_threshold, linestyle=':', color='red', alpha=0.5)
            else:
                ax.axhline(y=fc_threshold, linestyle=':', color='red', alpha=0.5)
        
        # Set title
        sig_count = plot_data[sig_col].sum()
        total_count = len(plot_data)
        ax.set_title(f'Tissue Enrichment in {region_name}\n'
                f'{sig_count} significant tissues (FDR < {fdr_threshold}, FC > {fc_threshold})'
                f'{" - showing top "+str(max_tissues) if max_tissues < total_count else ""}')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='red', alpha=0.7, label=f'Significant (FDR < {fdr_threshold}, FC > {fc_threshold})'),
            Patch(facecolor='gray', alpha=0.7, label='Not significant')
        ]
        ax.legend(handles=legend_elements, loc='best')
        
        # Adjust layout
        plt.tight_layout()
        
        # Save figure if requested
        if save_fig:
            if fig_filename is None:
                # Set target file name
                import datetime
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                orientation = "horizontal" if horizontal else "vertical"
                fig_filename = os.path.join(self.results_dir, 
                                        f"tissue_significance_{region_name}_{orientation}_{timestamp}.png")
            
            # Ensure directory exists
            os.makedirs(os.path.dirname(os.path.abspath(fig_filename)), exist_ok=True)
            
            try:
                plt.savefig(fig_filename, dpi=300, bbox_inches='tight')
                self.logger.info(f"Tissue significance figure saved: {fig_filename}")
            except Exception as e:
                self.logger.error(f"Figure save error: {e}")
        
        plt.show()
        return fig

    def save_tissue_significance_results(self, fdr_threshold=0.05, fc_threshold=1.5, filename=None):
        """
        Save tissue significance analysis results to CSV and display summary.
        
        Parameters:
        -----------
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        filename : str, optional
            Target file name. If None, generate automatically
            
        Returns:
        --------
        tuple
            (Path to saved CSV file, number of significant tissues)
        """
        # Refine analysis results
        results = self.refine_tissue_enrichment_results(
            fdr_threshold=fdr_threshold, 
            fc_threshold=fc_threshold,
            rename_columns=True
        )
        
        if results is None or len(results) == 0:
            self.logger.warning("No analysis results to save.")
            return None, 0
        
        # Extract significant tissues
        sig_col = 'Significant' if 'Significant' in results.columns else '有意'
        sig_results = results[results[sig_col]]
        
        # Display results summary
        print(f"=== Tissue Enrichment Analysis Summary ===")
        print(f"Total tissues: {len(results)}")
        print(f"Significant tissues: {len(sig_results)} (FDR < {fdr_threshold}, Fold Change > {fc_threshold})")
        
        if len(sig_results) > 0:
            print("\n--- Top 10 Significant Tissues ---")
            fdr_col = 'FDR' if 'FDR' in results.columns else ('adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value')
            tissue_col = 'Tissue Type' if 'Tissue Type' in results.columns else '組織型'
            fc_col = 'Fold Change' if 'Fold Change' in results.columns else 'フォールドチェンジ'
            rc_col = 'Region Count' if 'Region Count' in results.columns else 'リージョンカウント'
            tc_col = 'Total Count' if 'Total Count' in results.columns else '全体カウント'
            
            # Sort by FDR
            top_sig = sig_results.sort_values(fdr_col).head(10)
            
            for i, (_, row) in enumerate(top_sig.iterrows(), 1):
                tissue = row[tissue_col]
                fc = row[fc_col]
                fdr = row[fdr_col]
                rc = row[rc_col]
                tc = row[tc_col]
                
                print(f"{i}. {tissue}: {fc:.2f}x (FDR={fdr:.1e}, {rc}/{tc} samples)")
        
        # Save to CSV
        csv_path = self.save_enrichment_results_to_csv(results, filename)
        
        return csv_path, len(sig_results)

    def create_enrichment_summary_table(self, results=None, fdr_threshold=0.05, fc_threshold=1.5, 
                                    sort_by='significance', max_tissues=None):
        """
        Create a summary table of tissue enrichment analysis results.
        
        Parameters:
        -----------
        results : pandas.DataFrame, optional
            Analysis results. If None, use the last analysis results
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        sort_by : str, optional
            Sort method ('significance', 'fold_change', 'count')
        max_tissues : int, optional
            Maximum number of tissues to display
            
        Returns:
        --------
        pandas.DataFrame
            Formatted summary table
        """
        # Get or refine analysis results
        if results is None:
            results = self.refine_tissue_enrichment_results(
                fdr_threshold=fdr_threshold, 
                fc_threshold=fc_threshold,
                rename_columns=True
            )
            
        if results is None or len(results) == 0:
            self.logger.warning("No analysis results to display.")
            return None
        
        # Check column names
        tissue_col = 'Tissue Type' if 'Tissue Type' in results.columns else '組織型'
        fc_col = 'Fold Change' if 'Fold Change' in results.columns else 'フォールドチェンジ'
        fdr_col = 'FDR' if 'FDR' in results.columns else ('adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value')
        region_count_col = 'Region Count' if 'Region Count' in results.columns else 'リージョンカウント'
        total_count_col = 'Total Count' if 'Total Count' in results.columns else '全体カウント'
        sig_col = 'Significant' if 'Significant' in results.columns else '有意'
        
        # Add significance column if missing
        if sig_col not in results.columns:
            results[sig_col] = (results[fdr_col] < fdr_threshold) & (results[fc_col] > fc_threshold)
        
        # Sort
        if sort_by == 'significance':
            # Sort by significance + FDR (significant first, then by FDR)
            results['_sort_key'] = ~results[sig_col]  # Significant first
            sorted_results = results.sort_values(['_sort_key', fdr_col])
        elif sort_by == 'fold_change':
            # Sort by fold change (descending)
            sorted_results = results.sort_values(fc_col, ascending=False)
        elif sort_by == 'count':
            # Sort by count (descending)
            sorted_results = results.sort_values(region_count_col, ascending=False)
        else:
            self.logger.warning(f"Unknown sort method: {sort_by}. Sorting by significance.")
            results['_sort_key'] = ~results[sig_col]
            sorted_results = results.sort_values(['_sort_key', fdr_col])
        
        # Limit to maximum display count
        if max_tissues is not None and max_tissues > 0:
            display_results = sorted_results.head(max_tissues)
        else:
            display_results = sorted_results
        
        # Create summary table
        summary = pd.DataFrame({
            'Tissue Type': display_results[tissue_col],
            'Fold Change': display_results[fc_col].round(2),
            'FDR': display_results[fdr_col].apply(lambda x: f"{x:.1e}"),
            'Count': display_results[region_count_col].astype(int).astype(str) + '/' + 
                    display_results[total_count_col].astype(int).astype(str),
            'Significant': display_results[sig_col].map({True: '★', False: ''})
        })
        
        # Get region name
        if hasattr(self, 'region') and self.region is not None:
            region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        else:
            region_name = "Unknown Region"
        
        # Display summary information
        print(f"=== Tissue Enrichment Analysis: {region_name} ===")
        print(f"Thresholds: FDR < {fdr_threshold}, Fold Change > {fc_threshold}")
        print(f"Significant tissues: {display_results[sig_col].sum()} / {len(display_results)}")
        
        if '_sort_key' in summary.columns:
            summary = summary.drop(columns=['_sort_key'])
        
        return summary










    def plot_stacked_reads_simple(self, chrom=None, start=None, end=None,
                                title=None, xlabel="Genomic Position", ylabel="Read Stack",
                                color="blue", grid=True, linewidth=1.0,
                                save_dir=None, save_filename=None, 
                                save_svg=False, save_png=False, save_pdf=False, save_eps=False,
                                save_dpi=600, save_transparent=False, save_region_bed=None, 
                                show_plot=True, title_fontsize=18, xlabel_fontsize=16, 
                                ylabel_fontsize=16, tick_fontsize=12, annotation_fontsize=12, 
                                figsize=(15, 8)):
        """
        Display reads in the specified genomic region in a simple stacked layout.
        
        Parameters:
        -----------
        chrom : str, optional
            Chromosome name. If not specified, retrieved from self.region
        start : int, optional
            Start position. If not specified, retrieved from self.region
        end : int, optional
            End position. If not specified, retrieved from self.region
        title : str, optional
            Plot title
        xlabel : str, optional
            X-axis label
        ylabel : str, optional
            Y-axis label
        color : str, optional
            Color for all reads (default: blue)
        grid : bool, optional
            Whether to display grid lines
        linewidth : float, optional
            Width of read lines
        save_dir : str, optional
            Directory to save output files. If None, current directory is used.
        save_filename : str, optional
            Base filename for saved files (without extension). If None, auto-generated from region.
        save_svg : bool, optional
            Whether to save as SVG format (default: False)
        save_png : bool, optional
            Whether to save as PNG format (default: False)
        save_pdf : bool, optional
            Whether to save as PDF format (default: False)
        save_eps : bool, optional
            Whether to save as EPS format (default: False)
        save_dpi : int, optional
            DPI (resolution) for raster formats like PNG (default: 600)
        save_transparent : bool, optional
            Whether to save with transparent background (default: False)
        save_region_bed : str, optional
            Path to save reads in the region as a BED file
        show_plot : bool, optional
            Whether to display the plot
        title_fontsize : int, optional
            Font size for the plot title (default: 18)
        xlabel_fontsize : int, optional
            Font size for the x-axis label (default: 16)
        ylabel_fontsize : int, optional
            Font size for the y-axis label (default: 16)
        tick_fontsize : int, optional
            Font size for the axis tick labels (default: 12)
        annotation_fontsize : int, optional
            Font size for text annotations like sample count (default: 12)
        figsize : tuple, optional
            Figure size as (width, height) in inches (default: (15, 8))
        
        Returns:
        --------
        matplotlib.figure.Figure
            The created plot's Figure object
        """
        from matplotlib.ticker import ScalarFormatter
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import os
        
        # Check and retrieve position information
        if chrom is None or start is None or end is None:
            if self.region is not None:
                chrom = self.region.get('chr')
                start = self.region.get('start')
                end = self.region.get('end')
            else:
                self.logger.error("Region information not specified. Please specify chrom, start, end or set region first.")
                return None
        
        # Check if we have overlapping data
        if self.final_data is None or self.final_data.empty:
            self.logger.error("No overlapping data found. Please run extract_overlapping_se first.")
            return None
        
        # Filter data for the specified region
        region_data = self.final_data[(self.final_data['se_chr'] == chrom) & 
                                (self.final_data['se_start'] >= start) & 
                                (self.final_data['se_end'] <= end)]
        
        if region_data.empty:
            self.logger.warning(f"No reads found in the specified region: {chrom}:{start}-{end}")
            return None
        
        # Count the number of samples
        sample_count = region_data['cell_id'].nunique()
        
        # Display read count
        self.logger.info(f"Number of reads to display: {len(region_data)}")
        self.logger.info(f"Number of samples: {sample_count}")
        
        # Sort reads by start position
        region_data = region_data.sort_values('se_start')
        
        # Save region BED file if requested
        if save_region_bed:
            try:
                from pybedtools import BedTool
                # Create a BED format DataFrame
                bed_df = region_data[['se_chr', 'se_start', 'se_end', 'cell_id']].copy()
                bed_df.columns = ['chrom', 'start', 'end', 'name']
                region_bed = BedTool.from_dataframe(bed_df)
                region_bed.saveas(save_region_bed)
                self.logger.info(f"Reads in region saved as BED file: {save_region_bed}")
            except Exception as e:
                self.logger.error(f"Error saving BED file: {e}")
        
        # Convert DataFrame to list of read dictionaries
        reads = []
        for _, row in region_data.iterrows():
            reads.append({
                'chrom': row['se_chr'],
                'start': row['se_start'],
                'end': row['se_end'],
                'sample': row['cell_id']
            })
        
        # Sort reads by start position
        reads.sort(key=lambda r: r['start'])
        
        # Improved greedy interval scheduling algorithm for read stacking with margin
        margin = 5  # Small margin to prevent visual overlap in adjacent reads
        stacks = []
        for read in reads:
            placed = False
            for stack_idx, stack in enumerate(stacks):
                # Find the first stack where this read doesn't overlap with any existing read
                if stack and read['start'] > (stack[-1]['end'] + margin):
                    stack.append(read)
                    placed = True
                    break
            
            if not placed:
                # Create a new stack
                stacks.append([read])
        
        # Create the plot
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()  # Get current axes
        
        # Explicitly turn off all gridlines first
        ax.grid(False)
        
        # Display each read as a horizontal line
        for stack_idx, stack in enumerate(stacks):
            y_pos = stack_idx
            for read in stack:
                plt.hlines(y=y_pos, xmin=read['start'], xmax=read['end'], color=color, linewidth=linewidth)
        
        # Set title
        if title is None:
            title = f"Stacked read plot for {chrom}:{start:,}-{end:,}"
        plt.title(title, fontsize=title_fontsize)
        
        # Set axis labels
        plt.xlabel(xlabel, fontsize=xlabel_fontsize)
        plt.ylabel(ylabel, fontsize=ylabel_fontsize)
        
        # Set axis limits
        plt.xlim(start, end)
        plt.ylim(-1, len(stacks))
        
        # Configure x-axis for better display of large numbers
        plt.gca().xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        plt.ticklabel_format(style='plain', axis='x')
        
        # Set tick font sizes
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Add vertical grid lines only
        if grid:
            ax.xaxis.grid(True, linestyle='--', color='lightgray', alpha=0.7)
            ax.yaxis.grid(False)  # Explicitly disable horizontal grid lines
        
        # Add sample count to the plot
        plt.text(0.98, 0.98, f"Total: {sample_count} samples", 
            transform=ax.transAxes, ha='right', va='top',
            fontsize=annotation_fontsize,
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
        
        # Prepare for saving in multiple formats
        if save_svg or save_png or save_pdf or save_eps:
            # Create directory if needed
            if save_dir is None:
                save_dir = os.getcwd()  # Current directory
            elif not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
                self.logger.info(f"Created output directory: {save_dir}")
            
            # Generate base filename if not provided
            if save_filename is None:
                save_filename = f"stacked_reads_{chrom}_{start}_{end}"
            
            # Create full path without extension
            base_path = os.path.join(save_dir, save_filename)
            
            # Save in each requested format
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved plot in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        # Show plot
        if show_plot:
            plt.tight_layout()
            plt.show()
        else:
            plt.close(fig)
        
        return fig






    def plot_stacked_reads_with_peaks(self, chrom=None, start=None, end=None,
                                title=None, xlabel="Genomic position", ylabel="Read stack",
                                color="blue", grid=True, linewidth=1.0,
                                save_dir=None, save_filename=None, 
                                save_svg=False, save_png=False, save_pdf=False, save_eps=False,
                                save_dpi=600, save_transparent=False, save_region_bed=None, 
                                show_plot=True, title_fontsize=18, xlabel_fontsize=16, 
                                ylabel_fontsize=16, tick_fontsize=12, annotation_fontsize=12, 
                                figsize=(15, 8), detect_peaks=True, bin_size=50, peak_threshold=0.5, 
                                min_peak_distance=200, peak_width_factor=0.5,
                                peak_color='red', peak_alpha=0.2, annotate_peaks=True,
                                peak_annotation_fontsize=12, output_peak_bed=None,
                                smooth_window=3):
        """
        Display reads in the specified genomic region with peak detection.
        
        Parameters:
        -----------
        chrom : str, optional
            Chromosome name. If not specified, retrieved from self.region
        start : int, optional
            Start position. If not specified, retrieved from self.region
        end : int, optional
            End position. If not specified, retrieved from self.region
        title : str, optional
            Plot title
        xlabel : str, optional
            X-axis label
        ylabel : str, optional
            Y-axis label
        color : str, optional
            Color for all reads (default: blue)
        grid : bool, optional
            Whether to display grid lines
        linewidth : float, optional
            Width of read lines
        save_dir : str, optional
            Directory to save output files. If None, current directory is used
        save_filename : str, optional
            Base filename for saved files (without extension). If None, auto-generated from region
        save_svg : bool, optional
            Whether to save as SVG format (default: False)
        save_png : bool, optional
            Whether to save as PNG format (default: False)
        save_pdf : bool, optional
            Whether to save as PDF format (default: False)
        save_eps : bool, optional
            Whether to save as EPS format (default: False)
        save_dpi : int, optional
            DPI (resolution) for raster formats like PNG (default: 600)
        save_transparent : bool, optional
            Whether to save with transparent background (default: False)
        save_region_bed : str, optional
            Path to save reads in the region as a BED file
        show_plot : bool, optional
            Whether to display the plot
        title_fontsize : int, optional
            Font size for the plot title (default: 18)
        xlabel_fontsize : int, optional
            Font size for the x-axis label (default: 16)
        ylabel_fontsize : int, optional
            Font size for the y-axis label (default: 16)
        tick_fontsize : int, optional
            Font size for the axis tick labels (default: 12)
        annotation_fontsize : int, optional
            Font size for text annotations like sample count (default: 12)
        figsize : tuple, optional
            Figure size as (width, height) in inches (default: (15, 8))
        detect_peaks : bool, optional
            Whether to detect peaks in the region (default: True)
        bin_size : int, optional
            Size of bins for read counting and peak detection (default: 50)
        peak_threshold : float, optional
            Threshold for peak detection (default: 0.5)
        min_peak_distance : int, optional
            Minimum distance between peaks (default: 200)
        peak_width_factor : float, optional
            Factor to determine peak width (default: 0.5)
        peak_color : str, optional
            Color for peak highlighting (default: 'red')
        peak_alpha : float, optional
            Alpha value for peak highlighting (default: 0.2)
        annotate_peaks : bool, optional
            Whether to annotate peaks with labels (default: True)
        peak_annotation_fontsize : int, optional
            Font size for peak annotations (default: 12)
        output_peak_bed : str, optional
            Path to save detected peaks as a BED file
        smooth_window : int, optional
            Window size for smoothing bin counts (default: 3)
            
        Returns:
        --------
        tuple
            (matplotlib.figure.Figure, dict of detected peaks)
            
            The peak information dictionary contains:
            {
                'peak_1': {
                    'position': peak position (bp),
                    'start': peak left boundary (bp),
                    'end': peak right boundary (bp),
                    'width': peak width (bp),
                    'reads': {
                        'total': total read count in peak region,
                        'max': maximum read count in any bin within peak,
                        'average': average read count per bin within peak
                    },
                    'samples': {
                        'count': number of samples with reads in peak region,
                        'percentage': percentage of all samples,
                        'list': list of sample IDs with reads in peak region
                    },
                    'normalized_height': normalized peak height (0-1),
                    'density': read density (reads/bp)
                },
                'peak_2': {...},
                ...
            }
        """
        from matplotlib.ticker import ScalarFormatter
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import os
        from scipy import signal
        from scipy.ndimage import gaussian_filter1d
        from collections import defaultdict
        
        # Check and retrieve position information
        if chrom is None or start is None or end is None:
            if self.region is not None:
                chrom = self.region.get('chr')
                start = self.region.get('start')
                end = self.region.get('end')
            else:
                self.logger.error("Region information not specified. Please specify chrom, start, end or set region first.")
                return None, {}
        
        # Check if we have overlapping data
        if self.final_data is None or self.final_data.empty:
            self.logger.error("No overlapping data found. Please run extract_overlapping_se first.")
            return None, {}
        
        # Filter data for the specified region
        region_data = self.final_data[(self.final_data['se_chr'] == chrom) & 
                                (self.final_data['se_start'] >= start) & 
                                (self.final_data['se_end'] <= end)]
        
        if region_data.empty:
            self.logger.warning(f"No reads found in the specified region: {chrom}:{start}-{end}")
            return None, {}
        
        # Count the number of samples
        sample_count = region_data['cell_id'].nunique()
        
        # Display read count
        self.logger.info(f"Number of reads to display: {len(region_data)}")
        self.logger.info(f"Number of samples: {sample_count}")
        
        # Sort reads by start position
        region_data = region_data.sort_values('se_start')
        
        # Save region BED file if requested
        if save_region_bed:
            try:
                from pybedtools import BedTool
                # Create a BED format DataFrame
                bed_df = region_data[['se_chr', 'se_start', 'se_end', 'cell_id']].copy()
                bed_df.columns = ['chrom', 'start', 'end', 'name']
                region_bed = BedTool.from_dataframe(bed_df)
                region_bed.saveas(save_region_bed)
                self.logger.info(f"Reads in region saved as BED file: {save_region_bed}")
            except Exception as e:
                self.logger.error(f"Error saving BED file: {e}")
        
        # Convert DataFrame to list of read dictionaries
        reads = []
        for _, row in region_data.iterrows():
            reads.append({
                'chrom': row['se_chr'],
                'start': row['se_start'],
                'end': row['se_end'],
                'sample': row['cell_id']
            })
        
        # Sort reads by start position
        reads.sort(key=lambda r: r['start'])
        
        # Improved greedy interval scheduling algorithm for read stacking with margin
        margin = 5  # Small margin to prevent visual overlap in adjacent reads
        stacks = []
        for read in reads:
            placed = False
            for stack_idx, stack in enumerate(stacks):
                # Find the first stack where this read doesn't overlap with any existing read
                if stack and read['start'] > (stack[-1]['end'] + margin):
                    stack.append(read)
                    placed = True
                    break
            
            if not placed:
                # Create a new stack
                stacks.append([read])
        
        # Create the plot
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()  # Get current axes
        
        # Explicitly turn off all gridlines first
        ax.grid(False)
        
        # Display each read as a horizontal line
        for stack_idx, stack in enumerate(stacks):
            y_pos = stack_idx
            for read in stack:
                plt.hlines(y=y_pos, xmin=read['start'], xmax=read['end'], color=color, linewidth=linewidth)
        
        # ===== Peak detection part (with enhanced data structure) =====
        peak_data = {}
        peak_display_data = {}  # For figure display (keeps original format)
        
        if detect_peaks:
            # Create bins for the region
            region_range = end - start
            num_bins = int(region_range / bin_size) + 1
            bins = np.linspace(start, end, num_bins)
            bin_counts = np.zeros(num_bins - 1)
            
            # Track samples by bin
            sample_bins = defaultdict(set)
            bin_samples = [set() for _ in range(num_bins - 1)]
            
            # Count reads in each bin
            for read in reads:
                read_start_bin = max(0, int((read['start'] - start) / bin_size))
                read_end_bin = min(num_bins - 2, int((read['end'] - start) / bin_size))
                
                sample_name = read['sample']
                
                for bin_idx in range(read_start_bin, read_end_bin + 1):
                    bin_counts[bin_idx] += 1
                    bin_samples[bin_idx].add(sample_name)
                    sample_bins[sample_name].add(bin_idx)
            
            # Save original bin counts
            original_bin_counts = bin_counts.copy()
            
            # Smooth the bin counts
            if smooth_window > 1:
                smoothed_counts = gaussian_filter1d(bin_counts, sigma=smooth_window)
            else:
                smoothed_counts = bin_counts.copy()
            
            # Normalize counts for peak detection
            if smoothed_counts.max() > 0:
                normalized_counts = smoothed_counts / smoothed_counts.max()
            else:
                normalized_counts = smoothed_counts.copy()
            
            # Find peaks
            actual_threshold = peak_threshold * normalized_counts.max()
            peaks, peak_properties = signal.find_peaks(
                normalized_counts, 
                height=actual_threshold,
                distance=int(min_peak_distance / bin_size)
            )
            
            # Calculate peak positions in genomic coordinates
            peak_positions = start + peaks * bin_size
            peak_heights = peak_properties['peak_heights']
            
            # Calculate peak widths based on threshold
            peak_widths = signal.peak_widths(
                normalized_counts, peaks, rel_height=peak_width_factor
            )
            
            # Convert width indices to genomic coordinates
            peak_left_idx = peak_widths[2].astype(int)
            peak_right_idx = peak_widths[3].astype(int)
            peak_left_pos = start + peak_left_idx * bin_size
            peak_right_pos = start + peak_right_idx * bin_size
            
            # Total unique sample count
            total_samples = len(set(sample['sample'] for sample in reads))
            
            # Generate peak information (for display and return value)
            for i, peak_pos in enumerate(peak_positions):
                peak_id = f"peak_{i+1}"
                
                # Original data format for display (same as old code)
                peak_display_data[peak_id] = {
                    'position': int(peak_pos),
                    'height': float(peak_heights[i]),
                    'left': int(peak_left_pos[i]),
                    'right': int(peak_right_pos[i]),
                    'width': int(peak_right_pos[i] - peak_left_pos[i])
                }
                
                # Enhanced data format for return value
                
                # Get bin indices within peak region
                peak_bin_indices = list(range(
                    max(0, peak_left_idx[i]),
                    min(num_bins - 1, peak_right_idx[i] + 1)
                ))
                
                # Calculate read statistics in peak region
                if peak_bin_indices:
                    peak_region_counts = original_bin_counts[peak_bin_indices]
                    peak_total_reads = int(peak_region_counts.sum())
                    peak_max_reads = int(peak_region_counts.max())
                    peak_avg_reads = float(peak_region_counts.mean())
                else:
                    peak_total_reads = 0
                    peak_max_reads = 0
                    peak_avg_reads = 0
                
                # Calculate unique sample count in peak region
                peak_samples = set()
                for bin_idx in peak_bin_indices:
                    if 0 <= bin_idx < len(bin_samples):
                        peak_samples.update(bin_samples[bin_idx])
                
                peak_sample_count = len(peak_samples)
                peak_sample_percentage = (peak_sample_count / total_samples * 100) if total_samples > 0 else 0
                
                # Peak width
                peak_width = int(peak_right_pos[i] - peak_left_pos[i])
                
                # Read density (reads per bp)
                density = peak_total_reads / peak_width if peak_width > 0 else 0
                
                # Store enhanced peak data with sample list
                peak_data[peak_id] = {
                    'position': int(peak_pos),
                    'start': int(peak_left_pos[i]),
                    'end': int(peak_right_pos[i]),
                    'width': peak_width,
                    'reads': {
                        'total': peak_total_reads,
                        'max': peak_max_reads,
                        'average': round(peak_avg_reads, 2)
                    },
                    'samples': {
                        'count': peak_sample_count,
                        'percentage': round(peak_sample_percentage, 2),
                        'list': sorted(list(peak_samples))  # Added sample list
                    },
                    'normalized_height': float(peak_heights[i]),
                    'density': round(density, 3)
                }
            
            # Sort peaks by position for better labeling
            sorted_peaks = sorted(peak_display_data.items(), key=lambda x: x[1]['position'])
            
            # Calculate peak label positions
            vertical_positions = {}
            for i, (peak_id, peak) in enumerate(sorted_peaks):
                # Position labels inside the top part of the graph
                vertical_positions[peak_id] = len(stacks) * 0.9  # 90% height of the graph
                
                # Adjust position if too close to previous peak
                if i > 0:
                    prev_id, prev_peak = sorted_peaks[i-1]
                    if (peak['position'] - prev_peak['position']) / (end - start) < 0.04:
                        vertical_positions[peak_id] = vertical_positions[prev_id] - len(stacks) * 0.1
            
            # Highlight peaks with rectangles - use display data for visualization
            for peak_id, peak in peak_display_data.items():
                # Highlight peak region with colored rectangle
                rect = plt.Rectangle(
                    (peak['left'], -1),  # (x, y)
                    peak['width'],  # width
                    len(stacks) + 2,  # height
                    facecolor=peak_color,
                    alpha=peak_alpha,
                    zorder=-1
                )
                ax.add_patch(rect)
                
                # Add peak annotation with improved positioning
                if annotate_peaks:
                    plt.text(
                        peak['position'],
                        vertical_positions[peak_id],
                        f"{peak_id}\n{peak['position']:,}",
                        ha='center',
                        va='center',
                        fontsize=peak_annotation_fontsize,
                        bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.2')
                    )
            
            # Output peaks to BED file - use enhanced data format
            if output_peak_bed:
                try:
                    with open(output_peak_bed, 'w') as f:
                        f.write("#chrom\tstart\tend\tpeak_id\tsample_count\tsample_percentage\ttotal_reads\tdensity\tnormalized_height\n")
                        for peak_id, peak in peak_data.items():
                            f.write(f"{chrom}\t{peak['start']}\t{peak['end']}\t{peak_id}\t"
                                f"{peak['samples']['count']}\t{peak['samples']['percentage']}\t"
                                f"{peak['reads']['total']}\t{peak['density']}\t{peak['normalized_height']}\n")
                    self.logger.info(f"Detected peaks saved to BED file: {output_peak_bed}")
                except Exception as e:
                    self.logger.error(f"Error saving peak BED file: {e}")
        
        # Set title
        if title is None:
            title = f"Stacked read plot for {chrom}:{start:,}-{end:,}"
        plt.title(title, fontsize=title_fontsize)
        
        # Set axis labels
        plt.xlabel(xlabel, fontsize=xlabel_fontsize)
        plt.ylabel(ylabel, fontsize=ylabel_fontsize)
        
        # Set axis limits
        plt.xlim(start, end)
        plt.ylim(-1, len(stacks))
        
        # Configure x-axis to display large numbers properly
        plt.gca().xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        plt.ticklabel_format(style='plain', axis='x')
        
        # Set tick font sizes
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Add only vertical grid lines
        if grid:
            ax.xaxis.grid(True, linestyle='--', color='lightgray', alpha=0.7)
            ax.yaxis.grid(False)  # Explicitly disable horizontal grid lines
        
        # Add sample count and peak info in top right corner
        if peak_data:
            info_text = f"Total: {sample_count} samples, {len(peak_data)} peaks"
            plt.text(0.98, 0.98, info_text, 
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=annotation_fontsize,
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
        else:
            # Simple text if no peaks
            plt.text(0.98, 0.98, f"Total: {sample_count} samples", 
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=annotation_fontsize,
                    bbox=dict(facecolor='white', alpha=0.8, boxstyle='round'))
        
        # Prepare for saving in multiple formats
        if save_svg or save_png or save_pdf or save_eps:
            # Create directory if needed
            if save_dir is None:
                save_dir = os.getcwd()  # Current directory
            elif not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
                self.logger.info(f"Created output directory: {save_dir}")
            
            # Generate base filename if not provided
            if save_filename is None:
                save_filename = f"stacked_reads_{chrom}_{start}_{end}"
                if detect_peaks:
                    save_filename += "_with_peaks"
            
            # Create full path without extension
            base_path = os.path.join(save_dir, save_filename)
            
            # Save in each requested format
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved plot in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        # Show plot
        if show_plot:
            plt.tight_layout()
            plt.show()
        else:
            plt.close(fig)
        
        return fig, peak_data














    def analyze_region_comprehensive(self, chrom, start, end, region_name=None, 
                                output_dir=None, run_enrichment=True,
                                peak_threshold=0.2, enrichment_correction='fdr_bh'):
        """
        Perform comprehensive analysis of a genomic region in one step.
        
        This method automates the entire analysis workflow for a genomic region, 
        including overlap detection, distribution analysis, enrichment analysis,
        and read distribution visualization.
        
        Parameters:
        -----------
        chrom : str
            Chromosome name (e.g., 'chr7')
        start : int
            Start position
        end : int
            End position
        region_name : str, optional
            Region name (auto-generated if None)
        output_dir : str, optional
            Output directory (default: self.results_dir)
        run_enrichment : bool, optional
            Whether to run enrichment analysis (default: True)
        peak_threshold : float, optional
            Threshold for peak detection (default: 0.2)
        enrichment_correction : str, optional
            Multiple testing correction method for enrichment (default: 'fdr_bh')
                
        Returns:
        --------
        Dict[str, Any]
            Dictionary containing all analysis results
        """
        from datetime import datetime
        
        analysis_start_time = datetime.now()
        
        # Create results container
        results = {
            'region_info': {
                'chrom': chrom,
                'start': start,
                'end': end,
                'name': region_name or f"{chrom}:{start}-{end}",
                'size': end - start
            },
            'analysis_time': analysis_start_time.strftime("%Y-%m-%d %H:%M:%S"),
            'steps_completed': []
        }
        
        self.logger.info(f"Starting comprehensive analysis for region {chrom}:{start:,}-{end:,}")
        
        # Step 1: Create region of interest and analyze overlaps
        try:
            # Extract overlaps using SEdbRegionAnalyzer's extract_overlapping_se method
            self.extract_overlapping_se(chrom, start, end, region_name=region_name)
            results['interest_region'] = self.region
            results['steps_completed'].append('extract_overlapping_se')
            
            # Store overlapping samples information
            if hasattr(self, 'overlapping_se') and self.overlapping_se is not None:
                results['overlapping_results'] = self.overlapping_se
                results['steps_completed'].append('store_overlapping_results')
            
            # Get unique samples
            if hasattr(self, 'unique_sample_count'):
                results['overlapping_samples_count'] = self.unique_sample_count
                results['steps_completed'].append('calculate_sample_count')
            
            # Get unique samples list
            if hasattr(self, 'sample_summary'):
                sample_list = self.sample_summary['cell_id'].tolist()
                results['overlapping_samples'] = sample_list
                results['steps_completed'].append('store_sample_list')
            
            self.logger.info(f"Found {getattr(self, 'unique_sample_count', 0)} samples overlapping with region")
            
        except Exception as e:
            self.logger.error(f"Error in overlap detection: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            results['errors'] = results.get('errors', []) + [f"Overlap detection: {str(e)}"]
        
        # Step 2: Store tissue and biosample type distributions
        try:
            if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None:
                results['tissue_counts'] = self.region_tissue_counts
                results['steps_completed'].append('analyze_tissue_distribution')
            
            if hasattr(self, 'region_biosample_counts') and self.region_biosample_counts is not None:
                results['biosample_counts'] = self.region_biosample_counts
                results['steps_completed'].append('analyze_biosample_distribution')
        except Exception as e:
            self.logger.error(f"Error in distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            results['errors'] = results.get('errors', []) + [f"Distribution analysis: {str(e)}"]
        
        # Step 3: Run enrichment analysis if requested
        if run_enrichment:
            try:
                # Tissue type enrichment
                tissue_enrichment = self.analyze_region_tissue_enrichment(
                    correction=enrichment_correction,
                    plot=False
                )
                
                if tissue_enrichment is not None:
                    results['tissue_enrichment'] = tissue_enrichment
                    results['steps_completed'].append('analyze_region_tissue_enrichment')
            except Exception as e:
                self.logger.error(f"Error in tissue enrichment analysis: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
                results['errors'] = results.get('errors', []) + [f"Tissue enrichment analysis: {str(e)}"]
        
        # Step 4: Generate read distribution plots
        try:
            # Simple stacked plot result is just the figure, no data returned
            simple_plot = self.plot_stacked_reads_simple(
                chrom=chrom,
                start=start,
                end=end,
                title=f"SE Read Distribution in {results['region_info']['name']}",
                show_plot=False
            )
            results['steps_completed'].append('plot_stacked_reads_simple')
            
            # Peak detection plot returns both figure and peak data
            peak_plot, peak_data = self.plot_stacked_reads_with_peaks(
                chrom=chrom,
                start=start,
                end=end,
                title=f"SE Read Distribution with Peaks in {results['region_info']['name']}",
                peak_threshold=peak_threshold,
                show_plot=False
            )
            
            if peak_data:
                results['peak_data'] = peak_data
                results['steps_completed'].append('plot_stacked_reads_with_peaks')
        except Exception as e:
            self.logger.error(f"Error in read distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            results['errors'] = results.get('errors', []) + [f"Read distribution analysis: {str(e)}"]
        
        # Calculate total analysis time
        analysis_end_time = datetime.now()
        analysis_duration = (analysis_end_time - analysis_start_time).total_seconds()
        results['analysis_duration_seconds'] = analysis_duration
        
        self.logger.info(f"Comprehensive analysis completed in {analysis_duration:.2f} seconds")
        self.logger.info(f"Steps completed: {', '.join(results['steps_completed'])}")
        
        if 'errors' in results:
            self.logger.warning(f"Analysis completed with {len(results['errors'])} errors")
        
        # Store results in instance for later use
        self.comprehensive_results = results
        
        return results






    def generate_region_report(self, analysis_results=None, chrom=None, start=None, end=None, 
                            region_name=None, output_dir='results/region_report',
                            figure_formats=None, save_tables=True, dpi=300, 
                            fig_width=12, fig_height=8, pvalue_threshold=0.05, 
                            fc_threshold=1.5):
        """
        Generate an HTML report from comprehensive region analysis results.
        
        This method can either use previously computed analysis results or 
        automatically run the analysis if given region coordinates.
        
        Parameters:
        -----------
        analysis_results : dict, optional
            Results from analyze_region_comprehensive. If None, coordinates must be provided.
        chrom : str, optional
            Chromosome name (required if analysis_results is None)
        start : int, optional
            Start position (required if analysis_results is None)
        end : int, optional
            End position (required if analysis_results is None)
        region_name : str, optional
            Region name (auto-generated if None)
        output_dir : str, optional
            Base directory for all output files (default: 'results/region_report')
        figure_formats : list, optional
            List of figure formats to save (default: ['png', 'svg', 'pdf'])
        save_tables : bool, optional
            Whether to save distribution tables as CSV files (default: True)
        dpi : int, optional
            Resolution for raster images like PNG (default: 300)
        fig_width : int, optional
            Width of the figure in inches (default: 12)
        fig_height : int, optional
            Height of the figure in inches (default: 8)
        pvalue_threshold : float, optional
            P-value threshold for enrichment analysis display (default: 0.05)
        fc_threshold : float, optional
            Fold change threshold for enrichment analysis display (default: 1.5)
                
        Returns:
        --------
        str
            Path to the generated HTML report
        """
        import os
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        from datetime import datetime
        
        # Run analysis if results not provided
        if analysis_results is None:
            if chrom is None or start is None or end is None:
                self.logger.error("Must provide either analysis_results or region coordinates (chrom, start, end)")
                return None
            
            self.logger.info(f"No analysis results provided. Running comprehensive analysis for {chrom}:{start}-{end}")
            analysis_results = self.analyze_region_comprehensive(
                chrom=chrom, 
                start=start, 
                end=end, 
                region_name=region_name,
                output_dir=output_dir
            )
        
        # Extract region info
        region_info = analysis_results.get('region_info', {})
        if not region_info:
            self.logger.error("No region information in analysis results")
            return None
        
        chrom = region_info.get('chrom')
        start = region_info.get('start')
        end = region_info.get('end')
        region_name = region_info.get('name', f"{chrom}:{start}-{end}")
        
        # Set default figure formats
        if figure_formats is None:
            figure_formats = ['png', 'svg', 'pdf']
        
        # Create timestamp for filenames
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create directory structure with relative paths
        figures_dir = 'figures'  # Use relative paths for portable HTML
        tables_dir = 'tables'
        
        # Create full directory paths
        figures_full_path = os.path.join(output_dir, figures_dir)
        tables_full_path = os.path.join(output_dir, tables_dir)
        
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(figures_full_path, exist_ok=True)
        os.makedirs(tables_full_path, exist_ok=True)
        
        self.logger.info(f"Creating region analysis report in directory: {output_dir}")
        
        # Dictionary to track all generated files
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Prepare save options for figures
        save_opts = {}
        for fmt in figure_formats:
            save_opts[f'save_{fmt}'] = True
        
        # Set Matplotlib params for better quality figures
        plt.rcParams['figure.figsize'] = (fig_width, fig_height)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # 1. Generate tissue distribution visualizations if available
        tissue_dist_figures = []
        if 'region_tissue_counts' in analysis_results or hasattr(self, 'region_tissue_counts'):
            try:
                # Use results from analysis or class attribute
                tissue_data = analysis_results.get('region_tissue_counts', self.region_tissue_counts)
                
                # If tissue_data is a Series, convert to DataFrame for compatibility
                if isinstance(tissue_data, pd.Series):
                    tissue_data = pd.DataFrame({
                        'Tissue Type': tissue_data.index,
                        'Region Count': tissue_data.values
                    })
                
                # Bar chart - use self.visualize_region_tissues() if available
                tissue_bar_basename = f'region_{chrom}_{start}_{end}_tissue_bar_{timestamp}'
                try:
                    # If visualize_region_tissues exists, use it
                    if hasattr(self, 'visualize_region_tissues'):
                        tissue_bar_fig = self.visualize_region_tissues(
                            plot_type='horizontal',
                            figsize=(fig_width, fig_height),
                            save_dir=figures_full_path,
                            save_filename=tissue_bar_basename,
                            **save_opts
                        )

                        for fmt in figure_formats:
                            figure_filename = f'{tissue_bar_basename}.{fmt}'
                            figure_relative_path = os.path.join(figures_dir, figure_filename)
                            figure_full_path = os.path.join(figures_full_path, figure_filename)
                            
                            if os.path.exists(figure_full_path):
                                generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                                if fmt == 'png':
                                    tissue_dist_figures.append(('Tissue Distribution Bar Chart', figure_relative_path))


                    else:
                        # Fallback to direct matplotlib
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
                        sorted_tissues = tissue_data.sort_values('Region Count', ascending=False)
                        top_tissues = sorted_tissues.head(15)
                        bars = ax.barh(top_tissues['Tissue Type'], top_tissues['Region Count'], color='skyblue')
                        ax.set_xlabel('Sample Count')
                        ax.set_ylabel('Tissue Type')
                        ax.set_title(f'Tissue Type Distribution in {region_name}')
                        for i, bar in enumerate(bars):
                            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                                    str(int(bar.get_width())), va='center')
                        plt.tight_layout()
                        
                        # Save figure manually
                        for fmt in figure_formats:
                            figure_filename = f'{tissue_bar_basename}.{fmt}'
                            figure_path = os.path.join(figures_full_path, figure_filename)
                            plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                            
                            # Use relative path for HTML
                            figure_relative_path = os.path.join(figures_dir, figure_filename)
                            generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_dist_figures.append((f'Tissue Distribution Bar Chart', figure_relative_path))
                except Exception as e:
                    self.logger.warning(f"Error generating tissue bar chart: {e}")
                
                # Save tissue distribution table
                if save_tables:
                    tissue_table_filename = f'region_{chrom}_{start}_{end}_tissue_dist_{timestamp}.csv'
                    tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
                    tissue_data.to_csv(tissue_table_path, index=False)
                    generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating tissue distribution visualizations: {e}")
        
        # 2. Generate biosample type distribution visualizations if available
        biosample_dist_figures = []
        if 'region_biosample_counts' in analysis_results or hasattr(self, 'region_biosample_counts'):
            try:
                # Use results from analysis or class attribute
                biosample_data = analysis_results.get('region_biosample_counts', self.region_biosample_counts)
                
                # If biosample_data is a Series, convert to DataFrame for compatibility
                if isinstance(biosample_data, pd.Series):
                    biosample_data = pd.DataFrame({
                        'Biosample Type': biosample_data.index,
                        'Region Count': biosample_data.values
                    })
                
                # Bar chart
                biosample_bar_basename = f'region_{chrom}_{start}_{end}_biosample_bar_{timestamp}'
                try:
                    # If visualize_biosample_distribution exists, use it
                    if hasattr(self, 'visualize_biosample_distribution'):
                        biosample_bar_fig = self.visualize_biosample_distribution(
                            horizontal=True,
                            figsize=(fig_width, fig_height),
                            save_dir=figures_full_path,
                            save_filename=biosample_bar_basename,
                            **save_opts
                        )
                    else:
                        # Fallback to direct matplotlib
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
                        sorted_biosamples = biosample_data.sort_values('Region Count', ascending=False)
                        top_biosamples = sorted_biosamples.head(15)
                        bars = ax.barh(top_biosamples['Biosample Type'], top_biosamples['Region Count'], color='lightgreen')
                        ax.set_xlabel('Sample Count')
                        ax.set_ylabel('Biosample Type')
                        ax.set_title(f'Biosample Type Distribution in {region_name}')
                        for i, bar in enumerate(bars):
                            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, 
                                    str(int(bar.get_width())), va='center')
                        plt.tight_layout()
                        
                        # Save figure manually
                        for fmt in figure_formats:
                            figure_filename = f'{biosample_bar_basename}.{fmt}'
                            figure_path = os.path.join(figures_full_path, figure_filename)
                            plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                            
                            # Use relative path for HTML
                            figure_relative_path = os.path.join(figures_dir, figure_filename)
                            generated_files['figures'].append((f'Biosample Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                biosample_dist_figures.append((f'Biosample Distribution Bar Chart', figure_relative_path))
                except Exception as e:
                    self.logger.warning(f"Error generating biosample bar chart: {e}")
                
                # Save biosample distribution table
                if save_tables:
                    biosample_table_filename = f'region_{chrom}_{start}_{end}_biosample_dist_{timestamp}.csv'
                    biosample_table_path = os.path.join(tables_full_path, biosample_table_filename)
                    biosample_data.to_csv(biosample_table_path, index=False)
                    generated_files['tables'].append(('Biosample Distribution', os.path.join(tables_dir, biosample_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating biosample distribution visualizations: {e}")
        
        # 3. Generate tissue enrichment visualizations
        tissue_enrichment_figures = []
        if 'tissue_enrichment' in analysis_results and analysis_results['tissue_enrichment'] is not None:
            try:
                tissue_enrichment = analysis_results['tissue_enrichment']
                
                # Generate enrichment visualization
                tissue_enrich_basename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}'
                
                # Use visualize_tissue_enrichment method if available
                if hasattr(self, 'visualize_tissue_enrichment'):
                    tissue_enrich_fig = self.visualize_tissue_enrichment(
                        enrichment_results=tissue_enrichment,
                        title=f'Tissue Type Enrichment Analysis - {region_name}',
                        region_name=region_name,
                        pvalue_threshold=pvalue_threshold,
                        fc_threshold=fc_threshold,
                        plot_type='bar_fdr_horizontal',
                        save_dir=figures_full_path,
                        save_filename=tissue_enrich_basename,
                        **save_opts
                    )
                    
                    # Register visualization
                    for fmt in figure_formats:
                        figure_filename = f'{tissue_enrich_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Enrichment Analysis ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_enrichment_figures.append((f'Tissue Enrichment Analysis', figure_relative_path))
                else:
                    # Basic fallback visualization using matplotlib
                    self.logger.warning(f"visualize_tissue_enrichment method not available, using basic visualization")
                
                # Save enrichment table
                if save_tables:
                    tissue_enrich_filename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}.csv'
                    tissue_enrich_path = os.path.join(tables_full_path, tissue_enrich_filename)
                    tissue_enrichment.to_csv(tissue_enrich_path, index=False)
                    generated_files['tables'].append(('Tissue Enrichment Analysis', os.path.join(tables_dir, tissue_enrich_filename)))
            except Exception as e:
                self.logger.warning(f"Error generating tissue enrichment analysis: {e}")
        
        # 4. Generate stacked read plots
        stacked_plot_figures = []
        
        # 4.1 Simple stacked plot
        try:
            # Create a new plot for the report
            simple_stacked_basename = f'region_{chrom}_{start}_{end}_simple_stacked_{timestamp}'
            simple_stacked_fig = self.plot_stacked_reads_simple(
                chrom=chrom,
                start=start,
                end=end,
                title=f"SE Read Distribution in {region_name}",
                save_dir=figures_full_path,
                save_filename=simple_stacked_basename,
                show_plot=False,
                **save_opts
            )
            
            # Register files
            for fmt in figure_formats:
                figure_filename = f'{simple_stacked_basename}.{fmt}'
                figure_relative_path = os.path.join(figures_dir, figure_filename)
                figure_full_path = os.path.join(figures_full_path, figure_filename)
                
                if os.path.exists(figure_full_path):
                    generated_files['figures'].append((f'Simple Stacked Read Plot ({fmt.upper()})', figure_relative_path))
                    if fmt == 'png':
                        stacked_plot_figures.append((f'Simple Read Distribution', figure_relative_path))
        except Exception as e:
            self.logger.warning(f"Error generating simple stacked read plot: {e}")
        
        # 4.2 Stacked plot with peaks
        try:
            peak_stacked_basename = f'region_{chrom}_{start}_{end}_peak_stacked_{timestamp}'
            peak_stacked_fig, peak_data = self.plot_stacked_reads_with_peaks(
                chrom=chrom,
                start=start,
                end=end,
                title=f"SE Read Distribution with Peaks in {region_name}",
                peak_threshold=0.2,  # Use the default or from analysis_results
                save_dir=figures_full_path,
                save_filename=peak_stacked_basename,
                show_plot=False,
                **save_opts
            )
            
            # Register files
            for fmt in figure_formats:
                figure_filename = f'{peak_stacked_basename}.{fmt}'
                figure_relative_path = os.path.join(figures_dir, figure_filename)
                figure_full_path = os.path.join(figures_full_path, figure_filename)
                
                if os.path.exists(figure_full_path):
                    generated_files['figures'].append((f'Stacked Read Plot with Peaks ({fmt.upper()})', figure_relative_path))
                    if fmt == 'png':
                        stacked_plot_figures.append((f'Read Distribution with Peaks', figure_relative_path))
            
            # Get peak data from analysis results if available
            if peak_data is None and 'peak_data' in analysis_results:
                peak_data = analysis_results['peak_data']
            
            # Save peak data if available
            if peak_data and save_tables:
                peak_data_filename = f'region_{chrom}_{start}_{end}_peak_data_{timestamp}.csv'
                peak_data_path = os.path.join(tables_full_path, peak_data_filename)
                
                # Convert peak data to DataFrame
                peak_rows = []
                for peak_id, peak_info in peak_data.items():
                    # Flatten nested structure
                    row = {
                        'peak_id': peak_id,
                        'position': peak_info['position'],
                        'start': peak_info['start'],
                        'end': peak_info['end'],
                        'width': peak_info['width'],
                        'total_reads': peak_info['reads']['total'],
                        'max_reads': peak_info['reads']['max'],
                        'average_reads': peak_info['reads']['average'],
                        'sample_count': peak_info['samples']['count'],
                        'sample_percentage': peak_info['samples']['percentage'],
                        'normalized_height': peak_info['normalized_height'],
                        'density': peak_info['density']
                    }
                    peak_rows.append(row)
                
                if peak_rows:
                    pd.DataFrame(peak_rows).to_csv(peak_data_path, index=False)
                    generated_files['tables'].append(('Peak Analysis Data', os.path.join(tables_dir, peak_data_filename)))
                    
                    # Also save list of samples per peak
                    for peak_id, peak_info in peak_data.items():
                        if peak_info['samples']['list']:
                            peak_samples_filename = f'region_{chrom}_{start}_{end}_{peak_id}_samples_{timestamp}.txt'
                            peak_samples_path = os.path.join(tables_full_path, peak_samples_filename)
                            
                            with open(peak_samples_path, 'w') as f:
                                for sample in peak_info['samples']['list']:
                                    f.write(f"{sample}\n")
                            
                            generated_files['tables'].append((f'Samples in {peak_id}', os.path.join(tables_dir, peak_samples_filename)))
        except Exception as e:
            self.logger.warning(f"Error generating stacked read plot with peaks: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Generate HTML report
        html_filename = f"region_analysis_report_{chrom}_{start}_{end}_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        try:
            with open(html_path, 'w', encoding='utf-8') as f:
                # Start HTML document
                f.write("<!DOCTYPE html>\n")
                f.write("<html>\n")
                f.write("<head>\n")
                f.write("  <title>SEdb Region Analysis Report</title>\n")
                f.write("  <style>\n")
                f.write("    body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }\n")
                f.write("    h1 { color: #2c3e50; }\n")
                f.write("    h2 { color: #3498db; margin-top: 30px; border-bottom: 2px solid #3498db; padding-bottom: 5px; }\n")
                f.write("    h3 { color: #2980b9; margin-top: 20px; }\n")
                f.write("    h4 { color: #27ae60; }\n")
                f.write("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
                f.write("    th { background-color: #3498db; color: white; text-align: left; padding: 8px; }\n")
                f.write("    td { border: 1px solid #ddd; padding: 8px; }\n")
                f.write("    tr:nth-child(even) { background-color: #f2f2f2; }\n")
                f.write("    img { max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; box-shadow: 0 0 5px rgba(0,0,0,0.1); }\n")
                f.write("    .file-formats { margin: 10px 0; }\n")
                f.write("    .file-formats a { margin-right: 10px; color: #3498db; text-decoration: none; }\n")
                f.write("    .file-formats a:hover { text-decoration: underline; }\n")
                f.write("    .chart-section { display: flex; flex-wrap: wrap; justify-content: space-around; }\n")
                f.write("    .chart-container { margin: 15px; max-width: 48%; }\n")
                f.write("    .highlight { background-color: #f8f9fa; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }\n")
                f.write("    .explanation { background-color: #e8f4f8; padding: 15px; border-radius: 5px; margin: 10px 0 20px 0; }\n")
                f.write("    .notes { font-style: italic; color: #777; }\n")
                f.write("    .section-intro { margin-bottom: 20px; }\n")
                f.write("    .important { color: #e74c3c; font-weight: bold; }\n")
                f.write("    .enriched { color: #27ae60; }\n")
                f.write("    .depleted { color: #e74c3c; }\n")
                f.write("  </style>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                
                # Title and timestamp
                f.write(f"  <h1>SEdb Region Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Comprehensive report introduction
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This report provides a comprehensive analysis of super enhancer (SE) data for the selected genomic region. The analysis includes:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li><strong>Distribution Analysis:</strong> Tissue and biosample type distribution in the region</li>\n")
                f.write("      <li><strong>Enrichment Analysis:</strong> Statistical enrichment of tissues</li>\n")
                f.write("      <li><strong>Read Distribution Analysis:</strong> Visualization of SE read patterns</li>\n")
                f.write("      <li><strong>Peak Detection:</strong> Identification of high SE activity regions</li>\n")
                f.write("    </ul>\n")
                f.write("  </div>\n")
                
                # Section 1: Region Information
                f.write("  <h2>1. Region Information</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write(f"    <p><strong>Chromosome:</strong> {chrom}</p>\n")
                f.write(f"    <p><strong>Start Position:</strong> {start:,}</p>\n")
                f.write(f"    <p><strong>End Position:</strong> {end:,}</p>\n")
                f.write(f"    <p><strong>Region Name:</strong> {region_name}</p>\n")
                f.write(f"    <p><strong>Region Size:</strong> {region_info.get('size', end - start):,} bp</p>\n")
                f.write("  </div>\n")
                
                # Section 2: Sample Overview
                sample_stats = {
                    'total_samples': 0,  # Unknown without original metadata
                    'overlapping_samples': analysis_results.get('overlapping_samples_count', getattr(self, 'unique_sample_count', 0)),
                }
                
                f.write("  <h2>2. Sample Overview</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section summarizes the samples in the SEdb database that overlap with the selected genomic region.</p>\n")
                f.write("  </div>\n")
                
                f.write("  <table>\n")
                f.write("    <tr><th>Metric</th><th>Value</th></tr>\n")
                f.write(f"    <tr><td>Samples Overlapping with Region</td><td>{sample_stats['overlapping_samples']:,}</td></tr>\n")
                f.write("  </table>\n")
                
                # Section 3: Distribution Analysis
                f.write("  <h2>3. Tissue and Biosample Type Distribution</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section shows the distribution of tissue types and biosample types among samples that overlap with the selected genomic region.</p>\n")
                f.write("    <p>The distribution charts help identify which tissues and biosample types are most frequent in this region.</p>\n")
                f.write("  </div>\n")
                
                # 3.1 Tissue type distribution
                if tissue_dist_figures:
                    f.write("  <h3>3.1 Tissue Type Distribution</h3>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>The charts below show the distribution of tissue types in samples overlapping with this region.</p>\n")
                    f.write("    <p>This information is useful for understanding which tissues most commonly have super enhancers in this genomic region.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display tissue distribution charts
                    f.write("  <div class='chart-container' style='max-width: 90%; margin: 20px auto;'>\n")
                    for figure_name, figure_path in tissue_dist_figures:
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                    f.write("  </div>\n")
                    
                    # Add download links for other formats
                    f.write("  <div class='file-formats'>Download figures as: ")
                    for name, path in generated_files['figures']:
                        if 'Tissue Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                    f.write("  </div>\n")
                else:
                    f.write("  <h3>3.1 Tissue Type Distribution</h3>\n")
                    f.write("  <p class='notes'>No tissue type distribution data available for this region.</p>\n")
                
                # 3.2 Biosample type distribution
                if biosample_dist_figures:
                    f.write("  <h3>3.2 Biosample Type Distribution</h3>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>These charts illustrate the biosample type distribution in the overlapping samples.</p>\n")
                    f.write("    <p>Different biosample types (e.g., cell lines, primary cells, tissues) may have distinctive super enhancer patterns.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display biosample distribution charts
                    f.write("  <div class='chart-container' style='max-width: 90%; margin: 20px auto;'>\n")
                    for figure_name, figure_path in biosample_dist_figures:
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                    f.write("  </div>\n")
                    
                    # Add download links for other formats
                    f.write("  <div class='file-formats'>Download figures as: ")
                    for name, path in generated_files['figures']:
                        if 'Biosample Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                    f.write("  </div>\n")
                else:
                    f.write("  <h3>3.2 Biosample Type Distribution</h3>\n")
                    f.write("  <p class='notes'>No biosample type distribution data available for this region.</p>\n")
                
                # Section 4: Enrichment Analysis
                f.write("  <h2>4. Enrichment Analysis</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>Enrichment analysis identifies tissues that are statistically over-represented or under-represented in this region compared to the background distribution across all samples.</p>\n")
                f.write("    <p>This analysis helps to identify biological contexts where this genomic region might have specialized regulatory functions.</p>\n")
                f.write("  </div>\n")
                
                # 4.1 Tissue Enrichment Analysis
                f.write("  <h3>4.1 Tissue Type Enrichment</h3>\n")
                
                if tissue_enrichment_figures:
                    tissue_enrichment = analysis_results.get('tissue_enrichment')
                    if tissue_enrichment is not None:
                        # Determine which column contains p-values
                        p_col = 'adjusted p-value' if 'adjusted p-value' in tissue_enrichment.columns else 'p-value'
                        significant_tissues = sum(tissue_enrichment[p_col] < pvalue_threshold)
                        
                        f.write("  <div class='explanation'>\n")
                        f.write("    <p>This analysis uses Fisher's exact test to determine which tissue types are significantly enriched or depleted in the region compared to the genome-wide background distribution.</p>\n")
                        f.write("    <p><strong>How to interpret the results:</strong></p>\n")
                        f.write("    <ul>\n")
                        f.write("      <li><span class='enriched'>Enriched tissues (Fold Change > 1)</span> have more samples in this region than expected by chance</li>\n")
                        f.write("      <li><span class='depleted'>Depleted tissues (Fold Change < 1)</span> have fewer samples than expected by chance</li>\n")
                        f.write("      <li>The height of the bars indicates statistical significance (-log10 of p-value)</li>\n")
                        f.write("      <li>The dotted vertical line indicates the significance threshold (p-value = 0.05)</li>\n")
                        f.write("    </ul>\n")
                        f.write("    <p>Significant enrichment suggests that this genomic region may have tissue-specific super enhancer activity, potentially regulating genes important for that tissue's function.</p>\n")
                        f.write("  </div>\n")
                        
                        f.write(f"  <p><strong>Significant tissues:</strong> {significant_tissues} out of {len(tissue_enrichment)} tissue types (p < {pvalue_threshold}).</p>\n")
                        
                        # Top significant tissues table
                        if significant_tissues > 0:
                            sig_tissues = tissue_enrichment[tissue_enrichment[p_col] < pvalue_threshold].sort_values(p_col)
                            f.write("  <h4>Top Significant Tissues</h4>\n")
                            f.write("  <table>\n")
                            f.write("    <tr><th>Tissue Type</th><th>Fold Change</th><th>p-value</th><th>Adjusted p-value</th><th>Enrichment</th></tr>\n")
                            
                            for i, (_, row) in enumerate(sig_tissues.head(5).iterrows()):
                                # Handle both 'adjusted p-value' and FDR column formats
                                adj_p = row['adjusted p-value'] if 'adjusted p-value' in row else row.get(p_col, 'N/A')
                                
                                enrichment_type = "Enriched" if row['Fold Change'] > 1 else "Depleted"
                                enrichment_class = "enriched" if row['Fold Change'] > 1 else "depleted"
                                f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Fold Change']:.2f}x</td>" +
                                    f"<td>{row['p-value']:.2e}</td><td>{adj_p if isinstance(adj_p, str) else adj_p:.2e}</td>" +
                                    f"<td class='{enrichment_class}'>{enrichment_type}</td></tr>\n")
                            
                            f.write("  </table>\n")
                    
                    # Display enrichment visualization
                    for figure_name, figure_path in tissue_enrichment_figures:
                        f.write("  <div class='chart-container' style='max-width: 90%; margin: 20px auto;'>\n")
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("  </div>\n")
                    
                    # Add download links for other formats
                    f.write("  <div class='file-formats'>Download figures as: ")
                    for name, path in generated_files['figures']:
                        if 'Tissue Enrichment' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                    f.write("  </div>\n")
                else:
                    f.write("  <p class='notes'>No tissue enrichment analysis results available for this region.</p>\n")
                
                # Section 5: Read Distribution Analysis
                f.write("  <h2>5. Read Distribution Analysis</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section visualizes how super enhancer regions are distributed across the genomic region. Two complementary visualizations are provided:</p>\n")
                f.write("    <ol>\n")
                f.write("      <li><strong>Simple Read Distribution</strong>: Shows the raw distribution of super enhancers from all samples.</li>\n")
                f.write("      <li><strong>Read Distribution with Peaks</strong>: Identifies areas of high super enhancer activity (peaks).</li>\n")
                f.write("    </ol>\n")
                f.write("    <p>These visualizations help identify sub-regions with concentrated enhancer activity.</p>\n")
                f.write("  </div>\n")
                
                # 5.1 Simple Stacked Reads Plot
                f.write("  <h3>5.1 Simple Read Distribution</h3>\n")
                
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This plot shows super enhancer regions stacked horizontally across the genomic coordinate. Each horizontal line represents an enhancer from a sample that overlaps with this region. The visualization helps to:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li>Identify patterns in super enhancer coverage across the region</li>\n")
                f.write("      <li>Visualize the density of enhancers at different positions</li>\n")
                f.write("      <li>Observe the overall distribution of super enhancer activity</li>\n")
                f.write("    </ul>\n")
                f.write("    <p>Areas with many stacked reads indicate regions of high enhancer activity, suggesting important regulatory elements.</p>\n")
                f.write("  </div>\n")
                
                # Display simple stacked plot
                simple_plot_path = next((path for name, path in stacked_plot_figures if 'Simple Read Distribution' in name), None)
                if simple_plot_path:
                    f.write("  <div class='chart-container' style='max-width: 100%; margin: 20px auto;'>\n")
                    f.write(f"    <img src='{simple_plot_path}' alt='Simple Read Distribution'>\n")
                    f.write("  </div>\n")
                    
                    # Add download links for other formats
                    f.write("  <div class='file-formats'>Download figure as: ")
                    for name, path in generated_files['figures']:
                        if 'Simple Stacked Read Plot' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                    f.write("  </div>\n")
                else:
                    f.write("  <p class='notes'>No simple read distribution plot available for this region.</p>\n")
                
                # 5.2 Stacked Reads with Peaks
                f.write("  <h3>5.2 Read Distribution with Peaks</h3>\n")
                
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This enhanced visualization identifies peaks of super enhancer activity within the region. Peaks represent sub-regions with significantly higher enhancer density. Key features include:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li>Highlighted peak regions in red background</li>\n")
                f.write("      <li>Peak identifiers with genomic positions</li>\n")
                f.write("      <li>Quantification of enhancer density in each peak</li>\n")
                f.write("    </ul>\n")
                f.write("    <p>Peaks often represent the core of regulatory elements with high activity. These regions are prime candidates for further experimental investigation.</p>\n")
                f.write("  </div>\n")
                
                # Display peak plot
                peak_plot_path = next((path for name, path in stacked_plot_figures if 'Read Distribution with Peaks' in name), None)
                
                # Try alternative searches if not found
                if peak_plot_path is None:
                    for name, path in generated_files['figures']:
                        if 'peak_stacked' in path.lower() or 'stacked_peak' in path.lower() or 'Peaks' in name:
                            peak_plot_path = path
                            break
                
                # Check for peak data
                has_peak_data = 'peak_data' in analysis_results and analysis_results['peak_data']
                has_peak_plot = peak_plot_path is not None
                
                if has_peak_plot:
                    f.write("  <div class='chart-container' style='max-width: 100%; margin: 20px auto;'>\n")
                    f.write(f"    <img src='{peak_plot_path}' alt='Read Distribution with Peaks'>\n")
                    f.write("  </div>\n")
                    
                    # Add download links for other formats
                    f.write("  <div class='file-formats'>Download figure as: ")
                    for name, path in generated_files['figures']:
                        if 'Stacked Read Plot with Peaks' in name or 'peak' in name.lower():
                            format_name = name.split('(')[1].split(')')[0] if '(' in name else path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                    f.write("  </div>\n")
                elif not has_peak_plot and has_peak_data:
                    # Peak data exists but plot not found
                    f.write("  <p class='notes'>Peaks were detected, but the visualization is not available in this report.</p>\n")
                else:
                    f.write("  <p class='notes'>No peak detection plot available for this region.</p>\n")
                
                # Peak information if available
                peak_data = analysis_results.get('peak_data')
                if peak_data:
                    f.write("  <h4>Peak Analysis</h4>\n")
                    f.write("  <p>The following table summarizes the peaks of super enhancer activity detected in this region:</p>\n")
                    
                    f.write("  <table>\n")
                    f.write("    <tr><th>Peak ID</th><th>Position</th><th>Width</th><th>Total Reads</th><th>Sample Count</th><th>Sample %</th><th>Density</th></tr>\n")
                    
                    # Sort peaks by position
                    sorted_peaks = sorted(peak_data.items(), key=lambda x: x[1]['position'])
                    
                    for peak_id, peak_info in sorted_peaks:
                        f.write(f"    <tr><td>{peak_id}</td>" +
                            f"<td>{peak_info['position']:,}</td>" +
                            f"<td>{peak_info['width']:,} bp</td>" +
                            f"<td>{peak_info['reads']['total']}</td>" +
                            f"<td>{peak_info['samples']['count']}</td>" +
                            f"<td>{peak_info['samples']['percentage']:.1f}%</td>" +
                            f"<td>{peak_info['density']:.3f} reads/bp</td></tr>\n")
                    
                    f.write("  </table>\n")
                    
                    # Link to sample lists for each peak
                    peak_sample_files = []
                    for peak_id in peak_data.keys():
                        files = [path for name, path in generated_files['tables'] 
                            if peak_id in name and 'samples' in name.lower()]
                        if files:
                            peak_sample_files.append((peak_id, files[0], peak_data[peak_id]['samples']['count']))
                    
                    if peak_sample_files:
                        f.write("  <h4>Peak-Specific Sample Lists</h4>\n")
                        f.write("  <p>Download lists of samples contributing to each peak:</p>\n")
                        f.write("  <ul>\n")
                        for peak_id, file_path, count in peak_sample_files:
                            f.write(f"    <li><a href='{file_path}'>{peak_id} - {count} samples</a></li>\n")
                        f.write("  </ul>\n")
                
                # Section 6: Generated Files
                f.write("  <h2>6. Generated Files</h2>\n")
                
                # List of figures
                if generated_files['figures']:
                    f.write("  <h3>Figures</h3>\n")
                    f.write("  <ul>\n")
                    for name, path in generated_files['figures']:
                        f.write(f"    <li><a href='{path}'>{name}</a></li>\n")
                    f.write("  </ul>\n")
                
                # List of tables
                if generated_files['tables']:
                    f.write("  <h3>Tables</h3>\n")
                    f.write("  <ul>\n")
                    for name, path in generated_files['tables']:
                        f.write(f"    <li><a href='{path}'>{name}</a></li>\n")
                    f.write("  </ul>\n")
                
                # Section 7: Summary and Interpretation
                f.write("  <h2>7. Summary and Interpretation</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write("    <h3>Key Findings</h3>\n")
                f.write("    <ul>\n")
                
                # Sample overlap summary
                f.write(f"      <li><strong>Sample Coverage:</strong> {sample_stats['overlapping_samples']:,} samples have super enhancers in this region</li>\n")
                
                # Tissue findings
                tissue_enrichment = analysis_results.get('tissue_enrichment')
                if tissue_enrichment is not None:
                    p_col = 'adjusted p-value' if 'adjusted p-value' in tissue_enrichment.columns else 'p-value'
                    sig_tissues = sum(tissue_enrichment[p_col] < pvalue_threshold)
                    if sig_tissues > 0:
                        f.write(f"      <li><strong>Tissue Specificity:</strong> {sig_tissues} tissue types are significantly enriched or depleted (p < {pvalue_threshold})</li>\n")
                        
                        # Add top enriched tissue
                        top_enriched = tissue_enrichment[(tissue_enrichment['Fold Change'] > fc_threshold) & 
                                                (tissue_enrichment[p_col] < pvalue_threshold)].sort_values(p_col)
                        if not top_enriched.empty:
                            tissue = top_enriched.iloc[0]['Tissue Type']
                            fold_change = top_enriched.iloc[0]['Fold Change']
                            p_value = top_enriched.iloc[0][p_col]
                            f.write(f"      <li><strong>Most Enriched Tissue:</strong> <span class='enriched'>{tissue}</span> ({fold_change:.2f}x enrichment, p = {p_value:.2e})</li>\n")
                        
                        # Add top depleted tissue
                        top_depleted = tissue_enrichment[(tissue_enrichment['Fold Change'] < 1/fc_threshold) & 
                                                (tissue_enrichment[p_col] < pvalue_threshold)].sort_values(p_col)
                        if not top_depleted.empty:
                            tissue = top_depleted.iloc[0]['Tissue Type']
                            fold_change = top_depleted.iloc[0]['Fold Change']
                            p_value = top_depleted.iloc[0][p_col]
                            f.write(f"      <li><strong>Most Depleted Tissue:</strong> <span class='depleted'>{tissue}</span> ({fold_change:.2f}x depletion, p = {p_value:.2e})</li>\n")
                
                # Peak findings
                peak_data = analysis_results.get('peak_data')
                if peak_data:
                    f.write(f"      <li><strong>SE Activity Peaks:</strong> {len(peak_data)} distinct peaks of super enhancer activity identified</li>\n")
                    
                    # Highlight highest density peak
                    peak_densities = [(peak_id, peak_info['density'], peak_info['position']) for peak_id, peak_info in peak_data.items()]
                    if peak_densities:
                        highest_peak_id, highest_density, peak_pos = max(peak_densities, key=lambda x: x[1])
                        f.write(f"      <li><strong>Highest Activity Region:</strong> {highest_peak_id} at position {peak_pos:,} with {highest_density:.3f} reads/bp</li>\n")
                
                f.write("    </ul>\n")
                
                # Interpretation
                f.write("    <h3>Biological Interpretation</h3>\n")
                f.write("    <p>Based on the analysis results, this genomic region appears to")
                
                # Add interpretations based on results
                interpretations = []
                
                # Tissue specificity interpretation
                if tissue_enrichment is not None:
                    p_col = 'adjusted p-value' if 'adjusted p-value' in tissue_enrichment.columns else 'p-value'
                    significant_enriched = tissue_enrichment[(tissue_enrichment['Fold Change'] > fc_threshold) & 
                                                    (tissue_enrichment[p_col] < pvalue_threshold)]
                    if not significant_enriched.empty:
                        tissue = significant_enriched.iloc[0]['Tissue Type']
                        interpretations.append(f"have super enhancer activity that is particularly strong in <strong>{tissue}</strong> tissue")
                
                # Peak distribution interpretation
                if peak_data and len(peak_data) > 0:
                    if len(peak_data) == 1:
                        interpretations.append("contain a single well-defined super enhancer element")
                    elif len(peak_data) <= 3:
                        interpretations.append(f"contain {len(peak_data)} distinct super enhancer elements that may function independently or cooperatively")
                    else:
                        interpretations.append(f"represent a complex regulatory region with {len(peak_data)} distinct areas of super enhancer activity")
                
                # Default interpretation if no specific findings
                if not interpretations:
                    interpretations.append("show super enhancer activity across multiple tissues and cell types")
                
                # Join interpretations with proper conjunctions
                if len(interpretations) == 1:
                    f.write(f" {interpretations[0]}.")
                elif len(interpretations) == 2:
                    f.write(f" {interpretations[0]} and {interpretations[1]}.")
                else:
                    last = interpretations.pop()
                    f.write(f" {', '.join(interpretations)}, and {last}.")
                
                f.write(" The presence of super enhancers in this region suggests it may play an important role in gene regulation.</p>\n")
                f.write("  </div>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by SEdbRegionAnalyzer</p>\n")
                
                # End HTML document
                f.write("</body>\n")
                f.write("</html>\n")
            
            self.logger.info(f"Region analysis report generated: {html_path}")
            return html_path
            
        except Exception as e:
            self.logger.error(f"Error generating region analysis report: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None





    def analyze_and_report_region(self, chrom, start, end, region_name=None, output_dir=None, 
                                figure_formats=None, save_tables=True, peak_threshold=0.2,
                                pvalue_threshold=0.05, fdr_threshold=0.05, fc_threshold=1.5):
        """
        One-step function to analyze a genomic region and generate a comprehensive report.
        
        This is a convenience wrapper function that performs both the analysis and reporting.
        
        Parameters:
        -----------
        chrom : str
            Chromosome name (e.g., 'chr7')
        start : int
            Start position
        end : int
            End position
        region_name : str, optional
            Region name (auto-generated if None)
        output_dir : str, optional
            Output directory (default: 'results/[region_name]')
        figure_formats : list, optional
            List of figure formats to save (default: ['png', 'svg', 'pdf'])
        save_tables : bool, optional 
            Whether to save tables as CSV files (default: True)
        peak_threshold : float, optional
            Threshold for peak detection (default: 0.2)
        pvalue_threshold : float, optional
            P-value threshold for enrichment significance (default: 0.05)
        fdr_threshold : float, optional
            FDR threshold for enrichment significance (default: 0.05)
        fc_threshold : float, optional
            Fold change threshold for enrichment display (default: 1.5)
            
        Returns:
        --------
        dict
            Dictionary containing:
            - report_path: Path to the generated HTML report
            - analysis_results: Complete analysis results
            - has_overlaps: Boolean indicating if overlapping regions were found
            - overlap_sample_count: Number of samples overlapping with the region
            - overlap_read_count: Number of reads overlapping with the region
            - significant_tissues: List of statistically significant enriched tissues
                [{'tissue': tissue_name, 'fold_change': fc, 'pvalue': p, 'fdr': q}, ...]
        """
        # Set default region name if not provided
        if region_name is None:
            region_name = f"{chrom}_{start}_{end}"
        
        # Set default output directory if not provided
        if output_dir is None:
            output_dir = f"results/{region_name}"
        
        self.logger.info(f"Starting analysis and report generation for region {chrom}:{start:,}-{end:,}")
        
        # Step 1: Run comprehensive analysis
        analysis_results = self.analyze_region_comprehensive(
            chrom=chrom,
            start=start,
            end=end,
            region_name=region_name,
            output_dir=output_dir,
            peak_threshold=peak_threshold
        )
        
        # Step 2: Extract key information from analysis results
        
        # Check if overlapping regions were found
        has_overlaps = analysis_results.get('overlapping_samples_count', 0) > 0
        
        # Get sample count
        overlap_sample_count = analysis_results.get('overlapping_samples_count', 0)
        
        # Get read count (number of overlapping records)
        overlap_read_count = 0
        if 'overlapping_results' in analysis_results and analysis_results['overlapping_results'] is not None:
            overlap_read_count = len(analysis_results['overlapping_results'])
        
        # Extract significant tissues (positive enrichment, FDR significant)
        significant_tissues = []
        if 'tissue_enrichment' in analysis_results and analysis_results['tissue_enrichment'] is not None:
            tissue_df = analysis_results['tissue_enrichment']
            # Check if the adjusted p-value column exists, otherwise use p-value
            fdr_col = 'adjusted p-value' if 'adjusted p-value' in tissue_df.columns else 'p-value'
            
            # Filter for significant positively enriched tissues
            sig_tissues = tissue_df[(tissue_df[fdr_col] < fdr_threshold) & 
                                (tissue_df['Fold Change'] > fc_threshold)]
            
            # Sort by significance
            sig_tissues = sig_tissues.sort_values(fdr_col)
            
            # Convert to list of dictionaries
            for _, row in sig_tissues.iterrows():
                tissue_info = {
                    'tissue': row['Tissue Type'],
                    'fold_change': row['Fold Change'],
                    'pvalue': row['p-value'],
                    'fdr': row[fdr_col] if fdr_col in row else row['p-value']
                }
                significant_tissues.append(tissue_info)
        
        # Step 3: Generate report from analysis results
        report_path = self.generate_region_report(
            analysis_results=analysis_results,
            output_dir=output_dir,
            figure_formats=figure_formats,
            save_tables=save_tables,
            pvalue_threshold=pvalue_threshold,
            fc_threshold=fc_threshold
        )
        
        # Return comprehensive result dictionary
        return {
            'report_path': report_path,
            'analysis_results': analysis_results,
            'has_overlaps': has_overlaps,
            'overlap_sample_count': overlap_sample_count,
            'overlap_read_count': overlap_read_count,
            'significant_tissues': significant_tissues
        }

    def plot_simple_pie(self, data: pd.DataFrame, category: str, count_col: str = 'Count',
                    title: Optional[str] = None, limit: int = 10,
                    figsize: tuple = (12, 8),
                    colormap: str = 'tab20',
                    title_y: float = 1.1,  # タイトル位置を上に調整
                    save_dir: Optional[str] = None, 
                    save_filename: Optional[str] = None,
                    save_svg: bool = False, 
                    save_png: bool = False, 
                    save_pdf: bool = False, 
                    save_eps: bool = False,
                    save_dpi: int = 600, 
                    save_transparent: bool = False) -> plt.Figure:
        """
        Create a simple, clear pie chart
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame containing category and count columns
        category : str
            Category column name to display in pie chart
        count_col : str, default='Count'
            Column name for the count values
        title : str, optional
            Chart title
        limit : int, default=10
            Number of top categories to display
        figsize : tuple, default=(12, 8)
            Figure size
        colormap : str, default='tab20'
            Matplotlib colormap to use ('tab20', 'tab20b', 'tab20c', 'Set3', etc.)
        title_y : float, default=1.1
            Y position of the title (higher value moves title up)
        save_dir : str, optional
            Directory to save output files. If None, current directory is used
        save_filename : str, optional
            Base filename for saved files (without extension). If None, auto-generated
        save_svg : bool, optional
            Whether to save as SVG format (default: False)
        save_png : bool, optional
            Whether to save as PNG format (default: False)
        save_pdf : bool, optional
            Whether to save as PDF format (default: False)
        save_eps : bool, optional
            Whether to save as EPS format (default: False)
        save_dpi : int, optional
            DPI (resolution) for raster formats like PNG (default: 600)
        save_transparent : bool, optional
            Whether to save with transparent background (default: False)
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        # Check if required columns exist
        if category not in data.columns:
            self.logger.warning(f"Column '{category}' does not exist in the dataframe")
            return None
        
        if count_col not in data.columns:
            self.logger.warning(f"Column '{count_col}' does not exist in the dataframe")
            return None
        
        # Sort data by count (descending)
        sorted_data = data.sort_values(count_col, ascending=False)
        total_categories = len(sorted_data)
        
        # Set title
        if title is None:
            title = f'{category} Distribution'
            
        # Divide into top N and others
        if total_categories > limit:
            top_data = sorted_data.head(limit)
            others_count = sorted_data.iloc[limit:][count_col].sum()
            
            # Create a DataFrame with 'Others' row
            others_row = pd.DataFrame({category: ['Others'], count_col: [others_count]})
            plot_data = pd.concat([top_data, others_row], ignore_index=True)
            title_suffix = f" (Top {limit}, other {total_categories-limit} types grouped as 'Others')"
        else:
            plot_data = sorted_data
            title_suffix = ""
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # Get labels and values
        labels = plot_data[category]
        values = plot_data[count_col]
        
        # Generate distinct colors for pie segments
        cmap = plt.cm.get_cmap(colormap)
        num_colors = len(plot_data)
        
        # 重要な修正: Othersカテゴリと他のカテゴリの色が被らないようにする
        colors = []
        others_idx = None
        
        # まずOthersのインデックスを見つける（存在する場合）
        for i, label in enumerate(labels):
            if label == 'Others':
                others_idx = i
                break
        
        # カラーマップから色を生成
        for i in range(num_colors):
            if i == others_idx:
                # Othersには特別な色（灰色）を使用
                colors.append('grey')
            else:
                # カラーマップからインデックスを計算（Othersがある場合はずらす）
                if others_idx is not None and i > others_idx:
                    # Othersの後のインデックスは-1して詰める
                    colors.append(cmap((i-1) % cmap.N))
                else:
                    colors.append(cmap(i % cmap.N))
        
        # Create pie chart (maintaining original design)
        plt.pie(
            values, 
            labels=labels,
            autopct='%1.1f%%',
            startangle=90,
            shadow=False,
            wedgeprops={'edgecolor': 'none'},  # Remove edge lines
            textprops={'fontsize': 10},
            colors=colors  # 修正した色のリストを使用
        )
        
        # タイトル位置を上に調整し、余白を確保
        plt.title(f'{title}{title_suffix}', y=title_y)
        plt.axis('equal')
        
        # タイトル用に十分なスペースを確保
        plt.subplots_adjust(top=0.85)  # これによりタイトル用のスペースを確保
        plt.tight_layout(pad=1.5)  # パディングを少し増やす
        
        # Save figure in requested formats
        if save_svg or save_png or save_pdf or save_eps:
            import os
            
            # Create directory if needed
            if save_dir is None:
                save_dir = os.getcwd()  # Current directory
            elif not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
                self.logger.info(f"Created output directory: {save_dir}")
            
            # Generate base filename if not provided
            if save_filename is None:
                save_filename = f"pie_chart_{category.lower().replace(' ', '_')}"
            
            # Create full path without extension
            base_path = os.path.join(save_dir, save_filename)
            
            # Save in each requested format
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved pie chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        # Display detailed information
        total_count = values.sum()
        print(f"\n{category} distribution:")
        print(f"Total {total_categories} types, showing top {min(limit, total_categories)}")
        
        for i, (label, count) in enumerate(zip(labels, values)):
            if label != 'Others':
                print(f"  {label}: {int(count)} samples ({count/total_count*100:.1f}%)")
        
        if 'Others' in labels.values:
            others_idx = labels.tolist().index('Others')
            others_count = values.iloc[others_idx]
            print(f"  Others ({total_categories - limit} types): {int(others_count)} samples ({others_count/total_count*100:.1f}%)")
        
        return fig


    # 更新版 plot_sorted_bar 関数
    def plot_sorted_bar(self, data: pd.DataFrame, category: str, count_col: str = 'Count', 
                    title: Optional[str] = None, limit: int = 15,
                    horizontal: bool = True, color: str = 'skyblue', 
                    figsize: tuple = (12, 8),
                    save_dir: Optional[str] = None, 
                    save_filename: Optional[str] = None,
                    save_svg: bool = False, 
                    save_png: bool = False, 
                    save_pdf: bool = False, 
                    save_eps: bool = False,
                    save_dpi: int = 600, 
                    save_transparent: bool = False) -> plt.Figure:
        """
        Helper method to plot sorted bar chart for distribution visualization.
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame containing category and count columns
        category : str
            Column name for the category
        count_col : str, optional
            Column name for the count values (default: 'Count')
        title : str, optional
            Plot title
        limit : int, optional
            Number of top categories to display
        horizontal : bool, optional
            Whether to use horizontal bar chart (default: True)
        color : str, optional
            Bar color
        figsize : tuple, optional
            Figure size
        save_dir : str, optional
            Directory to save the figure
        save_filename : str, optional
            Base filename for saved figure (without extension)
        save_svg, save_png, save_pdf, save_eps : bool, optional
            Whether to save in respective formats
        save_dpi : int, optional
            DPI for saved images
        save_transparent : bool, optional
            Whether to save with transparent background
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        # Sort data and limit to top N categories
        if data is None or len(data) == 0:
            self.logger.warning(f"No data provided for plotting {category} distribution")
            return None
        
        sorted_data = data.sort_values(count_col, ascending=False)
        if limit > 0 and len(sorted_data) > limit:
            plot_data = sorted_data.head(limit)
        else:
            plot_data = sorted_data
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        if horizontal:
            # Horizontal bar chart (better for many categories with long names)
            # Reverse order for better display (highest at top)
            plot_data = plot_data.iloc[::-1].reset_index(drop=True)
            bars = ax.barh(plot_data[category], plot_data[count_col], color=color, alpha=0.8)
            
            # Add count values to the right of bars
            for i, bar in enumerate(bars):
                ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                    f"{int(bar.get_width())}", va='center')
            
            ax.set_xlabel('Count')
            ax.set_ylabel(category)
        else:
            # Vertical bar chart
            bars = ax.bar(plot_data[category], plot_data[count_col], color=color, alpha=0.8)
            
            # Add count values above bars
            for i, bar in enumerate(bars):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f"{int(bar.get_height())}", ha='center', va='bottom')
            
            ax.set_xlabel(category)
            ax.set_ylabel('Count')
            plt.xticks(rotation=45, ha='right')
        
        # Set title
        if title:
            ax.set_title(title)
        else:
            ax.set_title(f'{category} Distribution')
        
        plt.tight_layout()
        
        # Save figure if requested
        if save_dir and save_filename:
            import os
            os.makedirs(save_dir, exist_ok=True)
            base_path = os.path.join(save_dir, save_filename)
            
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved bar chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return fig



    def generate_sedb_report(self, output_dir='results/sedb_report',
                        top_n=20, figure_formats=None, save_tables=True,
                        dpi=300, fig_width=12, fig_height=8):
        """
        Generate a comprehensive HTML report of the SEdb metadata analysis.
        
        Parameters:
        -----------
        output_dir : str, optional
            Base directory for all output files (default: 'results/sedb_report')
        top_n : int, optional
            Number of top categories to display in visualizations (default: 20)
        figure_formats : list, optional
            List of figure formats to save (default: ['png', 'svg', 'pdf'])
        save_tables : bool, optional
            Whether to save distribution tables as CSV files (default: True)
        dpi : int, optional
            Resolution for raster images like PNG (default: 300)
        fig_width : int, optional
            Width of the figure in inches (default: 12)
        fig_height : int, optional
            Height of the figure in inches (default: 8)
            
        Returns:
        --------
        str
            Path to the generated HTML report
        """
        # Check if the necessary data is loaded
        if self.bed_df is None or self.sample_info is None:
            self.logger.warning("SE BED data or sample information not loaded. Call load_databases() first.")
            return None
        
        # Default figure formats if not provided
        if figure_formats is None:
            figure_formats = ['png', 'svg', 'pdf']
        
        # Create timestamp for filenames
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create directory structure with relative paths
        figures_dir = 'figures'  # Use relative paths for portable HTML
        tables_dir = 'tables'
        
        # Create full directory paths
        figures_full_path = os.path.join(output_dir, figures_dir)
        tables_full_path = os.path.join(output_dir, tables_dir)
        
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(figures_full_path, exist_ok=True)
        os.makedirs(tables_full_path, exist_ok=True)
        
        self.logger.info(f"Creating SEdb metadata report in directory: {output_dir}")
        
        # Dictionary to track all generated files
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Collect basic database stats
        db_stats = {
            'total_samples': len(self.sample_info) if self.sample_info is not None else 0,
            'total_se_records': len(self.bed_df) if self.bed_df is not None else 0,
            'bed_columns': list(self.bed_df.columns) if self.bed_df is not None else [],
            'sample_info_columns': list(self.sample_info.columns) if self.sample_info is not None else [],
            'timestamp': timestamp
        }
        
        # Calculate missing values in sample_info
        if self.sample_info is not None:
            missing_values = self.sample_info.isnull().sum()
            missing_values = missing_values[missing_values > 0]
            db_stats['missing_values'] = missing_values.to_dict()
        
        # Set Matplotlib params for better quality figures
        plt.rcParams['figure.figsize'] = (fig_width, fig_height)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # Step 1: Tissue Distribution Analysis
        self.logger.info("Calculating and visualizing tissue type distribution...")
        
        tissue_dist_figures = []
        tissue_csv_path = None
        
        try:
            # 1. 直接保存オプションを使用してtissue distributionを生成
            tissue_viz_basename = f'tissue_distribution_direct_{timestamp}'
            
            # 重要な修正：直接ファイル保存オプションを使用
            self.visualize_tissue_distribution(
                top_n=top_n, 
                horizontal=True, 
                store_all=True,
                figsize=(fig_width, fig_height),
                save_dir=figures_full_path,
                save_filename=tissue_viz_basename,
                save_png=True,
                save_svg=('svg' in figure_formats),
                save_pdf=('pdf' in figure_formats),
                save_eps=('eps' in figure_formats),
                save_dpi=dpi
            )
            
            # 保存されたファイルを登録
            for fmt in figure_formats:
                if fmt == 'png' or f'save_{fmt}' in locals():
                    figure_filename = f'{tissue_viz_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution (Direct) ({fmt.upper()})', figure_relative_path))
                        
                        if fmt == 'png':
                            tissue_dist_figures.append(('Tissue Distribution (Direct)', figure_relative_path))
            
            # Now create supplementary visualizations using the stored distribution data
            if hasattr(self, '_tissue_distribution') and self._tissue_distribution is not None:
                # Save tissue data as CSV
                if save_tables:
                    tissue_table_filename = f'tissue_distribution_{timestamp}.csv'
                    tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
                    self._tissue_distribution.to_csv(tissue_table_path, index=False)
                    generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
                    tissue_csv_path = os.path.join(tables_dir, tissue_table_filename)
                
                # 1. Bar chart with sorted_bar helper (directly save)
                tissue_bar_basename = f'tissue_distribution_bar_{timestamp}'
                self.plot_sorted_bar(
                    data=self._tissue_distribution,
                    category='Tissue Type',
                    count_col='Sample Count',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='skyblue',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=tissue_bar_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                
                # Register the figure
                for fmt in figure_formats:
                    if fmt == 'png' or f'save_{fmt}' in locals():
                        figure_filename = f'{tissue_bar_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_dist_figures.append(('Tissue Distribution Bar Chart', figure_relative_path))
                
                # 2. Pie chart with simple_pie helper (directly save)
                tissue_pie_basename = f'tissue_distribution_pie_{timestamp}'
                self.plot_simple_pie(
                    data=self._tissue_distribution,
                    category='Tissue Type',
                    count_col='Sample Count',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    figsize=(fig_width, fig_height),
                    colormap='tab20',
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                
                # Register pie chart
                for fmt in figure_formats:
                    if fmt == 'png' or f'save_{fmt}' in locals():
                        figure_filename = f'{tissue_pie_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_dist_figures.append(('Tissue Distribution Pie Chart', figure_relative_path))
                
                print(f"Successfully created tissue distribution visualizations")
            else:
                self.logger.warning("No tissue distribution data stored after visualization.")
                
        except Exception as e:
            self.logger.error(f"Error in tissue distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Step 2: Biosample Distribution Analysis
        self.logger.info("Calculating and visualizing biosample type distribution...")
        
        biosample_dist_figures = []
        biosample_csv_path = None
        
        try:
            # visualize_biosample_distributionメソッドが存在しないため、直接実装
            # バイオサンプルタイプの分布を計算
            if 'Biosample type' in self.sample_info.columns:
                biosample_counts = self.sample_info['Biosample type'].value_counts()
                total_samples = len(self.sample_info)
                
                # 計算結果をDataFrameに変換
                biosample_data = pd.DataFrame({
                    'Biosample Type': biosample_counts.index,
                    'Sample Count': biosample_counts.values,
                    'Percentage (%)': (biosample_counts / total_samples * 100).round(2)
                })
                
                # CSV出力
                if save_tables:
                    biosample_table_filename = f'biosample_distribution_{timestamp}.csv'
                    biosample_table_path = os.path.join(tables_full_path, biosample_table_filename)
                    biosample_data.to_csv(biosample_table_path, index=False)
                    generated_files['tables'].append(('Biosample Distribution', os.path.join(tables_dir, biosample_table_filename)))
                    biosample_csv_path = os.path.join(tables_dir, biosample_table_filename)
                
                # バイオサンプルタイプの分布を直接可視化
                plt.figure(figsize=(fig_width, fig_height))
                
                # トップNを表示用に選択
                if top_n > 0 and len(biosample_data) > top_n:
                    display_data = biosample_data.head(top_n)
                else:
                    display_data = biosample_data
                
                # 水平バーチャート
                plt.barh(display_data['Biosample Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='lightgreen')
                plt.xlabel('Sample Count')
                plt.ylabel('Biosample Type')
                plt.title(f'Biosample Type Distribution (Top {len(display_data)} of {len(biosample_data)})')
                plt.tight_layout()
                
                # 分布図を保存
                biosample_viz_basename = f'biosample_distribution_direct_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{biosample_viz_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # 相対パスで登録
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Biosample Distribution (Direct) ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        biosample_dist_figures.append(('Biosample Distribution (Direct)', figure_relative_path))
                
                plt.close()
                
                # 追加の可視化 - バーチャート
                biosample_bar_basename = f'biosample_distribution_bar_{timestamp}'
                self.plot_sorted_bar(
                    data=biosample_data,
                    category='Biosample Type',
                    count_col='Sample Count',
                    title=f'Biosample Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='lightgreen',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=biosample_bar_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                
                # Register the figure
                for fmt in figure_formats:
                    figure_filename = f'{biosample_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Biosample Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            biosample_dist_figures.append(('Biosample Distribution Bar Chart', figure_relative_path))
                
                # 追加の可視化 - 円グラフ
                biosample_pie_basename = f'biosample_distribution_pie_{timestamp}'
                self.plot_simple_pie(
                    data=biosample_data,
                    category='Biosample Type',
                    count_col='Sample Count',
                    title=f'Biosample Type Distribution (Top {top_n})',
                    limit=top_n,
                    figsize=(fig_width, fig_height),
                    colormap='Set3',
                    save_dir=figures_full_path,
                    save_filename=biosample_pie_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                
                # Register pie chart
                for fmt in figure_formats:
                    figure_filename = f'{biosample_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Biosample Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            biosample_dist_figures.append(('Biosample Distribution Pie Chart', figure_relative_path))
                
                print(f"Successfully created biosample distribution visualizations")
            else:
                self.logger.warning("'Biosample type' column not found in sample information")
        except Exception as e:
            self.logger.error(f"Error in biosample distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Step 3: Tissue-Biosample Relationship Analysis
        self.logger.info("Analyzing relationship between tissue types and biosample types...")
        
        relationship_figures = []
        relationship_csv_path = None
        
        try:
            # 手動でクロスタブとヒートマップを作成
            if 'Tissue type' in self.sample_info.columns and 'Biosample type' in self.sample_info.columns:
                # クロスタブを作成
                cross_tab = pd.crosstab(
                    self.sample_info['Tissue type'], 
                    self.sample_info['Biosample type']
                )
                
                # クロスタブデータを保存
                if save_tables:
                    cross_tab_filename = f'tissue_biosample_crosstab_{timestamp}.csv'
                    cross_tab_path = os.path.join(tables_full_path, cross_tab_filename)
                    cross_tab.to_csv(cross_tab_path)
                    generated_files['tables'].append(('Tissue-Biosample Cross-tabulation', os.path.join(tables_dir, cross_tab_filename)))
                    relationship_csv_path = os.path.join(tables_dir, cross_tab_filename)
                
                # 上位の組織とバイオサンプルタイプを取得
                row_sums = cross_tab.sum(axis=1)
                col_sums = cross_tab.sum(axis=0)
                
                row_order = row_sums.sort_values(ascending=False).index[:min(top_n, len(row_sums))]
                col_order = col_sums.sort_values(ascending=False).index[:min(top_n, len(col_sums))]
                
                # 表示用にデータをフィルタリング
                filtered_data = cross_tab.loc[row_order, col_order]
                
                # ヒートマップ図作成のためseabornをインポート
                import seaborn as sns
                
                # 新しい図を作成
                plt.figure(figsize=(fig_width, fig_height * 1.2))
                
                # ヒートマップ作成 - mask パラメータを使用して0の値をマスク
                cmap = sns.diverging_palette(220, 10, as_cmap=True)
                mask = filtered_data == 0  # 値が0のセルをマスク
                
                # ヒートマップを描画
                heatmap = sns.heatmap(
                    filtered_data, 
                    annot=True,  # セルに値を表示
                    fmt='d',     # 整数フォーマット
                    cmap=cmap,   # カラーマップ
                    linewidths=0.5, 
                    mask=mask,   # マスク適用
                    cbar_kws={'label': 'Sample Count'}
                )
                
                plt.title('Relationship Between Tissue Type and Biosample Type')
                plt.xlabel('Biosample Type')
                plt.ylabel('Tissue Type')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()  # レイアウトを調整
                
                # ヒートマップを保存
                heatmap_basename = f'tissue_biosample_heatmap_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{heatmap_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # 相対パスで登録
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Tissue-Biosample Heatmap ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        relationship_figures.append(('Tissue-Biosample Heatmap', figure_relative_path))
                
                plt.close()
                print(f"Successfully created tissue-biosample heatmap visualization")
            else:
                self.logger.warning("'Tissue type' or 'Biosample type' column not found in sample information")
        except Exception as e:
            self.logger.error(f"Error in tissue-biosample relationship analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Step 4: Cell ID Occurrence Analysis (specific to SEdb)
        self.logger.info("Analyzing cell_id occurrences...")
        
        cell_id_figures = []
        cell_id_csv_path = None
        
        try:
            # Count cell_id occurrences if not already done
            if not hasattr(self, '_cell_id_counts') or self._cell_id_counts is None:
                self.count_cell_id_occurrences()
            
            if hasattr(self, '_cell_id_counts') and self._cell_id_counts is not None:
                # Save cell_id counts data
                if save_tables:
                    cellid_table_filename = f'cell_id_counts_{timestamp}.csv'
                    cellid_table_path = os.path.join(tables_full_path, cellid_table_filename)
                    self._cell_id_counts.to_csv(cellid_table_path, index=False)
                    generated_files['tables'].append(('Cell ID Counts', os.path.join(tables_dir, cellid_table_filename)))
                    cell_id_csv_path = os.path.join(tables_dir, cellid_table_filename)
                
                # Generate histogram of cell_id occurrences
                plt.figure(figsize=(fig_width, fig_height))
                
                # Create histogram
                plt.hist(self._cell_id_counts['count'], bins=30, color='skyblue', edgecolor='black')
                plt.title('Distribution of cell_id Occurrences')
                plt.xlabel('Super Enhancers per Sample')
                plt.ylabel('Number of Samples')
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
                
                # Save histogram
                cellid_hist_basename = f'cell_id_histogram_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{cellid_hist_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # Register with relative path
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Cell ID Occurrence Histogram ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        cell_id_figures.append(('Cell ID Occurrence Histogram', figure_relative_path))
                
                plt.close()
                
                # Generate summary statistics for display in report
                count_stats = self._cell_id_counts['count'].describe().to_dict()
                db_stats['cell_id_stats'] = count_stats
                
                # Display counts per sample using helper method
                cellid_bars_basename = f'cell_id_top_counts_{timestamp}'
                
                # Create a DataFrame of top samples by SE count
                top_samples = self._cell_id_counts.sort_values('count', ascending=False).head(top_n)
                
                # Create bar chart of top samples
                self.plot_sorted_bar(
                    data=top_samples,
                    category='cell_id',
                    count_col='count',
                    title=f'Top {top_n} Samples by SE Count',
                    limit=top_n,
                    horizontal=True,
                    color='lightcoral',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=cellid_bars_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                
                # Register the figure
                for fmt in figure_formats:
                    figure_filename = f'{cellid_bars_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Top Samples by SE Count ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_id_figures.append(('Top Samples by SE Count', figure_relative_path))
                
                print(f"Successfully created cell_id occurrence visualizations")
        except Exception as e:
            self.logger.error(f"Error in cell_id occurrence analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Generate HTML Report
        self.logger.info("Generating HTML report...")
        
        html_filename = f"sedb_metadata_report_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        try:
            with open(html_path, 'w', encoding='utf-8') as f:
                # Start HTML document
                f.write("<!DOCTYPE html>\n")
                f.write("<html>\n")
                f.write("<head>\n")
                f.write("  <title>SEdb Metadata Report</title>\n")
                f.write("  <style>\n")
                f.write("    body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }\n")
                f.write("    h1 { color: #2c3e50; }\n")
                f.write("    h2 { color: #3498db; margin-top: 30px; border-bottom: 2px solid #3498db; padding-bottom: 5px; }\n")
                f.write("    h3 { color: #2980b9; margin-top: 20px; }\n")
                f.write("    h4 { color: #27ae60; margin-top: 15px; }\n")
                f.write("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
                f.write("    th { background-color: #3498db; color: white; text-align: left; padding: 8px; }\n")
                f.write("    td { border: 1px solid #ddd; padding: 8px; }\n")
                f.write("    tr:nth-child(even) { background-color: #f2f2f2; }\n")
                f.write("    img { max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; box-shadow: 0 0 5px rgba(0,0,0,0.1); }\n")
                f.write("    .file-formats { margin: 10px 0; }\n")
                f.write("    .file-formats a { margin-right: 10px; color: #3498db; text-decoration: none; }\n")
                f.write("    .file-formats a:hover { text-decoration: underline; }\n")
                f.write("    .chart-section { display: flex; flex-wrap: wrap; justify-content: space-around; }\n")
                f.write("    .chart-container { margin: 15px; max-width: 48%; }\n")
                f.write("    .highlight { background-color: #f8f9fa; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }\n")
                f.write("    .explanation { background-color: #e8f4f8; padding: 15px; border-radius: 5px; margin: 10px 0; }\n")
                f.write("    .data-link { background-color: #f0f7ff; padding: 8px; margin: 5px 0; border-radius: 3px; }\n")
                f.write("    .data-link a { color: #0066cc; text-decoration: none; }\n")
                f.write("    .data-link a:hover { text-decoration: underline; }\n")
                f.write("  </style>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                
                # Title and timestamp
                f.write(f"  <h1>SEdb Metadata Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Database overview
                f.write("  <h2>1. Database Overview</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This section provides an overview of the SEdb database.</p>\n")
                f.write("  </div>\n")
                
                f.write("  <table>\n")
                f.write("    <tr><th>Metric</th><th>Value</th></tr>\n")
                f.write(f"    <tr><td>Total samples</td><td>{db_stats['total_samples']:,}</td></tr>\n")
                f.write(f"    <tr><td>Total SE records</td><td>{db_stats['total_se_records']:,}</td></tr>\n")
                
                # Add cell_id stats if available
                if 'cell_id_stats' in db_stats:
                    stats = db_stats['cell_id_stats']
                    f.write(f"    <tr><td>Average SE per sample</td><td>{stats['mean']:.2f}</td></tr>\n")
                    f.write(f"    <tr><td>Maximum SE per sample</td><td>{int(stats['max'])}</td></tr>\n")
                    f.write(f"    <tr><td>Minimum SE per sample</td><td>{int(stats['min'])}</td></tr>\n")
                
                f.write("  </table>\n")
                
                # Data schema section
                f.write("  <h3>Data Schema</h3>\n")
                
                # BED file columns
                f.write("  <h4>BED File Columns</h4>\n")
                f.write("  <ul>\n")
                for col in db_stats['bed_columns']:
                    f.write(f"    <li>{col}</li>\n")
                f.write("  </ul>\n")
                
                # Sample info columns
                f.write("  <h4>Sample Info Columns</h4>\n")
                f.write("  <ul>\n")
                for col in db_stats['sample_info_columns']:
                    f.write(f"    <li>{col}</li>\n")
                f.write("  </ul>\n")
                
                # Missing values section
                if 'missing_values' in db_stats and db_stats['missing_values']:
                    f.write("  <h3>Missing Values</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Column</th><th>Missing Count</th><th>Percentage</th></tr>\n")
                    
                    for col, count in db_stats['missing_values'].items():
                        percentage = (count / db_stats['total_samples']) * 100
                        f.write(f"    <tr><td>{col}</td><td>{count}</td><td>{percentage:.2f}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                
                # Tissue Distribution section
                f.write("  <h2>2. Tissue Type Distribution</h2>\n")
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This section shows the distribution of tissue types among samples in the SEdb database.</p>\n")
                f.write("    <p>The visualizations below help identify which tissues are most represented in the database.</p>\n")
                f.write("  </div>\n")
                
                # Direct visualization from visualize_tissue_distribution
                tissue_direct_png = [path for name, path in tissue_dist_figures if 'Direct' in name]
                if tissue_direct_png:
                    f.write("  <h3>Tissue Type Distribution (Direct Visualization)</h3>\n")
                    f.write(f"  <img src='{tissue_direct_png[0]}' alt='Tissue Distribution Direct Visualization'>\n")
                    
                    # Add data link for CSV if available
                    if tissue_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{tissue_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    tissue_direct_formats = [tup for tup in generated_files['figures'] if 'Tissue Distribution (Direct)' in tup[0]]
                    if len(tissue_direct_formats) > 1:
                        f.write("  <div class='file-formats'>Download in other formats: ")
                        for name, path in tissue_direct_formats:
                            if 'PNG' not in name:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[-1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Tissue distribution table
                if hasattr(self, '_tissue_distribution') and self._tissue_distribution is not None:
                    # Table of top tissues
                    f.write(f"  <h3>Top {min(top_n, len(self._tissue_distribution))} Tissue Types</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Tissue Type</th><th>Sample Count</th><th>Percentage</th></tr>\n")
                    
                    top_tissues = self._tissue_distribution.head(top_n)
                    for _, row in top_tissues.iterrows():
                        f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Sample Count']}</td><td>{row['Percentage (%)']}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                    
                    # Additional visualizations
                    f.write("  <h3>Alternative Visualizations</h3>\n")
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    tissue_bar_png = [path for name, path in tissue_dist_figures if 'Bar Chart' in name]
                    if tissue_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{tissue_bar_png[0]}' alt='Tissue Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    tissue_pie_png = [path for name, path in tissue_dist_figures if 'Pie Chart' in name]
                    if tissue_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{tissue_pie_png[0]}' alt='Tissue Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    tissue_figures_all = [tup for tup in generated_files['figures'] if 'Tissue Distribution' in tup[0] and 'Direct' not in tup[0]]
                    if len(tissue_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        # Group by format
                        format_groups = {}
                        for name, path in tissue_figures_all:
                            format_name = name.split('(')[1].split(')')[0]
                            if format_name not in format_groups:
                                format_groups[format_name] = []
                            format_groups[format_name].append((name, path))
                        
                        # Add links by format
                        for format_name, items in sorted(format_groups.items()):
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Biosample Distribution section
                f.write("  <h2>3. Biosample Type Distribution</h2>\n")
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This section shows the distribution of biosample types (cell lines, primary cells, tissues, etc.) among samples in the SEdb database.</p>\n")
                f.write("    <p>These visualizations help understand the types of biological samples represented in the database.</p>\n")
                f.write("  </div>\n")
                
                # Direct visualization from biosample distribution
                biosample_direct_png = [path for name, path in biosample_dist_figures if 'Direct' in name]
                if biosample_direct_png:
                    f.write("  <h3>Biosample Type Distribution (Direct Visualization)</h3>\n")
                    f.write(f"  <img src='{biosample_direct_png[0]}' alt='Biosample Distribution Direct Visualization'>\n")
                    
                    # Add data link for CSV if available
                    if biosample_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{biosample_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    biosample_direct_formats = [tup for tup in generated_files['figures'] if 'Biosample Distribution (Direct)' in tup[0]]
                    if len(biosample_direct_formats) > 1:
                        f.write("  <div class='file-formats'>Download in other formats: ")
                        for name, path in biosample_direct_formats:
                            if 'PNG' not in name:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[-1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Additional biosample visualizations
                other_biosample_png = [path for name, path in biosample_dist_figures if 'Direct' not in name]
                if other_biosample_png:
                    f.write("  <h3>Alternative Visualizations</h3>\n")
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    biosample_bar_png = [path for name, path in biosample_dist_figures if 'Bar Chart' in name]
                    if biosample_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{biosample_bar_png[0]}' alt='Biosample Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    biosample_pie_png = [path for name, path in biosample_dist_figures if 'Pie Chart' in name]
                    if biosample_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{biosample_pie_png[0]}' alt='Biosample Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    biosample_figures_all = [tup for tup in generated_files['figures'] if 'Biosample Distribution' in tup[0] and 'Direct' not in tup[0]]
                    if len(biosample_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        # Group by format
                        format_groups = {}
                        for name, path in biosample_figures_all:
                            format_name = name.split('(')[1].split(')')[0]
                            if format_name not in format_groups:
                                format_groups[format_name] = []
                            format_groups[format_name].append((name, path))
                        
                        # Add links by format
                        for format_name, items in sorted(format_groups.items()):
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Tissue-Biosample Relationship section
                if relationship_figures:
                    f.write("  <h2>4. Tissue-Biosample Relationship</h2>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>This section shows the relationship between tissue types and biosample types in the SEdb database.</p>\n")
                    f.write("    <p>The heatmap indicates the number of samples with both the corresponding tissue type and biosample type.</p>\n")
                    f.write("    <p>Areas with darker colors represent more samples with that specific tissue and biosample combination.</p>\n")
                    f.write("  </div>\n")
                    
                    # Find PNG relationship image
                    relationship_png = relationship_figures[0][1] if relationship_figures else None
                    if relationship_png:
                        f.write(f"  <img src='{relationship_png}' alt='Tissue-Biosample Relationship'>\n")
                        
                        # Add data link for CSV if available
                        if relationship_csv_path:
                            f.write(f"  <div class='data-link'>📊 <a href='{relationship_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    relationship_figures_all = [tup for tup in generated_files['figures'] if ('Relationship' in tup[0] or 'Heatmap' in tup[0])]
                    if len(relationship_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figure as: ")
                        for name, path in relationship_figures_all:
                            if '.png' not in path:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Cell ID Occurrence Analysis section (SEdb specific)
                if cell_id_figures and 'cell_id_stats' in db_stats:
                    f.write("  <h2>5. Cell ID Occurrence Analysis</h2>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>This section analyzes the distribution of super enhancer occurrences per sample (cell_id).</p>\n")
                    f.write("    <p>It shows how many super enhancers (SEs) are found in each sample, helping to identify")
                    f.write("       samples with unusual SE counts that might be of special interest.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display statistics
                    stats = db_stats['cell_id_stats']
                    f.write("  <h3>Statistics</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Metric</th><th>Value</th></tr>\n")
                    f.write(f"    <tr><td>Mean</td><td>{stats['mean']:.2f} SEs per sample</td></tr>\n")
                    f.write(f"    <tr><td>Standard Deviation</td><td>{stats['std']:.2f}</td></tr>\n")
                    f.write(f"    <tr><td>Minimum</td><td>{int(stats['min'])} SEs</td></tr>\n")
                    f.write(f"    <tr><td>25th Percentile</td><td>{stats['25%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>Median (50th Percentile)</td><td>{stats['50%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>75th Percentile</td><td>{stats['75%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>Maximum</td><td>{int(stats['max'])} SEs</td></tr>\n")
                    f.write("  </table>\n")
                    
                    # Add data link for CSV if available
                    if cell_id_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{cell_id_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Find PNG histogram image
                    histogram_png = [path for name, path in cell_id_figures if 'Histogram' in name]
                    if histogram_png:
                        f.write("  <h3>Distribution</h3>\n")
                        f.write(f"  <img src='{histogram_png[0]}' alt='Cell ID Occurrence Histogram'>\n")
                    
                    # Top samples bar chart
                    top_samples_png = [path for name, path in cell_id_figures if 'Top Samples' in name]
                    if top_samples_png:
                        f.write("  <h3>Top Samples by Super Enhancer Count</h3>\n")
                        f.write(f"  <img src='{top_samples_png[0]}' alt='Top Samples by SE Count'>\n")
                    
                    # Add links to other formats
                    if len(cell_id_figures) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        for fmt in ['SVG', 'PDF', 'EPS']:  # Using common formats
                            format_files = [tup for tup in generated_files['figures'] if 'Cell ID' in tup[0] and fmt in tup[0]]
                            for name, path in format_files:
                                chart_type = "Histogram" if "Histogram" in name else "Bar Chart"
                                f.write(f"<a href='{path}'>{fmt} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Generated Files section
                f.write("  <h2>6. Generated Files</h2>\n")
                
                # List of all data files (tables)
                if generated_files['tables']:
                    f.write("  <h3>Data Files</h3>\n")
                    f.write("  <ul>\n")
                    for name, path in generated_files['tables']:
                        f.write(f"    <li><a href='{path}'>{name} (CSV)</a></li>\n")
                    f.write("  </ul>\n")
                
                # List of figures by category
                if generated_files['figures']:
                    f.write("  <h3>Figures</h3>\n")
                    
                    # Group figures by category
                    figure_groups = {}
                    for name, path in generated_files['figures']:
                        if 'Tissue' in name and 'Biosample' not in name:
                            key = 'Tissue Type Distribution'
                        elif 'Biosample' in name:
                            key = 'Biosample Type Distribution'
                        elif 'Relationship' in name or 'Heatmap' in name:
                            key = 'Tissue-Biosample Relationship'
                        elif 'Cell ID' in name or 'Samples by SE Count' in name:
                            key = 'Cell ID Analysis'
                        else:
                            key = 'Other Figures'
                        
                        if key not in figure_groups:
                            figure_groups[key] = []
                        figure_groups[key].append((name, path))
                    
                    # Display figures by group
                    for group_name, figures in figure_groups.items():
                        f.write(f"  <h4>{group_name}</h4>\n")
                        f.write("  <ul>\n")
                        for name, path in figures:
                            f.write(f"    <li><a href='{path}'>{name}</a></li>\n")
                        f.write("  </ul>\n")
                
                # Summary
                f.write("  <h2>7. Summary</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This report provides a comprehensive overview of the SEdb database metadata.</p>\n")
                f.write("    <p>Key findings include:</p>\n")
                f.write("    <ul>\n")
                f.write(f"      <li>The database contains {db_stats['total_samples']:,} samples with a total of {db_stats['total_se_records']:,} super enhancer records.</li>\n")
                
                if hasattr(self, '_tissue_distribution') and self._tissue_distribution is not None and len(self._tissue_distribution) > 0:
                    top_tissue = self._tissue_distribution.iloc[0]
                    f.write(f"      <li>The most common tissue type is '{top_tissue['Tissue Type']}' with {top_tissue['Sample Count']} samples ({top_tissue['Percentage (%)']}%).</li>\n")
                
                if 'cell_id_stats' in db_stats:
                    f.write(f"      <li>On average, each sample has {db_stats['cell_id_stats']['mean']:.2f} super enhancers.</li>\n")
                
                f.write("    </ul>\n")
                f.write("    <p>The visualizations show the distribution of samples across different biological categories using direct SEdb methods and additional representations.</p>\n")
                f.write("  </div>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by SEdbRegionAnalyzer</p>\n")
                
                # End HTML document
                f.write("</body>\n")
                f.write("</html>\n")
            
            self.logger.info(f"SEdb metadata report generated: {html_path}")
            return html_path
            
        except Exception as e:
            self.logger.error(f"Error generating HTML report: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None




    # Properties
    @property
    def cell_id_counts(self) -> Optional[pd.DataFrame]:
        """
        Access the internal results of cell_id occurrence counts.
        If count_cell_id_occurrences() has not been called, returns None with a warning.
        
        Returns:
            pandas.DataFrame: DataFrame created by count_cell_id_occurrences() or None
        """
        if hasattr(self, '_cell_id_counts'):
            return self._cell_id_counts
        else:
            self.logger.warning("cell_id counts have not been computed yet. Call count_cell_id_occurrences() first.")
            return None

    @property
    def bed_data(self):
        """Access the BED dataframe"""
        if self.bed_df is None:
            self.logger.warning("BED data not loaded yet. Call load_databases first.")
            return None
        return self.bed_df

    @property
    def sample_data(self):
        """Access the sample information dataframe"""
        if self.sample_info is None:
            self.logger.warning("Sample information not loaded yet. Call load_databases first.")
            return None
        return self.sample_info

    @property
    def bed_columns(self):
        """Get the original BED file column names"""
        if not hasattr(self, '_original_bed_columns'):
            self.logger.warning("BED columns not available. Load BED file first.")
            return None
        return self._original_bed_columns
        
    @property
    def original_bed_data(self):
        """Access the original, unfiltered BED dataframe."""
        if self._original_bed_df is None:
            self.logger.warning("Original BED data not available. Load data or run filter_human_chromosomes first.")
            return None
        return self._original_bed_df

    @property
    def removed_bed_data(self):
        """Access the removed BED dataframe."""
        if not hasattr(self,"_removed_bed_df"):
            self.logger.warning("Removed BED data not available. Run filter_human_chromosomes first.")
            return None
        return self._removed_bed_df

    @property
    def tissue_distribution(self):
        """
        Access tissue type distribution data.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing tissue type, sample count, and percentage.
            Returns None if data has not been calculated.
        """
        if hasattr(self, '_tissue_distribution'):
            return self._tissue_distribution
        else:
            self.logger.warning("Tissue distribution data has not been calculated. Run visualize_tissue_distribution() first.")
            return None

    @property
    def overlapping_ses(self):
        """
        Access the dataframe of overlapping SEs.
        
        Returns:
        pandas.DataFrame: DataFrame of extracted SEs
        """
        if self.overlapping_se is None:
            self.logger.warning("Overlapping SE data has not been extracted yet. Run extract_overlapping_se() first.")
            return None
        return self.overlapping_se

    @property
    def region_data(self):
        """
        Access current region information.
        
        Returns:
        dict: Dictionary containing chromosome, start, end, name, etc.
        """
        if self.region is None:
            self.logger.warning("Region data has not been set yet. Run extract_overlapping_se() first.")
            return None
        return self.region

    @property
    def region_tissues(self):
        """
        Access tissue type counts in the region.
        
        Returns:
        pandas.Series: Tissue types and their counts
        """
        if self.region_tissue_counts is None:
            self.logger.warning("Tissue type counts have not been calculated yet. Run extract_overlapping_se() first.")
            return None
        return self.region_tissue_counts

    @property
    def region_biosamples(self):
        """
        Access biosample type counts in the region.
        
        Returns:
        pandas.Series: Biosample types and their counts
        """
        if self.region_biosample_counts is None:
            self.logger.warning("Biosample type counts have not been calculated yet. Run extract_overlapping_se() first.")
            return None
        return self.region_biosample_counts

    @property
    def region_genes(self):
        """
        Access gene counts in the region.
        Gene counts represent the frequency of genes associated with SEs.
        The same gene may be annotated to multiple SEs, even from the same sample.
        
        Returns:
        pandas.Series: Genes and their counts
        """
        if self.region_gene_counts is None:
            self.logger.warning("Gene counts have not been calculated yet. Run extract_overlapping_se() first.")
            return None
        return self.region_gene_counts

    @property
    def region_enrichment(self):
        """
        Access tissue enrichment analysis results for the region.
        
        Returns:
        pandas.DataFrame: Tissue enrichment analysis results
        """
        if not hasattr(self, 'region_tissue_enrichment') or self.region_tissue_enrichment is None:
            self.logger.warning("Tissue enrichment analysis results have not been calculated yet. Run analyze_region_tissue_enrichment() first.")
            return None
        return self.region_tissue_enrichment

    @property
    def region_final_data(self):
        """
        Access the final combined data for the region analysis.
        Contains overlapping SE data joined with sample information.
        
        Returns:
        pandas.DataFrame: Combined final dataframe
        """
        if self.final_data is None:
            self.logger.warning("Final data has not been prepared yet. Run extract_overlapping_se() first.")
            return None
        return self.final_data