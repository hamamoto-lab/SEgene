import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
import os
import pybedtools
from io import StringIO
from scipy import stats
from scipy import signal
from scipy.ndimage import gaussian_filter1d
from statsmodels.stats.multitest import multipletests
import logging
import warnings
import time
from typing import Optional, Dict, List, Tuple, Any
import japanize_matplotlib
import networkx as nx
from datetime import datetime
from typing import Self

from matplotlib.ticker import ScalarFormatter

from collections import defaultdict



# # データローダー分割追加
# from SEgene_analyzer.data_loader.loaders import SEdbDataLoader

# # プロセッサー分割
# from SEgene_analyzer.data_loader.processors import SEdbDataProcessor

# #アナライザー分割
# from SEgene_analyzer.analysis.region import SEdbRegionOverlapAnalyzer

# # 統計解析分割
# from SEgene_analyzer.analysis.statistics import SEdbStatisticsAnalyzer

# # 基本プロット分割:
# from SEgene_analyzer.visualization.plots import SEdbPlotGenerator

# from SEgene_analyzer.visualization.distribution_plots import SEdbDistributionPlotter


# # レポートジェネレータ
# from SEgene_analyzer.reporting.report_generator import SEdbReportGenerator


# # バッチツール
# from SEgene_analyzer.batch.batch_analyzer import SEdbBatchAnalyzer


# # peak_detector実装
# from SEgene_analyzer.analysis.peak_detector import SEdbPeakDetector

# # bed-tsv対応の為に臨時追加
# from SEgene_analyzer.batch.batch_analyzer_extended import SEdbBatchAnalyzerExtended

# 相対インポートに変更
from .data_loader.loaders import SEdbDataLoader
from .data_loader.processors import SEdbDataProcessor
from .analysis.region import SEdbRegionOverlapAnalyzer
from .analysis.statistics import SEdbStatisticsAnalyzer
from .visualization.plots import SEdbPlotGenerator
from .visualization.distribution_plots import SEdbDistributionPlotter
from .reporting.report_generator import SEdbReportGenerator
from .batch.batch_analyzer import SEdbBatchAnalyzer
from .analysis.peak_detector import SEdbPeakDetector
from .batch.batch_analyzer_extended import SEdbBatchAnalyzerExtended


class SEdbRegionAnalyzer:
    """
    Super-enhancer (SE) genomic region analyzer for SEdb data.
    
    This class loads SEdb database files and provides methods to analyze
    specific genomic regions for overlapping super enhancers.
    """
    
    def __init__(self, results_dir='results', logger=None, 
                se_bed_file=None, sample_info_file=None,
                enable_peak_analysis=False):  # peakdetecterオフの為に新しいパラメータを追加):
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

        # peak_detector設定フラグを保存
        self.enable_peak_analysis = enable_peak_analysis

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


        # 各種コンポーネントの初期化
        self.data_loader = SEdbDataLoader(logger=self.logger)
        self.data_processor = SEdbDataProcessor(logger=self.logger)
        self.region_analyzer = SEdbRegionOverlapAnalyzer(logger=self.logger)
        self.statistics_analyzer = SEdbStatisticsAnalyzer(logger=self.logger)
        self.plot_generator = SEdbPlotGenerator(logger=self.logger)
        self.distribution_plotter = SEdbDistributionPlotter(logger=self.logger)

        self.report_generator = SEdbReportGenerator(self)

        self.peak_detector = SEdbPeakDetector(logger=self.logger)

        self.logger.info("SEgeneRegionAnalyzer initialization completed successfully")


  # 元のメソッドを修正して新しいローダーを使用
    def _read_bed_file(self, file_path):
        """データローダーを使用してBEDファイルを読み込む内部メソッド"""
        se_df, self._original_bed_columns = self.data_loader.read_bed_file(file_path)
        return se_df
            
    def _read_sample_info(self, file_path):
        """データローダーを使用してサンプル情報を読み込む内部メソッド"""
        return self.data_loader.read_sample_info(file_path)

    def _remove_na_sample_ids(self) -> None:
        """
        Remove rows from sample_info where 'Sample ID' is NaN.
        This method is intended for internal use (private).
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Nothing to remove.")
            return

        # データプロセッサを使用してNaN Sample IDを除去
        self.sample_info = self.data_processor.remove_na_sample_ids(self.sample_info)




    def _add_cell_id_to_sample_info(self) -> None:
        """
        Add 'cell_id' column to sample_info, derived from 'Sample ID', and make it the first column.
        Internal method.
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot add cell_id.")
            return

        # データプロセッサを使用してcell_idカラムを追加
        self.sample_info = self.data_processor.add_cell_id_to_sample_info(self.sample_info)


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


    def filter_human_chromosomes(self, chromosomes: Optional[List[str]] = None) -> 'SEdbRegionAnalyzer':
        """
        Filter the BED data to include only specified chromosomes (default: chr1-22, chrX, chrY).
        Keeps the original data in _original_bed_df, and filtered data in bed_df.

        Parameters:
            chromosomes (list of str, optional): List of chromosomes to keep.
                Defaults to standard human chromosomes (chr1-chr22, chrX, chrY).

        Returns:
            SEdbRegionAnalyzer: self, for method chaining
        """
        if self.bed_df is None:
            self.logger.warning("BED data is not loaded. Nothing to filter.")
            return self

        # データプロセッサを使用して染色体フィルタリング
        filtered_df, original_df, removed_df = self.data_processor.filter_human_chromosomes(
            self.bed_df, chromosomes
        )
        
        # 結果を保存
        self._original_bed_df = original_df
        self.bed_df = filtered_df
        self._removed_bed_df = removed_df
        
        # 削除したデータが存在する場合のみログ出力
        if removed_df is not None and not removed_df.empty:
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

        # データプロセッサを使用してcell_idの出現回数をカウント
        self._cell_id_counts = self.data_processor.count_cell_id_occurrences(
            self.bed_df, self.sample_info
        )
        
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
        
        # 分布プロッターに処理を委譲
        tissue_distribution = self.distribution_plotter.visualize_tissue_distribution(
            self.sample_info,
            top_n=top_n,
            horizontal=horizontal,
            store_all=store_all,
            figsize=figsize,
            save_dir=save_dir,
            save_filename=save_filename,
            save_png=save_png,
            save_svg=save_svg,
            save_pdf=save_pdf,
            save_eps=save_eps,
            save_dpi=save_dpi,
            save_transparent=save_transparent
        )
        
        # 全データを保存（オプション）
        if store_all and tissue_distribution is not None:
            self._tissue_distribution = tissue_distribution
        
        return tissue_distribution







    # analyze_tissue_biosample_relationship メソッドを更新
    def analyze_tissue_biosample_relationship(self, figsize=(15, 12), min_count=5):
        """
        組織タイプとバイオサンプルタイプの関係を分析します。
        
        Parameters:
        -----------
        figsize : tuple, optional
            図のサイズ
        min_count : int, optional
            ヒートマップに含める最小カウント数
        
        Returns:
        --------
        pandas.DataFrame
            組織タイプとバイオサンプルタイプのクロスタブ
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Please run load_databases().")
            return None
        
        # 分布プロッターに可視化を委譲
        cross_tab = self.distribution_plotter.visualize_tissue_biosample_relationship(
            self.sample_info,
            figsize=figsize,
            min_count=min_count
        )
        
        return cross_tab





    # analyze_data_sources メソッドを更新
    def analyze_data_sources(self, figsize=(12, 8)):
        """
        データソースごとの分析を行います。
        
        Parameters:
        -----------
        figsize : tuple, optional
            図のサイズ
        """
        if self.sample_info is None:
            self.logger.warning("Sample information is not loaded. Please run load_databases().")
            return
        
        # 分布プロッターに可視化を委譲
        _, summary = self.distribution_plotter.visualize_data_sources(
            self.sample_info,
            figsize=figsize
        )
        
        if summary:
            # サマリー情報を表示
            source_counts = summary['source_counts']
            
            print("\n=== Sample Count by Data Source ===")
            source_df = pd.DataFrame({
                'Data Source': source_counts.index,
                'Sample Count': source_counts.values,
                'Percentage (%)': (source_counts / source_counts.sum() * 100).round(2).values
            })
            print(source_df)
            
            if summary.get('tissue_source') is not None:
                print("\n=== Cross-tabulation of Tissue Type and Data Source ===")
                print(summary['tissue_source'])
            
            if summary.get('biosample_source') is not None:
                print("\n=== Cross-tabulation of Biosample Type and Data Source ===")
                print(summary['biosample_source'])




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
        
        # 統計アナライザを使用して組織タイプのエンリッチメントをテスト
        return self.statistics_analyzer.test_tissue_enrichment(
            region_tissue_counts,
            self._tissue_distribution,
            method=method,
            correction=correction
        )



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
        # 統計アナライザを使用して複数領域を比較
        fig, comparison_data = self.statistics_analyzer.compare_multiple_regions(
            region_results_dict,
            plot_type=plot_type,
            figsize=figsize,
            p_threshold=p_threshold
        )
        
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
        SEdbRegionAnalyzer: self for method chaining
        """
        if self.bed_df is None:
            self.logger.warning("BED data not loaded. Call load_databases first.")
            return self
        
        self.logger.info(f"Extracting SEs overlapping with region: {region_chr}:{region_start}-{region_end}")
        
        # 領域分析クラスを使用して重複するSEを抽出
        results = self.region_analyzer.extract_overlapping_se(
            self.bed_df,
            self.sample_info,
            region_chr,
            region_start,
            region_end,
            region_name,
            extraction_method
        )
        
        # 結果を保存
        if results:
            self.region = results['region']
            self.extraction_method = results['extraction_method']
            self.overlapping_se = results['overlapping_se']
            self.unique_se_data = results['unique_se_data']
            self.final_data = results['final_data']
            self.unique_sample_count = results['unique_sample_count']
            self.region_tissue_counts = results['region_tissue_counts']
            self.region_biosample_counts = results['region_biosample_counts']
            self.region_gene_counts = results['region_gene_counts']
            
            # サンプル概要があれば保存
            if 'sample_summary' in results:
                self.sample_summary = results['sample_summary']
        
        return self


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
        
        # 統計アナライザを使用してサンプル分布を分析
        results = self.statistics_analyzer.analyze_samples_distribution(
            self.sample_summary,
            figsize=figsize
        )
        
        if results is None:
            return
        
        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        
        # 結果を表示
        print(f"=== Sample Distribution Analysis for Region: {region_name} ===")
        print(f"Total unique samples with overlapping SEs: {results['sample_count']}")
        
        # SE count per sample statistics
        se_count_stats = results['se_count_stats']
        print("\nSE Count per Sample Statistics:")
        for stat, value in se_count_stats.items():
            print(f"  {stat}: {value:.2f}" if isinstance(value, float) else f"  {stat}: {value}")
        
        # Figure object displays the plots
        if 'figure' in results:
            plt = results['figure']
        
        # Print additional statistics
        if results['tissue_counts'] is not None:
            print("\nSample Count by Tissue Type (each sample counted once):")
            print(results['tissue_counts'])
        
        if results['biosample_counts'] is not None:
            print("\nSample Count by Biosample Type (each sample counted once):")
            print(results['biosample_counts'])
        
        return results



    def analyze_region_tissue_enrichment(self, method='fisher', correction='fdr_bh', plot=True, **kwargs):
        """
        Analyze tissue enrichment in the overlapping SEs compared to the background distribution.
        Counts are based on unique samples (each sample is counted only once).
        
        Parameters:
        -----------
        method : str, optional
            Statistical test method ('fisher' or 'chi2')
        correction : str, optional
            Multiple testing correction method
        plot : bool, optional
            Whether to plot the results
        **kwargs: Additional parameters passed to visualize_tissue_enrichment
        
        Returns:
        --------
        pandas.DataFrame
            Results of the enrichment analysis
        """
        if self.region_tissue_counts is None:
            self.logger.warning("No region tissue counts available. Call extract_overlapping_se first.")
            return None
        
        region_name = self.region.get('name', f"{self.region['chr']}:{self.region['start']}-{self.region['end']}")
        
        # 統計アナライザを使用して組織タイプのエンリッチメント分析を実行
        enrichment_results = self.statistics_analyzer.analyze_region_tissue_enrichment(
            self.region_tissue_counts,
            self._tissue_distribution,
            region_name=region_name,
            method=method,
            correction=correction
        )
        
        # 結果を保存
        self.region_tissue_enrichment = enrichment_results
        
        # リクエストに応じて結果を可視化
        if plot and enrichment_results is not None:
            self.visualize_tissue_enrichment(
                enrichment_results=enrichment_results,
                region_name=region_name,
                **kwargs
            )
        
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
        
        # 領域分析クラスを使用してユニークなSEデータを取得
        unique_se_data = self.region_analyzer.get_unique_se_data(self.overlapping_se)
        
        # 結果を保存
        self.unique_se_data = unique_se_data
        
        return unique_se_data




    def display_cell_id_summary(self, top_n=10, bottom_n=10, show_histogram=True, bins=30, figsize=(15, 10)):
        """
        cell_idの出現回数に関するサマリー情報を表示します。
        
        Parameters:
        -----------
        top_n : int, optional
            表示する最頻出cell_idの数
        bottom_n : int, optional
            表示する最少出現cell_idの数
        show_histogram : bool, optional
            ヒストグラムを表示するかどうか
        bins : int, optional
            ヒストグラムのビン数
        figsize : tuple, optional
            図のサイズ
        """
        if not hasattr(self, '_cell_id_counts') or self._cell_id_counts is None:
            self.logger.warning("cell_id counts not available. Running count_cell_id_occurrences()...")
            self.count_cell_id_occurrences()
            if not hasattr(self, '_cell_id_counts') or self._cell_id_counts is None:
                return
        
        # 分布プロッターに可視化を委譲
        _, summary = self.distribution_plotter.visualize_cell_id_distribution(
            self._cell_id_counts,
            top_n=top_n,
            bottom_n=bottom_n,
            show_histogram=show_histogram,
            bins=bins,
            figsize=figsize
        )
        
        if summary:
            # 基本統計情報を表示
            stats = summary['stats']
            print("=== Basic Statistics of cell_id Occurrences ===")
            print(pd.Series(stats))
            
            # ゼロカウントの表示
            zero_counts = summary['zero_counts']
            if zero_counts > 0:
                print(f"\ncell_ids with zero occurrences: {zero_counts} ({summary['zero_percent']:.2f}%)")
            
            # 最頻出cell_ids
            print(f"\n=== Most Frequent cell_ids (Top {top_n}) ===")
            print(summary['top_counts'])
            
            # 最低頻出cell_ids
            if len(summary['bottom_counts']) > 0:
                print(f"\n=== Least Frequent cell_ids (Bottom {bottom_n}, excluding zeros) ===")
                print(summary['bottom_counts'])



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
        # 解析結果を取得
        if enrichment_results is None:
            if hasattr(self, 'region_tissue_enrichment'):
                enrichment_results = self.region_tissue_enrichment
            else:
                self.logger.warning("No enrichment analysis results available. Please run analyze_region_tissue_enrichment() first.")
                return None
        
        # 統計アナライザを使用してエンリッチメント結果を精緻化
        return self.statistics_analyzer.refine_tissue_enrichment_results(
            enrichment_results,
            fdr_threshold=fdr_threshold,
            fc_threshold=fc_threshold,
            top_n=top_n,
            rename_columns=rename_columns
        )


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
        

        # 緊急修理
        # # Filter data for the specified region
        # region_data = self.final_data[(self.final_data['se_chr'] == chrom) & 
        #                         (self.final_data['se_start'] >= start) & 
        #                         (self.final_data['se_end'] <= end)]

        region_data = self.final_data[(self.final_data['se_chr'] == chrom) & 
                                (self.final_data['se_start'] <= end) & 
                                (self.final_data['se_end'] >= start)]


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
        

        # ピーク検出（オプション）
        peak_data = {}
        if detect_peaks:
            # ピーク検出を専用クラスに委譲
            peak_data, _, _ = self.peak_detector.detect_peaks(
                reads, start, end, bin_size, peak_threshold,
                min_peak_distance, peak_width_factor, smooth_window
            )
            
            # 表示用データの取得
            peak_display_data = self.peak_detector.get_peak_display_data(peak_data)
        

            
            # ピークハイライトを専用メソッドに委譲
            self._highlight_peaks(
                ax=ax, 
                peak_data=peak_display_data,
                stacks=stacks,
                peak_color=peak_color,
                peak_alpha=peak_alpha,
                annotate_peaks=annotate_peaks,
                peak_annotation_fontsize=peak_annotation_fontsize,
                start=start,
                end=end
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
        
        This is a simplified wrapper around plot_stacked_reads_with_peaks
        that doesn't perform peak detection.
        
        Parameters:
        [Parameter documentation preserved from original method]
        
        Returns:
        matplotlib.figure.Figure: The created plot's Figure object
        """
        # Simply call plot_stacked_reads_with_peaks with detect_peaks=False
        fig, _ = self.plot_stacked_reads_with_peaks(
            chrom=chrom, 
            start=start, 
            end=end,
            title=title, 
            xlabel=xlabel, 
            ylabel=ylabel,
            color=color, 
            grid=grid, 
            linewidth=linewidth,
            save_dir=save_dir, 
            save_filename=save_filename,
            save_svg=save_svg, 
            save_png=save_png, 
            save_pdf=save_pdf, 
            save_eps=save_eps,
            save_dpi=save_dpi, 
            save_transparent=save_transparent, 
            save_region_bed=save_region_bed,
            show_plot=show_plot, 
            title_fontsize=title_fontsize, 
            xlabel_fontsize=xlabel_fontsize,
            ylabel_fontsize=ylabel_fontsize, 
            tick_fontsize=tick_fontsize, 
            annotation_fontsize=annotation_fontsize,
            figsize=figsize,
            # Peak detection related parameters - turn off
            detect_peaks=False  # Turn off peak detection
        )
        
        return fig







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
            
            # Peak detection plot returns both figure and peak data (only if enabled)
            if self.enable_peak_analysis:
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
        tissue_dist_tables = []
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
                    tissue_table_relative_path = os.path.join(tables_dir, tissue_table_filename)
                    generated_files['tables'].append(('Tissue Distribution', tissue_table_relative_path))
                    tissue_dist_tables.append(('Tissue Distribution Data', tissue_table_relative_path))
            
            except Exception as e:
                self.logger.warning(f"Error generating tissue distribution visualizations: {e}")
        
        # 2. Generate biosample type distribution visualizations if available
        biosample_dist_figures = []
        biosample_dist_tables = []
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
                    biosample_table_relative_path = os.path.join(tables_dir, biosample_table_filename)
                    generated_files['tables'].append(('Biosample Distribution', biosample_table_relative_path))
                    biosample_dist_tables.append(('Biosample Distribution Data', biosample_table_relative_path))
            
            except Exception as e:
                self.logger.warning(f"Error generating biosample distribution visualizations: {e}")
        
        # 3. Generate tissue enrichment visualizations
        tissue_enrichment_figures = []
        tissue_enrichment_tables = []
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
                    tissue_enrich_relative_path = os.path.join(tables_dir, tissue_enrich_filename)
                    generated_files['tables'].append(('Tissue Enrichment Analysis', tissue_enrich_relative_path))
                    tissue_enrichment_tables.append(('Tissue Enrichment Data', tissue_enrich_relative_path))
            except Exception as e:
                self.logger.warning(f"Error generating tissue enrichment analysis: {e}")
        
        # 4. Generate stacked read plots
        stacked_plot_figures = []
        stacked_plot_tables = []
        
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
        
        # # 4.2 Stacked plot with peaks
        # peak_data_tables = []

        # ▼▼▼【修正箇所 1/3】ここからif文を追加 ▼▼▼
        # 4.2 Stacked plot with peaks (条件付きで実行)
        peak_data_tables = []
        if self.enable_peak_analysis: # フラグをチェック
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
                        peak_data_relative_path = os.path.join(tables_dir, peak_data_filename)
                        generated_files['tables'].append(('Peak Analysis Data', peak_data_relative_path))
                        peak_data_tables.append(('Peak Analysis Data', peak_data_relative_path))
                        
                        # Also save list of samples per peak
                        peak_sample_files = []
                        for peak_id, peak_info in peak_data.items():
                            if peak_info['samples']['list']:
                                peak_samples_filename = f'region_{chrom}_{start}_{end}_{peak_id}_samples_{timestamp}.txt'
                                peak_samples_path = os.path.join(tables_full_path, peak_samples_filename)
                                
                                with open(peak_samples_path, 'w') as f:
                                    for sample in peak_info['samples']['list']:
                                        f.write(f"{sample}\n")
                                
                                peak_sample_relative_path = os.path.join(tables_dir, peak_samples_filename)
                                generated_files['tables'].append((f'Samples in {peak_id}', peak_sample_relative_path))
                                peak_data_tables.append((f'Samples in {peak_id}', peak_sample_relative_path))
            except Exception as e:
                self.logger.warning(f"Error generating stacked read plot with peaks: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
        # ▲▲▲【修正箇所 1/3】ここまでif文 ▲▲▲




        
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
                f.write("    .download-group { display: flex; flex-wrap: wrap; margin: 15px 0; background-color: #f9f9f9; padding: 10px; border-radius: 5px; }\n")
                f.write("    .download-section { margin-right: 20px; }\n")
                f.write("    .download-section-title { font-weight: bold; margin-bottom: 5px; color: #555; }\n")
                f.write("    .download-links a { margin-right: 8px; }\n")
                f.write("    .download-image-links a { color: #3498db; }\n")
                f.write("    .download-data-links a { color: #e67e22; }\n")
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
                    
                    # Add combined download links (figures and data together)
                    f.write("  <div class='download-group'>\n")
                    
                    # Image download section
                    f.write("    <div class='download-section'>\n")
                    f.write("      <div class='download-section-title'>Download figures:</div>\n")
                    f.write("      <div class='download-links download-image-links'>\n")
                    for name, path in generated_files['figures']:
                        if 'Tissue Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"        <a href='{path}'>{format_name}</a>\n")
                    f.write("      </div>\n")
                    f.write("    </div>\n")
                    
                    # Data download section
                    if tissue_dist_tables:
                        f.write("    <div class='download-section'>\n")
                        f.write("      <div class='download-section-title'>Download data:</div>\n")
                        f.write("      <div class='download-links download-data-links'>\n")
                        for name, path in tissue_dist_tables:
                            f.write(f"        <a href='{path}'>CSV</a>\n")
                        f.write("      </div>\n")
                        f.write("    </div>\n")
                    
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
                    
                    # Add combined download links (figures and data together)
                    f.write("  <div class='download-group'>\n")
                    
                    # Image download section
                    f.write("    <div class='download-section'>\n")
                    f.write("      <div class='download-section-title'>Download figures:</div>\n")
                    f.write("      <div class='download-links download-image-links'>\n")
                    for name, path in generated_files['figures']:
                        if 'Biosample Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"        <a href='{path}'>{format_name}</a>\n")
                    f.write("      </div>\n")
                    f.write("    </div>\n")
                    
                    # Data download section
                    if biosample_dist_tables:
                        f.write("    <div class='download-section'>\n")
                        f.write("      <div class='download-section-title'>Download data:</div>\n")
                        f.write("      <div class='download-links download-data-links'>\n")
                        for name, path in biosample_dist_tables:
                            f.write(f"        <a href='{path}'>CSV</a>\n")
                        f.write("      </div>\n")
                        f.write("    </div>\n")
                    
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
                    
                    # Add combined download links (figures and data together)
                    f.write("  <div class='download-group'>\n")
                    
                    # Image download section
                    f.write("    <div class='download-section'>\n")
                    f.write("      <div class='download-section-title'>Download figures:</div>\n")
                    f.write("      <div class='download-links download-image-links'>\n")
                    for name, path in generated_files['figures']:
                        if 'Tissue Enrichment' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"        <a href='{path}'>{format_name}</a>\n")
                    f.write("      </div>\n")
                    f.write("    </div>\n")
                    
                    # Data download section
                    if tissue_enrichment_tables:
                        f.write("    <div class='download-section'>\n")
                        f.write("      <div class='download-section-title'>Download data:</div>\n")
                        f.write("      <div class='download-links download-data-links'>\n")
                        for name, path in tissue_enrichment_tables:
                            f.write(f"        <a href='{path}'>CSV</a>\n")
                        f.write("      </div>\n")
                        f.write("    </div>\n")
                    
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
                    
                    # Add combined download links
                    f.write("  <div class='download-group'>\n")
                    
                    # Image download section
                    f.write("    <div class='download-section'>\n")
                    f.write("      <div class='download-section-title'>Download figures:</div>\n")
                    f.write("      <div class='download-links download-image-links'>\n")
                    for name, path in generated_files['figures']:
                        if 'Simple Stacked Read Plot' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"        <a href='{path}'>{format_name}</a>\n")
                    f.write("      </div>\n")
                    f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                else:
                    f.write("  <p class='notes'>No simple read distribution plot available for this region.</p>\n")
                
                # ▼▼▼【修正箇所 2/3】ここからif文を追加 ▼▼▼
                # 5.2 Stacked Reads with Peaks
                if self.enable_peak_analysis: # フラグをチェック


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
                        
                        # Add combined download links for peaks (images and data)
                        f.write("  <div class='download-group'>\n")
                        
                        # Image download section
                        f.write("    <div class='download-section'>\n")
                        f.write("      <div class='download-section-title'>Download figures:</div>\n")
                        f.write("      <div class='download-links download-image-links'>\n")
                        for name, path in generated_files['figures']:
                            if 'Stacked Read Plot with Peaks' in name or 'peak' in name.lower():
                                format_name = name.split('(')[1].split(')')[0] if '(' in name else path.split('.')[-1].upper()
                                f.write(f"        <a href='{path}'>{format_name}</a>\n")
                        f.write("      </div>\n")
                        f.write("    </div>\n")
                        
                        # Data download section for peak data
                        if peak_data_tables:
                            f.write("    <div class='download-section'>\n")
                            f.write("      <div class='download-section-title'>Download peak data:</div>\n")
                            f.write("      <div class='download-links download-data-links'>\n")
                            # For main peak data CSV
                            for name, path in peak_data_tables:
                                if name == 'Peak Analysis Data':
                                    f.write(f"        <a href='{path}'>CSV</a>\n")
                            f.write("      </div>\n")
                            f.write("    </div>\n")
                            
                            # Add section for peak sample lists
                            has_sample_lists = False
                            for name, path in peak_data_tables:
                                if name != 'Peak Analysis Data':
                                    has_sample_lists = True
                                    break
                            
                            if has_sample_lists:
                                f.write("    <div class='download-section'>\n")
                                f.write("      <div class='download-section-title'>Download peak sample lists:</div>\n")
                                f.write("      <div class='download-links download-data-links'>\n")
                                for name, path in peak_data_tables:
                                    if name != 'Peak Analysis Data':
                                        peak_id = name.replace('Samples in ', '')
                                        f.write(f"        <a href='{path}'>{peak_id}</a>\n")
                                f.write("      </div>\n")
                                f.write("    </div>\n")
                        
                        f.write("  </div>\n")
                    elif not has_peak_plot and has_peak_data:
                        # Peak data exists but plot not found
                        f.write("  <p class='notes'>Peaks were detected, but the visualization is not available in this report.</p>\n")
                    else:
                        f.write("  <p class='notes'>No peak detection plot available for this region.</p>\n")
                # ▲▲▲【修正箇所 2/3】ここまでif文 ▲▲▲


                # ▼▼▼【修正箇所 3/3】ここからif文を追加 ▼▼▼
                peak_data = analysis_results.get('peak_data')


                # ▼▼▼【修正箇所 3/3】ここからif文を追加 ▼▼▼

                # Peak information if available
                peak_data = analysis_results.get('peak_data')
                if self.enable_peak_analysis and peak_data: # フラグとデータの存在をチェック
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

                # ▲▲▲【修正箇所 3/3】ここまでif文 ▲▲▲


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
                        title_y: float = 1.1,
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
        
        この関数は plot_generator に委譲
        """
        return self.plot_generator.plot_simple_pie(
            data=data,
            category=category, 
            count_col=count_col,
            title=title, 
            limit=limit,
            figsize=figsize,
            colormap=colormap,
            title_y=title_y,
            save_dir=save_dir, 
            save_filename=save_filename,
            save_svg=save_svg, 
            save_png=save_png, 
            save_pdf=save_pdf, 
            save_eps=save_eps,
            save_dpi=save_dpi, 
            save_transparent=save_transparent
        )

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
        
        この関数は plot_generator に委譲
        """
        return self.plot_generator.plot_sorted_bar(
            data=data,
            category=category, 
            count_col=count_col,
            title=title, 
            limit=limit,
            horizontal=horizontal, 
            color=color,
            figsize=figsize,
            save_dir=save_dir, 
            save_filename=save_filename,
            save_svg=save_svg, 
            save_png=save_png, 
            save_pdf=save_pdf, 
            save_eps=save_eps,
            save_dpi=save_dpi, 
            save_transparent=save_transparent
        )

    def generate_sedb_report(self, output_dir='results/sedb_report',
                        top_n=20, figure_formats=None, save_tables=True,
                        dpi=600, fig_width=12, fig_height=8):
        """
        SEdbメタデータの総合HTMLレポートを生成します
        
        Parameters:
        -----------
        output_dir : str, optional
            出力ディレクトリ（デフォルト: 'results/sedb_report'）
        top_n : int, optional
            可視化で表示するトップカテゴリの数（デフォルト: 20）
        figure_formats : list, optional
            保存する図のフォーマット（デフォルト: ['png', 'svg', 'pdf']）
        save_tables : bool, optional
            分布テーブルをCSVファイルとして保存するかどうか（デフォルト: True）
        dpi : int, optional
            ラスター画像のDPI（デフォルト: 300）
        fig_width : int, optional
            図の幅（インチ単位）（デフォルト: 12）
        fig_height : int, optional
            図の高さ（インチ単位）（デフォルト: 8）
            
        Returns:
        --------
        str
            生成されたHTMLレポートのパス
        """
        # データの存在確認
        if self.bed_df is None or self.sample_info is None:
            self.logger.warning("SE BED data or sample information not loaded. Call load_databases() first.")
            return None
        
        # レポート生成処理をreport_generatorに委譲
        return self.report_generator.generate_sedb_report(
            output_dir=output_dir,
            top_n=top_n,
            figure_formats=figure_formats,
            save_tables=save_tables,
            dpi=dpi,
            fig_width=fig_width,
            fig_height=fig_height
        )



    def batch_analyze_regions_from_tsv(self, 
                                    bed_file: str,  # ディレクトリではなく単一ファイルに変更
                                    metadata_file: str,
                                    data_output_dir: str,
                                    report_output_dir: str, 
                                    regions_tsv_file: str,
                                    filter_chromosomes: bool = True,  # 追加
                                    count_cell_id: bool = True,       # 追加
                                    max_rows: int = None,
                                    **kwargs):
        """
        TSVファイルから複数のゲノム領域を一括解析し、結果を新たなTSVとして出力します。
        
        SEdbBatchAnalyzer クラスに委譲します。詳細はそちらのドキュメントを参照してください。
        
        Parameters:
        -----------
        bed_file : str
            BEDファイルのパス
        metadata_file : str
            SEdbメタデータファイルのパス
        data_output_dir : str
            解析データの出力先ディレクトリ
        report_output_dir : str
            レポートファイルの出力先ディレクトリ
        regions_tsv_file : str
            解析対象のゲノム領域が記載されたTSVファイルのパス
        filter_chromosomes : bool, optional
            人間の標準染色体でフィルタリングするかどうか（デフォルト: True）
        count_cell_id : bool, optional
            cell_idの出現回数をカウントするかどうか（デフォルト: True）
        max_rows : int, optional
            TSVファイルから読み込む最大行数。Noneの場合は全行を読み込みます
        **kwargs : dict
            その他のパラメータ（analyze_and_report_regionに渡されます）
            
        Returns:
        --------
        str
            生成された結果TSVファイルのパス
        """
        batch_analyzer = SEdbBatchAnalyzer(self)
        return batch_analyzer.batch_analyze_regions_from_tsv(
            bed_file=bed_file,  # 変更：bed_directory → bed_file
            metadata_file=metadata_file,
            data_output_dir=data_output_dir,
            report_output_dir=report_output_dir,
            regions_tsv_file=regions_tsv_file,
            filter_chromosomes=filter_chromosomes,  # 追加
            count_cell_id=count_cell_id,            # 追加
            max_rows=max_rows,
            **kwargs
        )



    # ver3
    def generate_region_report_extended(self, chrom, start, end, region_name=None, output_dir=None, 
                                    figure_formats=None, save_tables=True, peak_threshold=0.2,
                                    pvalue_threshold=0.05, fdr_threshold=0.05, fc_threshold=1.5):
        """
        Extended report generation - wrapper that returns the expected dictionary format
        
        This method calls analyze_and_report_region to ensure compatibility with batch processing.
        
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
        self.logger.info("Extended report generation: calling analyze_and_report_region for compatibility")
        
        try:
            # Call analyze_and_report_region which returns the expected dictionary format
            result_dict = self.analyze_and_report_region(
                chrom=chrom,
                start=start,
                end=end,
                region_name=region_name,
                output_dir=output_dir,
                figure_formats=figure_formats,
                save_tables=save_tables,
                peak_threshold=peak_threshold,
                pvalue_threshold=pvalue_threshold,
                fdr_threshold=fdr_threshold,
                fc_threshold=fc_threshold
            )
            
            self.logger.info(f"Extended report generation completed successfully")
            return result_dict
            
        except Exception as e:
            self.logger.error(f"Error in extended report generation: {str(e)}")
            import traceback
            self.logger.error(traceback.format_exc())
            
            # Return a minimal valid dictionary in case of error
            error_result = {
                'report_path': None,
                'analysis_results': None,
                'has_overlaps': False,
                'overlap_sample_count': 0,
                'overlap_read_count': 0,
                'significant_tissues': []
            }
            return error_result


    #bed-tsv対応の為に臨時追加
    # def batch_analyze_regions_from_file(self, **kwargs):
    #     """
    #     BED/TSVファイル対応の一括解析（拡張版）
    #     """
        
    #     extended_batch_analyzer = SEdbBatchAnalyzerExtended(self)
    #     return extended_batch_analyzer.batch_analyze_regions_from_file(**kwargs)

    def batch_analyze_regions_from_file(self, 
                                    bed_file: str,
                                    metadata_file: str,
                                    data_output_dir: str,
                                    report_output_dir: str, 
                                    regions_tsv_file: str = None,  # TSVファイル（後方互換性のため）
                                    regions_file: str = None,     # TSVまたはBEDファイル
                                    **kwargs):
        """
        BED/TSVファイル対応の一括解析（拡張版）
        
        Parameters:
        -----------
        bed_file : str
            BEDファイルのパス
        metadata_file : str
            SEdbメタデータファイルのパス
        data_output_dir : str
            解析データの出力先ディレクトリ
        report_output_dir : str
            レポートファイルの出力先ディレクトリ
        regions_tsv_file : str, optional
            解析対象のゲノム領域が記載されたTSVファイルのパス（後方互換性のため）
        regions_file : str, optional
            解析対象のゲノム領域ファイル（TSVまたはBED形式）
        **kwargs : dict
            その他のパラメータ
            
        Returns:
        --------
        str
            生成された結果TSVファイルのパス
        """
        # regions_file と regions_tsv_file の両方が指定された場合、regions_file を優先
        if regions_file is None:
            if regions_tsv_file is not None:
                regions_file = regions_tsv_file
            else:
                raise ValueError("Either 'regions_file' or 'regions_tsv_file' must be specified")
        
        extended_batch_analyzer = SEdbBatchAnalyzerExtended(self)
        return extended_batch_analyzer.batch_analyze_regions_from_file(
            bed_file=bed_file,
            metadata_file=metadata_file,
            data_output_dir=data_output_dir,
            report_output_dir=report_output_dir,
            regions_file=regions_file,
            **kwargs
        )


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
    



#####移動予定

    def _highlight_peaks(self, ax, peak_data, stacks, peak_color='red', peak_alpha=0.2,
                        annotate_peaks=True, peak_annotation_fontsize=12, start=None, end=None):
        """
        ピークをハイライト表示する内部メソッド
        
        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            描画対象の軸オブジェクト
        peak_data : dict
            ピーク情報を含む辞書
        stacks : list
            読み込みスタックのリスト
        peak_color : str, optional
            ピークハイライトの色 (デフォルト: 'red')
        peak_alpha : float, optional
            ピークハイライトの透明度 (デフォルト: 0.2)
        annotate_peaks : bool, optional
            ピークにラベルを付けるかどうか (デフォルト: True)
        peak_annotation_fontsize : int, optional
            ピークラベルのフォントサイズ (デフォルト: 12)
        start : int, optional
            領域の開始位置
        end : int, optional
            領域の終了位置
        """
        if not peak_data:
            return
        
        # Sort peaks by position
        sorted_peaks = sorted(peak_data.items(), key=lambda x: x[1]['position'])
        
        # Calculate peak label positions
        vertical_positions = {}
        for i, (peak_id, peak) in enumerate(sorted_peaks):
            # Position labels inside the top part of the graph
            vertical_positions[peak_id] = len(stacks) * 0.9
            
            # Adjust position if too close to previous peak
            if i > 0:
                prev_id, prev_peak = sorted_peaks[i-1]
                region_width = end - start if (start is not None and end is not None) else 1
                if (peak['position'] - prev_peak['position']) / region_width < 0.04:
                    vertical_positions[peak_id] = vertical_positions[prev_id] - len(stacks) * 0.1
        
        # Highlight peaks with rectangles - use display data for visualization
        for peak_id, peak in peak_data.items():
            # Determine left position and width based on available keys
            if 'left' in peak and 'width' in peak:
                # Use original structure
                left_pos = peak['left']
                width = peak['width']
            elif 'start' in peak and 'width' in peak:
                # Use new structure
                left_pos = peak['start']
                width = peak['width'] 
            elif 'start' in peak and 'end' in peak:
                # Calculate from start/end
                left_pos = peak['start']
                width = peak['end'] - peak['start']
            else:
                # Skip if required data not available
                continue
                
            # Highlight peak region with colored rectangle
            rect = plt.Rectangle(
                (left_pos, -1),  # (x, y)
                width,  # width
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