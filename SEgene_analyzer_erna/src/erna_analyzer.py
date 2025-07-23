import os
import pandas as pd
import numpy as np

# Set non-interactive backend before ANY matplotlib import
import matplotlib
matplotlib.use('Agg', force=True)  # Force non-interactive backend for WSL/headless environments
import matplotlib.pyplot as plt

# Disable matplotlib warnings for non-interactive backend
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='matplotlib')
warnings.filterwarnings('ignore', message='.*non-interactive.*cannot be shown.*')
from pybedtools import BedTool
import seaborn as sns
from collections import Counter
import logging
import glob
import time
from typing import List, Dict, Tuple, Any, Optional
from collections import defaultdict
from natsort import natsorted
from datetime import datetime

import os
import pandas as pd
import numpy as np
# matplotlib already imported above with Agg backend
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
from collections import defaultdict
from pybedtools import BedTool
from scipy import stats
from statsmodels.stats.multitest import multipletests
import re

# Import SEdbRegionAnalyzer
try:
    from .sedb_analyzer import SEdbRegionAnalyzer
except ImportError:
    # Fallback for direct execution
    import sys
    import os
    current_dir = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, current_dir)
    from sedb_analyzer import SEdbRegionAnalyzer

class ERNARegionAnalyzer(SEdbRegionAnalyzer):
    """
    eRNAbase genomic region analyzer.
    
    This class extends the SEdbRegionAnalyzer to analyze eRNAbase BED files 
    and their overlaps with specific genomic regions.
    """
    
    def __init__(self, results_dir='erna_results', logger=None):
        """
        Initialize the eRNA analyzer
        
        Parameters:
        results_dir (str): Directory to save results (default: 'erna_results')
        logger (Logger): Optional logger instance to use. If None, creates a basic console logger.
        """
        # 親クラスのコンストラクタを呼び出し
        super().__init__(results_dir=results_dir, logger=logger)
        
        self.logger.info("Initializing ERNARegionAnalyzer instance")
        
        # eRNAbase固有の属性を初期化
        # 内部属性には _ をプレフィックスとしてつける
        self._erna_bed_files = None
        self._erna_metadata = None
        self._interest_region = None
        self._interest_bed = None
        self._overlapping_results = None
        self._overlapping_samples = None
        self._overlap_sample_count = 0
        
        # 分布・エンリッチメント結果用の内部変数
        self._metadata_tissue_dist = None  # 全メタデータにおける組織型の分布
        self._metadata_cell_type_dist = None  # 全メタデータにおける細胞型の分布  
        self._region_tissue_counts_data = None  # 重複領域における組織型カウント
        self._region_cell_type_counts_data = None  # 重複領域における細胞型カウント
        self._tissue_enrichment_results_data = None  # 組織型のエンリッチメント解析結果
        self._cell_type_enrichment_results_data = None  # 細胞型のエンリッチメント解析結果
    
        
        # Default column definitions for BED files
        self.default_bed_columns = [
            "chrom", "start", "end", "peak_id", "fow_start2", "fow_end2", "size", 
            "risk_snp", "common_snp", "crisps", "eqtl", "tfbs", "enhancer", "super_enhancer"
        ]
        
        # Plot settings
        self._setup_plot_style()
        
        self.logger.info("ERNARegionAnalyzer initialization completed successfully")

    def _setup_plot_style(self):
        """Setup common style for plots"""
        plt.rcParams['figure.figsize'] = (14, 8)
        sns.set_style('whitegrid')
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
    
    def _normalize_chromosome(self, chrom: str) -> str:
        """
        Normalize chromosome name (e.g., convert 1 to chr1, Chr1 to chr1)
        
        Parameters:
        -----------
        chrom : str
            Input chromosome name
            
        Returns:
        --------
        str
            Normalized chromosome name
        """
        chrom = str(chrom).strip()
        
        # Handle chr prefix
        if not chrom.lower().startswith('chr'):
            if chrom.upper() in ['X', 'Y', 'M', 'MT'] or chrom.isdigit():
                chrom = f'chr{chrom}'
        
        # Standardize case (chr is lowercase, X/Y/M/MT are uppercase)
        if chrom.lower().startswith('chr'):
            prefix = 'chr'
            suffix = chrom[3:].upper() if chrom[3:].upper() in ['X', 'Y', 'M', 'MT'] else chrom[3:]
            chrom = prefix + suffix
        
        return chrom
        


    def load_erna_metadata(self, metadata_file: str, species: str = 'human', 
                        encoding: str = 'utf-8') -> 'ERNARegionAnalyzer':
        """
        Load eRNAbase metadata file (supports both CSV and Parquet formats)
        
        Parameters:
        -----------
        metadata_file : str
            Path to the metadata file (.csv or .parquet)
        species : str, optional
            Species filter option:
            - 'human': Only Homo sapiens data (default)
            - 'mouse': Only Mus musculus data
            - 'all': All species data
            - custom string: Exact match for Species column
        encoding : str, optional
            Encoding for CSV files (default: 'utf-8')
            
        Returns:
        --------
        ERNARegionAnalyzer
            self for method chaining
        """
        self.logger.info(f"Loading eRNAbase metadata from {metadata_file}")
        
        if not os.path.exists(metadata_file):
            self.logger.error(f"Metadata file not found: {metadata_file}")
            raise FileNotFoundError(f"File not found: {metadata_file}")
        
        try:
            # Determine file format based on extension
            file_ext = os.path.splitext(metadata_file)[1].lower()
            
            if file_ext == '.parquet':
                self.logger.info("Loading Parquet file...")
                metadata_df = pd.read_parquet(metadata_file)
            elif file_ext == '.csv':
                self.logger.info(f"Loading CSV file with encoding {encoding}...")
                metadata_df = pd.read_csv(metadata_file, encoding=encoding)
            else:
                self.logger.warning(f"Unrecognized file extension: {file_ext}. Attempting to load as CSV...")
                metadata_df = pd.read_csv(metadata_file, encoding=encoding)
            
            self.logger.info(f"Loaded metadata with {len(metadata_df)} samples and {len(metadata_df.columns)} columns")
            
            # Log column information
            self.logger.info(f"Metadata columns: {', '.join(metadata_df.columns)}")
            
            # Filter by species if requested and Species column exists
            original_count = len(metadata_df)
            if species and species != 'all' and 'Species' in metadata_df.columns:
                # Map common species names to their scientific names
                species_map = {
                    'human': 'Homo sapiens',
                    'mouse': 'Mus musculus',
                    # Additional mappings can be added here
                }
                
                # Get the scientific name or use the provided string directly
                scientific_name = species_map.get(species.lower(), species)
                
                # Apply filter
                metadata_df = metadata_df[metadata_df['Species'] == scientific_name]
                filtered_count = len(metadata_df)
                removed_count = original_count - filtered_count
                
                self.logger.info(f"Filtered for {scientific_name}: retained {filtered_count} samples " 
                            f"({filtered_count/original_count*100:.1f}%), "
                            f"removed {removed_count} samples")
                
                if filtered_count == 0:
                    self.logger.warning(f"No samples found for species '{scientific_name}' after filtering!")
            
            # Store the data
            self.erna_metadata = metadata_df
            
            # Log some basic statistics
            if not metadata_df.empty:
                if 'Species' in metadata_df.columns:
                    species_counts = metadata_df['Species'].value_counts()
                    for sp, count in species_counts.items():
                        self.logger.info(f"Species '{sp}': {count} samples ({count/len(metadata_df)*100:.1f}%)")
                    
                if 'Tissue Type' in metadata_df.columns:
                    tissue_count = metadata_df['Tissue Type'].nunique()
                    self.logger.info(f"Data contains {tissue_count} unique tissue types")
                    
                if 'Biosample Type' in metadata_df.columns:
                    biosample_count = metadata_df['Biosample Type'].nunique()
                    self.logger.info(f"Data contains {biosample_count} unique biosample types")
                    
            return self
            
        except Exception as e:
            self.logger.error(f"Error loading eRNAbase metadata: {e}")
            raise


    def get_bed_files(self, directory: str, pattern: str = "*.bed", 
                   show_samples: bool = True, sample_count: int = 5) -> List[str]:
        """
        Get a list of BED files from the specified directory
        
        Parameters:
        -----------
        directory : str
            Directory containing BED files
        pattern : str, optional
            File search pattern (default: "*.bed")
        show_samples : bool, optional
            Whether to display sample filenames (default: True)
        sample_count : int, optional
            Number of sample files to display (default: 5)
        
        Returns:
        --------
        List[str]
            List of BED file paths
        """
        # Get paths of BED files
        bed_files = glob.glob(os.path.join(directory, pattern))
        
        self.logger.info(f"Found {len(bed_files)} BED files in {directory}")

        # Natural sort the file paths using natsort
        bed_files = natsorted(bed_files)

        
        # Display sample filenames
        if show_samples and bed_files:
            self.logger.info("Sample files:")
            for file in bed_files[:sample_count]:
                self.logger.info(f"  - {os.path.basename(file)}")
        
        # Save to attribute
        self.erna_bed_files = bed_files
        
        return bed_files

    def create_interest_region(self, chrom: str, start: int, end: int, 
                             region_name: str | None = None, 
                             output_bed: str | None = None) -> Tuple[Dict, BedTool, str]:
        """
        Define a region of interest and save it as a BED file
        
        Parameters:
        -----------
        chrom : str
            Chromosome name (will be normalized)
        start : int
            Start position
        end : int
            End position
        region_name : str, optional
            Name for the region. Default is None, which generates automatically
        output_bed : str, optional
            Path for the output BED file. Default is None, which generates
            '{chrom}_{start}_{end}.bed'
        
        Returns:
        --------
        tuple
            (interest_region dict, BedTool object, output file path)
        """
        # Normalize chromosome name
        chrom = self._normalize_chromosome(chrom)
        
        # Validate coordinates
        if start < 0:
            raise ValueError(f"Start position must be non-negative, got {start}")
        if end <= start:
            raise ValueError(f"End position ({end}) must be greater than start ({start})")
        
        # Region size warning
        region_size = end - start
        if region_size > 10_000_000:  # 10Mb
            self.logger.warning(f"Large region size: {region_size:,} bp. This may affect performance.")
        
        # Define region of interest as a dictionary
        interest_region = {
            'chrom': chrom,
            'start': start,
            'end': end
        }
        
        # Generate region name if not specified
        if region_name is None:
            region_name = f"Region_{chrom}_{start}_{end}"
        
        # Generate output filename if not specified
        if output_bed is None:
            output_bed = os.path.join(self.results_dir, f"{chrom}_{start}_{end}.bed")
        
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_bed), exist_ok=True)
        
        # Save as BED file
        with open(output_bed, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\t{region_name}\n")
        
        # Create BedTool object
        interest_bed = BedTool(output_bed)
        
        self.logger.info(f"Created interest region: {chrom}:{start:,}-{end:,}")
        self.logger.info(f"Saved BED file: {output_bed}")
        
        # Save to attributes
        self.interest_region = interest_region
        self.interest_bed = interest_bed
        
        return interest_region, interest_bed, output_bed


    def analyze_file_overlap(self, bed_file: str, columns: List[str] | None = None) -> pd.DataFrame:
        """
        Analyze overlaps between a BED file and the region of interest
        
        Parameters:
        -----------
        bed_file : str
            BED file to analyze
        columns : List[str], optional
            Column names for the BED file. Uses default if None
        
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing overlapping records, or empty DataFrame if no overlaps
        """
        if self.interest_bed is None:
            self.logger.error("Interest region not defined. Call create_interest_region first.")
            return pd.DataFrame()
        
        # Set column names
        if columns is None:
            columns = self.default_bed_columns
        
        try:
            # Load the file
            file_bed = BedTool(bed_file)
            
            # Detect overlaps
            overlaps = file_bed.intersect(self.interest_bed, wa=True)
            
            # Process results if overlaps exist
            if len(overlaps) > 0:
                overlap_df = pd.read_table(overlaps.fn, header=None, names=columns)
                
                # Double-check chromosome consistency (defensive programming)
                if hasattr(self, 'interest_region') and self.interest_region:
                    expected_chrom = self.interest_region['chrom']
                    if 'chrom' in overlap_df.columns:
                        unexpected_chroms = overlap_df[overlap_df['chrom'] != expected_chrom]
                        if len(unexpected_chroms) > 0:
                            self.logger.warning(
                                f"Unexpected chromosomes found in {os.path.basename(bed_file)}: "
                                f"{unexpected_chroms['chrom'].unique()}. Filtering to {expected_chrom} only."
                            )
                            overlap_df = overlap_df[overlap_df['chrom'] == expected_chrom]
                
                overlap_df['source_file'] = os.path.splitext(os.path.basename(bed_file))[0]
                return overlap_df
            else:
                return pd.DataFrame(columns=columns + ['source_file'])
        except Exception as e:
            self.logger.error(f"Error ({os.path.basename(bed_file)}): {e}")
            return pd.DataFrame(columns=columns + ['source_file'])

    def merge_all_overlaps(self, columns: List[str] | None = None, 
                         max_files: int | None = None) -> Tuple[pd.DataFrame, List[str], int]:
        """
        Merge overlaps from all BED files
        
        Parameters:
        -----------
        columns : List[str], optional
            Column names for the BED files. Uses default if None
        max_files : int, optional
            Maximum number of files to analyze (for debugging)
        
        Returns:
        --------
        Tuple[pd.DataFrame, List[str], int]
            (merged overlapping records DataFrame, list of overlapping sample names, count of overlapping samples)
        """
        if self.erna_bed_files is None:
            self.logger.error("BED files not loaded. Call get_bed_files first.")
            return pd.DataFrame(), [], 0
            
        if self.interest_bed is None:
            self.logger.error("Interest region not defined. Call create_interest_region first.")
            return pd.DataFrame(), [], 0
        
        # Set column names
        if columns is None:
            columns = self.default_bed_columns
        
        results = []
        overlapping_samples = []  # List to store names of samples with overlaps
        
        # Limit number of files to analyze (for debugging)
        if max_files:
            files_to_analyze = self.erna_bed_files[:max_files]
        else:
            files_to_analyze = self.erna_bed_files
        
        total_files = len(files_to_analyze)
        
        # Start time for processing
        start_time = time.time()
        
        for i, file in enumerate(files_to_analyze):
            if i % 50 == 0:
                elapsed = time.time() - start_time
                self.logger.info(f"Processing: {i}/{total_files} files (elapsed: {elapsed:.2f}s)")
            
            result = self.analyze_file_overlap(file, columns)
            if not result.empty:
                results.append(result)
                
                # Extract sample name from file path
                sample_name = os.path.basename(file).replace('.bed', '')
                overlapping_samples.append(sample_name)
        
        # Show summary information
        overlap_sample_count = len(overlapping_samples)
        elapsed_total = time.time() - start_time
        self.logger.info(f"Completed analysis in {elapsed_total:.2f}s")
        self.logger.info(f"Overlapping samples found: {overlap_sample_count}/{total_files}")
        
        # Merge results
        if results:
            combined = pd.concat(results, ignore_index=True)
            
            # Save results to attributes
            self.overlapping_results = combined
            self.overlapping_samples = overlapping_samples
            self.overlap_sample_count = overlap_sample_count
            
            return combined, overlapping_samples, overlap_sample_count
        else:
            # Return empty DataFrame with same column definition
            empty_df = pd.DataFrame(columns=columns + ['source_file'])
            
            # Save empty results to attributes
            self.overlapping_results = empty_df
            self.overlapping_samples = []
            self.overlap_sample_count = 0
            
            return empty_df, [], 0

    def save_overlap_results(self, output_dir: Optional[str] = None, 
                          display_samples: bool = True, sample_count: int = 5) -> Dict[str, str]:
        """
        Save overlap analysis results and display sample data
        
        Parameters:
        -----------
        output_dir : str, optional
            Directory to save output files. Uses results_dir if None
        display_samples : bool, optional
            Whether to display sample data (default: True)
        sample_count : int, optional
            Number of samples to display (default: 5)
        
        Returns:
        --------
        Dict[str, str]
            Dictionary containing paths of generated files
        """
        if self.overlapping_results is None or self.interest_region is None:
            self.logger.error("No overlap results or interest region. Run merge_all_overlaps first.")
            return {}
        
        # Reference results to local variables
        merged_overlaps = self.overlapping_results
        overlapping_samples = self.overlapping_samples
        interest_region = self.interest_region
        
        # Skip if empty DataFrame
        if merged_overlaps.empty:
            self.logger.warning("No overlaps found. Skipping result saving.")
            self.logger.info("Try a different interest region.")
            return {}
        
        # Set output directory
        if output_dir is None:
            output_dir = self.results_dir  # Default is results directory
        
        # Create directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            self.logger.info(f"Created output directory: {output_dir}")
        
        # Filename prefix
        prefix = f"{interest_region['chrom']}_{interest_region['start']}_{interest_region['end']}"
        
        # Display results summary
        self.logger.info(f"\nTotal {len(merged_overlaps)} overlapping records found")
        self.logger.info(f"Overlapping samples: {len(overlapping_samples)}")
        
        # Dictionary to store file paths
        saved_files = {}
        
        # Save overlapping sample names to text file
        overlap_samples_path = os.path.join(output_dir, f"{prefix}_overlapping_samples.txt")
        with open(overlap_samples_path, 'w') as f:
            for sample_name in overlapping_samples:
                f.write(f"{sample_name}\n")
        self.logger.info(f"Saved overlapping sample names to {overlap_samples_path}")
        saved_files['samples_txt'] = overlap_samples_path
        
        # Also save sample names as CSV
        samples_df = pd.DataFrame({'sample_name': overlapping_samples})
        samples_csv = os.path.join(output_dir, f"{prefix}_overlapping_samples.csv")
        samples_df.to_csv(samples_csv, index=False)
        self.logger.info(f"Saved overlapping sample names as CSV to {samples_csv}")
        saved_files['samples_csv'] = samples_csv
        
        # Save overlapping records to CSV
        output_file = os.path.join(output_dir, f"{prefix}_overlaps.csv")
        merged_overlaps.to_csv(output_file, index=False)
        self.logger.info(f"Saved overlapping records to {output_file}")
        saved_files['overlaps_csv'] = output_file
        
        # Display samples
        if display_samples and not merged_overlaps.empty:
            # Show examples of overlapping records
            print(f"\nOverlapping records examples (top {sample_count}):")
            print(merged_overlaps.head(sample_count))
            
            # Also show examples of overlapping samples
            if overlapping_samples:
                print(f"\nOverlapping samples examples (first {sample_count}):")
                for sample in overlapping_samples[:sample_count]:
                    print(f"  - {sample}")
        
        return saved_files

    def analyze_overlapping_samples_with_metadata(self) -> pd.DataFrame:
        """
        Analyze metadata for samples with overlaps
        
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing metadata for overlapping samples
        """
        if self.overlapping_samples is None or self.erna_metadata is None:
            self.logger.error("Overlapping samples or metadata not loaded. Run merge_all_overlaps and load_erna_metadata first.")
            return pd.DataFrame()
        
        # Check sample ID column in metadata
        sample_id_col = 'Sample ID'  # Default value
        
        # Check if sample ID column exists
        if sample_id_col not in self.erna_metadata.columns:
            # Check alternative column names
            possible_columns = ['sample_id', 'SampleID', 'Sample', 'sample_name', 'SampleName']
            for col in possible_columns:
                if col in self.erna_metadata.columns:
                    sample_id_col = col
                    self.logger.info(f"Using '{col}' as sample ID column")
                    break
            else:
                self.logger.error(f"Error: '{sample_id_col}' column not found in metadata")
                return pd.DataFrame()
        
        # Extract metadata for overlapping samples
        overlapping_df = self.erna_metadata[self.erna_metadata[sample_id_col].isin(self.overlapping_samples)].copy()
        
        self.logger.info(f"Found metadata for {len(overlapping_df)} out of {len(self.overlapping_samples)} overlapping samples")
        
        if len(overlapping_df) == 0:
            self.logger.warning("No metadata found for overlapping samples. Analysis stopped.")
            return overlapping_df
        
        # Display metadata summary
        print("\n========== Overlapping Samples Metadata Analysis ==========")
        
        # Save results
        self.overlapping_metadata = overlapping_df
        
        return overlapping_df

    def analyze_erna_complete(self, output_plots=True) -> pd.DataFrame:
        """
        Perform complete analysis of eRNAbase data
        
        Parameters:
        -----------
        output_plots : bool, optional
            Whether to output plots (default: True)
            
        Returns:
        --------
        pandas.DataFrame
            DataFrame used for analysis
        """
        if self.overlapping_metadata is None:
            if self.overlapping_samples is not None and self.erna_metadata is not None:
                self.analyze_overlapping_samples_with_metadata()
                
            if self.overlapping_metadata is None:
                self.logger.error("No overlapping metadata available for analysis")
                return pd.DataFrame()
        
        # DataFrame for analysis
        df = self.overlapping_metadata
        
        print("========== eRNAbase Metadata Analysis ==========")
        
        # 1. Basic information about the data
        print(f"Data size: {df.shape[0]} rows x {df.shape[1]} columns")
        print(f"Columns: {', '.join(df.columns)}")
        
        # Check species distribution
        if 'Species' in df.columns:
            species_counts = df['Species'].value_counts()
            print("\nSpecies distribution:")
            for species, count in species_counts.items():
                print(f"  {species}: {count} samples ({count/len(df)*100:.1f}%)")
        
        # Check for missing values
        missing_values = df.isnull().sum()
        if missing_values.sum() > 0:
            print("\nMissing values:")
            for col, count in missing_values.items():
                if count > 0:
                    print(f"  {col}: {count} values ({count/len(df)*100:.1f}%)")
        else:
            print("\nNo missing values")
        
        # 2. Visualize distributions with pie charts
        if output_plots:
            print("\n========== Category Distribution Pie Charts ==========")
            categories = ['Tissue Type', 'Biosample Type', 'Cell Type', 'Experiment Type']
            
            for category in categories:
                if category in df.columns:
                    self.plot_simple_pie(df, category)
        
        # 3. Analyze relationships with cross-tabulation
        if output_plots:
            print("\n========== Category Relationship Analysis ==========")
            combinations = [
                ('Tissue Type', 'Biosample Type'),
                ('Tissue Type', 'Cell Type'),
                ('Biosample Type', 'Cell Type')
            ]
            
            for cat1, cat2 in combinations:
                if cat1 in df.columns and cat2 in df.columns:
                    self.plot_cross_heatmap(df, cat1, cat2)
        
        # 4. Analysis summary and suggestions
        print("\n========== Analysis Summary and Suggestions ==========")
        print("Additional analyses you can perform:")
        print("1. Focus on specific tissue type")
        print("   Example: analyzer.analyze_focus('Tissue Type', 'Bone marrow')")
        print("2. Focus on specific cell type")
        print("   Example: analyzer.analyze_focus('Cell Type', 'HeLa')")
        print("3. Focus on specific biosample type")
        print("   Example: analyzer.analyze_focus('Biosample Type', 'Cell line')")
        
        print("\nAnalysis completed. Run the suggested commands for more detailed analysis.")
        
        return df


    def plot_simple_pie(self, data: pd.DataFrame, category: str, 
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
            DataFrame to analyze
        category : str
            Category column name to display in pie chart
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
        # Count category distribution
        if category not in data.columns:
            self.logger.warning(f"Column '{category}' does not exist in the dataframe")
            return None
            
        value_counts = data[category].value_counts()
        total_categories = len(value_counts)
        
        # Set title
        if title is None:
            title = f'{category} Distribution'
            
        # Divide into top N and others
        if total_categories > limit:
            top_counts = value_counts.head(limit)
            others_count = value_counts.iloc[limit:].sum()
            plot_data = pd.Series({**top_counts.to_dict(), 'Others': others_count})
            title_suffix = f" (Top {limit}, other {total_categories-limit} types grouped as 'Others')"
        else:
            plot_data = value_counts
            title_suffix = ""
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # Generate distinct colors for pie segments
        cmap = plt.cm.get_cmap(colormap)
        num_colors = len(plot_data)
        colors = [cmap(i % cmap.N) for i in range(num_colors)]
        
        # Create pie chart (maintaining original design)
        plt.pie(
            plot_data, 
            labels=plot_data.index,
            autopct='%1.1f%%',
            startangle=90,
            shadow=False,
            wedgeprops={'edgecolor': 'none'},  # Maintain original edge style
            textprops={'fontsize': 10},
            colors=colors  # Use our custom color list
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
        
        plt.show()
        
        # Display detailed information
        print(f"\n{category} distribution:")
        print(f"Total {total_categories} types, showing top {min(limit, total_categories)}")
        
        for i, (value, count) in enumerate(value_counts.head(limit).items()):
            print(f"  {value}: {count} samples ({count/len(data)*100:.1f}%)")
        
        if total_categories > limit:
            remaining_count = value_counts.iloc[limit:].sum()
            print(f"  Others ({total_categories - limit} types): {remaining_count} samples ({remaining_count/len(data)*100:.1f}%)")
        
        return fig


    def plot_simple_bar(self, data: pd.DataFrame, category: str, 
                       title: Optional[str] = None, limit: int = 20) -> None:
        """
        Create a simple bar plot
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame to analyze
        category : str
            Category column name to display in bar plot
        title : str, optional
            Chart title
        limit : int, default=20
            Number of top categories to display
        """
        if category not in data.columns:
            self.logger.warning(f"Column '{category}' does not exist in the dataframe")
            return
        
        # Count category distribution
        value_counts = data[category].value_counts().head(limit)
        
        # Set title
        if title is None:
            title = f'Top {limit} {category}'
        
        # Create bar plot
        plt.figure(figsize=(14, 8))
        sns.barplot(x=value_counts.values, y=value_counts.index)
        plt.title(title)
        plt.xlabel('Number of Samples')
        plt.tight_layout()
        plt.show()

    def plot_cross_heatmap(self, data: pd.DataFrame, cat1: str, cat2: str, 
                         title: Optional[str] = None, limit: int = 15) -> None:
        """
        Display cross-tabulation of two categories as a heatmap
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame to analyze
        cat1 : str
            Column name for category to use for rows
        cat2 : str
            Column name for category to use for columns
        title : str, optional
            Chart title
        limit : int, default=15
            Maximum number of top categories to display for each axis
        """
        if cat1 not in data.columns or cat2 not in data.columns:
            self.logger.warning(f"Column '{cat1}' or '{cat2}' does not exist in the dataframe")
            return
        
        # Cross-tabulation
        cross_tab = pd.crosstab(data[cat1], data[cat2])
        
        # For large cross-tabulation tables, display only top entries
        if cross_tab.shape[0] > limit or cross_tab.shape[1] > limit:
            # Calculate row and column totals
            row_sums = cross_tab.sum(axis=1)
            col_sums = cross_tab.sum(axis=0)
            
            # Extract top rows and columns
            top_rows = row_sums.nlargest(limit).index
            top_cols = col_sums.nlargest(limit).index
            
            # Display only top rows/columns
            heat_data = cross_tab.loc[top_rows, top_cols]
            title_suffix = f" (Top {limit}x{limit} only)"
        else:
            heat_data = cross_tab
            title_suffix = ""
        
        # Set title
        if title is None:
            title = f'Relationship between {cat1} and {cat2}'
        
        # Create heatmap
        plt.figure(figsize=(16, 12))
        sns.heatmap(heat_data, annot=True, cmap='YlGnBu', fmt='d', cbar_kws={'label': 'Number of Samples'})
        plt.title(f'{title}{title_suffix}')
        plt.tight_layout()
        plt.show()

    def analyze_focus(self, category: str, value: str) -> None:
        """
        Perform detailed analysis focused on a specific category value
        
        Parameters:
        -----------
        category : str
            Category column name to filter on
        value : str
            Category value to focus on
        """
        if self.overlapping_metadata is None:
            self.logger.error("No overlapping metadata available. Run analyze_overlapping_samples_with_metadata first.")
            return
        
        data = self.overlapping_metadata
        
        if category not in data.columns:
            self.logger.warning(f"Column '{category}' does not exist in the dataframe")
            return
        
        # Extract samples matching the category
        subset = data[data[category] == value]
        
        if len(subset) == 0:
            self.logger.warning(f"No samples found with '{category}' = '{value}'")
            return
        
        print(f"\n===== {category} = '{value}' Detailed Analysis ({len(subset)} samples) =====")
        
        # Analyze distribution of other categories
        other_categories = ['Tissue Type', 'Biosample Type', 'Cell Type', 'Experiment Type']
        other_categories = [cat for cat in other_categories if cat != category and cat in data.columns]
        
        # Display each category in pie charts
        for cat in other_categories:
            self.plot_simple_pie(subset, cat, title=f'{cat} Distribution for {value}')

    def analyze_region(self, chrom: str, start: int, end: int, region_name: str | None = None,
                     metadata_file: str | None = None, bed_directory: str | None = None,
                     bed_pattern: str = "*.bed", output_dir: str | None = None,
                     full_analysis: bool = True) -> Dict[str, Any]:
        """
        Perform complete region analysis workflow
        
        Parameters:
        -----------
        chrom : str
            Chromosome name (e.g., 'chr7')
        start : int
            Start position
        end : int
            End position
        region_name : str, optional
            Region name (can be auto-generated)
        metadata_file : str, optional
            Path to metadata file (if not already loaded)
        bed_directory : str, optional
            Directory with BED files (if not already loaded)
        bed_pattern : str, optional
            BED file pattern (default: "*.bed")
        output_dir : str, optional
            Output directory (default: self.results_dir)
        full_analysis : bool, optional
            Whether to perform full analysis (default: True)
            
        Returns:
        --------
        Dict[str, Any]
            Dictionary of analysis results
        """
        results = {}
        
        # Load metadata if needed
        if self.erna_metadata is None and metadata_file is not None:
            self.load_erna_metadata(metadata_file)
        
        # Load BED files if needed
        if self.erna_bed_files is None and bed_directory is not None:
            self.get_bed_files(bed_directory, pattern=bed_pattern)
        
        # Verify settings
        if self.erna_bed_files is None:
            self.logger.error("BED files not loaded. Provide bed_directory or call get_bed_files first.")
            return results
        
        # Create region of interest
        interest_region, interest_bed, output_bed = self.create_interest_region(
            chrom, start, end, region_name=region_name
        )
        results['interest_region'] = interest_region
        results['output_bed'] = output_bed
        
        # Analysis start time
        start_time = time.time()
        self.logger.info(f"Starting analysis for region {chrom}:{start}-{end}")
        
        # Analyze overlaps in all files
        merged_overlaps, overlapping_samples, overlap_sample_count = self.merge_all_overlaps()
        results['overlapping_samples_count'] = overlap_sample_count
        results['overlapping_samples'] = overlapping_samples
        
        # Save results
        if not merged_overlaps.empty:
            saved_files = self.save_overlap_results(output_dir=output_dir)
            results['saved_files'] = saved_files
            
            # Analyze metadata (if available and samples overlap)
            if self.erna_metadata is not None and overlapping_samples:
                overlapping_metadata = self.analyze_overlapping_samples_with_metadata()
                
                # Perform full analysis
                if full_analysis and not overlapping_metadata.empty:
                    analysis_df = self.analyze_erna_complete()
                    results['analysis_df'] = analysis_df
        
        # Elapsed time
        elapsed_time = time.time() - start_time
        self.logger.info(f"Analysis completed in {elapsed_time:.2f} seconds")
        
        # Display main analysis results
        print("\n========== Analysis Summary ==========")
        print(f"Region: {chrom}:{start:,}-{end:,}")
        print(f"Total BED files analyzed: {len(self.erna_bed_files)}")
        print(f"Overlapping samples found: {overlap_sample_count}")
        
        if self.overlapping_metadata is not None:
            print(f"Samples with metadata: {len(self.overlapping_metadata)}")
        
        print(f"Total analysis time: {elapsed_time:.2f} seconds")
        
        # Also save results to attribute
        self.analysis_results = results
        
        return results

    # Example of using SEdbRegionAnalyzer methods
    def compare_erna_sedb_enrichment(self, chrom: str, start: int, end: int, 
                                   region_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Compare enrichment results between eRNAbase and SEdb for the same region
        
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
            
        Returns:
        --------
        Dict[str, Any]
            Dictionary of comparison results
        """
        # eRNAbase analysis results (assuming already run)
        if self.overlapping_metadata is None:
            self.logger.error("eRNAbase analysis results not available. Run analyze_region first.")
            return {}
            
        # SEdb analysis (using parent class methods)
        print("\n========== Analyzing the same region with SEdb ==========")
        
        # Check if SEdb data is loaded
        if self.bed_df is None:
            self.logger.warning("SEdb data not loaded. Load databases before comparison.")
            return {}
            
        # Analyze the same region using parent class methods
        print(f"Extracting SEs from region {chrom}:{start}-{end}")
        self.extract_overlapping_se(chrom, start, end, region_name=region_name)
        
        # Tissue enrichment analysis
        if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None:
            print("\nPerforming tissue enrichment analysis with SEdb data")
            sedb_enrichment = self.analyze_region_tissue_enrichment(plot=True)
            
            # Save results to attribute
            self.sedb_enrichment_results = sedb_enrichment
            
            # Display comparison of eRNAbase and SEdb results
            print("\n========== Comparison of eRNAbase and SEdb Results ==========")
            print(f"eRNAbase samples with region overlap: {len(self.overlapping_samples)}")
            print(f"SEdb samples with region overlap: {getattr(self, 'unique_sample_count', 0)}")
            
            # Create results dictionary
            comparison_results = {
                'region': {
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'name': region_name or f"{chrom}:{start}-{end}"
                },
                'erna_samples_count': len(self.overlapping_samples),
                'sedb_samples_count': getattr(self, 'unique_sample_count', 0),
                'erna_tissue_data': self.overlapping_metadata,
                'sedb_tissue_data': self.sample_summary if hasattr(self, 'sample_summary') else None,
                'sedb_enrichment': sedb_enrichment
            }
            
            return comparison_results
        else:
            self.logger.warning("SEdb analysis did not produce tissue counts. Cannot compare enrichment.")
            return {}

    # Additional utility methods
    def export_results_to_excel(self, filename: Optional[str] = None) -> str:
        """
        Export all analysis results to an Excel file
        
        Parameters:
        -----------
        filename : str, optional
            Output filename (default: "erna_analysis_results_{timestamp}.xlsx")
            
        Returns:
        --------
        str
            Path to the exported file
        """
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            self.logger.warning("No analysis results available to export.")
            return ""
            
        # Generate filename
        if filename is None:
            import datetime
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = os.path.join(self.results_dir, f"erna_analysis_results_{timestamp}.xlsx")
        
        # Ensure output directory exists
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        
        # Write to Excel file
        with pd.ExcelWriter(filename) as writer:
            # Basic information sheet
            info_df = pd.DataFrame([{
                'Analysis Date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'Region': f"{self.interest_region['chrom']}:{self.interest_region['start']}-{self.interest_region['end']}",
                'Overlapping Samples': self.overlap_sample_count,
                'Total BED Files': len(self.erna_bed_files) if self.erna_bed_files else 0
            }])
            info_df.to_excel(writer, sheet_name='Analysis_Info', index=False)
            
            # Overlapping samples sheet
            if self.overlapping_samples:
                samples_df = pd.DataFrame({'Sample': self.overlapping_samples})
                samples_df.to_excel(writer, sheet_name='Overlapping_Samples', index=False)
            
            # Overlapping records sheet
            if self.overlapping_results is not None and not self.overlapping_results.empty:
                self.overlapping_results.to_excel(writer, sheet_name='Overlapping_Records', index=False)
            
            # Metadata sheet
            if self.overlapping_metadata is not None and not self.overlapping_metadata.empty:
                self.overlapping_metadata.to_excel(writer, sheet_name='Sample_Metadata', index=False)
        
        self.logger.info(f"Analysis results exported to {filename}")
        return filename
    


    def calculate_metadata_distributions(self):
        """
        Calculate the distribution of Tissue Type and Cell Type in the entire metadata.
        Sets the tissue_distribution and cell_type_distribution properties.
        
        Returns:
        --------
        bool
            True if calculation was successful
        """
        if self.erna_metadata is None:
            self.logger.error("Metadata not loaded. Call load_erna_metadata first.")
            return False
        
        # Calculate Tissue Type distribution
        if 'Tissue Type' in self.erna_metadata.columns:
            tissue_counts = self.erna_metadata['Tissue Type'].value_counts()
            tissue_percentages = (tissue_counts / len(self.erna_metadata) * 100).round(2)
            
            # 修正: 内部変数に保存
            self._metadata_tissue_dist = pd.DataFrame({
                'Tissue Type': tissue_counts.index,
                'Count': tissue_counts.values,
                'Percentage (%)': tissue_percentages.values
            })
            
            self.logger.info(f"Calculated distribution for {len(tissue_counts)} tissue types")
        else:
            self.logger.warning("'Tissue Type' column not found in metadata")
        
        # Calculate Cell Type distribution
        if 'Cell Type' in self.erna_metadata.columns:
            cell_type_counts = self.erna_metadata['Cell Type'].value_counts()
            cell_type_percentages = (cell_type_counts / len(self.erna_metadata) * 100).round(2)
            
            # 修正: 内部変数に保存
            self._metadata_cell_type_dist = pd.DataFrame({
                'Cell Type': cell_type_counts.index,
                'Count': cell_type_counts.values,
                'Percentage (%)': cell_type_percentages.values
            })
            
            self.logger.info(f"Calculated distribution for {len(cell_type_counts)} cell types")
        else:
            self.logger.warning("'Cell Type' column not found in metadata")
        
        return True



    def analyze_overlapping_sample_distributions(self):
        """
        Analyze the distribution of Tissue Type and Cell Type in overlapping samples.
        Sets the region_tissue_counts and region_cell_type_counts properties.
        
        Returns:
        --------
        bool
            True if analysis was successful
        """
        if self.overlapping_metadata is None:
            if self.overlapping_samples is not None and self.erna_metadata is not None:
                self.analyze_overlapping_samples_with_metadata()
                
            if self.overlapping_metadata is None:
                self.logger.error("No overlapping metadata available for analysis")
                return False
        
        # 重要: df 変数を定義（これが欠けていたことが原因）
        df = self.overlapping_metadata
        
        # Calculate Tissue Type counts in overlapping samples
        if 'Tissue Type' in df.columns:
            tissue_counts = df['Tissue Type'].value_counts()
            tissue_percentages = (tissue_counts / len(df) * 100).round(2)
            
            # region_tissue_counts に直接設定
            self.region_tissue_counts = pd.DataFrame({
                'Tissue Type': tissue_counts.index,
                'Region Count': tissue_counts.values,
                'Region Percentage (%)': tissue_percentages.values
            })
            
            # 内部変数にもコピーしておく
            self._region_tissue_counts_data = self.region_tissue_counts.copy()
            
            self.logger.info(f"Found {len(tissue_counts)} tissue types in overlapping samples")
        else:
            self.logger.warning("'Tissue Type' column not found in overlapping metadata")
        
        # Calculate Cell Type counts in overlapping samples
        if 'Cell Type' in df.columns:
            cell_type_counts = df['Cell Type'].value_counts()
            cell_type_percentages = (cell_type_counts / len(df) * 100).round(2)
            
            self.region_cell_type_counts = pd.DataFrame({
                'Cell Type': cell_type_counts.index,
                'Region Count': cell_type_counts.values,
                'Region Percentage (%)': cell_type_percentages.values
            })
            
            # 内部変数にもコピー
            self._region_cell_type_counts_data = self.region_cell_type_counts.copy()
            
            self.logger.info(f"Found {len(cell_type_counts)} cell types in overlapping samples")
        else:
            self.logger.warning("'Cell Type' column not found in overlapping metadata")
        
        return True




    def analyze_tissue_enrichment(self, correction='fdr_bh', plot=True, **kwargs):
        """
        Perform tissue type enrichment analysis using Fisher's exact test.
        
        Parameters:
        -----------
        correction : str, optional
            Multiple testing correction method (default: 'fdr_bh')
        plot : bool, optional
            Whether to plot the results (default: True)
        **kwargs : dict
            Additional parameters for the visualization function
            
        Returns:
        --------
        pandas.DataFrame
            Enrichment analysis results
        """
        # Ensure we have the necessary data
        if not hasattr(self, '_metadata_tissue_dist') or self._metadata_tissue_dist is None:
            self.logger.info("Tissue distribution not calculated. Calculating now...")
            self.calculate_metadata_distributions()
            
        if self.region_tissue_counts is None:
            self.logger.info("Region tissue counts not calculated. Calculating now...")
            self.analyze_overlapping_sample_distributions()
        
        # Check that both the background and region distributions are available
        if self._metadata_tissue_dist is None or self.region_tissue_counts is None:
            self.logger.error("Could not calculate necessary distributions for enrichment analysis")
            return None
        
        # Get background distribution - use internal variable
        background = self._metadata_tissue_dist.set_index('Tissue Type')['Count']
        total_background = background.sum()
        
        # Get region distribution - use parent class property
        region_counts = self.region_tissue_counts.set_index('Tissue Type')['Region Count']
        total_region = region_counts.sum()
            
        # Prepare results container
        results = []
        
        # Perform Fisher's exact test for each tissue type
        from scipy import stats
        from statsmodels.stats.multitest import multipletests
        
        for tissue in background.index:
            # Count in the region (0 if not found)
            a = region_counts.get(tissue, 0)
            
            # Region total minus tissue count
            b = total_region - a
            
            # Background count for this tissue
            c = background[tissue]
            
            # Background total minus this tissue count
            d = total_background - c
            
            # Create 2x2 contingency table
            contingency = np.array([[a, b], [c, d]])
            
            # Perform Fisher's exact test
            odds_ratio, p_value = stats.fisher_exact(contingency)
            
            # Calculate expected value and fold change
            expected = (c / total_background) * total_region
            fold_change = (a / total_region) / (c / total_background) if c > 0 else float('nan')
            
            # Store results
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
        
        # Convert to DataFrame
        result_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        if len(result_df) > 1:
            reject, p_corr, _, _ = multipletests(result_df['p-value'].values, method=correction)
            result_df['adjusted p-value'] = p_corr
            result_df['Significant'] = reject
        
        # Sort by p-value
        result_df = result_df.sort_values('p-value')
        
        # Store results
        self.region_tissue_enrichment = result_df
        
        # Print summary of enrichment results
        sig_count = sum(result_df['p-value'] < 0.05)
        self.logger.info(f"Tissue enrichment analysis: {sig_count}/{len(result_df)} significant at p<0.05")
        
        # Top significant tissues
        if sig_count > 0:
            top_sig = result_df[result_df['p-value'] < 0.05].sort_values('p-value').head(5)
            self.logger.info("Top enriched tissues:")
            for _, row in top_sig.iterrows():
                self.logger.info(f"  {row['Tissue Type']}: {row['Fold Change']:.2f}x (p={row['p-value']:.1e})")
        
        # Visualize if requested
        if plot:
            self.visualize_tissue_enrichment(
                enrichment_results=result_df,
                plot_type=kwargs.get('plot_type', 'volcano'),
                title=kwargs.get('title', None),
                figsize=kwargs.get('figsize', (14, 8)),
                top_n=kwargs.get('top_n', 15),
                pvalue_threshold=kwargs.get('pvalue_threshold', 0.05),
                fc_threshold=kwargs.get('fc_threshold', 1.5),
                region_name=kwargs.get('region_name', None)
            )
        
        return result_df



    def analyze_cell_type_enrichment(self, correction='fdr_bh', plot=True, **kwargs):
        """
        Perform cell type enrichment analysis using Fisher's exact test.
        
        Parameters:
        -----------
        correction : str, optional
            Multiple testing correction method (default: 'fdr_bh')
        plot : bool, optional
            Whether to plot the results (default: True)
        **kwargs : dict
            Additional parameters for the visualization function
            
        Returns:
        --------
        pandas.DataFrame
            Enrichment analysis results
        """
        # Ensure we have the necessary data
        if self.metadata_cell_type_distribution is None:
            self.logger.info("Cell type distribution not calculated. Calculating now...")
            self.calculate_metadata_distributions()
            
        if self.region_cell_type_counts is None:
            self.logger.info("Region cell type counts not calculated. Calculating now...")
            self.analyze_overlapping_sample_distributions()
        
        if self.metadata_cell_type_distribution is None or self.region_cell_type_counts is None:
            self.logger.error("Could not calculate necessary distributions for enrichment analysis")
            return None
        
        # Get background distribution
        background = self.metadata_cell_type_distribution.set_index('Cell Type')['Count']
        total_background = background.sum()
        
        # Get region distribution
        region_counts = self.region_cell_type_counts.set_index('Cell Type')['Region Count']
        total_region = region_counts.sum()
        
        # Prepare results container
        results = []
        
        # Perform Fisher's exact test for each cell type
        from scipy import stats
        from statsmodels.stats.multitest import multipletests
        
        for cell_type in background.index:
            # Count in the region (0 if not found)
            a = region_counts.get(cell_type, 0)
            
            # Region total minus cell type count
            b = total_region - a
            
            # Background count for this cell type
            c = background[cell_type]
            
            # Background total minus this cell type count
            d = total_background - c
            
            # Create 2x2 contingency table
            contingency = np.array([[a, b], [c, d]])
            
            # Perform Fisher's exact test
            odds_ratio, p_value = stats.fisher_exact(contingency)
            
            # Calculate expected value and fold change
            expected = (c / total_background) * total_region
            fold_change = (a / total_region) / (c / total_background) if c > 0 else float('nan')
            
            # Store results
            results.append({
                'Cell Type': cell_type,
                'Region Count': a,
                'Region Percentage': (a / total_region * 100) if total_region > 0 else 0,
                'Total Count': c,
                'Total Percentage': (c / total_background * 100),
                'Expected Value': expected,
                'Odds Ratio': odds_ratio,
                'Fold Change': fold_change,
                'p-value': p_value
            })
        
        # Convert to DataFrame
        result_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        if len(result_df) > 1:
            reject, p_corr, _, _ = multipletests(result_df['p-value'].values, method=correction)
            result_df['adjusted p-value'] = p_corr
            result_df['Significant'] = reject
        
        # Sort by p-value
        result_df = result_df.sort_values('p-value')
        
        # Store results
        self.region_cell_type_enrichment = result_df
        
        # Print summary of enrichment results
        sig_count = sum(result_df['p-value'] < 0.05)
        self.logger.info(f"Cell type enrichment analysis: {sig_count}/{len(result_df)} significant at p<0.05")
        
        # Top significant cell types
        if sig_count > 0:
            top_sig = result_df[result_df['p-value'] < 0.05].sort_values('p-value').head(5)
            self.logger.info("Top enriched cell types:")
            for _, row in top_sig.iterrows():
                self.logger.info(f"  {row['Cell Type']}: {row['Fold Change']:.2f}x (p={row['p-value']:.1e})")
        
        # Visualize if requested
        if plot:
            # Use the same visualization method as for tissues, with adjusted parameters
            self.visualize_tissue_enrichment(
                enrichment_results=result_df,
                plot_type=kwargs.get('plot_type', 'volcano'),
                title=kwargs.get('title', "Cell Type Enrichment Analysis"),
                figsize=kwargs.get('figsize', (14, 8)),
                top_n=kwargs.get('top_n', 15),
                pvalue_threshold=kwargs.get('pvalue_threshold', 0.05),
                fc_threshold=kwargs.get('fc_threshold', 1.5),
                region_name=kwargs.get('region_name', None)
            )
        
        return result_df



    def visualize_tissue_enrichment(self, enrichment_results=None, **kwargs):
        """
        Visualize tissue or cell type enrichment results.
        Delegates to the parent class visualization method if available.
        
        Parameters:
        -----------
        enrichment_results : pandas.DataFrame, optional
            Enrichment results to visualize. If None, uses stored results.
        **kwargs : dict
            Additional parameters for visualization
        
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        # Use results from parameter, or stored tissue results if not provided
        if enrichment_results is None:
            if self.region_tissue_enrichment is not None:
                enrichment_results = self.region_tissue_enrichment
            elif self.region_cell_type_enrichment is not None:
                enrichment_results = self.region_cell_type_enrichment
            else:
                self.logger.error("No enrichment results available. Run analyze_tissue_enrichment or analyze_cell_type_enrichment first.")
                return None
        
        # If parent class has visualization method, use it
        if hasattr(super(), 'visualize_tissue_enrichment'):
            return super().visualize_tissue_enrichment(
                enrichment_results=enrichment_results,
                **kwargs
            )
        else:
            self.logger.warning("Parent class does not have visualize_tissue_enrichment method. Visualization skipped.")
            return None

    def visualize_fdr_horizontal_bar(self, enrichment_results=None, 
                                    title=None, figsize=(14, 10), top_n=15, 
                                    pvalue_threshold=0.05, fc_threshold=1.5,
                                    region_name=None, 
                                    enriched_color='#e66101', depleted_color='#5e3c99', 
                                    nonsig_color='#b2b2b2',
                                    # 追加オプション
                                    save_dir=None, save_filename=None,
                                    save_formats=None, save_dpi=300):
        """
        エンリッチメント解析結果をFDRに基づく横棒グラフで可視化します。
        
        Parameters:
        -----------
        enrichment_results : pandas.DataFrame, optional
            エンリッチメント解析結果。指定しない場合は、保存されている結果を使用
        title : str, optional
            グラフのタイトル。指定しない場合は自動生成
        figsize : tuple, optional
            グラフのサイズ（デフォルト: (14, 10)）
        top_n : int, optional
            表示する上位の組織/細胞タイプの数（デフォルト: 15）
        pvalue_threshold : float, optional
            有意性の閾値（デフォルト: 0.05）
        fc_threshold : float, optional
            フォールドチェンジの閾値（デフォルト: 1.5）
        region_name : str, optional
            領域名。タイトルに使用
        enriched_color : str, optional
            エンリッチされた組織/細胞のバーの色（デフォルト: オレンジ）
        depleted_color : str, optional
            デプリートされた組織/細胞のバーの色（デフォルト: 紫）
        nonsig_color : str, optional
            有意でない組織/細胞のバーの色（デフォルト: グレー）
        save_dir : str, optional
            図を保存するディレクトリ
        save_filename : str, optional
            保存する図のファイル名（拡張子なし）
        save_formats : list, optional
            保存する図の形式のリスト（例: ['png', 'svg', 'pdf']）
        save_dpi : int, optional
            保存する図の解像度（デフォルト: 300）
            
        Returns:
        --------
        matplotlib.figure.Figure
            作成されたグラフのFigureオブジェクト
        """
        # 結果の取得
        if enrichment_results is None:
            if hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None:
                enrichment_results = self.region_tissue_enrichment
            elif hasattr(self, 'region_cell_type_enrichment') and self.region_cell_type_enrichment is not None:
                enrichment_results = self.region_cell_type_enrichment
            else:
                self.logger.error("エンリッチメント解析結果が見つかりません。analyze_tissue_enrichmentまたはanalyze_cell_type_enrichmentを実行してください。")
                return None
        
        # 結果が空の場合
        if enrichment_results is None or enrichment_results.empty:
            self.logger.warning("エンリッチメント解析結果が空です。")
            return None
        
        # カラム名の確認と取得
        is_tissue = 'Tissue Type' in enrichment_results.columns
        is_cell = 'Cell Type' in enrichment_results.columns
        
        if is_tissue:
            type_col = 'Tissue Type'
            type_name = 'Tissue'
        elif is_cell:
            type_col = 'Cell Type'
            type_name = 'Cell Type'
        else:
            self.logger.error("組織型または細胞型のカラムが見つかりません")
            return None
        
        # FDRカラムの確認
        if 'adjusted p-value' in enrichment_results.columns:
            fdr_col = 'adjusted p-value'
        elif 'FDR' in enrichment_results.columns:
            fdr_col = 'FDR'
        else:
            fdr_col = 'p-value'
            self.logger.warning(f"FDRカラムが見つからないため、p値カラム'{fdr_col}'を使用します")
        
        # その他のカラムを確認
        fc_col = 'Fold Change' if 'Fold Change' in enrichment_results.columns else 'フォールドチェンジ'
        count_col = 'Region Count' if 'Region Count' in enrichment_results.columns else 'リージョンカウント'
        
        # -log10(FDR)の計算
        results_df = enrichment_results.copy()
        results_df['-log10(FDR)'] = -np.log10(results_df[fdr_col])
        
        # FDRで並べ替えて上位を抽出
        sorted_df = results_df.sort_values('-log10(FDR)', ascending=False)
        if top_n and len(sorted_df) > top_n:
            plot_data = sorted_df.head(top_n)
        else:
            plot_data = sorted_df
        
        # 横棒グラフ用に逆順にする（有意性の高いものが上に来るように）
        plot_data = plot_data.iloc[::-1].reset_index(drop=True)
        
        # グラフの作成
        fig, ax = plt.subplots(figsize=figsize)
        
        # 横棒グラフを作成
        bars = ax.barh(
            y=plot_data[type_col],
            width=plot_data['-log10(FDR)'],
            color=[
                enriched_color if row[fc_col] >= fc_threshold else
                depleted_color if row[fc_col] <= 1/fc_threshold else
                nonsig_color
                for _, row in plot_data.iterrows()
            ],
            alpha=0.7
        )
        
        # 有意性閾値の線を追加
        significance_line = -np.log10(pvalue_threshold)
        ax.axvline(x=significance_line, linestyle='--', color='black', alpha=0.7)
        
        # バーにラベルを追加
        for i, bar in enumerate(bars):
            row = plot_data.iloc[i]
            fc = row[fc_col]
            count = row[count_col] if count_col in row else 0
            fdr_value = row[fdr_col]
            bar_width = bar.get_width()
            y_pos = bar.get_y() + bar.get_height() * 0.5
            
            # テキストの位置を決定（バーが有意性の線を超える場合は線の外に）
            if bar_width >= significance_line:
                x_pos = max(bar_width + 0.1, significance_line + 0.3)
            else:
                x_pos = bar_width + 0.1
            
            # すべての情報を含むコンパクトなラベルを作成
            label_text = f"FDR={fdr_value:.1e}, n={int(count)}, FC={fc:.2f}"
            
            ax.text(
                x=x_pos,
                y=y_pos,
                s=label_text,
                va='center',
                ha='left',
                fontsize=9,
                color='dimgray'
            )
        
        # グラフの書式設定
        ax.set_xlabel('-log10(FDR)')
        ax.set_ylabel(f'{type_name}')
        
        # タイトルの設定
        if title is None:
            title_prefix = f"{type_name} Enrichment Analysis"
            if region_name:
                title_prefix += f" - {region_name}"
            title = f'{title_prefix} ({len(plot_data)} {type_name.lower()}s)'
        
        ax.set_title(f'{title} - Statistical Significance (Horizontal Bar Plot)')
        
        # 凡例の追加
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        
        # デプリート閾値を計算
        depleted_threshold = round(1/fc_threshold, 2)
        
        legend_elements = [
            Patch(facecolor=enriched_color, alpha=0.7, label=f'Enriched (FC >= {fc_threshold})'),
            Patch(facecolor=depleted_color, alpha=0.7, label=f'Depleted (FC <= {depleted_threshold})'),
            Patch(facecolor=nonsig_color, alpha=0.7, label='Not enriched/depleted'),
            Line2D([0], [0], color='black', linestyle='--', 
                label=f'FDR = {pvalue_threshold} threshold')
        ]
        ax.legend(handles=legend_elements, loc='lower right', bbox_to_anchor=(0.98, 0.02))
        
        # x軸の調整（すべてのバーとラベルが見えるように）
        max_width = max(plot_data['-log10(FDR)'])
        ax.set_xlim(0, max_width * 1.3)  # ラベル用に30%の余白を追加
        
        plt.tight_layout()
        
        # 図を保存するオプションを追加
        if save_dir and save_filename:
            # デフォルトの保存形式
            if save_formats is None:
                save_formats = ['png', 'svg', 'pdf']
            
            # 保存ディレクトリが存在しない場合は作成
            import os
            if not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
                self.logger.info(f"Created output directory: {save_dir}")
            
            # 各形式で保存
            saved_files = []
            for fmt in save_formats:
                file_path = os.path.join(save_dir, f"{save_filename}.{fmt}")
                fig.savefig(file_path, format=fmt, dpi=save_dpi, bbox_inches='tight')
                saved_files.append(f"{fmt.upper()}: {file_path}")
            
            if saved_files:
                self.logger.info(f"Saved FDR horizontal bar chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        # 表示しない（後でレポート生成側で処理するため）
        # plt.show()
        
        return fig









    def visualize_region_reads(self, chrom=None, start=None, end=None, 
                            figsize=(15, 8), max_samples=50, sample_colors=None,
                            highlight_interest_region=True, title=None,
                            show_sample_names=True, group_by=None):
        """
        抽出されたリージョン内のeRNAリードの分布をゲノム座標上に可視化します。
        
        Parameters:
        -----------
        chrom : str, optional
            染色体名。指定しない場合はself.interest_regionから取得
        start : int, optional
            開始位置。指定しない場合はself.interest_regionから取得
        end : int, optional
            終了位置。指定しない場合はself.interest_regionから取得
        figsize : tuple, optional
            グラフのサイズ
        max_samples : int, optional
            表示するサンプルの最大数（多すぎる場合はランダムにサンプリング）
        sample_colors : dict, optional
            サンプル名：色のマッピング辞書。指定しない場合は自動生成
        highlight_interest_region : bool, optional
            関心領域を強調表示するかどうか
        title : str, optional
            グラフのタイトル
        show_sample_names : bool, optional
            サンプル名を表示するかどうか
        group_by : str, optional
            サンプルをグループ化する列名（'Tissue Type'や'Cell Type'など）
            
        Returns:
        --------
        matplotlib.figure.Figure
            作成されたグラフのFigureオブジェクト
        """
        # 位置情報の確認と取得
        if chrom is None or start is None or end is None:
            if self.interest_region is not None:
                chrom = self.interest_region.get('chrom')
                start = self.interest_region.get('start')
                end = self.interest_region.get('end')
            else:
                self.logger.error("領域情報が指定されていません。chrom, start, endを指定するか、create_interest_regionを先に実行してください。")
                return None
        
        # オーバーラップ結果の確認
        if self.overlapping_results is None or self.overlapping_results.empty:
            self.logger.error("オーバーラップ結果が見つかりません。merge_all_overlapsを先に実行してください。")
            return None
        
        # データの準備
        reads_df = self.overlapping_results.copy()
        
        # 列名の確認と取得
        required_cols = ['chrom', 'start', 'end', 'source_file']
        if not all(col in reads_df.columns for col in required_cols):
            missing = [col for col in required_cols if col not in reads_df.columns]
            self.logger.error(f"必要な列が見つかりません: {', '.join(missing)}")
            return None
        
        # 指定された染色体のデータだけを抽出
        chrom_df = reads_df[reads_df['chrom'] == chrom].copy()
        if chrom_df.empty:
            self.logger.warning(f"染色体 {chrom} のデータが見つかりません。")
            return None
        
        # サンプルのリストを取得
        samples = chrom_df['source_file'].unique()
        n_samples = len(samples)
        
        # サンプル数が多い場合のサンプリング
        if n_samples > max_samples:
            self.logger.warning(f"サンプル数 ({n_samples}) が最大表示数 ({max_samples}) を超えています。ランダムにサンプリングします。")
            samples = np.random.choice(samples, max_samples, replace=False)
            chrom_df = chrom_df[chrom_df['source_file'].isin(samples)]
        
        # メタデータをマージ（利用可能な場合）
        if self.overlapping_metadata is not None and group_by is not None:
            if group_by in self.overlapping_metadata.columns:
                # メタデータとサンプルをマージするためのマッピングを作成
                sample_id_col = 'Sample ID' if 'Sample ID' in self.overlapping_metadata.columns else 'sample_id'
                sample_map = {}
                
                for _, row in self.overlapping_metadata.iterrows():
                    sample_id = row[sample_id_col]
                    # ファイル名はサンプルIDから.bedを取り除いたもの
                    if sample_id in samples:
                        sample_map[sample_id] = row[group_by]
                
                # サンプルにグループ情報を追加
                chrom_df['group'] = chrom_df['source_file'].map(sample_map)
                
                # グループでソート
                chrom_df = chrom_df.sort_values('group')
                
                # グループごとの色を設定
                if sample_colors is None:
                    groups = chrom_df['group'].unique()
                    cmap = plt.cm.get_cmap('tab20', len(groups))
                    group_colors = {group: cmap(i) for i, group in enumerate(groups)}
                    sample_colors = {sample: group_colors[chrom_df[chrom_df['source_file'] == sample]['group'].iloc[0]] 
                                for sample in chrom_df['source_file'].unique()}
        
        # 色の設定
        if sample_colors is None:
            unique_samples = chrom_df['source_file'].unique()
            cmap = plt.cm.get_cmap('tab20', len(unique_samples))
            sample_colors = {sample: cmap(i) for i, sample in enumerate(unique_samples)}
        
        # グラフの作成
        fig, ax = plt.subplots(figsize=figsize)
        
        # サンプルのインデックスを追跡
        sample_indices = {}
        current_index = 0
        
        # グループのインデックス（グループごとに区切り線を引く用）
        group_boundaries = []
        current_group = None
        
        # 各サンプルについてリードを描画
        for i, (sample, sample_df) in enumerate(chrom_df.groupby('source_file')):
            # サンプルのインデックスを記録
            sample_indices[sample] = current_index
            
            # グループの変更を検出（グループでソートされている前提）
            if 'group' in chrom_df.columns:
                group = sample_df['group'].iloc[0]
                if current_group is None:
                    current_group = group
                elif current_group != group:
                    group_boundaries.append(current_index - 0.5)
                    current_group = group
            
            # 各リードを描画
            for _, read in sample_df.iterrows():
                read_start = int(read['start'])
                read_end = int(read['end'])
                
                # リードを描画（水平な長方形）
                rect = plt.Rectangle(
                    (read_start, current_index - 0.4),
                    read_end - read_start,
                    0.8,
                    facecolor=sample_colors.get(sample, 'gray'),
                    alpha=0.7,
                    edgecolor='black',
                    linewidth=0.5
                )
                ax.add_patch(rect)
            
            # 次のサンプルへ
            current_index += 1
        
        # 関心領域の強調表示
        if highlight_interest_region:
            # 関心領域の背景をハイライト
            highlight = plt.Rectangle(
                (start, -1),
                end - start,
                current_index + 1,
                facecolor='yellow',
                alpha=0.1,
                zorder=-1
            )
            ax.add_patch(highlight)
            
            # 境界線
            ax.axvline(x=start, color='red', linestyle='--', alpha=0.5)
            ax.axvline(x=end, color='red', linestyle='--', alpha=0.5)
        
        # グループの区切り線を描画
        for boundary in group_boundaries:
            ax.axhline(y=boundary, color='black', linestyle='-', alpha=0.3)
        
        # 軸の設定
        ax.set_xlim(start - (end - start) * 0.05, end + (end - start) * 0.05)  # 余白を追加
        ax.set_ylim(-1, current_index)
        
        # y軸のティックとラベル
        if show_sample_names:
            # サンプル名を表示
            yticks = list(sample_indices.values())
            yticklabels = list(sample_indices.keys())
            
            # グループラベルを追加（利用可能な場合）
            if 'group' in chrom_df.columns and group_boundaries:
                # グループラベルの位置を計算
                group_positions = []
                group_labels = []
                
                current_start = 0
                for boundary in group_boundaries + [current_index]:
                    position = (current_start + boundary) / 2
                    group_positions.append(position)
                    
                    # このグループのサンプルを取得
                    group_samples = list(sample_indices.keys())[current_start:int(boundary + 0.5)]
                    if group_samples:
                        group_name = chrom_df[chrom_df['source_file'] == group_samples[0]]['group'].iloc[0]
                        group_labels.append(f"{group_name} ({len(group_samples)})")
                    
                    current_start = int(boundary + 0.5)
                
                # グループラベルを別の軸に表示
                ax2 = ax.twinx()
                ax2.set_yticks(group_positions)
                ax2.set_yticklabels(group_labels, fontsize=12, fontweight='bold')
                ax2.set_ylim(ax.get_ylim())
                ax2.spines['right'].set_visible(False)
                ax2.tick_params(right=False)
        else:
            # サンプル名なしの場合、単純なインデックスを表示
            yticks = range(current_index)
            yticklabels = range(1, current_index + 1)
        
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        
        # x軸の形式設定（大きな数値の場合に見やすく）
        def format_bp(x, pos):
            if x >= 1e6:
                return f'{x/1e6:.2f} Mb'
            elif x >= 1e3:
                return f'{x/1e3:.0f} kb'
            else:
                return f'{x:.0f} bp'
        
        from matplotlib.ticker import FuncFormatter
        ax.xaxis.set_major_formatter(FuncFormatter(format_bp))
        
        # グラフタイトルと軸ラベル
        if title is None:
            title = f"eRNA Reads Distribution in {chrom}:{start:,}-{end:,}"
        
        ax.set_title(title)
        ax.set_xlabel('Genomic Position')
        ax.set_ylabel('Samples')
        
        # 凡例（グループごとに）
        if 'group' in chrom_df.columns:
            from matplotlib.patches import Patch
            
            groups = chrom_df['group'].unique()
            legend_elements = [
                Patch(facecolor=group_colors[group], edgecolor='black', alpha=0.7, label=group)
                for group in groups
            ]
            ax.legend(handles=legend_elements, loc='upper right', title='Groups')
        
        plt.tight_layout()
        plt.show()
        
        return fig



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
        Display eRNA reads in the specified genomic region in a simple stacked layout,
        similar to traditional genomic read displays.
        Style matches plot_stacked_reads_with_peaks.
        
        Parameters:
        -----------
        chrom : str, optional
            Chromosome name. If not specified, retrieved from self.interest_region
        start : int, optional
            Start position. If not specified, retrieved from self.interest_region
        end : int, optional
            End position. If not specified, retrieved from self.interest_region
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
            if self.interest_region is not None:
                chrom = self.interest_region.get('chrom')
                start = self.interest_region.get('start')
                end = self.interest_region.get('end')
            else:
                self.logger.error("Region information not specified. Please specify chrom, start, end or run create_interest_region first.")
                return None
        
        # Check overlap results
        if self.overlapping_results is None or self.overlapping_results.empty:
            self.logger.error("No overlap results found. Please run merge_all_overlaps first.")
            return None
        
        # #緊急バグフィックス
        # # Filter data for the specified region
        # reads_df = self.overlapping_results
        # region_data = reads_df[(reads_df['chrom'] == chrom) & 
        #                     (reads_df['start'] >= start) & 
        #                     (reads_df['end'] <= end)]

       # Filter data for the specified region
        reads_df = self.overlapping_results
        region_data = reads_df[(reads_df['chrom'] == chrom) & 
                            (reads_df['start'] <= start) & 
                            (reads_df['end'] >= end)]
 

        
        if region_data.empty:
            self.logger.warning(f"No reads found in the specified region: {chrom}:{start}-{end}")
            return None
        
        # Count the number of samples
        sample_count = region_data['source_file'].nunique()
        
        # Display read count
        self.logger.info(f"Number of reads to display: {len(region_data)}")
        self.logger.info(f"Number of samples: {sample_count}")
        
        # Sort reads by start position
        region_data = region_data.sort_values('start')
        
        # Save region BED file if requested
        if save_region_bed:
            try:
                from pybedtools import BedTool
                region_bed = BedTool.from_dataframe(region_data)
                region_bed.saveas(save_region_bed)
                self.logger.info(f"Reads in region saved as BED file: {save_region_bed}")
            except Exception as e:
                self.logger.error(f"Error saving BED file: {e}")
        
        # Convert DataFrame to list of read dictionaries
        reads = []
        for _, row in region_data.iterrows():
            reads.append({
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['end'],
                'sample': row['source_file']
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
        
        # Create the plot - matching plot_stacked_reads_with_peaks style
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
        
        # Configure x-axis to match the style in plot_stacked_reads_with_peaks
        plt.gca().xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        plt.ticklabel_format(style='plain', axis='x')
        
        # Set tick font sizes
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Add ONLY vertical grid lines to match plot_stacked_reads_with_peaks
        if grid:
            # Draw only vertical grid lines
            ax.xaxis.grid(True, linestyle='--', color='lightgray', alpha=0.7)
            ax.yaxis.grid(False)  # Explicitly disable horizontal grid lines
        
        # Add sample count to the plot - match the style of plot_stacked_reads_with_peaks
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
        Display eRNA reads in the specified genomic region with peak detection.
        Grid style matches region_visualize_tools.py: vertical grid lines only by default.
        
        Parameters:
        -----------
        chrom : str, optional
            Chromosome name. If not specified, retrieved from self.interest_region
        start : int, optional
            Start position. If not specified, retrieved from self.interest_region
        end : int, optional
            End position. If not specified, retrieved from self.interest_region
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
            if self.interest_region is not None:
                chrom = self.interest_region.get('chrom')
                start = self.interest_region.get('start')
                end = self.interest_region.get('end')
            else:
                self.logger.error("Region information not specified. Please specify chrom, start, end or run create_interest_region first.")
                return None, {}
        
        # Check overlap results
        if self.overlapping_results is None or self.overlapping_results.empty:
            self.logger.error("No overlap results found. Please run merge_all_overlaps first.")
            return None, {}

        #緊急バグフィックス        
        # # Filter data for the specified region
        # reads_df = self.overlapping_results
        # region_data = reads_df[(reads_df['chrom'] == chrom) & 
        #                     (reads_df['start'] >= start) & 
        #                     (reads_df['end'] <= end)]

        # Filter data for the specified region
        reads_df = self.overlapping_results
        region_data = reads_df[(reads_df['chrom'] == chrom) & 
                            (reads_df['start'] <= start) & 
                            (reads_df['end'] >= end)]


        if region_data.empty:
            self.logger.warning(f"No reads found in the specified region: {chrom}:{start}-{end}")
            return None, {}
        
        # Count the number of samples
        sample_count = region_data['source_file'].nunique()
        
        # Display read count
        self.logger.info(f"Number of reads to display: {len(region_data)}")
        self.logger.info(f"Number of samples: {sample_count}")
        
        # Sort reads by start position
        region_data = region_data.sort_values('start')
        
        # Save region BED file if requested
        if save_region_bed:
            try:
                from pybedtools import BedTool
                region_bed = BedTool.from_dataframe(region_data)
                region_bed.saveas(save_region_bed)
                self.logger.info(f"Reads in region saved as BED file: {save_region_bed}")
            except Exception as e:
                self.logger.error(f"Error saving BED file: {e}")
        
        # Convert DataFrame to list of read dictionaries
        reads = []
        for _, row in region_data.iterrows():
            reads.append({
                'chrom': row['chrom'],
                'start': row['start'],
                'end': row['end'],
                'sample': row['source_file']
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
        
        # Create the plot - matching region_visualize_tools.py style
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
            
            # Find peaks - same as original code
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
            
            # Sort peaks by position for better labeling - same as original code
            sorted_peaks = sorted(peak_display_data.items(), key=lambda x: x[1]['position'])
            
            # Calculate peak label positions - same as original code
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
                
                # Add peak annotation with improved positioning - same as original code
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
        
        # Configure x-axis to match the style in region_visualize_tools.py
        plt.gca().xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        plt.ticklabel_format(style='plain', axis='x')
        
        # Set tick font sizes
        plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
        
        # Add ONLY vertical grid lines - explicitly control grid lines to match image b
        if grid:
            # Draw only vertical grid lines
            ax.xaxis.grid(True, linestyle='--', color='lightgray', alpha=0.7)
            ax.yaxis.grid(False)  # Explicitly disable horizontal grid lines
        
        # Add sample count and peak info in top right corner
        if peak_data:
            info_text = f"Total: {sample_count} samples, {len(peak_data)} peaks"
            plt.text(0.98, 0.98, info_text, 
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



    def plot_sorted_bar(self, data: pd.DataFrame, category: str, count_col: str = 'Count',
                    title: Optional[str] = None, limit: int = 15, 
                    figsize: tuple = (10, 8), horizontal: bool = True,
                    color: str = 'skyblue', grid: bool = True,
                    title_fontsize: int = 14, label_fontsize: int = 12, tick_fontsize: int = 10,
                    save_dir: Optional[str] = None, save_filename: Optional[str] = None,
                    save_svg: bool = False, save_png: bool = False, 
                    save_pdf: bool = False, save_eps: bool = False,
                    save_dpi: int = 600, save_transparent: bool = False) -> plt.Figure:
        """
        Create a sorted bar chart showing frequencies in descending order
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame to analyze
        category : str
            Category column name (e.g., 'Tissue Type', 'Cell Type')
        count_col : str
            Column name containing count values (default: 'Count')
        title : str, optional
            Chart title
        limit : int, default=15
            Number of top categories to display
        figsize : tuple, default=(10, 8)
            Figure size
        horizontal : bool, default=True
            Whether to create horizontal bar chart (True) or vertical (False)
        color : str, default='skyblue'
            Bar color
        grid : bool, default=True
            Whether to display grid lines
        title_fontsize : int, default=14
            Font size for the plot title
        label_fontsize : int, default=12
            Font size for axis labels
        tick_fontsize : int, default=10
            Font size for tick labels
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
        # Check if data exists
        if data is None or data.empty:
            self.logger.warning(f"Empty data provided for plotting")
            return None
            
        # Ensure the columns exist
        if category not in data.columns or count_col not in data.columns:
            missing_cols = []
            if category not in data.columns:
                missing_cols.append(category)
            if count_col not in data.columns:
                missing_cols.append(count_col)
            self.logger.warning(f"Columns not found in data: {', '.join(missing_cols)}")
            return None
        
        # Sort by count in descending order and limit to top entries
        plot_data = data.sort_values(count_col, ascending=False).head(limit)
        
        # Set default title if not provided
        if title is None:
            title = f"Top {min(limit, len(plot_data))} {category} by {count_col}"
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create horizontal or vertical bar chart
        if horizontal:
            # For horizontal bars, we need to reverse the order so highest is at top
            plot_data = plot_data.iloc[::-1]
            bars = ax.barh(plot_data[category], plot_data[count_col], color=color)
            ax.set_xlabel(count_col, fontsize=label_fontsize)
            ax.set_ylabel(category, fontsize=label_fontsize)
            
            # Add value labels to bars
            for bar in bars:
                width = bar.get_width()
                ax.text(width + (width * 0.02), 
                        bar.get_y() + bar.get_height()/2, 
                        f'{width:,.0f}', 
                        va='center', 
                        fontsize=tick_fontsize)
            
            # Add grid only on x-axis for horizontal bars
            if grid:
                ax.xaxis.grid(True, linestyle='--', alpha=0.6)
                ax.yaxis.grid(False)
        else:
            # For vertical bars
            bars = ax.bar(plot_data[category], plot_data[count_col], color=color)
            ax.set_xlabel(category, fontsize=label_fontsize)
            ax.set_ylabel(count_col, fontsize=label_fontsize)
            
            # Rotate x-tick labels for better readability
            plt.xticks(rotation=45, ha='right', fontsize=tick_fontsize)
            
            # Add value labels to bars
            for bar in bars:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2, 
                        height + (height * 0.02), 
                        f'{height:,.0f}', 
                        ha='center', 
                        fontsize=tick_fontsize)
            
            # Add grid only on y-axis for vertical bars
            if grid:
                ax.yaxis.grid(True, linestyle='--', alpha=0.6)
                ax.xaxis.grid(False)
        
        # Set general plot attributes
        ax.set_title(title, fontsize=title_fontsize)
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        # Tight layout
        plt.tight_layout()
        
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
                plot_type = "horizontal_bar" if horizontal else "vertical_bar"
                save_filename = f"{plot_type}_{category.lower().replace(' ', '_')}"
            
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
                self.logger.info(f"Saved bar chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        plt.show()
        return fig




    def generate_erna_db_report(self, metadata_file=None, output_dir='results/report',
                            species='human', top_n=15, figure_formats=None, save_tables=True,
                            dpi=300, fig_width=10, fig_height=6):
        """
        Generate a simple HTML report of the eRNA metadata analysis with relative links and multiple figure formats.
        
        Parameters:
        -----------
        metadata_file : str, optional
            Path to the metadata file (.csv or .parquet). If None, uses already loaded metadata.
        output_dir : str, optional
            Base directory for all output files (default: 'results/report')
        species : str, optional
            Species filter option (default: 'human')
        top_n : int, optional
            Number of top categories to display in visualizations (default: 15)
        figure_formats : list, optional
            List of figure formats to save (default: ['png', 'svg', 'pdf', 'eps'])
        save_tables : bool, optional
            Whether to save distribution tables as CSV files (default: True)
        dpi : int, optional
            Resolution for raster images like PNG (default: 300)
        fig_width : int, optional
            Width of the figure in inches (default: 10)
        fig_height : int, optional
            Height of the figure in inches (default: 6)
            
        Returns:
        --------
        str
            Path to the generated HTML report
        """
        
        # Default figure formats if not provided
        if figure_formats is None:
            figure_formats = ['png', 'svg', 'pdf', 'eps']
        
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
        
        self.logger.info(f"Creating simple report in directory: {output_dir}")
        
        # Dictionary to track all generated files
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Step 1: Load metadata if needed
        if metadata_file is not None or self.erna_metadata is None:
            if metadata_file is None:
                self.logger.error("No metadata loaded and no metadata file provided")
                return None
            
            self.load_erna_metadata(metadata_file, species=species)
        
        if self.erna_metadata is None:
            self.logger.error("Failed to load metadata")
            return None
        
        # Collect basic metadata stats
        metadata_stats = {
            'total_samples': len(self.erna_metadata),
            'columns': list(self.erna_metadata.columns),
            'timestamp': timestamp,
            'species': species
        }
        
        # Calculate missing values
        missing_values = self.erna_metadata.isnull().sum()
        missing_values = missing_values[missing_values > 0]
        metadata_stats['missing_values'] = missing_values.to_dict()
        
        # Step 2: Calculate distributions if not already done
        self.calculate_metadata_distributions()
        
        # Set Matplotlib params for better quality figures
        plt.rcParams['figure.figsize'] = (fig_width, fig_height)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # Step 3: Generate visualizations in multiple formats
        if self.metadata_tissue_distribution is not None:
            # Prepare save options
            save_opts = {}
            for fmt in ['png', 'svg', 'pdf', 'eps']:
                save_opts[f'save_{fmt}'] = fmt in figure_formats
            
            try:
                # Save tissue bar chart
                tissue_basename = f'tissue_distribution_{timestamp}'
                fig = self.plot_sorted_bar(
                    data=self.metadata_tissue_distribution,
                    category='Tissue Type',
                    count_col='Count',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='skyblue',
                    save_dir=figures_full_path,
                    save_filename=tissue_basename,
                    **save_opts  # Unpacks into save_png=True/False, save_svg=True/False, etc.
                )
                
                # Register all generated figure files with relative paths
                for fmt in figure_formats:
                    figure_filename = f'{tissue_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)  # Relative path for HTML
                    figure_full_path = os.path.join(figures_full_path, figure_filename)  # Full path to check if file exists
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                
                # Additional: Create and save tissue pie chart
                tissue_pie_basename = f'tissue_pie_chart_{timestamp}'
                pie_fig = self.plot_simple_pie(
                    data=self.metadata_tissue_distribution,
                    category='Tissue Type',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{tissue_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                
            except Exception as e:
                self.logger.warning(f"Error generating tissue distribution chart: {e}")
        
        if self.metadata_cell_type_distribution is not None:
            try:
                # Save cell type bar chart
                cell_basename = f'cell_type_distribution_{timestamp}'
                fig = self.plot_sorted_bar(
                    data=self.metadata_cell_type_distribution,
                    category='Cell Type',
                    count_col='Count',
                    title=f'Cell Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='lightgreen',
                    save_dir=figures_full_path,
                    save_filename=cell_basename,
                    **save_opts  # Unpacks into save_png=True/False, save_svg=True/False, etc.
                )
                
                # Register all generated figure files with relative paths
                for fmt in figure_formats:
                    figure_filename = f'{cell_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)  # Relative path for HTML
                    figure_full_path = os.path.join(figures_full_path, figure_filename)  # Full path to check if file exists
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                
                # Additional: Create and save cell type pie chart
                cell_pie_basename = f'cell_type_pie_chart_{timestamp}'
                pie_fig = self.plot_simple_pie(
                    data=self.metadata_cell_type_distribution,
                    category='Cell Type',
                    title=f'Cell Type Distribution (Top {top_n})',
                    limit=top_n,
                    save_dir=figures_full_path,
                    save_filename=cell_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{cell_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        
            except Exception as e:
                self.logger.warning(f"Error generating cell type distribution chart: {e}")
        
        # Step 4: Generate heatmap if possible
        if 'Tissue Type' in self.erna_metadata.columns and 'Cell Type' in self.erna_metadata.columns:
            try:
                # Create cross-tabulation
                tissue_cell_cross = pd.crosstab(
                    self.erna_metadata['Tissue Type'],
                    self.erna_metadata['Cell Type']
                )
                
                # Save cross-tabulation data
                if save_tables:
                    cross_table_filename = f'tissue_cell_crosstab_{timestamp}.csv'
                    cross_table_path = os.path.join(tables_full_path, cross_table_filename)
                    tissue_cell_cross.to_csv(cross_table_path)
                    generated_files['tables'].append(('Tissue-Cell Cross-tabulation', os.path.join(tables_dir, cross_table_filename)))
                
                # Generate heatmap visualization
                try:
                    import seaborn as sns
                    
                    # Limit to top categories to avoid too large heatmap
                    top_tissues = tissue_cell_cross.sum(axis=1).nlargest(top_n).index
                    top_cells = tissue_cell_cross.sum(axis=0).nlargest(top_n).index
                    
                    # Create subset of the cross-tab for display
                    display_data = tissue_cell_cross.loc[top_tissues, top_cells]
                    
                    # Create heatmap figure
                    plt.figure(figsize=(fig_width, fig_height))
                    sns.heatmap(display_data, annot=True, cmap='YlGnBu', fmt='d', 
                            cbar_kws={'label': 'Number of Samples'})
                    plt.title(f'Relationship between Tissue Type and Cell Type (Top {min(top_n, len(tissue_cell_cross))} each)')
                    plt.tight_layout()
                    
                    # Save heatmap in multiple formats
                    heatmap_basename = f'tissue_cell_heatmap_{timestamp}'
                    
                    for fmt in figure_formats:
                        figure_filename = f'{heatmap_basename}.{fmt}'
                        figure_path = os.path.join(figures_full_path, figure_filename)
                        plt.savefig(figure_path, format=fmt, bbox_inches='tight', dpi=dpi)
                        
                        # Register with relative path
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        generated_files['figures'].append((f'Tissue-Cell Heatmap ({fmt.upper()})', figure_relative_path))
                    
                    plt.close()
                except Exception as e:
                    self.logger.warning(f"Could not generate heatmap visualization: {e}")
            except Exception as e:
                self.logger.warning(f"Could not generate cross-tabulation data: {e}")
        
        # Step 5: Save distribution tables
        if save_tables and self.metadata_tissue_distribution is not None:
            tissue_table_filename = f'tissue_distribution_{timestamp}.csv'
            tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
            self.metadata_tissue_distribution.to_csv(tissue_table_path, index=False)
            generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
        
        if save_tables and self.metadata_cell_type_distribution is not None:
            cell_table_filename = f'cell_type_distribution_{timestamp}.csv'
            cell_table_path = os.path.join(tables_full_path, cell_table_filename)
            self.metadata_cell_type_distribution.to_csv(cell_table_path, index=False)
            generated_files['tables'].append(('Cell Type Distribution', os.path.join(tables_dir, cell_table_filename)))
        
        # Step 6: Generate simple HTML report
        html_filename = f"simple_report_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        try:
            with open(html_path, 'w', encoding='utf-8') as f:
                # Start HTML document
                f.write("<!DOCTYPE html>\n")
                f.write("<html>\n")
                f.write("<head>\n")
                f.write("  <title>eRNAbase Metadata Report</title>\n")
                f.write("  <style>\n")
                f.write("    body { font-family: Arial, sans-serif; margin: 20px; }\n")
                f.write("    h1 { color: #2c3e50; }\n")
                f.write("    h2 { color: #3498db; margin-top: 30px; }\n")
                f.write("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
                f.write("    th { background-color: #3498db; color: white; text-align: left; padding: 8px; }\n")
                f.write("    td { border: 1px solid #ddd; padding: 8px; }\n")
                f.write("    tr:nth-child(even) { background-color: #f2f2f2; }\n")
                f.write("    img { max-width: 100%; height: auto; margin: 20px 0; }\n")
                f.write("    .file-formats { margin: 10px 0; }\n")
                f.write("    .file-formats a { margin-right: 10px; }\n")
                f.write("    .chart-section { display: flex; flex-wrap: wrap; justify-content: space-around; }\n")
                f.write("    .chart-container { margin: 15px; max-width: 48%; }\n")
                f.write("  </style>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                
                # Title and timestamp
                f.write(f"  <h1>eRNAbase Metadata Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Overview section
                f.write("  <h2>1. Overview</h2>\n")
                f.write("  <ul>\n")
                f.write(f"    <li><strong>Total samples</strong>: {metadata_stats['total_samples']}</li>\n")
                f.write(f"    <li><strong>Species filter</strong>: {metadata_stats['species']}</li>\n")
                f.write(f"    <li><strong>Columns in dataset</strong>: {', '.join(metadata_stats['columns'])}</li>\n")
                f.write("  </ul>\n")
                
                # Missing values section
                if metadata_stats['missing_values']:
                    f.write("  <h3>Missing Values</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Column</th><th>Missing Count</th><th>Percentage</th></tr>\n")
                    
                    for col, count in metadata_stats['missing_values'].items():
                        percentage = (count / metadata_stats['total_samples']) * 100
                        f.write(f"    <tr><td>{col}</td><td>{count}</td><td>{percentage:.2f}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                
                # Tissue Distribution section
                if self.metadata_tissue_distribution is not None:
                    f.write("  <h2>2. Tissue Type Distribution</h2>\n")
                    
                    # Table of top tissues
                    f.write(f"  <h3>Top {min(top_n, len(self.metadata_tissue_distribution))} Tissue Types</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Tissue Type</th><th>Count</th><th>Percentage</th></tr>\n")
                    
                    top_tissues = self.metadata_tissue_distribution.head(top_n)
                    for _, row in top_tissues.iterrows():
                        f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Count']}</td><td>{row['Percentage (%)']}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                    
                    # Tissue distribution charts - both bar and pie chart
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    tissue_bar_png = [path for name, path in generated_files['figures'] 
                                    if 'Tissue Distribution Bar Chart' in name and path.endswith('.png')]
                    if tissue_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{tissue_bar_png[0]}' alt='Tissue Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    tissue_pie_png = [path for name, path in generated_files['figures'] 
                                    if 'Tissue Distribution Pie Chart' in name and path.endswith('.png')]
                    if tissue_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{tissue_pie_png[0]}' alt='Tissue Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    tissue_figures_all = [tup for tup in generated_files['figures'] if 'Tissue Distribution' in tup[0]]
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
                        for format_name, items in format_groups.items():
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Cell Type Distribution section
                if self.metadata_cell_type_distribution is not None:
                    f.write("  <h2>3. Cell Type Distribution</h2>\n")
                    
                    # Table of top cell types
                    f.write(f"  <h3>Top {min(top_n, len(self.metadata_cell_type_distribution))} Cell Types</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Cell Type</th><th>Count</th><th>Percentage</th></tr>\n")
                    
                    top_cell_types = self.metadata_cell_type_distribution.head(top_n)
                    for _, row in top_cell_types.iterrows():
                        f.write(f"    <tr><td>{row['Cell Type']}</td><td>{row['Count']}</td><td>{row['Percentage (%)']}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                    
                    # Cell type distribution charts - both bar and pie chart
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    cell_bar_png = [path for name, path in generated_files['figures'] 
                                if 'Cell Type Distribution Bar Chart' in name and path.endswith('.png')]
                    if cell_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{cell_bar_png[0]}' alt='Cell Type Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    cell_pie_png = [path for name, path in generated_files['figures'] 
                                if 'Cell Type Distribution Pie Chart' in name and path.endswith('.png')]
                    if cell_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{cell_pie_png[0]}' alt='Cell Type Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    cell_figures_all = [tup for tup in generated_files['figures'] if 'Cell Type Distribution' in tup[0]]
                    if len(cell_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        # Group by format
                        format_groups = {}
                        for name, path in cell_figures_all:
                            format_name = name.split('(')[1].split(')')[0]
                            if format_name not in format_groups:
                                format_groups[format_name] = []
                            format_groups[format_name].append((name, path))
                        
                        # Add links by format
                        for format_name, items in format_groups.items():
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Heatmap section
                heatmap_figures_png = [path for name, path in generated_files['figures'] 
                                    if 'Tissue-Cell Heatmap' in name and path.endswith('.png')]
                if heatmap_figures_png:
                    f.write("  <h2>4. Tissue Type and Cell Type Relationship</h2>\n")
                    f.write("  <p>This heatmap shows the relationship between tissue types and cell types in the dataset. The numbers in each cell represent the count of samples with both the corresponding tissue type and cell type.</p>\n")
                    f.write(f"  <img src='{heatmap_figures_png[0]}' alt='Tissue-Cell Relationship'>\n")
                    
                    # Add links to other formats
                    heatmap_figures_all = [tup for tup in generated_files['figures'] if 'Tissue-Cell Heatmap' in tup[0]]
                    if len(heatmap_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figure as: ")
                        for name, path in heatmap_figures_all:
                            format_name = name.split('(')[1].split(')')[0]  # Extract format from "Heatmap (FORMAT)"
                            f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Generated Files section
                f.write("  <h2>5. Generated Files</h2>\n")
                
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
                
                # Summary
                f.write("  <h2>6. Summary</h2>\n")
                f.write("  <p>This report provides an overview of the eRNAbase metadata distribution by tissue type and cell type.</p>\n")
                f.write("  <p>The visualizations show the distribution of samples across different biological categories using both bar charts and pie charts.</p>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by ERNARegionAnalyzer</p>\n")
                
                # End HTML document
                f.write("</body>\n")
                f.write("</html>\n")
            
            self.logger.info(f"Simple HTML report generated: {html_path}")
            return html_path
            
        except Exception as e:
            self.logger.error(f"Error generating simple HTML report: {e}")
            return None






    def generate_region_analysis_report(self, region_info=None, 
                                    output_dir='results/region_report',
                                    figure_formats=None, save_tables=True,
                                    dpi=300, fig_width=10, fig_height=6,
                                    run_enrichment=True, run_stacked_plots=True,
                                    peak_threshold=0.2, max_samples=100,
                                    pvalue_threshold=0.05, fc_threshold=1.5):
        """
        Generate an HTML report of the genomic region analysis results with figures and tables.
        
        Parameters:
        -----------
        region_info : dict, optional
            Dictionary containing region information (chrom, start, end, name).
            If None, uses the current interest_region.
        output_dir : str, optional
            Base directory for all output files (default: 'results/region_report')
        figure_formats : list, optional
            List of figure formats to save (default: ['png', 'svg', 'pdf'])
        save_tables : bool, optional
            Whether to save distribution tables as CSV files (default: True)
        dpi : int, optional
            Resolution for raster images like PNG (default: 300)
        fig_width : int, optional
            Width of the figure in inches (default: 10)
        fig_height : int, optional
            Height of the figure in inches (default: 6)
        run_enrichment : bool, optional
            Whether to run enrichment analysis if not already done (default: True)
        run_stacked_plots : bool, optional
            Whether to generate stacked read plots (default: True)
        peak_threshold : float, optional
            Threshold for peak detection in stacked plots (default: 0.2)
        max_samples : int, optional
            Maximum number of samples to display in stacked plots (default: 100)
        pvalue_threshold : float, optional
            P-value threshold for enrichment analysis (default: 0.05)
        fc_threshold : float, optional
            Fold change threshold for enrichment analysis (default: 1.5)
                
        Returns:
        --------
        str
            Path to the generated HTML report
        """
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
        
        self.logger.info(f"Creating region analysis report in directory: {output_dir}")
        
        # Dictionary to track all generated files
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Step 1: Get region information
        if region_info is None:
            if self.interest_region is None:
                self.logger.error("No region information available. Please provide region_info or run create_interest_region first.")
                return None
            region_info = self.interest_region
        
        chrom = region_info.get('chrom')
        start = region_info.get('start')
        end = region_info.get('end')
        region_name = region_info.get('name', f"{chrom}:{start}-{end}")
        
        # Basic region stats
        region_stats = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'region_name': region_name,
            'region_size': end - start,
            'timestamp': timestamp
        }
        
        # Step 2: Check if overlapping samples were analyzed
        if self.overlapping_samples is None or self.overlapping_results is None:
            self.logger.error("No overlapping samples analysis found. Run analyze_region or merge_all_overlaps first.")
            return None
        
        # Get sample stats
        sample_stats = {
            'total_samples': len(self.erna_bed_files) if self.erna_bed_files else 0,
            'overlapping_samples': len(self.overlapping_samples),
            'overlap_percentage': (len(self.overlapping_samples) / len(self.erna_bed_files) * 100) if self.erna_bed_files else 0
        }
        
        # Step 3: Analyze overlapping sample distributions if not already done
        if not hasattr(self, 'region_tissue_counts') or self.region_tissue_counts is None:
            self.analyze_overlapping_sample_distributions()
        
        # Set Matplotlib params for better quality figures
        plt.rcParams['figure.figsize'] = (fig_width, fig_height)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # Step 4: Generate visualizations for tissue and cell type distributions
        # Prepare save options
        save_opts = {}
        for fmt in figure_formats:
            save_opts[f'save_{fmt}'] = True
        
        # Step 4.1: Generate tissue distribution visualizations if available
        tissue_dist_figures = []
        if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None:
            try:
                # Bar chart
                tissue_bar_basename = f'region_{chrom}_{start}_{end}_tissue_bar_{timestamp}'
                tissue_bar_fig = self.plot_sorted_bar(
                    data=self.region_tissue_counts,
                    category='Tissue Type',
                    count_col='Region Count',
                    title=f'Tissue Type Distribution in {region_name}',
                    limit=15,
                    horizontal=True,
                    color='skyblue',
                    save_dir=figures_full_path,
                    save_filename=tissue_bar_basename,
                    **save_opts
                )
                
                # Register bar chart files
                for fmt in figure_formats:
                    figure_filename = f'{tissue_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            tissue_dist_figures.append((f'Tissue Distribution Bar Chart', figure_relative_path))
                
                # Pie chart
                tissue_pie_basename = f'region_{chrom}_{start}_{end}_tissue_pie_{timestamp}'
                tissue_pie_fig = self.plot_simple_pie(
                    data=self.region_tissue_counts,
                    category='Tissue Type',
                    title=f'Tissue Type Distribution in {region_name}',
                    limit=10,
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{tissue_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            tissue_dist_figures.append((f'Tissue Distribution Pie Chart', figure_relative_path))
                
                # Save tissue distribution table
                if save_tables:
                    tissue_table_filename = f'region_{chrom}_{start}_{end}_tissue_dist_{timestamp}.csv'
                    tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
                    self.region_tissue_counts.to_csv(tissue_table_path, index=False)
                    generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating tissue distribution visualizations: {e}")
        
        # Step 4.2: Generate cell type distribution visualizations if available
        cell_dist_figures = []
        if hasattr(self, 'region_cell_type_counts') and self.region_cell_type_counts is not None:
            try:
                # Bar chart
                cell_bar_basename = f'region_{chrom}_{start}_{end}_cell_bar_{timestamp}'
                cell_bar_fig = self.plot_sorted_bar(
                    data=self.region_cell_type_counts,
                    category='Cell Type',
                    count_col='Region Count',
                    title=f'Cell Type Distribution in {region_name}',
                    limit=15,
                    horizontal=True,
                    color='lightgreen',
                    save_dir=figures_full_path,
                    save_filename=cell_bar_basename,
                    **save_opts
                )
                
                # Register bar chart files
                for fmt in figure_formats:
                    figure_filename = f'{cell_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_dist_figures.append((f'Cell Type Distribution Bar Chart', figure_relative_path))
                
                # Pie chart
                cell_pie_basename = f'region_{chrom}_{start}_{end}_cell_pie_{timestamp}'
                cell_pie_fig = self.plot_simple_pie(
                    data=self.region_cell_type_counts,
                    category='Cell Type',
                    title=f'Cell Type Distribution in {region_name}',
                    limit=10,
                    save_dir=figures_full_path,
                    save_filename=cell_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{cell_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_dist_figures.append((f'Cell Type Distribution Pie Chart', figure_relative_path))
                
                # Save cell type distribution table
                if save_tables:
                    cell_table_filename = f'region_{chrom}_{start}_{end}_cell_dist_{timestamp}.csv'
                    cell_table_path = os.path.join(tables_full_path, cell_table_filename)
                    self.region_cell_type_counts.to_csv(cell_table_path, index=False)
                    generated_files['tables'].append(('Cell Type Distribution', os.path.join(tables_dir, cell_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating cell type distribution visualizations: {e}")
        
        # Step 5: Run enrichment analysis if requested or get existing results
        tissue_enrichment_figures = []
        cell_enrichment_figures = []
        
        # Step 5.1: Tissue enrichment analysis
        if run_enrichment or hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None:
            try:
                # Run analysis if needed and not already done
                if run_enrichment and (not hasattr(self, 'region_tissue_enrichment') or self.region_tissue_enrichment is None):
                    self.analyze_tissue_enrichment(plot=False, correction='fdr_bh')
                
                if hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None:
                    # Generate horizontal bar chart
                    tissue_enrich_basename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}'
                    tissue_enrich_fig = self.visualize_fdr_horizontal_bar(
                        enrichment_results=self.region_tissue_enrichment,
                        title=f'Tissue Type Enrichment Analysis - {region_name}',
                        region_name=region_name,
                        pvalue_threshold=pvalue_threshold,
                        fc_threshold=fc_threshold,
                        save_dir=figures_full_path,
                        save_filename=tissue_enrich_basename,
                        **save_opts
                    )
                    
                    # Register files
                    for fmt in figure_formats:
                        figure_filename = f'{tissue_enrich_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Enrichment Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_enrichment_figures.append((f'Tissue Enrichment Analysis', figure_relative_path))
                    
                    # Save enrichment table
                    if save_tables:
                        tissue_enrich_filename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}.csv'
                        tissue_enrich_path = os.path.join(tables_full_path, tissue_enrich_filename)
                        self.region_tissue_enrichment.to_csv(tissue_enrich_path, index=False)
                        generated_files['tables'].append(('Tissue Enrichment Analysis', os.path.join(tables_dir, tissue_enrich_filename)))
            except Exception as e:
                self.logger.warning(f"Error generating tissue enrichment analysis: {e}")
        
        # Step 5.2: Cell type enrichment analysis
        if run_enrichment or hasattr(self, 'region_cell_type_enrichment') and self.region_cell_type_enrichment is not None:
            try:
                # Run analysis if needed and not already done
                if run_enrichment and (not hasattr(self, 'region_cell_type_enrichment') or self.region_cell_type_enrichment is None):
                    self.analyze_cell_type_enrichment(plot=False, correction='fdr_bh')
                
                if hasattr(self, 'region_cell_type_enrichment') and self.region_cell_type_enrichment is not None:
                    # Generate horizontal bar chart
                    cell_enrich_basename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}'
                    cell_enrich_fig = self.visualize_fdr_horizontal_bar(
                        enrichment_results=self.region_cell_type_enrichment,
                        title=f'Cell Type Enrichment Analysis - {region_name}',
                        region_name=region_name,
                        pvalue_threshold=pvalue_threshold,
                        fc_threshold=fc_threshold,
                        save_dir=figures_full_path,
                        save_filename=cell_enrich_basename,
                        **save_opts
                    )
                    
                    # Register files
                    for fmt in figure_formats:
                        figure_filename = f'{cell_enrich_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Cell Type Enrichment Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                cell_enrichment_figures.append((f'Cell Type Enrichment Analysis', figure_relative_path))
                    
                    # Save enrichment table
                    if save_tables:
                        cell_enrich_filename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}.csv'
                        cell_enrich_path = os.path.join(tables_full_path, cell_enrich_filename)
                        self.region_cell_type_enrichment.to_csv(cell_enrich_path, index=False)
                        generated_files['tables'].append(('Cell Type Enrichment Analysis', os.path.join(tables_dir, cell_enrich_filename)))
            except Exception as e:
                self.logger.warning(f"Error generating cell type enrichment analysis: {e}")
        
        # Step 6: Generate stacked read plots if requested
        stacked_plot_figures = []
        peak_data = None
        
        if run_stacked_plots:
            try:
                # Simple stacked plot
                simple_stacked_basename = f'region_{chrom}_{start}_{end}_simple_stacked_{timestamp}'
                simple_stacked_fig = self.plot_stacked_reads_simple(
                    chrom=chrom,
                    start=start,
                    end=end,
                    title=f"eRNA Read Distribution in {region_name}",
                    save_dir=figures_full_path,
                    save_filename=simple_stacked_basename,
                    show_plot=True,
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
                
                # Stacked plot with peaks
                peak_stacked_basename = f'region_{chrom}_{start}_{end}_peak_stacked_{timestamp}'
                peak_stacked_fig, peak_data = self.plot_stacked_reads_with_peaks(
                    chrom=chrom,
                    start=start,
                    end=end,
                    title=f"eRNA Read Distribution with Peaks in {region_name}",
                    peak_threshold=peak_threshold,
                    max_samples=max_samples,
                    save_dir=figures_full_path,
                    save_filename=peak_stacked_basename,
                    show_plot=True,
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
                self.logger.warning(f"Error generating stacked read plots: {e}")
        
        # Step 7: Generate HTML report
        html_filename = f"region_analysis_report_{chrom}_{start}_{end}_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        try:
            with open(html_path, 'w', encoding='utf-8') as f:
                # Start HTML document
                f.write("<!DOCTYPE html>\n")
                f.write("<html>\n")
                f.write("<head>\n")
                f.write("  <title>eRNAbase Region Analysis Report</title>\n")
                f.write("  <style>\n")
                f.write("    body { font-family: Arial, sans-serif; margin: 20px; }\n")
                f.write("    h1 { color: #2c3e50; }\n")
                f.write("    h2 { color: #3498db; margin-top: 30px; }\n")
                f.write("    h3 { color: #2980b9; margin-top: 20px; }\n")
                f.write("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
                f.write("    th { background-color: #3498db; color: white; text-align: left; padding: 8px; }\n")
                f.write("    td { border: 1px solid #ddd; padding: 8px; }\n")
                f.write("    tr:nth-child(even) { background-color: #f2f2f2; }\n")
                f.write("    img { max-width: 100%; height: auto; margin: 20px 0; }\n")
                f.write("    .file-formats { margin: 10px 0; }\n")
                f.write("    .file-formats a { margin-right: 10px; }\n")
                f.write("    .chart-section { display: flex; flex-wrap: wrap; justify-content: space-around; }\n")
                f.write("    .chart-container { margin: 15px; max-width: 48%; }\n")
                f.write("    .highlight { background-color: #ffffcc; padding: 10px; border-left: 4px solid #ffd700; }\n")
                f.write("    .notes { font-style: italic; color: #777; }\n")
                f.write("  </style>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                
                # Title and timestamp
                f.write(f"  <h1>eRNAbase Region Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Section 1: Region Information
                f.write("  <h2>1. Region Information</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write(f"    <p><strong>Chromosome:</strong> {chrom}</p>\n")
                f.write(f"    <p><strong>Start Position:</strong> {start:,}</p>\n")
                f.write(f"    <p><strong>End Position:</strong> {end:,}</p>\n")
                f.write(f"    <p><strong>Region Name:</strong> {region_name}</p>\n")
                f.write(f"    <p><strong>Region Size:</strong> {region_stats['region_size']:,} bp</p>\n")
                f.write("  </div>\n")
                
                # Section 2: Sample Overview
                f.write("  <h2>2. Sample Overview</h2>\n")
                f.write("  <table>\n")
                f.write("    <tr><th>Metric</th><th>Value</th><th>Percentage</th></tr>\n")
                f.write(f"    <tr><td>Total Samples Analyzed</td><td>{sample_stats['total_samples']:,}</td><td>100%</td></tr>\n")
                f.write(f"    <tr><td>Samples Overlapping with Region</td><td>{sample_stats['overlapping_samples']:,}</td><td>{sample_stats['overlap_percentage']:.2f}%</td></tr>\n")
                f.write("  </table>\n")
                
                # Sample list preview
                if self.overlapping_samples:
                    f.write("  <h3>Sample Preview</h3>\n")
                    f.write("  <p>Examples of samples overlapping with the region:</p>\n")
                    f.write("  <ul>\n")
                    for sample in self.overlapping_samples[:min(5, len(self.overlapping_samples))]:
                        f.write(f"    <li>{sample}</li>\n")
                    f.write("  </ul>\n")
                    
                    if len(self.overlapping_samples) > 5:
                        f.write(f"  <p class='notes'>... and {len(self.overlapping_samples) - 5} more samples.</p>\n")
                    
                    # Add link to full sample list
                    sample_files = [path for name, path in generated_files['tables'] if 'Samples' in name and not 'Peak' in name]
                    if sample_files:
                        f.write(f"  <p><a href='{sample_files[0]}'>Download complete sample list</a></p>\n")
                
                # Section 3: Distribution Analysis
                if hasattr(self, 'region_tissue_counts') or hasattr(self, 'region_cell_type_counts'):
                    f.write("  <h2>3. Tissue and Cell Type Distribution</h2>\n")
                    
                    # Tissue type distribution
                    if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None and not self.region_tissue_counts.empty:
                        f.write("  <h3>3.1 Tissue Type Distribution</h3>\n")
                        
                        # Table of top tissues
                        f.write(f"  <h4>Top Tissue Types in Overlapping Samples</h4>\n")
                        f.write("  <table>\n")
                        f.write("    <tr><th>Tissue Type</th><th>Count</th><th>Percentage</th></tr>\n")
                        
                        top_tissues = self.region_tissue_counts.sort_values('Region Count', ascending=False).head(10)
                        for _, row in top_tissues.iterrows():
                            percentage = row['Region Percentage'] if 'Region Percentage' in row else (row['Region Count'] / self.overlap_sample_count * 100) if hasattr(self, 'overlap_sample_count') else 0
                            f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Region Count']}</td><td>{percentage:.2f}%</td></tr>\n")
                        
                        f.write("  </table>\n")
                        
                        # Display tissue distribution charts
                        if tissue_dist_figures:
                            f.write("  <div class='chart-section'>\n")
                            for figure_name, figure_path in tissue_dist_figures:
                                f.write("    <div class='chart-container'>\n")
                                f.write(f"      <h4>{figure_name}</h4>\n")
                                f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                                f.write("    </div>\n")
                            f.write("  </div>\n")
                            

                            # # 差し替え
                            # # Add download links for other formats
                            # f.write("  <div class='file-formats'>Download figures as: ")
                            # for name, path in generated_files['figures']:
                            #     if 'Tissue Distribution' in name:
                            #         format_name = name.split('(')[1].split(')')[0]
                            #         chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                            #         f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                            # f.write("  </div>\n")
                    

                            # データテーブルファイルへのリンクを取得
                            tissue_data_files = [path for name, path in generated_files['tables'] 
                                            if 'Tissue Distribution' in name or 'tissue_dist' in path.lower()]

                            # 画像リンクとデータリンクを組み合わせて表示
                            f.write("  <div class='file-formats'>Download: ")

                            # 画像形式ごとのリンク
                            f.write("<strong>Figures:</strong> ")
                            for name, path in generated_files['figures']:
                                if 'Tissue Distribution' in name:
                                    format_name = name.split('(')[1].split(')')[0]
                                    chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                    f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")

                            # データファイルへのリンク
                            if tissue_data_files:
                                f.write("<strong>Data:</strong> ")
                                for path in tissue_data_files:
                                    file_ext = path.split('.')[-1].upper()
                                    f.write(f"<a href='{path}'>{file_ext}</a> ")

                            f.write("  </div>\n")




                    # Cell type distribution
                    if hasattr(self, 'region_cell_type_counts') and self.region_cell_type_counts is not None and not self.region_cell_type_counts.empty:
                        f.write("  <h3>3.2 Cell Type Distribution</h3>\n")
                        
                        # Table of top cell types
                        f.write(f"  <h4>Top Cell Types in Overlapping Samples</h4>\n")
                        f.write("  <table>\n")
                        f.write("    <tr><th>Cell Type</th><th>Count</th><th>Percentage</th></tr>\n")
                        
                        top_cells = self.region_cell_type_counts.sort_values('Region Count', ascending=False).head(10)
                        for _, row in top_cells.iterrows():
                            percentage = row['Region Percentage'] if 'Region Percentage' in row else (row['Region Count'] / self.overlap_sample_count * 100) if hasattr(self, 'overlap_sample_count') else 0
                            f.write(f"    <tr><td>{row['Cell Type']}</td><td>{row['Region Count']}</td><td>{percentage:.2f}%</td></tr>\n")
                        
                        f.write("  </table>\n")
                        
                        # Display cell type distribution charts
                        if cell_dist_figures:
                            f.write("  <div class='chart-section'>\n")
                            for figure_name, figure_path in cell_dist_figures:
                                f.write("    <div class='chart-container'>\n")
                                f.write(f"      <h4>{figure_name}</h4>\n")
                                f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                                f.write("    </div>\n")
                            f.write("  </div>\n")
                            


                            # #　差し替え
                            # # Add download links for other formats
                            # f.write("  <div class='file-formats'>Download figures as: ")
                            # for name, path in generated_files['figures']:
                            #     if 'Cell Type Distribution' in name:
                            #         format_name = name.split('(')[1].split(')')[0]
                            #         chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                            #         f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                            # f.write("  </div>\n")
                


                            # データテーブルファイルへのリンクを取得
                            cell_data_files = [path for name, path in generated_files['tables'] 
                                            if 'Cell Type Distribution' in name or 'cell_dist' in path.lower()]

                            # 画像リンクとデータリンクを組み合わせて表示
                            f.write("  <div class='file-formats'>Download: ")

                            # 画像形式ごとのリンク
                            f.write("<strong>Figures:</strong> ")
                            for name, path in generated_files['figures']:
                                if 'Cell Type Distribution' in name:
                                    format_name = name.split('(')[1].split(')')[0]
                                    chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                    f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")

                            # データファイルへのリンク
                            if cell_data_files:
                                f.write("<strong>Data:</strong> ")
                                for path in cell_data_files:
                                    file_ext = path.split('.')[-1].upper()
                                    f.write(f"<a href='{path}'>{file_ext}</a> ")

                            f.write("  </div>\n")


                # Section 4: Enrichment Analysis
                if hasattr(self, 'region_tissue_enrichment') or hasattr(self, 'region_cell_type_enrichment'):
                    f.write("  <h2>4. Enrichment Analysis</h2>\n")
                    
                    # Tissue enrichment analysis
                    if hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None and not self.region_tissue_enrichment.empty:
                        f.write("  <h3>4.1 Tissue Type Enrichment</h3>\n")
                        
                        # Enrichment summary
                        significant_tissues = sum(self.region_tissue_enrichment['p-value'] < pvalue_threshold)
                        f.write(f"  <p>This analysis shows which tissue types are significantly enriched or depleted in the region compared to the background distribution.</p>\n")
                        f.write(f"  <p><strong>Significant tissues:</strong> {significant_tissues} out of {len(self.region_tissue_enrichment)} tissue types (p < {pvalue_threshold}).</p>\n")
                        
                        # Top significant tissues
                        if significant_tissues > 0:
                            sig_tissues = self.region_tissue_enrichment[self.region_tissue_enrichment['p-value'] < pvalue_threshold].sort_values('p-value')
                            f.write("  <h4>Top Significant Tissues</h4>\n")
                            f.write("  <table>\n")
                            f.write("    <tr><th>Tissue Type</th><th>Fold Change</th><th>p-value</th><th>Adjusted p-value</th><th>Enrichment</th></tr>\n")
                            
                            for i, (_, row) in enumerate(sig_tissues.head(5).iterrows()):
                                enrichment_type = "Enriched" if row['Fold Change'] > 1 else "Depleted"
                                f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Fold Change']:.2f}x</td>" +
                                    f"<td>{row['p-value']:.2e}</td><td>{row['adjusted p-value']:.2e}</td>" +
                                    f"<td>{enrichment_type}</td></tr>\n")
                            
                            f.write("  </table>\n")
                        
                        # Display enrichment visualizations
                        if tissue_enrichment_figures:
                            f.write("  <div class='chart-section'>\n")
                            for figure_name, figure_path in tissue_enrichment_figures:
                                f.write("    <div class='chart-container' style='max-width: 90%;'>\n")
                                f.write(f"      <h4>{figure_name}</h4>\n")
                                f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                                f.write("    </div>\n")
                            f.write("  </div>\n")
                            


                            # # 差し替え
                            # # Add download links for other formats
                            # f.write("  <div class='file-formats'>Download figures as: ")
                            # for name, path in generated_files['figures']:
                            #     if 'Tissue Enrichment' in name:
                            #         format_name = name.split('(')[1].split(')')[0]
                            #         f.write(f"<a href='{path}'>{format_name}</a> ")
                            # f.write("  </div>\n")



                            # データテーブルファイルへのリンクを取得
                            tissue_enrich_files = [path for name, path in generated_files['tables'] 
                                                if 'Tissue Enrichment' in name or 'tissue_enrichment' in path.lower()]

                            # 画像リンクとデータリンクを組み合わせて表示
                            f.write("  <div class='file-formats'>Download: ")

                            # 画像形式ごとのリンク
                            f.write("<strong>Figures:</strong> ")
                            for name, path in generated_files['figures']:
                                if 'Tissue Enrichment' in name:
                                    format_name = name.split('(')[1].split(')')[0]
                                    f.write(f"<a href='{path}'>{format_name}</a> ")

                            # データファイルへのリンク
                            if tissue_enrich_files:
                                f.write("<strong>Data:</strong> ")
                                for path in tissue_enrich_files:
                                    file_ext = path.split('.')[-1].upper()
                                    f.write(f"<a href='{path}'>{file_ext}</a> ")

                            f.write("  </div>\n")


                    # Cell type enrichment analysis
                    if hasattr(self, 'region_cell_type_enrichment') and self.region_cell_type_enrichment is not None and not self.region_cell_type_enrichment.empty:
                        f.write("  <h3>4.2 Cell Type Enrichment</h3>\n")
                        
                        # Enrichment summary
                        significant_cells = sum(self.region_cell_type_enrichment['p-value'] < pvalue_threshold)
                        f.write(f"  <p>This analysis shows which cell types are significantly enriched or depleted in the region compared to the background distribution.</p>\n")
                        f.write(f"  <p><strong>Significant cell types:</strong> {significant_cells} out of {len(self.region_cell_type_enrichment)} cell types (p < {pvalue_threshold}).</p>\n")
                        
                        # Top significant cell types
                        if significant_cells > 0:
                            sig_cells = self.region_cell_type_enrichment[self.region_cell_type_enrichment['p-value'] < pvalue_threshold].sort_values('p-value')
                            f.write("  <h4>Top Significant Cell Types</h4>\n")
                            f.write("  <table>\n")
                            f.write("    <tr><th>Cell Type</th><th>Fold Change</th><th>p-value</th><th>Adjusted p-value</th><th>Enrichment</th></tr>\n")
                            
                            for i, (_, row) in enumerate(sig_cells.head(5).iterrows()):
                                enrichment_type = "Enriched" if row['Fold Change'] > 1 else "Depleted"
                                f.write(f"    <tr><td>{row['Cell Type']}</td><td>{row['Fold Change']:.2f}x</td>" +
                                    f"<td>{row['p-value']:.2e}</td><td>{row['adjusted p-value']:.2e}</td>" +
                                    f"<td>{enrichment_type}</td></tr>\n")
                            
                            f.write("  </table>\n")
                        
                        # Display enrichment visualizations
                        if cell_enrichment_figures:
                            f.write("  <div class='chart-section'>\n")
                            for figure_name, figure_path in cell_enrichment_figures:
                                f.write("    <div class='chart-container' style='max-width: 90%;'>\n")
                                f.write(f"      <h4>{figure_name}</h4>\n")
                                f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                                f.write("    </div>\n")
                            f.write("  </div>\n")
                            


                            # # 差し替え
                            # # Add download links for other formats
                            # f.write("  <div class='file-formats'>Download figures as: ")
                            # for name, path in generated_files['figures']:
                            #     if 'Cell Type Enrichment' in name:
                            #         format_name = name.split('(')[1].split(')')[0]
                            #         f.write(f"<a href='{path}'>{format_name}</a> ")
                            # f.write("  </div>\n")
                
                            # データテーブルファイルへのリンクを取得
                            cell_enrich_files = [path for name, path in generated_files['tables'] 
                                            if 'Cell Type Enrichment' in name or 'cell_enrichment' in path.lower()]

                            # 画像リンクとデータリンクを組み合わせて表示
                            f.write("  <div class='file-formats'>Download: ")

                            # 画像形式ごとのリンク
                            f.write("<strong>Figures:</strong> ")
                            for name, path in generated_files['figures']:
                                if 'Cell Type Enrichment' in name:
                                    format_name = name.split('(')[1].split(')')[0]
                                    f.write(f"<a href='{path}'>{format_name}</a> ")

                            # データファイルへのリンク
                            if cell_enrich_files:
                                f.write("<strong>Data:</strong> ")
                                for path in cell_enrich_files:
                                    file_ext = path.split('.')[-1].upper()
                                    f.write(f"<a href='{path}'>{file_ext}</a> ")

                            f.write("  </div>\n")


                # Section 5: Read Distribution Analysis
                if stacked_plot_figures:
                    f.write("  <h2>5. Read Distribution Analysis</h2>\n")
                    f.write("  <p>These plots show the distribution of eRNA reads across the genomic region.</p>\n")
                    
                    # Display each stacked plot
                    for figure_name, figure_path in stacked_plot_figures:
                        f.write("  <div class='chart-container' style='max-width: 100%;'>\n")
                        f.write(f"    <h3>{figure_name}</h3>\n")
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("  </div>\n")
                        
                        # Add download links for other formats
                        f.write("  <div class='file-formats'>Download figure as: ")
                        for name, path in generated_files['figures']:
                            if (figure_name == 'Simple Read Distribution' and 'Simple Stacked' in name) or \
                            (figure_name == 'Read Distribution with Peaks' and 'Stacked Read Plot with Peaks' in name):
                                format_name = name.split('(')[1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                    
                    # Peak information if available
                    peak_data_files = [path for name, path in generated_files['tables'] if 'Peak Analysis Data' in name]
                    if peak_data_files:
                        f.write("  <h3>5.1 Peak Analysis</h3>\n")
                        f.write(f"  <p>Detected peaks in the read distribution. <a href='{peak_data_files[0]}'>Download complete peak data</a>.</p>\n")
                        
                        # If peak data is available, show summary
                        if peak_data:
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
                            f.write("  <h4>Samples in Each Peak</h4>\n")
                            f.write("  <ul>\n")
                            for peak_id in peak_data.keys():
                                peak_sample_files = [path for name, path in generated_files['tables'] 
                                                if peak_id in name and 'samples' in name.lower()]
                                if peak_sample_files:
                                    f.write(f"    <li><a href='{peak_sample_files[0]}'>{peak_id} - {peak_data[peak_id]['samples']['count']} samples</a></li>\n")
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
                
                # Summary
                f.write("  <h2>7. Summary</h2>\n")
                f.write("  <p>This report provides an analysis of eRNAbase data for the genomic region " +
                    f"{chrom}:{start:,}-{end:,} ({region_name}).</p>\n")
                
                # Key findings
                f.write("  <h3>Key Findings</h3>\n")
                f.write("  <ul>\n")
                f.write(f"    <li>{sample_stats['overlapping_samples']} samples ({sample_stats['overlap_percentage']:.2f}%) " +
                    f"from eRNAbase overlap with this region.</li>\n")
                
                # Add tissue findings if available
                if hasattr(self, 'region_tissue_enrichment') and self.region_tissue_enrichment is not None:
                    sig_tissues = sum(self.region_tissue_enrichment['p-value'] < pvalue_threshold)
                    if sig_tissues > 0:
                        f.write(f"    <li>{sig_tissues} tissue types are significantly enriched or depleted in this region.</li>\n")
                        
                        # Add top enriched tissue
                        top_tissue = self.region_tissue_enrichment[self.region_tissue_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                        if not top_tissue.empty:
                            tissue = top_tissue.iloc[0]['Tissue Type']
                            fold_change = top_tissue.iloc[0]['Fold Change']
                            f.write(f"    <li>Most significantly enriched tissue: {tissue} ({fold_change:.2f}x enrichment)</li>\n")
                
                # Add cell type findings if available
                if hasattr(self, 'region_cell_type_enrichment') and self.region_cell_type_enrichment is not None:
                    sig_cells = sum(self.region_cell_type_enrichment['p-value'] < pvalue_threshold)
                    if sig_cells > 0:
                        f.write(f"    <li>{sig_cells} cell types are significantly enriched or depleted in this region.</li>\n")
                        
                        # Add top enriched cell type
                        top_cell = self.region_cell_type_enrichment[self.region_cell_type_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                        if not top_cell.empty:
                            cell = top_cell.iloc[0]['Cell Type']
                            fold_change = top_cell.iloc[0]['Fold Change']
                            f.write(f"    <li>Most significantly enriched cell type: {cell} ({fold_change:.2f}x enrichment)</li>\n")
                
                # Add peak findings if available
                if peak_data:
                    f.write(f"    <li>{len(peak_data)} peaks of eRNA activity were detected in this region.</li>\n")
                    
                    # Highlight highest density peak
                    peak_densities = [(peak_id, peak_info['density']) for peak_id, peak_info in peak_data.items()]
                    if peak_densities:
                        highest_peak_id, highest_density = max(peak_densities, key=lambda x: x[1])
                        peak_pos = peak_data[highest_peak_id]['position']
                        f.write(f"    <li>Highest eRNA activity density at {highest_peak_id} (position {peak_pos:,}) " +
                            f"with {highest_density:.3f} reads/bp.</li>\n")
                    
                f.write("  </ul>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by ERNARegionAnalyzer</p>\n")
                
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





    def analyze_region_comprehensive(self, chrom, start, end, region_name=None, 
                                max_files=None, output_dir=None,
                                run_enrichment=True, peak_threshold=0.2, 
                                enrichment_correction='fdr_bh'):
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
        max_files : int, optional
            Maximum number of files to analyze (for debugging)
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
            interest_region, interest_bed, output_bed = self.create_interest_region(
                chrom, start, end, region_name=region_name
            )
            results['interest_region'] = interest_region
            results['steps_completed'].append('create_interest_region')
            
            merged_overlaps, overlapping_samples, overlap_sample_count = self.merge_all_overlaps(max_files=max_files)
            results['overlapping_results'] = merged_overlaps
            results['overlapping_samples'] = overlapping_samples
            results['overlapping_samples_count'] = overlap_sample_count
            results['steps_completed'].append('merge_all_overlaps')
            
            self.logger.info(f"Found {overlap_sample_count} samples overlapping with region")
            
            # Save overlap results
            if not merged_overlaps.empty:
                saved_files = self.save_overlap_results(output_dir=output_dir, display_samples=False)
                results['saved_files'] = saved_files
                results['steps_completed'].append('save_overlap_results')
        except Exception as e:
            self.logger.error(f"Error in overlap detection: {e}")
            results['errors'] = results.get('errors', []) + [f"Overlap detection: {str(e)}"]
        
        # Step 2: Analyze metadata for overlapping samples
        try:
            if self.erna_metadata is not None and overlapping_samples:
                overlapping_metadata = self.analyze_overlapping_samples_with_metadata()
                results['overlapping_metadata'] = overlapping_metadata
                results['steps_completed'].append('analyze_overlapping_samples_with_metadata')
                
                # Analyze sample distributions
                self.analyze_overlapping_sample_distributions()
                
                if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None:
                    results['tissue_counts'] = self.region_tissue_counts
                    results['steps_completed'].append('analyze_tissue_distribution')
                
                if hasattr(self, 'region_cell_type_counts') and self.region_cell_type_counts is not None:
                    results['cell_type_counts'] = self.region_cell_type_counts
                    results['steps_completed'].append('analyze_cell_type_distribution')
        except Exception as e:
            self.logger.error(f"Error in metadata analysis: {e}")
            results['errors'] = results.get('errors', []) + [f"Metadata analysis: {str(e)}"]
        
        # Step 3: Run enrichment analysis if requested
        if run_enrichment:
            try:
                # Make sure distributions are calculated
                if not hasattr(self, '_metadata_tissue_dist') or self._metadata_tissue_dist is None:
                    self.calculate_metadata_distributions()
                
                # Tissue type enrichment
                tissue_enrichment = self.analyze_tissue_enrichment(
                    correction=enrichment_correction,
                    plot=False
                )
                
                if tissue_enrichment is not None:
                    results['tissue_enrichment'] = tissue_enrichment
                    results['steps_completed'].append('analyze_tissue_enrichment')
                
                # Cell type enrichment
                cell_enrichment = self.analyze_cell_type_enrichment(
                    correction=enrichment_correction,
                    plot=False
                )
                
                if cell_enrichment is not None:
                    results['cell_enrichment'] = cell_enrichment
                    results['steps_completed'].append('analyze_cell_type_enrichment')
            except Exception as e:
                self.logger.error(f"Error in enrichment analysis: {e}")
                results['errors'] = results.get('errors', []) + [f"Enrichment analysis: {str(e)}"]
        
        # Step 4: Generate read distribution plots
        try:
            # Simple stacked plot result is just the figure, no data returned
            simple_plot = self.plot_stacked_reads_simple(
                chrom=chrom,
                start=start,
                end=end,
                title=f"eRNA Read Distribution in {results['region_info']['name']}",
                show_plot=False
            )
            results['steps_completed'].append('plot_stacked_reads_simple')
            
            # Peak detection plot returns both figure and peak data
            # 修正：max_samplesパラメータを削除
            peak_plot, peak_data = self.plot_stacked_reads_with_peaks(
                chrom=chrom,
                start=start,
                end=end,
                title=f"eRNA Read Distribution with Peaks in {results['region_info']['name']}",
                peak_threshold=peak_threshold,
                show_plot=False
            )
            
            if peak_data:
                results['peak_data'] = peak_data
                results['steps_completed'].append('plot_stacked_reads_with_peaks')
        except Exception as e:
            self.logger.error(f"Error in read distribution analysis: {e}")
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
        if 'tissue_counts' in analysis_results and analysis_results['tissue_counts'] is not None:
            try:
                tissue_data = analysis_results['tissue_counts']
                
                # Bar chart
                tissue_bar_basename = f'region_{chrom}_{start}_{end}_tissue_bar_{timestamp}'
                tissue_bar_fig = self.plot_sorted_bar(
                    data=tissue_data,
                    category='Tissue Type',
                    count_col='Region Count',
                    title=f'Tissue Type Distribution in {region_name}',
                    limit=15,
                    horizontal=True,
                    color='skyblue',
                    save_dir=figures_full_path,
                    save_filename=tissue_bar_basename,
                    **save_opts
                )
                
                # Register bar chart files
                for fmt in figure_formats:
                    figure_filename = f'{tissue_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            tissue_dist_figures.append((f'Tissue Distribution Bar Chart', figure_relative_path))
                
                # Pie chart
                tissue_pie_basename = f'region_{chrom}_{start}_{end}_tissue_pie_{timestamp}'
                tissue_pie_fig = self.plot_simple_pie(
                    data=tissue_data,
                    category='Tissue Type',
                    title=f'Tissue Type Distribution in {region_name}',
                    limit=10,
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{tissue_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            tissue_dist_figures.append((f'Tissue Distribution Pie Chart', figure_relative_path))
                
                # Save tissue distribution table
                if save_tables:
                    tissue_table_filename = f'region_{chrom}_{start}_{end}_tissue_dist_{timestamp}.csv'
                    tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
                    tissue_data.to_csv(tissue_table_path, index=False)
                    generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating tissue distribution visualizations: {e}")
        
        # 2. Generate cell type distribution visualizations if available
        cell_dist_figures = []
        if 'cell_type_counts' in analysis_results and analysis_results['cell_type_counts'] is not None:
            try:
                cell_data = analysis_results['cell_type_counts']
                
                # Bar chart
                cell_bar_basename = f'region_{chrom}_{start}_{end}_cell_bar_{timestamp}'
                cell_bar_fig = self.plot_sorted_bar(
                    data=cell_data,
                    category='Cell Type',
                    count_col='Region Count',
                    title=f'Cell Type Distribution in {region_name}',
                    limit=15,
                    horizontal=True,
                    color='lightgreen',
                    save_dir=figures_full_path,
                    save_filename=cell_bar_basename,
                    **save_opts
                )
                
                # Register bar chart files
                for fmt in figure_formats:
                    figure_filename = f'{cell_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_dist_figures.append((f'Cell Type Distribution Bar Chart', figure_relative_path))
                
                # Pie chart
                cell_pie_basename = f'region_{chrom}_{start}_{end}_cell_pie_{timestamp}'
                cell_pie_fig = self.plot_simple_pie(
                    data=cell_data,
                    category='Cell Type',
                    title=f'Cell Type Distribution in {region_name}',
                    limit=10,
                    save_dir=figures_full_path,
                    save_filename=cell_pie_basename,
                    **save_opts
                )
                
                # Register pie chart files
                for fmt in figure_formats:
                    figure_filename = f'{cell_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_dist_figures.append((f'Cell Type Distribution Pie Chart', figure_relative_path))
                
                # Save cell type distribution table
                if save_tables:
                    cell_table_filename = f'region_{chrom}_{start}_{end}_cell_dist_{timestamp}.csv'
                    cell_table_path = os.path.join(tables_full_path, cell_table_filename)
                    cell_data.to_csv(cell_table_path, index=False)
                    generated_files['tables'].append(('Cell Type Distribution', os.path.join(tables_dir, cell_table_filename)))
            
            except Exception as e:
                self.logger.warning(f"Error generating cell type distribution visualizations: {e}")
        

        # 3. Generate tissue enrichment visualizations
        tissue_enrichment_figures = []
        if 'tissue_enrichment' in analysis_results and analysis_results['tissue_enrichment'] is not None:
            try:
                tissue_enrichment = analysis_results['tissue_enrichment']
                
                # Generate horizontal bar chart
                tissue_enrich_basename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}'
                
                # 修正: visualize_fdr_horizontal_barに保存オプションを追加
                tissue_enrich_fig = self.visualize_fdr_horizontal_bar(
                    enrichment_results=tissue_enrichment,
                    title=f'Tissue Type Enrichment Analysis - {region_name}',
                    region_name=region_name,
                    pvalue_threshold=pvalue_threshold,
                    fc_threshold=fc_threshold,
                    # 追加: 保存オプション
                    save_dir=figures_full_path,
                    save_filename=tissue_enrich_basename,
                    save_formats=figure_formats,
                    save_dpi=dpi
                )
                
                # 画像ファイルをgenerated_filesに登録（保存は関数内で済んでいる）
                for fmt in figure_formats:
                    figure_filename = f'{tissue_enrich_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Enrichment Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            tissue_enrichment_figures.append((f'Tissue Enrichment Analysis', figure_relative_path))
                
                # Save enrichment table
                if save_tables:
                    tissue_enrich_filename = f'region_{chrom}_{start}_{end}_tissue_enrichment_{timestamp}.csv'
                    tissue_enrich_path = os.path.join(tables_full_path, tissue_enrich_filename)
                    tissue_enrichment.to_csv(tissue_enrich_path, index=False)
                    generated_files['tables'].append(('Tissue Enrichment Analysis', os.path.join(tables_dir, tissue_enrich_filename)))
            except Exception as e:
                self.logger.warning(f"Error generating tissue enrichment analysis: {e}")

        # 4. Generate cell type enrichment visualizations
        cell_enrichment_figures = []
        if 'cell_enrichment' in analysis_results and analysis_results['cell_enrichment'] is not None:
            try:
                cell_enrichment = analysis_results['cell_enrichment']
                
                # Generate horizontal bar chart
                cell_enrich_basename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}'
                
                # 修正: visualize_fdr_horizontal_barに保存オプションを追加
                cell_enrich_fig = self.visualize_fdr_horizontal_bar(
                    enrichment_results=cell_enrichment,
                    title=f'Cell Type Enrichment Analysis - {region_name}',
                    region_name=region_name,
                    pvalue_threshold=pvalue_threshold,
                    fc_threshold=fc_threshold,
                    # 追加: 保存オプション
                    save_dir=figures_full_path,
                    save_filename=cell_enrich_basename,
                    save_formats=figure_formats,
                    save_dpi=dpi
                )
                
                # 画像ファイルをgenerated_filesに登録（保存は関数内で済んでいる）
                for fmt in figure_formats:
                    figure_filename = f'{cell_enrich_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Cell Type Enrichment Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_enrichment_figures.append((f'Cell Type Enrichment Analysis', figure_relative_path))
                
                # Save enrichment table
                if save_tables:
                    cell_enrich_filename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}.csv'
                    cell_enrich_path = os.path.join(tables_full_path, cell_enrich_filename)
                    cell_enrichment.to_csv(cell_enrich_path, index=False)
                    generated_files['tables'].append(('Cell Type Enrichment Analysis', os.path.join(tables_dir, cell_enrich_filename)))
            except Exception as e:
                self.logger.warning(f"Error generating cell type enrichment analysis: {e}")












        
        # # 4. Generate cell type enrichment visualizations
        # cell_enrichment_figures = []
        # if 'cell_enrichment' in analysis_results and analysis_results['cell_enrichment'] is not None:
        #     try:
        #         cell_enrichment = analysis_results['cell_enrichment']
                
        #         # Generate horizontal bar chart
        #         cell_enrich_basename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}'
                
        #         # 修正: 既存のメソッドに合わせて引数を調整
        #         cell_enrich_fig = self.visualize_fdr_horizontal_bar(
        #             enrichment_results=cell_enrichment,
        #             title=f'Cell Type Enrichment Analysis - {region_name}',
        #             region_name=region_name,
        #             pvalue_threshold=pvalue_threshold,
        #             fc_threshold=fc_threshold
        #             # save_dir と save_filenameパラメータを削除
        #         )
                
        #         # 手動で図を保存
        #         for fmt in figure_formats:
        #             figure_filename = f'{cell_enrich_basename}.{fmt}'
        #             figure_path = os.path.join(figures_full_path, figure_filename)
        #             figure_relative_path = os.path.join(figures_dir, figure_filename)
                    
        #             # 各形式で保存
        #             cell_enrich_fig.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
        #             # 登録
        #             generated_files['figures'].append((f'Cell Type Enrichment Bar Chart ({fmt.upper()})', figure_relative_path))
        #             if fmt == 'png':
        #                 cell_enrichment_figures.append((f'Cell Type Enrichment Analysis', figure_relative_path))
                
        #         # Save enrichment table
        #         if save_tables:
        #             cell_enrich_filename = f'region_{chrom}_{start}_{end}_cell_enrichment_{timestamp}.csv'
        #             cell_enrich_path = os.path.join(tables_full_path, cell_enrich_filename)
        #             cell_enrichment.to_csv(cell_enrich_path, index=False)
        #             generated_files['tables'].append(('Cell Type Enrichment Analysis', os.path.join(tables_dir, cell_enrich_filename)))
        #     except Exception as e:
        #         self.logger.warning(f"Error generating cell type enrichment analysis: {e}")
        
        # 5. Generate stacked read plots
        stacked_plot_figures = []
        
        # 5.1 Simple stacked plot
        try:
            # Create a new plot for the report
            simple_stacked_basename = f'region_{chrom}_{start}_{end}_simple_stacked_{timestamp}'
            simple_stacked_fig = self.plot_stacked_reads_simple(
                chrom=chrom,
                start=start,
                end=end,
                title=f"eRNA Read Distribution in {region_name}",
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
        
        # 5.2 Stacked plot with peaks
        try:
            peak_stacked_basename = f'region_{chrom}_{start}_{end}_peak_stacked_{timestamp}'
            # 修正：max_samplesパラメータを削除
            peak_stacked_fig, peak_data = self.plot_stacked_reads_with_peaks(
                chrom=chrom,
                start=start,
                end=end,
                title=f"eRNA Read Distribution with Peaks in {region_name}",
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
                f.write("  <title>eRNAbase Region Analysis Report</title>\n")
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
                f.write(f"  <h1>eRNAbase Region Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Comprehensive report introduction
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This report provides a comprehensive analysis of enhancer RNA (eRNA) data for the selected genomic region. The analysis includes:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li><strong>Distribution Analysis:</strong> Tissue and cell type distribution in the region</li>\n")
                f.write("      <li><strong>Enrichment Analysis:</strong> Statistical enrichment of tissues and cell types</li>\n")
                f.write("      <li><strong>Read Distribution Analysis:</strong> Visualization of eRNA read patterns</li>\n")
                f.write("      <li><strong>Peak Detection:</strong> Identification of high eRNA activity regions</li>\n")
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
                    'total_samples': len(self.erna_bed_files) if hasattr(self, 'erna_bed_files') and self.erna_bed_files else 0,
                    'overlapping_samples': len(analysis_results.get('overlapping_samples', [])),
                }
                sample_stats['overlap_percentage'] = (sample_stats['overlapping_samples'] / sample_stats['total_samples'] * 100) if sample_stats['total_samples'] else 0
                
                f.write("  <h2>2. Sample Overview</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section summarizes the samples in the eRNAbase database that overlap with the selected genomic region.</p>\n")
                f.write("  </div>\n")
                
                f.write("  <table>\n")
                f.write("    <tr><th>Metric</th><th>Value</th><th>Percentage</th></tr>\n")
                f.write(f"    <tr><td>Total Samples Analyzed</td><td>{sample_stats['total_samples']:,}</td><td>100%</td></tr>\n")
                f.write(f"    <tr><td>Samples Overlapping with Region</td><td>{sample_stats['overlapping_samples']:,}</td><td>{sample_stats['overlap_percentage']:.2f}%</td></tr>\n")
                f.write("  </table>\n")
                
                # Section 3: Distribution Analysis
                f.write("  <h2>3. Tissue and Cell Type Distribution</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section shows the distribution of tissue types and cell types among samples that overlap with the selected genomic region.</p>\n")
                f.write("    <p>The distribution charts help identify which tissues and cell types are most frequent in this region.</p>\n")
                f.write("  </div>\n")
                
                # 3.1 Tissue type distribution
                if tissue_dist_figures:
                    f.write("  <h3>3.1 Tissue Type Distribution</h3>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>The charts below show the distribution of tissue types in samples overlapping with this region. Both bar charts and pie charts provide different views of the same distribution data:</p>\n")
                    f.write("    <ul>\n")
                    f.write("      <li><strong>Bar Charts</strong> are better for comparing the exact counts/percentages across different tissue types</li>\n")
                    f.write("      <li><strong>Pie Charts</strong> are better for visualizing the proportion of each tissue type relative to the whole</li>\n")
                    f.write("    </ul>\n")
                    f.write("    <p>This information is useful for understanding which tissues most commonly express eRNAs in this genomic region.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display tissue distribution charts
                    f.write("  <div class='chart-section'>\n")
                    for figure_name, figure_path in tissue_dist_figures:
                        f.write("    <div class='chart-container'>\n")
                        f.write(f"      <h4>{figure_name}</h4>\n")
                        f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("    </div>\n")
                    f.write("  </div>\n")
                    

                    # # 差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figures as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Tissue Distribution' in name:
                    #         format_name = name.split('(')[1].split(')')[0]
                    #         chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                    #         f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                    # f.write("  </div>\n")

                    # データテーブルファイルへのリンクを取得
                    tissue_data_files = [path for name, path in generated_files['tables'] 
                                    if 'Tissue Distribution' in name or 'tissue_dist' in path.lower()]

                    # 画像リンクとデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Tissue Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                            f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")

                    # データファイルへのリンク
                    if tissue_data_files:
                        f.write("<strong>Data:</strong> ")
                        for path in tissue_data_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")

                    f.write("  </div>\n")





                else:
                    f.write("  <h3>3.1 Tissue Type Distribution</h3>\n")
                    f.write("  <p class='notes'>No tissue type distribution data available for this region.</p>\n")
                
                # 3.2 Cell type distribution
                if cell_dist_figures:
                    f.write("  <h3>3.2 Cell Type Distribution</h3>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>These charts illustrate the cell type distribution in the overlapping samples. Specific cell types may have distinctive eRNA expression patterns due to:</p>\n")
                    f.write("    <ul>\n")
                    f.write("      <li>Cell-type specific transcription factors and regulatory mechanisms</li>\n")
                    f.write("      <li>Different chromatin organization and accessibility</li>\n")
                    f.write("      <li>Specialized cellular functions that require specific enhancer activity</li>\n")
                    f.write("    </ul>\n")
                    f.write("    <p>The predominant cell types in this region may indicate its regulatory importance in specific cellular contexts.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display cell distribution charts
                    f.write("  <div class='chart-section'>\n")
                    for figure_name, figure_path in cell_dist_figures:
                        f.write("    <div class='chart-container'>\n")
                        f.write(f"      <h4>{figure_name}</h4>\n")
                        f.write(f"      <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("    </div>\n")
                    f.write("  </div>\n")
                    

                    # # 差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figures as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Cell Type Distribution' in name:
                    #         format_name = name.split('(')[1].split(')')[0]
                    #         chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                    #         f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                    # f.write("  </div>\n")


                    # データテーブルファイルへのリンクを取得
                    cell_data_files = [path for name, path in generated_files['tables'] 
                                    if 'Cell Type Distribution' in name or 'cell_dist' in path.lower()]

                    # 画像リンクとデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Cell Type Distribution' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                            f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")

                    # データファイルへのリンク
                    if cell_data_files:
                        f.write("<strong>Data:</strong> ")
                        for path in cell_data_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")

                    f.write("  </div>\n")


                else:
                    f.write("  <h3>3.2 Cell Type Distribution</h3>\n")
                    f.write("  <p class='notes'>No cell type distribution data available for this region.</p>\n")
                
                # Section 4: Enrichment Analysis
                f.write("  <h2>4. Enrichment Analysis</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>Enrichment analysis identifies tissues and cell types that are statistically over-represented or under-represented in this region compared to the background distribution across all samples.</p>\n")
                f.write("    <p>This analysis helps to identify biological contexts where this genomic region might have specialized regulatory functions.</p>\n")
                f.write("  </div>\n")
                
                # 4.1 Tissue Enrichment Analysis
                f.write("  <h3>4.1 Tissue Type Enrichment</h3>\n")
                
                if tissue_enrichment_figures:
                    tissue_enrichment = analysis_results.get('tissue_enrichment')
                    if tissue_enrichment is not None:
                        significant_tissues = sum(tissue_enrichment['p-value'] < pvalue_threshold)
                        
                        f.write("  <div class='explanation'>\n")
                        f.write("    <p>This analysis uses Fisher's exact test to determine which tissue types are significantly enriched or depleted in the region compared to the genome-wide background distribution.</p>\n")
                        f.write("    <p><strong>How to interpret the results:</strong></p>\n")
                        f.write("    <ul>\n")
                        f.write("      <li><span class='enriched'>Enriched tissues (Fold Change > 1)</span> have more samples in this region than expected by chance</li>\n")
                        f.write("      <li><span class='depleted'>Depleted tissues (Fold Change < 1)</span> have fewer samples than expected by chance</li>\n")
                        f.write("      <li>The height of the bars indicates statistical significance (-log10 of p-value)</li>\n")
                        f.write("      <li>The dotted vertical line indicates the significance threshold (p-value = 0.05)</li>\n")
                        f.write("    </ul>\n")
                        f.write("    <p>Significant enrichment suggests that this genomic region may have tissue-specific enhancer activity, potentially regulating genes important for that tissue's function.</p>\n")
                        f.write("  </div>\n")
                        
                        f.write(f"  <p><strong>Significant tissues:</strong> {significant_tissues} out of {len(tissue_enrichment)} tissue types (p < {pvalue_threshold}).</p>\n")
                        
                        # Top significant tissues table
                        if significant_tissues > 0:
                            sig_tissues = tissue_enrichment[tissue_enrichment['p-value'] < pvalue_threshold].sort_values('p-value')
                            f.write("  <h4>Top Significant Tissues</h4>\n")
                            f.write("  <table>\n")
                            f.write("    <tr><th>Tissue Type</th><th>Fold Change</th><th>p-value</th><th>Adjusted p-value</th><th>Enrichment</th></tr>\n")
                            
                            for i, (_, row) in enumerate(sig_tissues.head(5).iterrows()):
                                enrichment_type = "Enriched" if row['Fold Change'] > 1 else "Depleted"
                                enrichment_class = "enriched" if row['Fold Change'] > 1 else "depleted"
                                f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Fold Change']:.2f}x</td>" +
                                    f"<td>{row['p-value']:.2e}</td><td>{row['adjusted p-value']:.2e}</td>" +
                                    f"<td class='{enrichment_class}'>{enrichment_type}</td></tr>\n")
                            
                            f.write("  </table>\n")
                    
                    # Display enrichment visualization
                    for figure_name, figure_path in tissue_enrichment_figures:
                        f.write("  <div class='chart-container' style='max-width: 90%; margin: 20px auto;'>\n")
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("  </div>\n")
                    

                    # # 差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figures as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Tissue Enrichment' in name:
                    #         format_name = name.split('(')[1].split(')')[0]
                    #         f.write(f"<a href='{path}'>{format_name}</a> ")
                    # f.write("  </div>\n")

                    # データテーブルファイルへのリンクを取得
                    tissue_enrich_files = [path for name, path in generated_files['tables'] 
                                        if 'Tissue Enrichment' in name or 'tissue_enrichment' in path.lower()]

                    # 画像リンクとデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Tissue Enrichment' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")

                    # データファイルへのリンク
                    if tissue_enrich_files:
                        f.write("<strong>Data:</strong> ")
                        for path in tissue_enrich_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")

                    f.write("  </div>\n")

                else:
                    f.write("  <p class='notes'>No tissue enrichment analysis results available for this region.</p>\n")
                
                # 4.2 Cell Type Enrichment Analysis
                f.write("  <h3>4.2 Cell Type Enrichment</h3>\n")
                
                if cell_enrichment_figures:
                    cell_enrichment = analysis_results.get('cell_enrichment')
                    if cell_enrichment is not None:
                        significant_cells = sum(cell_enrichment['p-value'] < pvalue_threshold)
                        
                        f.write("  <div class='explanation'>\n")
                        f.write("    <p>This analysis identifies cell types that are significantly over- or under-represented in this genomic region. Cell type enrichment provides insights into:</p>\n")
                        f.write("    <ul>\n")
                        f.write("      <li>Cell-specific enhancer activity that might regulate lineage-specific genes</li>\n")
                        f.write("      <li>Regulatory elements that function primarily in specific cellular contexts</li>\n")
                        f.write("      <li>Potential relationships between the enhancer region and cell type-specific functions</li>\n")
                        f.write("    </ul>\n")
                        f.write("    <p>Strong enrichment of specific cell types can help prioritize experimental validation in the most relevant cellular models.</p>\n")
                        f.write("  </div>\n")
                        
                        f.write(f"  <p><strong>Significant cell types:</strong> {significant_cells} out of {len(cell_enrichment)} cell types (p < {pvalue_threshold}).</p>\n")
                        
                        # Top significant cell types table
                        if significant_cells > 0:
                            sig_cells = cell_enrichment[cell_enrichment['p-value'] < pvalue_threshold].sort_values('p-value')
                            f.write("  <h4>Top Significant Cell Types</h4>\n")
                            f.write("  <table>\n")
                            f.write("    <tr><th>Cell Type</th><th>Fold Change</th><th>p-value</th><th>Adjusted p-value</th><th>Enrichment</th></tr>\n")
                            
                            for i, (_, row) in enumerate(sig_cells.head(5).iterrows()):
                                enrichment_type = "Enriched" if row['Fold Change'] > 1 else "Depleted"
                                enrichment_class = "enriched" if row['Fold Change'] > 1 else "depleted"
                                f.write(f"    <tr><td>{row['Cell Type']}</td><td>{row['Fold Change']:.2f}x</td>" +
                                    f"<td>{row['p-value']:.2e}</td><td>{row['adjusted p-value']:.2e}</td>" +
                                    f"<td class='{enrichment_class}'>{enrichment_type}</td></tr>\n")
                            
                            f.write("  </table>\n")
                    
                    # Display enrichment visualization
                    for figure_name, figure_path in cell_enrichment_figures:
                        f.write("  <div class='chart-container' style='max-width: 90%; margin: 20px auto;'>\n")
                        f.write(f"    <img src='{figure_path}' alt='{figure_name}'>\n")
                        f.write("  </div>\n")
                    
                    # #差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figures as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Cell Type Enrichment' in name:
                    #         format_name = name.split('(')[1].split(')')[0]
                    #         f.write(f"<a href='{path}'>{format_name}</a> ")
                    # f.write("  </div>\n")

                    # データテーブルファイルへのリンクを取得
                    cell_enrich_files = [path for name, path in generated_files['tables'] 
                                    if 'Cell Type Enrichment' in name or 'cell_enrichment' in path.lower()]

                    # 画像リンクとデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Cell Type Enrichment' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")

                    # データファイルへのリンク
                    if cell_enrich_files:
                        f.write("<strong>Data:</strong> ")
                        for path in cell_enrich_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")

                    f.write("  </div>\n")

                else:
                    f.write("  <p class='notes'>No cell type enrichment analysis results available for this region.</p>\n")
                
                # Section 5: Read Distribution Analysis
                f.write("  <h2>5. Read Distribution Analysis</h2>\n")
                f.write("  <div class='section-intro'>\n")
                f.write("    <p>This section visualizes how eRNA reads are distributed across the genomic region. Two complementary visualizations are provided:</p>\n")
                f.write("    <ol>\n")
                f.write("      <li><strong>Simple Read Distribution</strong>: Shows the raw distribution of reads from all samples.</li>\n")
                f.write("      <li><strong>Read Distribution with Peaks</strong>: Identifies areas of high eRNA activity (peaks).</li>\n")
                f.write("    </ol>\n")
                f.write("    <p>These visualizations help identify sub-regions with concentrated enhancer activity.</p>\n")
                f.write("  </div>\n")
                
                # 5.1 Simple Stacked Reads Plot
                f.write("  <h3>5.1 Simple Read Distribution</h3>\n")
                
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This plot shows eRNA reads stacked horizontally across the genomic coordinate. Each horizontal line represents a read from a sample that overlaps with this region. The visualization helps to:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li>Identify patterns in read coverage across the region</li>\n")
                f.write("      <li>Visualize the density of reads at different positions</li>\n")
                f.write("      <li>Observe the overall distribution of eRNA activity</li>\n")
                f.write("    </ul>\n")
                f.write("    <p>Areas with many stacked reads indicate regions of high enhancer RNA activity, suggesting active regulatory elements.</p>\n")
                f.write("  </div>\n")
                
                # Display simple stacked plot
                simple_plot_path = next((path for name, path in stacked_plot_figures if 'Simple Read Distribution' in name), None)
                if simple_plot_path:
                    f.write("  <div class='chart-container' style='max-width: 100%; margin: 20px auto;'>\n")
                    f.write(f"    <img src='{simple_plot_path}' alt='Simple Read Distribution'>\n")
                    f.write("  </div>\n")
                    


                    # # 差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figure as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Simple Stacked Read Plot' in name:
                    #         format_name = name.split('(')[1].split(')')[0]
                    #         f.write(f"<a href='{path}'>{format_name}</a> ")
                    # f.write("  </div>\n")

                    # 関連データファイルがあれば取得
                    simple_data_files = [path for name, path in generated_files['tables'] 
                                    if 'simple_stacked' in path.lower() or 'simple_read' in path.lower()]

                    # 画像リンクとデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Simple Stacked Read Plot' in name:
                            format_name = name.split('(')[1].split(')')[0]
                            f.write(f"<a href='{path}'>{format_name}</a> ")

                    # データファイルへのリンク
                    if simple_data_files:
                        f.write("<strong>Data:</strong> ")
                        for path in simple_data_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")

                    f.write("  </div>\n")


                else:
                    f.write("  <p class='notes'>No simple read distribution plot available for this region.</p>\n")
                
                # 5.2 Stacked Reads with Peaks
                f.write("  <h3>5.2 Read Distribution with Peaks</h3>\n")
                
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This enhanced visualization identifies peaks of eRNA activity within the region. Peaks represent sub-regions with significantly higher read density. Key features include:</p>\n")
                f.write("    <ul>\n")
                f.write("      <li>Highlighted peak regions in red background</li>\n")
                f.write("      <li>Peak identifiers with genomic positions</li>\n")
                f.write("      <li>Quantification of read density in each peak</li>\n")
                f.write("    </ul>\n")
                f.write("    <p>Peaks often represent the core of enhancer elements with high transcriptional activity. These regions are prime candidates for further experimental investigation.</p>\n")
                f.write("  </div>\n")
                
                # いずれかの条件でピークプロットを探す（複数条件で検索）
                peak_plot_path = None
                for name, path in stacked_plot_figures:
                    if 'Read Distribution with Peaks' in name or 'Stacked Read Plot with Peaks' in name or 'peak' in name.lower():
                        peak_plot_path = path
                        break
                
                # ファイル名でも検索（複数の異なる検索を試す）
                if peak_plot_path is None:
                    for name, path in generated_files['figures']:
                        if 'peak_stacked' in path.lower() or 'stacked_peak' in path.lower():
                            peak_plot_path = path
                            break
                
                # いずれのファイルも見つからなかった場合、リスト内に含まれているかもしれないファイル名を検査
                if peak_plot_path is None:
                    peak_related_files = [path for name, path in generated_files['figures'] 
                                    if 'peak' in path.lower() or 'peak' in name.lower()]
                    if peak_related_files:
                        peak_plot_path = peak_related_files[0]
                
                # ピーク関連のデータやプロットが見つかった場合に表示
                has_peak_data = 'peak_data' in analysis_results and analysis_results['peak_data']
                has_peak_plot = peak_plot_path is not None
                
                if has_peak_plot:
                    f.write("  <div class='chart-container' style='max-width: 100%; margin: 20px auto;'>\n")
                    f.write(f"    <img src='{peak_plot_path}' alt='Read Distribution with Peaks'>\n")
                    f.write("  </div>\n")
                    


                    # # 差し替え
                    # # Add download links for other formats
                    # f.write("  <div class='file-formats'>Download figure as: ")
                    # for name, path in generated_files['figures']:
                    #     if 'Stacked Read Plot with Peaks' in name or 'peak' in name.lower():
                    #         format_name = name.split('(')[1].split(')')[0] if '(' in name else path.split('.')[-1].upper()
                    #         f.write(f"<a href='{path}'>{format_name}</a> ")
                    # f.write("  </div>\n")


                    # ピークデータファイルを取得
                    peak_data_files = [path for name, path in generated_files['tables'] 
                                    if 'Peak Analysis Data' in name or 'peak_data' in path.lower()]

                    # 画像リンクとピークデータリンクを組み合わせて表示
                    f.write("  <div class='file-formats'>Download: ")

                    # 画像形式ごとのリンク
                    f.write("<strong>Figures:</strong> ")
                    for name, path in generated_files['figures']:
                        if 'Stacked Read Plot with Peaks' in name or 'peak' in name.lower():
                            format_name = name.split('(')[1].split(')')[0] if '(' in name else path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{format_name}</a> ")

                    # ピークデータファイルへのリンク
                    if peak_data_files:
                        f.write("<strong>Data:</strong> ")
                        for path in peak_data_files:
                            file_ext = path.split('.')[-1].upper()
                            f.write(f"<a href='{path}'>{file_ext}</a> ")
                            

                        # # 差し替え
                        # # ピークサンプルリストへのリンクもあれば追加
                        # peak_sample_files = [path for name, path in generated_files['tables'] 
                        #                 if 'peak' in path.lower() and 'samples' in path.lower()]
                        # if peak_sample_files:
                        #     f.write("<strong>Peak Samples:</strong> ")
                        #     f.write(f"<a href='#peak-samples'>View Sample Lists</a> ")

                        # ピークサンプルリストへのリンクもあれば追加
                        peak_sample_files = []
                        for peak_id, peak_info in peak_data.items():
                            # 各ピークIDのサンプルファイルを検索
                            files = [path for name, path in generated_files['tables'] 
                                    if peak_id in name.lower() and 'samples' in name.lower()]
                            if files:
                                peak_sample_files.append((peak_id, files[0]))

                        if peak_sample_files:
                            f.write("<strong>Peak Samples:</strong> ")
                            # 直接ファイルへのリンクを追加
                            for peak_id, file_path in peak_sample_files:
                                f.write(f"<a href='{file_path}'>{peak_id}</a> ")


                    f.write("  </div>\n")

                elif not has_peak_plot and has_peak_data:
                    # ピークデータは存在するがプロットが見つからない場合
                    f.write("  <p class='notes'>Peaks were detected, but the visualization is not available in this report.</p>\n")
                else:
                    f.write("  <p class='notes'>No peak detection plot available for this region.</p>\n")
                
                # Peak information if available
                peak_data = analysis_results.get('peak_data', None)
                if peak_data:
                    f.write("  <h4>Peak Analysis</h4>\n")
                    f.write("  <p>The following table summarizes the peaks of eRNA activity detected in this region:</p>\n")
                    
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
                f.write(f"      <li><strong>Sample Coverage:</strong> {sample_stats['overlapping_samples']:,} samples ({sample_stats['overlap_percentage']:.2f}%) from eRNAbase overlap with this region</li>\n")
                
                # Tissue findings
                tissue_enrichment = analysis_results.get('tissue_enrichment')
                if tissue_enrichment is not None:
                    sig_tissues = sum(tissue_enrichment['p-value'] < pvalue_threshold)
                    if sig_tissues > 0:
                        f.write(f"      <li><strong>Tissue Specificity:</strong> {sig_tissues} tissue types are significantly enriched or depleted (p < {pvalue_threshold})</li>\n")
                        
                        # Add top enriched tissue
                        top_enriched = tissue_enrichment[tissue_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                        if not top_enriched.empty:
                            tissue = top_enriched.iloc[0]['Tissue Type']
                            fold_change = top_enriched.iloc[0]['Fold Change']
                            p_value = top_enriched.iloc[0]['p-value']
                            f.write(f"      <li><strong>Most Enriched Tissue:</strong> <span class='enriched'>{tissue}</span> ({fold_change:.2f}x enrichment, p = {p_value:.2e})</li>\n")
                        
                        # Add top depleted tissue
                        top_depleted = tissue_enrichment[tissue_enrichment['Fold Change'] < 1/fc_threshold].sort_values('p-value')
                        if not top_depleted.empty:
                            tissue = top_depleted.iloc[0]['Tissue Type']
                            fold_change = top_depleted.iloc[0]['Fold Change']
                            p_value = top_depleted.iloc[0]['p-value']
                            f.write(f"      <li><strong>Most Depleted Tissue:</strong> <span class='depleted'>{tissue}</span> ({fold_change:.2f}x depletion, p = {p_value:.2e})</li>\n")
                
                # Cell type findings
                cell_enrichment = analysis_results.get('cell_enrichment')
                if cell_enrichment is not None:
                    sig_cells = sum(cell_enrichment['p-value'] < pvalue_threshold)
                    if sig_cells > 0:
                        f.write(f"      <li><strong>Cell Type Specificity:</strong> {sig_cells} cell types are significantly enriched or depleted (p < {pvalue_threshold})</li>\n")
                        
                        # Add top enriched cell type
                        top_cell = cell_enrichment[cell_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                        if not top_cell.empty:
                            cell = top_cell.iloc[0]['Cell Type']
                            fold_change = top_cell.iloc[0]['Fold Change']
                            p_value = top_cell.iloc[0]['p-value']
                            f.write(f"      <li><strong>Most Enriched Cell Type:</strong> <span class='enriched'>{cell}</span> ({fold_change:.2f}x enrichment, p = {p_value:.2e})</li>\n")
                
                # Peak findings
                peak_data = analysis_results.get('peak_data')
                if peak_data:
                    f.write(f"      <li><strong>eRNA Activity Peaks:</strong> {len(peak_data)} distinct peaks of enhancer RNA activity identified</li>\n")
                    
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
                if tissue_enrichment is not None and sum(tissue_enrichment['p-value'] < pvalue_threshold) > 0:
                    top_tissue = tissue_enrichment[tissue_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                    if not top_tissue.empty:
                        tissue = top_tissue.iloc[0]['Tissue Type']
                        interpretations.append(f"have enhancer activity that is particularly strong in <strong>{tissue}</strong> tissue")
                
                # Cell type specificity interpretation
                if cell_enrichment is not None and sum(cell_enrichment['p-value'] < pvalue_threshold) > 0:
                    top_cell = cell_enrichment[cell_enrichment['Fold Change'] > fc_threshold].sort_values('p-value')
                    if not top_cell.empty:
                        cell = top_cell.iloc[0]['Cell Type']
                        interpretations.append(f"show cell type-specific regulation in <strong>{cell}</strong> cells")
                
                # Peak distribution interpretation
                if peak_data and len(peak_data) > 0:
                    if len(peak_data) == 1:
                        interpretations.append("contain a single well-defined enhancer element")
                    elif len(peak_data) <= 3:
                        interpretations.append(f"contain {len(peak_data)} distinct enhancer elements that may function independently or cooperatively")
                    else:
                        interpretations.append(f"represent a complex regulatory region with {len(peak_data)} distinct areas of enhancer activity")
                
                # Default interpretation if no specific findings
                if not interpretations:
                    interpretations.append("show enhancer RNA activity across multiple tissues and cell types")
                
                # Join interpretations with proper conjunctions
                if len(interpretations) == 1:
                    f.write(f" {interpretations[0]}.")
                elif len(interpretations) == 2:
                    f.write(f" {interpretations[0]} and {interpretations[1]}.")
                else:
                    last = interpretations.pop()
                    f.write(f" {', '.join(interpretations)}, and {last}.")
                
                # # Additional biological context
                # if tissue_enrichment is not None and cell_enrichment is not None:
                #     sig_tissues = sum(tissue_enrichment['p-value'] < pvalue_threshold)
                #     sig_cells = sum(cell_enrichment['p-value'] < pvalue_threshold) 
                #     if sig_tissues > 0 and sig_cells > 0:
                #         f.write(" The tissue and cell type enrichment patterns suggest this may be an important regulatory element for lineage-specific gene expression.</p>\n")
                #     elif peak_data and len(peak_data) > 2:
                #         f.write(" The multiple peaks of activity indicate this could be a super-enhancer region with broad regulatory impact.</p>\n")
                #     else:
                #         f.write(" Further experimental validation would be valuable to confirm the regulatory function of this region.</p>\n")
                # else:
                #     f.write(" Additional experimental studies would help clarify the specific gene targets and functional impact of this enhancer region.</p>\n")
                
                f.write("  </div>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by ERNARegionAnalyzer</p>\n")
                
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
            - significant_cell_types: List of statistically significant enriched cell types
                [{'cell_type': cell_type_name, 'fold_change': fc, 'pvalue': p, 'fdr': q}, ...]
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
        
        # Extract significant cell types (positive enrichment, FDR significant)
        significant_cell_types = []
        if 'cell_enrichment' in analysis_results and analysis_results['cell_enrichment'] is not None:
            cell_df = analysis_results['cell_enrichment']
            # Check if the adjusted p-value column exists, otherwise use p-value
            fdr_col = 'adjusted p-value' if 'adjusted p-value' in cell_df.columns else 'p-value'
            
            # Filter for significant positively enriched cell types
            sig_cells = cell_df[(cell_df[fdr_col] < fdr_threshold) & 
                            (cell_df['Fold Change'] > fc_threshold)]
            
            # Sort by significance
            sig_cells = sig_cells.sort_values(fdr_col)
            
            # Convert to list of dictionaries
            for _, row in sig_cells.iterrows():
                cell_info = {
                    'cell_type': row['Cell Type'],
                    'fold_change': row['Fold Change'],
                    'pvalue': row['p-value'],
                    'fdr': row[fdr_col] if fdr_col in row else row['p-value']
                }
                significant_cell_types.append(cell_info)
        
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
            'analysis_results': analysis_results,
            'report_path': report_path,
            'has_overlaps': has_overlaps,
            'overlap_sample_count': overlap_sample_count,
            'overlap_read_count': overlap_read_count,
            'significant_tissues': significant_tissues,
            'significant_cell_types': significant_cell_types
        }


    def batch_analyze_regions_from_tsv(self, 
                                    bed_directory: str,
                                    metadata_file: str,
                                    data_output_dir: str,
                                    report_output_dir: str, 
                                    regions_tsv_file: str,
                                    max_rows: Optional[int] = None,
                                    species: str = 'human',
                                    **kwargs):
        """
        ユーザー補助メソッド：TSVファイルから複数のゲノム領域を一括解析し、結果を新たなTSVとして出力します。
        
        このメソッドの動作：
        1. generate_erna_db_report を実行して全体のメタデータ評価を行います
        2. 指定されたTSVファイルから読み込んだ各ゲノム領域に対して analyze_and_report_region を実行します
        3. 各領域の解析結果を元のTSVに追加カラムとして付加した新たなTSVファイルを生成します
        
        Parameters:
        -----------
        bed_directory : str
            BEDファイルが格納されているディレクトリパス
        metadata_file : str
            eRNAbaseメタデータファイルのパス
        data_output_dir : str
            解析データの出力先ディレクトリ
        report_output_dir : str
            レポートファイルの出力先ディレクトリ
        regions_tsv_file : str
            解析対象のゲノム領域が記載されたTSVファイルのパス
        max_rows : int, optional
            TSVファイルから読み込む最大行数。Noneの場合は全行を読み込みます
        species : str, optional
            対象種（default: 'human'）
        **kwargs : dict
            その他のパラメータ（analyze_and_report_regionに渡されます）
        
        Returns:
        --------
        str
            生成された結果TSVファイルのパス
            
        Notes:
        ------
        入力TSVファイルの形式：
        入力TSVファイルには少なくとも1列目にゲノム領域を示す情報（例：chr1_1109435_1174178）が
        含まれている必要があります。
        
        出力TSVファイルには以下のカラムが追加されます：
        - has_overlaps: 重複領域があるかどうか（True/False）
        - overlap_sample_count: 重複サンプル数
        - overlap_read_count: 重複リード数
        - significant_tissues: FDR有意なプラス方向に濃縮された組織（カンマ区切り）
        - significant_cell_types: FDR有意なプラス方向に濃縮されたセルライン（カンマ区切り）
        """

        
        # 出力ディレクトリが存在しない場合は作成
        os.makedirs(data_output_dir, exist_ok=True)
        os.makedirs(report_output_dir, exist_ok=True)
        
        # タイムスタンプを作成（ファイル名に使用）
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # BEDファイルとメタデータを読み込み
        self.logger.info(f"Loading BED files from {bed_directory}")
        self.get_bed_files(bed_directory)
        
        self.logger.info(f"Loading metadata from {metadata_file}")
        self.load_erna_metadata(metadata_file, species=species)
        
        # eRNAデータベース全体の解析を実行
        self.logger.info("Generating eRNAbase metadata report")
        db_report_path = self.generate_erna_db_report(
            output_dir=os.path.join(report_output_dir, "erna_db_report"),
            species=species
        )
        self.logger.info(f"eRNAbase metadata report generated: {db_report_path}")
        
        # TSVファイルを読み込み
        self.logger.info(f"Reading regions from TSV file: {regions_tsv_file}")
        regions_df = pd.read_csv(regions_tsv_file, sep='\t')
        
        # 行数を制限（指定がある場合）
        if max_rows is not None and max_rows > 0:
            regions_df = regions_df.head(max_rows)
            self.logger.info(f"Processing first {max_rows} rows from TSV file")
        
        # 結果を格納するための新しいカラムを追加
        regions_df['has_overlaps'] = False
        regions_df['overlap_sample_count'] = 0
        regions_df['overlap_read_count'] = 0
        regions_df['significant_tissues'] = ''
        regions_df['significant_cell_types'] = ''
        
        # 各リージョンに対して解析を実行
        for idx, row in regions_df.iterrows():
            # 最初の列をリージョン識別子として取得
            region_id = str(row.iloc[0])
            
            # リージョン文字列をパース（例: "chr1_123456_234567"）
            match = re.match(r'(chr\w+)_(\d+)_(\d+)', region_id)
            if match:
                chrom, start, end = match.groups()
                start = int(start)
                end = int(end)
                
                # リージョンの整形された名前
                region_name = f"Region_{idx+1}_{region_id}"
                region_dir = os.path.join(report_output_dir, region_name)
                
                self.logger.info(f"Analyzing region {idx+1}/{len(regions_df)}: {region_id} ({chrom}:{start}-{end})")
                
                try:
                    # リージョンの解析を実行
                    results = self.analyze_and_report_region(
                        chrom=chrom,
                        start=start,
                        end=end,
                        region_name=region_name,
                        output_dir=region_dir,
                        **kwargs
                    )
                    
                    # 結果をDataFrameに追加
                    regions_df.at[idx, 'has_overlaps'] = results['has_overlaps']
                    regions_df.at[idx, 'overlap_sample_count'] = results['overlap_sample_count']
                    regions_df.at[idx, 'overlap_read_count'] = results['overlap_read_count']
                    
                    # 有意な組織のリストをカンマ区切りの文字列に変換
                    if results['significant_tissues']:
                        tissues = [item['tissue'] for item in results['significant_tissues']]
                        regions_df.at[idx, 'significant_tissues'] = ','.join(tissues)
                    
                    # 有意なセルラインのリストをカンマ区切りの文字列に変換
                    if results['significant_cell_types']:
                        cell_types = [item['cell_type'] for item in results['significant_cell_types']]
                        regions_df.at[idx, 'significant_cell_types'] = ','.join(cell_types)
                    
                    self.logger.info(f"Analysis completed for region {region_id}")
                    self.logger.info(f"Report saved to {region_dir}")
                    
                except Exception as e:
                    self.logger.error(f"Error analyzing region {region_id}: {str(e)}")
                    import traceback
                    self.logger.error(traceback.format_exc())
            else:
                self.logger.warning(f"Row {idx+1}: Could not parse region ID '{region_id}'. Expected format: chr1_123456_234567. Skipping.")
        
        # 結果TSVファイルを生成
        output_tsv_path = os.path.join(data_output_dir, f"regions_analysis_results_{timestamp}.tsv")
        regions_df.to_csv(output_tsv_path, sep='\t', index=False)
        
        self.logger.info(f"Analysis completed for all regions")
        self.logger.info(f"Results saved to: {output_tsv_path}")
        
        return output_tsv_path


    # Additional property decorators for accessing the data

    @property
    def metadata_tissue_distribution(self):
        """
        Get tissue type distribution in the entire metadata.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing tissue type, count, and percentage.
        """
        if not hasattr(self, '_metadata_tissue_dist') or self._metadata_tissue_dist is None:
            self.logger.warning("Tissue distribution has not been calculated. Run calculate_metadata_distributions() first.")
        return self._metadata_tissue_dist


    @property
    def overlapping_tissue_counts(self):
        """
        Get tissue type counts in overlapping samples.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing tissue type and counts in overlapping samples.
        """
        if self.region_tissue_counts is None:
            self.logger.warning("Region tissue counts have not been calculated. Run analyze_overlapping_sample_distributions() first.")
        return self.region_tissue_counts

    @property
    def overlapping_cell_type_counts(self):
        """
        Get cell type counts in overlapping samples.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing cell type and counts in overlapping samples.
        """
        if self.region_cell_type_counts is None:
            self.logger.warning("Region cell type counts have not been calculated. Run analyze_overlapping_sample_distributions() first.")
        return self.region_cell_type_counts

    @property
    def tissue_enrichment_results(self):
        """
        Get tissue enrichment analysis results.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing enrichment analysis results for tissues.
        """
        if self.region_tissue_enrichment is None:
            self.logger.warning("Tissue enrichment analysis has not been performed. Run analyze_tissue_enrichment() first.")
        return self.region_tissue_enrichment

    @property
    def cell_type_enrichment_results(self):
        """
        Get cell type enrichment analysis results.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing enrichment analysis results for cell types.
        """
        if self.region_cell_type_enrichment is None:
            self.logger.warning("Cell type enrichment analysis has not been performed. Run analyze_cell_type_enrichment() first.")
        return self.region_cell_type_enrichment
    



    @property
    def metadata_cell_type_distribution(self):
        """
        Get cell type distribution in the entire metadata.
        
        Returns:
        --------
        pandas.DataFrame or None
            DataFrame containing cell type, count, and percentage.
        """
        if not hasattr(self, '_metadata_cell_type_dist') or self._metadata_cell_type_dist is None:
            self.logger.warning("Cell type distribution has not been calculated. Run calculate_metadata_distributions() first.")
        return self._metadata_cell_type_dist
    



    @property
    def erna_metadata(self):
        """
        Get the eRNAbase metadata
        
        Returns:
        --------
        pandas.DataFrame or None
            The loaded eRNAbase metadata
        """
        return self._erna_metadata

    @erna_metadata.setter
    def erna_metadata(self, value):
        """
        Set the eRNAbase metadata
        
        Parameters:
        -----------
        value : pandas.DataFrame
            The eRNAbase metadata to set
        """
        self._erna_metadata = value


    @property
    def erna_bed_files(self):
        """
        Get the list of eRNAbase BED files
        
        Returns:
        --------
        List[str] or None
            The list of eRNAbase BED file paths
        """
        return self._erna_bed_files

    @erna_bed_files.setter
    def erna_bed_files(self, value):
        """
        Set the list of eRNAbase BED files
        
        Parameters:
        -----------
        value : List[str]
            The list of eRNAbase BED file paths to set
        """
        self._erna_bed_files = value
    
    def validate_overlapping_results(self) -> bool:
        """
        Validate that overlapping results contain only expected chromosome
        
        Returns:
        --------
        bool
            True if validation passes, False otherwise
        """
        if self.overlapping_results is None:
            self.logger.warning("No overlapping results to validate")
            return False
            
        if not hasattr(self, 'interest_region') or self.interest_region is None:
            self.logger.warning("No interest region defined for validation")
            return False
        
        expected_chrom = self.interest_region['chrom']
        
        # Check if chrom column exists
        if 'chrom' not in self.overlapping_results.columns:
            self.logger.error("'chrom' column not found in overlapping results")
            return False
        
        # Get unique chromosomes
        unique_chroms = self.overlapping_results['chrom'].unique()
        
        # Check for multiple chromosomes
        if len(unique_chroms) > 1:
            self.logger.warning(
                f"Multiple chromosomes found in overlapping results: {unique_chroms}. "
                f"Expected only: {expected_chrom}"
            )
            # Filter to expected chromosome only
            self.overlapping_results = self.overlapping_results[
                self.overlapping_results['chrom'] == expected_chrom
            ]
            self.logger.info(f"Filtered results to {expected_chrom} only")
            return False
        
        # Check for chromosome mismatch
        if len(unique_chroms) == 1 and unique_chroms[0] != expected_chrom:
            self.logger.error(
                f"Chromosome mismatch: expected {expected_chrom}, got {unique_chroms[0]}"
            )
            return False
        
        # Validation passed
        self.logger.debug(f"Validation passed: all results are from {expected_chrom}")
        return True