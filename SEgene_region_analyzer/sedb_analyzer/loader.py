"""
sedb_analyzer.loader - Module for loading SEdb extracted data

This module provides functionality to load SEdb extracted data packages and convert them into appropriate data structures for analysis.
"""

import os
import pandas as pd
import json
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

class AnalysisDataLoader:
    """
    Loader class for SEdb extracted data packages.
    """
    
    def __init__(self, logger=None):
        """
        Initialization.
        
        Parameters
        ----------
        logger : Logger, optional
            Logger instance to use.
        """
        # ロガーの設定
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing AnalysisDataLoader")
    
    def _setup_default_logger(self):
        """Set up the default logger."""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        
        # 既存のハンドラを確認
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        
        return logger
    
    def load_extraction_package(self, package_dir: str) -> Dict[str, Any]:
        """
        Load the SEdb extraction package.
        
        Parameters
        ----------
        package_dir : str
            Directory path of the extraction package.
            
        Returns
        -------
        Dict[str, Any]
            The loaded package data:
            {
                'metadata': dict,            # Metadata
                'region': dict,              # Region information
                'overlapping_se': DataFrame, # Overlapping SE
                'sample_info': DataFrame,    # Sample information
                'summary': dict,             # Statistical summary
                'distributions': dict        # Distribution data
            }
            
        Raises
        ------
        ValueError
            If file reading fails.
        FileNotFoundError
            If a required file is not found.
        """
        package_dir = Path(package_dir)
        self.logger.info(f"Loading extraction package: {package_dir}")
        
        if not package_dir.exists() or not package_dir.is_dir():
            raise FileNotFoundError(f"Package directory not found: {package_dir}")
        
        result = {}
        
        # 1. メタデータの読み込み
        metadata_path = package_dir / "metadata.json"
        if metadata_path.exists():
            try:
                with open(metadata_path, 'r', encoding='utf-8') as f:
                    result['metadata'] = json.load(f)
                self.logger.debug("Metadata loaded")
            except Exception as e:
                self.logger.warning(f"Failed to load metadata: {e}")
                result['metadata'] = {}
        else:
            self.logger.warning("Metadata file not found")
            result['metadata'] = {}
        
        # 2. 領域情報の読み込み
        region_path = package_dir / "region_info.json"
        if not region_path.exists():
            raise FileNotFoundError(f"Region information file not found: {region_path}")
        
        try:
            with open(region_path, 'r', encoding='utf-8') as f:
                region_data = json.load(f)
            result['region'] = region_data
            chrom = region_data.get('chrom', region_data.get('chr', 'unknown'))
            start = region_data.get('start', 0)
            end = region_data.get('end', 0)
            self.logger.debug(f"Region information loaded: {chrom}:{start}-{end}")
        except Exception as e:
            raise ValueError(f"Failed to load region information: {e}")
        
        # 3. 重複SEデータの読み込み
        bed_path = package_dir / "bed_data.parquet"
        if bed_path.exists():
            try:
                result['overlapping_se'] = pd.read_parquet(bed_path)
                self.logger.debug(f"Loaded overlapping SE data: {len(result['overlapping_se'])} rows")
            except Exception as e:
                self.logger.warning(f"Failed to load overlapping SE data: {e}")
                result['overlapping_se'] = pd.DataFrame()
        else:
            self.logger.warning("Overlapping SE data file not found")
            result['overlapping_se'] = pd.DataFrame()
        
        # 4. サンプル情報の読み込み
        sample_path = package_dir / "sample_info.parquet"
        if sample_path.exists():
            try:
                result['sample_info'] = pd.read_parquet(sample_path)
                self.logger.debug(f"Loaded sample information: {len(result['sample_info'])} rows")
            except Exception as e:
                self.logger.warning(f"Failed to load sample information: {e}")
                result['sample_info'] = pd.DataFrame()
        else:
            self.logger.warning("Sample information file not found")
            result['sample_info'] = pd.DataFrame()
        
        # 5. 統計サマリーの読み込み
        summary_path = package_dir / "summary.json"
        if summary_path.exists():
            try:
                with open(summary_path, 'r', encoding='utf-8') as f:
                    result['summary'] = json.load(f)
                self.logger.debug("Loaded statistical summary")
            except Exception as e:
                self.logger.warning(f"Failed to load statistical summary: {e}")
                result['summary'] = {}
        else:
            self.logger.warning("Statistical summary file not found")
            result['summary'] = {}
        
        # 6. 分布データの読み込み
        result['distributions'] = self._load_distributions(package_dir)
        
        self.logger.info(f"Extraction package loading complete: {package_dir}")
        return result
    
    # def _load_distributions(self, package_dir: Path) -> Dict[str, pd.DataFrame]:
    #     """
    #     分布データディレクトリからデータを読み込む内部メソッド
    #        
    #     Parameters
    #     ----------
    #     package_dir : Path
    #         パッケージディレクトリ
    #            
    #     Returns
    #     -------
    #     Dict[str, pd.DataFrame]
    #         分布データの辞書
    #     """
    #     distributions = {}
    #     dist_dir = package_dir / "distributions"
    #        
    #     if not dist_dir.exists() or not dist_dir.is_dir():
    #         self.logger.warning("分布データディレクトリが見つかりません")
    #         return distributions
    #        
    #     # 各分布ファイルを読み込み
    #     for dist_type in ['tissue_counts', 'biosample_counts', 'gene_counts']:
    #         dist_path = dist_dir / f"{dist_type}.parquet"
    #         if dist_path.exists():
    #             try:
    #                 distributions[dist_type] = pd.read_parquet(dist_path)
    #                 self.logger.debug(f"{dist_type}を読み込みました")
    #             except Exception as e:
    #                 self.logger.error(f"{dist_type}の読み込みに失敗しました: {e}")
    #        
    #     return distributions

    def _load_distributions(self, package_dir: Path) -> Dict[str, Union[pd.DataFrame, pd.Series]]:
        """
        Internal method to load data from the distribution data directory.
        
        Parameters
        ----------
        package_dir : Path
            Package directory.
            
        Returns
        -------
        Dict[str, Union[pd.DataFrame, pd.Series]]
            Dictionary of distribution data (including those converted to Series).
        """
        distributions = {}
        dist_dir = package_dir / "distributions"
        
        if not dist_dir.exists() or not dist_dir.is_dir():
            self.logger.warning("Distribution data directory not found")
            return distributions
        
        # Load each distribution file
        for dist_type in ['tissue_counts', 'biosample_counts', 'gene_counts']:
            dist_path = dist_dir / f"{dist_type}.parquet"
            if dist_path.exists():
                try:
                    # First, read as DataFrame
                    df = pd.read_parquet(dist_path)
                    self.logger.info(f"Loaded {dist_type} as DataFrame: {len(df)} rows, columns={list(df.columns)}")
                    
                    # Processing for tissue distribution
                    if dist_type == 'tissue_counts':
                        if 'tissue_type' in df.columns and 'count' in df.columns:
                            # Standard format [tissue_type, count]
                            try:
                                series = df.set_index('tissue_type')['count']
                                self.logger.info(f"Converted {dist_type} from [tissue_type, count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        elif 'Tissue Type' in df.columns and 'Count' in df.columns:
                            # Alternative format [Tissue Type, Count]
                            try:
                                series = df.set_index('Tissue Type')['Count']
                                self.logger.info(f"Converted {dist_type} from [Tissue Type, Count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        elif 'term' in df.columns and 'count' in df.columns:
                            # Alternative format [term, count]
                            try:
                                series = df.set_index('term')['count']
                                self.logger.info(f"Converted {dist_type} from [term, count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        else:
                            # Generic handling for differing column names
                            self.logger.warning(f"Unrecognized column pattern for {dist_type}: {list(df.columns)}")
                            
                            # Attempt to guess count and type columns
                            count_cols = [col for col in df.columns if 'count' in col.lower() or 'freq' in col.lower()]
                            type_cols = [col for col in df.columns if 'type' in col.lower() or 'tissue' in col.lower() or 'term' in col.lower()]
                            
                            if len(count_cols) == 1 and len(type_cols) == 1:
                                count_col = count_cols[0]
                                type_col = type_cols[0]
                                try:
                                    series = df.set_index(type_col)[count_col]
                                    self.logger.info(f"Converted {dist_type} from [{type_col}, {count_col}] format to Series: {len(series)} items")
                                    distributions[dist_type] = series
                                except Exception as e:
                                    self.logger.error(f"Failed to convert using guessed column names: {e}")
                                    distributions[dist_type] = df
                            else:
                                self.logger.warning("Could not infer appropriate column names; using original DataFrame")
                                distributions[dist_type] = df
                    
                    # Processing for biosample distribution
                    elif dist_type == 'biosample_counts':
                        if 'biosample_type' in df.columns and 'count' in df.columns:
                            # Standard format [biosample_type, count]
                            try:
                                series = df.set_index('biosample_type')['count']
                                self.logger.info(f"Converted {dist_type} from [biosample_type, count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        elif 'Biosample Type' in df.columns and 'Count' in df.columns:
                            # Alternative format [Biosample Type, Count]
                            try:
                                series = df.set_index('Biosample Type')['Count']
                                self.logger.info(f"Converted {dist_type} from [Biosample Type, Count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        else:
                            # Generic handling for differing column names
                            self.logger.warning(f"Unrecognized column pattern for {dist_type}: {list(df.columns)}")
                            
                            # Attempt to guess count and biosample columns
                            count_cols = [col for col in df.columns if 'count' in col.lower() or 'freq' in col.lower()]
                            type_cols = [col for col in df.columns if 'biosample' in col.lower() or 'sample' in col.lower() or 'type' in col.lower()]
                            
                            if len(count_cols) == 1 and len(type_cols) == 1:
                                count_col = count_cols[0]
                                type_col = type_cols[0]
                                try:
                                    series = df.set_index(type_col)[count_col]
                                    self.logger.info(f"Converted {dist_type} from [{type_col}, {count_col}] format to Series: {len(series)} items")
                                    distributions[dist_type] = series
                                except Exception as e:
                                    self.logger.error(f"Failed to convert using guessed column names: {e}")
                                    distributions[dist_type] = df
                            else:
                                self.logger.warning("Could not infer appropriate column names; using original DataFrame")
                                distributions[dist_type] = df
                    
                    # Processing for gene counts distribution
                    elif dist_type == 'gene_counts':
                        if 'gene' in df.columns and 'count' in df.columns:
                            # Standard format [gene, count]
                            try:
                                series = df.set_index('gene')['count']
                                self.logger.info(f"Converted {dist_type} from [gene, count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        elif 'Gene' in df.columns and 'Count' in df.columns:
                            # Alternative format [Gene, Count]
                            try:
                                series = df.set_index('Gene')['Count']
                                self.logger.info(f"Converted {dist_type} from [Gene, Count] format to Series: {len(series)} items")
                                distributions[dist_type] = series
                            except Exception as e:
                                self.logger.error(f"Failed to convert {dist_type} to Series: {e}")
                                distributions[dist_type] = df
                        
                        else:
                            # Generic handling for differing column names
                            self.logger.warning(f"Unrecognized column pattern for {dist_type}: {list(df.columns)}")
                            
                            # Attempt to guess count and gene columns
                            count_cols = [col for col in df.columns if 'count' in col.lower() or 'freq' in col.lower()]
                            gene_cols = [col for col in df.columns if 'gene' in col.lower() or 'symbol' in col.lower()]
                            
                            if len(count_cols) == 1 and len(gene_cols) == 1:
                                count_col = count_cols[0]
                                gene_col = gene_cols[0]
                                try:
                                    series = df.set_index(gene_col)[count_col]
                                    self.logger.info(f"Converted {dist_type} from [{gene_col}, {count_col}] format to Series: {len(series)} items")
                                    distributions[dist_type] = series
                                except Exception as e:
                                    self.logger.error(f"Failed to convert using guessed column names: {e}")
                                    distributions[dist_type] = df
                            else:
                                self.logger.warning("Could not infer appropriate column names; using original DataFrame")
                                distributions[dist_type] = df
                    
                    # For other data types, keep as DataFrame
                    else:
                        self.logger.info(f"{dist_type} will be retained as DataFrame")
                        distributions[dist_type] = df
                    
                    # Log the final data type (for debugging)
                    data_type = type(distributions[dist_type]).__name__
                    if isinstance(distributions[dist_type], pd.Series):
                        self.logger.info(f"Final data type for {dist_type}: Series (dtype={distributions[dist_type].dtype}, items={len(distributions[dist_type])})")
                    else:
                        self.logger.info(f"Final data type for {dist_type}: {data_type} (rows={len(distributions[dist_type])})")
                    
                except Exception as e:
                    self.logger.error(f"Failed to load and convert {dist_type}: {e}")
        
        return distributions



    def load_standard_bed_file(self, file_path: str) -> Tuple[pd.DataFrame, List[str]]:
        """
        Load a standard BED file.
        
        Parameters
        ----------
        file_path : str
            Path to the BED file to be loaded.
            
        Returns
        -------
        Tuple[pd.DataFrame, List[str]]
            (DataFrame, original list of column names)
        """
        self.logger.info(f"Loading BED file: {file_path}")
        
        try:
            # ヘッダーの有無を自動検出
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                has_header = first_line.startswith('#') or not any(c.isdigit() for c in first_line.split('\t')[0])
            
            # BEDファイルの読み込み
            if has_header:
                df = pd.read_csv(file_path, sep='\t')
                original_columns = list(df.columns)
            else:
                # デフォルト列名を設定
                default_names = ['chr', 'start', 'end', 'name', 'score', 'strand', 
                                'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
                                'blockSizes', 'blockStarts']
                df = pd.read_csv(file_path, sep='\t', header=None, names=default_names)
                original_columns = default_names
            
            return df, original_columns
            
        except Exception as e:
            raise ValueError(f"Failed to load BED file: {file_path}, error: {e}")

    def load_sample_info_file(self, file_path: str) -> pd.DataFrame:
        """
        Load a sample information file.
        
        Parameters
        ----------
        file_path : str
            Path to the sample information file.
            
        Returns
        -------
        pd.DataFrame
            Loaded sample information.
        """
        self.logger.info(f"Loading sample information file: {file_path}")
        
        try:
            # ファイル形式によって読み込み方法を変更
            if file_path.endswith('.parquet'):
                return pd.read_parquet(file_path)
            
            # タブ区切りファイルの読み込み（複数エンコーディング対応）
            encodings = ['utf-8', 'utf-8-sig', 'utf-16', 'utf-16-le']
            
            for encoding in encodings:
                try:
                    df = pd.read_csv(file_path, sep='\t', encoding=encoding)
                    self.logger.debug(f"Successfully loaded sample information with {encoding} encoding")
                    return df
                except Exception as e:
                    self.logger.debug(f"Failed to load with {encoding} encoding: {str(e)}")
                    continue
            
            # すべてのエンコーディングが失敗した場合
            raise ValueError(f"Could not load the sample information file with any encoding: {file_path}")
            
        except Exception as e:
            raise ValueError(f"Failed to load sample information file: {file_path}, error: {e}")

    def load_preprocessed_package(self, package_dir: str) -> Dict[str, Any]:
        """
        Load the SEdb preprocessed package (for obtaining overall distribution data).
        
        Parameters
        ----------
        package_dir : str
            Directory path of the preprocessed package.
            
        Returns
        -------
        Dict[str, Any]
            The loaded overall data:
            {
                'metadata': dict,               # Metadata
                'statistics': dict,             # Statistical information
                'distributions': dict,          # Distribution data
                'tissue_distribution': Series,  # Tissue type distribution
                'biosample_distribution': Series # Biosample type distribution
            }
        """
        package_dir = Path(package_dir)
        self.logger.info(f"Loading preprocessed package: {package_dir}")
        
        if not package_dir.exists() or not package_dir.is_dir():
            raise FileNotFoundError(f"Package directory not found: {package_dir}")
        
        result = {}
        
        # 1. メタデータの読み込み
        metadata_path = package_dir / "metadata.json"
        if metadata_path.exists():
            try:
                with open(metadata_path, 'r', encoding='utf-8') as f:
                    result['metadata'] = json.load(f)
                self.logger.debug("Metadata loaded")
            except Exception as e:
                self.logger.warning(f"Failed to load metadata: {e}")
                result['metadata'] = {}
        else:
            self.logger.warning("Metadata file not found")
            result['metadata'] = {}
        
        # 2. 統計情報の読み込み
        stats_dir = package_dir / "statistics"
        if stats_dir.exists() and stats_dir.is_dir():
            # 統計ファイルの読み込み
            stats_file = stats_dir / "statistics.json"
            if stats_file.exists():
                try:
                    with open(stats_file, 'r', encoding='utf-8') as f:
                        result['statistics'] = json.load(f)
                    self.logger.debug("Statistical information loaded")
                except Exception as e:
                    self.logger.warning(f"Failed to load statistical information: {e}")
                    result['statistics'] = {}
            
            # 組織分布の読み込み
            tissue_dist_file = stats_dir / "tissue_distribution.json"
            if tissue_dist_file.exists():
                try:
                    with open(tissue_dist_file, 'r', encoding='utf-8') as f:
                        tissue_data = json.load(f)
                    
                    # Convert JSON dictionary to Series
                    tissue_counts = {}
                    for tissue, data in tissue_data.items():
                        if isinstance(data, dict) and 'count' in data:
                            tissue_counts[tissue] = data['count']
                    
                    result['tissue_distribution'] = pd.Series(tissue_counts)
                    self.logger.debug(f"Loaded tissue distribution: {len(tissue_counts)} types")
                except Exception as e:
                    self.logger.warning(f"Failed to load tissue distribution: {e}")
        else:
            self.logger.warning("Statistics directory not found")
        
        # 3. 分布データファイルの読み込み
        dist_dir = package_dir / "data"
        result['distributions'] = {}
        
        if dist_dir.exists() and dist_dir.is_dir():
            # tissue_distribution.parquetを探す
            tissue_file = dist_dir / "tissue_distribution.parquet"
            if tissue_file.exists():
                try:
                    df = pd.read_parquet(tissue_file)
                    # 形式に応じた変換を試みる
                    if 'tissue_type' in df.columns and 'count' in df.columns:
                        tissue_series = pd.Series(df['count'].values, index=df['tissue_type'])
                        result['tissue_distribution'] = tissue_series
                    elif 'Tissue Type' in df.columns and 'Sample Count' in df.columns:
                        tissue_series = pd.Series(df['Sample Count'].values, index=df['Tissue Type'])
                        result['tissue_distribution'] = tissue_series
                    else:
                        # そのままDataFrameとして保存
                        result['distributions']['tissue'] = df
                    
                    self.logger.debug("Loaded tissue distribution data")
                except Exception as e:
                    self.logger.warning(f"Failed to load tissue distribution data: {e}")
            
            # biosample_distribution.parquetを探す
            biosample_file = dist_dir / "biosample_distribution.parquet"
            if biosample_file.exists():
                try:
                    df = pd.read_parquet(biosample_file)
                    # 形式に応じた変換を試みる
                    if 'biosample_type' in df.columns and 'count' in df.columns:
                        biosample_series = pd.Series(df['count'].values, index=df['biosample_type'])
                        result['biosample_distribution'] = biosample_series
                    elif 'Biosample Type' in df.columns and 'Sample Count' in df.columns:
                        biosample_series = pd.Series(df['Sample Count'].values, index=df['Biosample Type'])
                        result['biosample_distribution'] = biosample_series
                    else:
                        # そのままDataFrameとして保存
                        result['distributions']['biosample'] = df
                    
                    self.logger.debug("Loaded biosample distribution data")
                except Exception as e:
                    self.logger.warning(f"Failed to load biosample distribution data: {e}")
        
        # 4. 代替手段: cell_id_counts.parquetから間接的に分布を構築
        cell_id_file = dist_dir / "cell_id_counts.parquet"
        if cell_id_file.exists() and 'tissue_distribution' not in result:
            try:
                # Attempt to construct distribution data from cell_id counts
                self.logger.debug("Attempting to construct distribution data from cell_id counts")
                
                # Read cell_id counts
                cell_counts = pd.read_parquet(cell_id_file)
                
                # If sample information is available, merge and calculate distributions
                sample_file = dist_dir / "sample_info.parquet"
                if sample_file.exists():
                    sample_info = pd.read_parquet(sample_file)
                    if 'cell_id' in sample_info.columns:
                        # Identify tissue type column
                        tissue_col = None
                        for col in ['Tissue type', 'Tissue']:
                            if col in sample_info.columns:
                                tissue_col = col
                                break
                        
                        # Identify biosample type column
                        biosample_col = None
                        for col in ['Biosample type', 'Biosample']:
                            if col in sample_info.columns:
                                biosample_col = col
                                break
                        
                        # Merge on cell_id and compute distributions
                        if tissue_col:
                            merged = pd.merge(cell_counts, sample_info[['cell_id', tissue_col]], on='cell_id', how='left')
                            tissue_counts = merged.groupby(tissue_col)['count'].sum()
                            result['tissue_distribution'] = tissue_counts
                            self.logger.debug(f"Constructed tissue distribution from cell_id counts: {len(tissue_counts)} types")
                        
                        if biosample_col:
                            merged = pd.merge(cell_counts, sample_info[['cell_id', biosample_col]], on='cell_id', how='left')
                            biosample_counts = merged.groupby(biosample_col)['count'].sum()
                            result['biosample_distribution'] = biosample_counts
                            self.logger.debug(f"Constructed biosample distribution from cell_id counts: {len(biosample_counts)} types")
            except Exception as e:
                self.logger.warning(f"Failed to construct distribution from cell_id counts: {e}")
        
        self.logger.info(f"Preprocessed package loading complete: {package_dir}")
        return result
