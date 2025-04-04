# sedb_preprocessor/processor.py
import pandas as pd
import numpy as np
import logging
import warnings
from collections import Counter

from sedb_common.logging import get_module_logger

# モジュールレベルのロガー
logger = get_module_logger(__name__)

class SEdbDataProcessor:
    """SEdb Data Preprocessing Class"""
    
    def __init__(self, logger=None):
        self.logger = logger or get_module_logger(f"{__name__}.loader")
        
    def _setup_default_logger(self):
        """Set up the default logger"""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        return logger
    
    def remove_na_sample_ids(self, sample_info):
        """Remove rows with NaN in 'Sample ID' from sample information"""
        if sample_info is None:
            self.logger.warning("Sample information is not loaded. Nothing to remove.")
            return None

        original_count = len(sample_info)
        filtered_info = sample_info.dropna(subset=['Sample ID']).copy()
        filtered_count = len(filtered_info)
        removed_count = original_count - filtered_count

        if removed_count > 0:
            self.logger.info(f"Removed {removed_count} rows with NaN 'Sample ID' from sample_info")
            
        return filtered_info
    
    def add_cell_id_to_sample_info(self, sample_info):
        """Add a 'cell_id' column to sample information derived from 'Sample ID'"""
        if sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot add cell_id.")
            return None

        try:
            sample_info_copy = sample_info.copy()
            sample_info_copy['cell_id'] = sample_info_copy['Sample ID'].str.replace('Sample_', 'SE_', regex=False)
            
            # 'cell_id'カラムを最初に配置
            cols = ['cell_id'] + [col for col in sample_info_copy.columns if col != 'cell_id']
            sample_info_copy = sample_info_copy[cols]
            
            self.logger.info("'cell_id' column added to sample_info and moved to the first position")
            return sample_info_copy
        except Exception as e:
            self.logger.error(f"Error creating 'cell_id' column: {e}")
            return sample_info
    
    def filter_human_chromosomes(self, bed_df, chromosomes=None):
        """Filter BED data to include only specified chromosomes (default: chr1-22, chrX, chrY)"""
        if bed_df is None:
            self.logger.warning("BED data is not loaded. Nothing to filter.")
            return None, None, None

        # 指定されていない場合は標準のヒト染色体をデフォルトとする
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

        self.logger.info(f"Filtering BED data for {len(chromosomes)} chromosomes")

        # 元のデータを保存
        original_df = bed_df.copy()

        # データをフィルタリング
        filtered_df = bed_df.loc[bed_df['se_chr'].isin(chromosomes)].copy()
        removed_df = bed_df.loc[~bed_df['se_chr'].isin(chromosomes)].copy()

        original_size = len(original_df)
        filtered_size = len(filtered_df)
        removed_size = len(removed_df)

        self.logger.info(f"Filtered {removed_size} rows, keeping {filtered_size} rows out of {original_size}")

        return filtered_df, original_df, removed_df
    
    def count_cell_id_occurrences(self, bed_df, sample_info):
        """Count occurrences of 'cell_id' in bed_df"""
        if sample_info is None or bed_df is None:
            self.logger.warning("Sample information or BED data is not loaded.")
            return None

        # bed_dfにおける'cell_id'の出現回数を取得
        bed_counts = bed_df['cell_id'].value_counts()

        result = []
        # sample_infoのcell_idに基づいてカウント
        for cell_id in sample_info['cell_id']:
            count = bed_counts.get(cell_id, 0)
            result.append({'cell_id': cell_id, 'count': count})
            if count == 0:
                self.logger.warning(f"cell_id '{cell_id}' not found in BED data")

        # 結果をDataFrameに変換
        cell_id_counts = pd.DataFrame(result)
        self.logger.info("Cell ID occurrence counting completed")
        return cell_id_counts
    
    def calculate_basic_statistics(self, bed_df, sample_info):
        """Calculate basic statistics"""
        stats = {}
        
        # サンプル統計
        if sample_info is not None:
            # 組織タイプの分布
            if 'Tissue type' in sample_info.columns:
                tissue_counts = sample_info['Tissue type'].value_counts()
                tissue_percent = (tissue_counts / len(sample_info) * 100).round(2)
                tissue_stats = pd.DataFrame({
                    'count': tissue_counts,
                    'percent': tissue_percent
                })
                stats['tissue_distribution'] = tissue_stats.to_dict(orient='index')
            
            # バイオサンプルタイプの分布
            if 'Biosample type' in sample_info.columns:
                biosample_counts = sample_info['Biosample type'].value_counts()
                biosample_percent = (biosample_counts / len(sample_info) * 100).round(2)
                biosample_stats = pd.DataFrame({
                    'count': biosample_counts,
                    'percent': biosample_percent
                })
                stats['biosample_distribution'] = biosample_stats.to_dict(orient='index')
        
        # BEDデータ統計
        if bed_df is not None:
            # 染色体分布
            chrom_counts = bed_df['se_chr'].value_counts()
            chrom_stats = pd.DataFrame({
                'count': chrom_counts,
                'percent': (chrom_counts / len(bed_df) * 100).round(2)
            })
            stats['chromosome_distribution'] = chrom_stats.to_dict(orient='index')
            
            # SE長さの統計
            bed_df['se_length'] = bed_df['se_end'] - bed_df['se_start']
            length_stats = {
                'min': int(bed_df['se_length'].min()),
                'max': int(bed_df['se_length'].max()),
                'mean': float(bed_df['se_length'].mean()),
                'median': float(bed_df['se_length'].median()),
                'std': float(bed_df['se_length'].std()),
                'total': int(bed_df['se_length'].sum())
            }
            stats['se_length_statistics'] = length_stats
            
            # サンプルあたりのSE数統計
            sample_se_counts = bed_df['cell_id'].value_counts()
            se_count_stats = {
                'min': int(sample_se_counts.min()),
                'max': int(sample_se_counts.max()),
                'mean': float(sample_se_counts.mean()),
                'median': float(sample_se_counts.median()),
                'std': float(sample_se_counts.std())
            }
            stats['se_per_sample_statistics'] = se_count_stats
        
        # 総合統計
        if bed_df is not None and sample_info is not None:
            stats['total_samples'] = len(sample_info)
            stats['total_se_regions'] = len(bed_df)
            stats['unique_cell_ids'] = bed_df['cell_id'].nunique()
        
        self.logger.info("Basic statistics calculation completed")
        return stats
