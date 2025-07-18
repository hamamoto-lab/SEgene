# sedb_tools/data_loader/processors.py
import pandas as pd
import numpy as np
import logging
import warnings
from collections import Counter

class SEdbDataProcessor:
    """Super-enhancer (SE) データのプロセッサークラス"""
    
    def __init__(self, logger=None):
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbDataProcessor")
    
    def _setup_default_logger(self):
        """デフォルトロガーのセットアップ"""
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
        """
        サンプル情報から'Sample ID'がNaNの行を削除
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            処理するサンプル情報DataFrame
            
        Returns:
        --------
        pandas.DataFrame
            処理後のDataFrame
        """
        if sample_info is None:
            self.logger.warning("Sample information is not loaded. Nothing to remove.")
            return None

        original_count = len(sample_info)
        # .copy()を使用して元のDataFrameを変更しないようにする
        filtered_info = sample_info.dropna(subset=['Sample ID']).copy()
        filtered_count = len(filtered_info)
        removed_count = original_count - filtered_count

        if removed_count > 0:
            self.logger.info(f"Removed rows with NaN 'Sample ID' from sample_info:")
            self.logger.info(f"  Original row count: {original_count}")
            self.logger.info(f"  Removed row count: {removed_count}")
            self.logger.info(f"  Filtered row count: {filtered_count}")
        else:
            self.logger.info("No rows with NaN 'Sample ID' found in sample_info.")
            
        return filtered_info

    def add_cell_id_to_sample_info(self, sample_info):
        """
        サンプル情報に'cell_id'カラムを追加し、'Sample ID'から派生させる
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            処理するサンプル情報DataFrame
            
        Returns:
        --------
        pandas.DataFrame
            処理後のDataFrame
        """
        if sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot add cell_id.")
            return None

        try:
            # 出力タイプのチェック
            self.logger.debug(f"Type of sample_info['Sample ID']: {type(sample_info['Sample ID'])}")
            
            # 単純な文字列置換を使用（より安全な方法）
            sample_info_copy = sample_info.copy()
            sample_info_copy['cell_id'] = sample_info_copy['Sample ID'].str.replace('Sample_', 'SE_', regex=False)
            
            # 'cell_id'カラムが正しく作成されたか確認
            if 'cell_id' not in sample_info_copy.columns:
                raise ValueError("Failed to create 'cell_id' column.")

        except (AttributeError, ValueError) as e:
            self.logger.error(f"Error creating 'cell_id' column: {e}.  Check 'Sample ID' format.")
            return sample_info

        try:
            # カラムを並べ替えて'cell_id'を最初に配置
            cols = ['cell_id'] + [col for col in sample_info_copy.columns if col != 'cell_id']
            sample_info_copy = sample_info_copy[cols]
        except KeyError as e:
            self.logger.error(f"Error reordering columns: {e}.  'cell_id' might not be present.")
            return sample_info
        except Exception as e:
            self.logger.error(f"An unexpected error occurred during column reordering")
            return sample_info

        self.logger.info("'cell_id' column added to sample_info and moved to the first position.")
        return sample_info_copy

    def filter_human_chromosomes(self, bed_df, chromosomes=None):
        """
        BEDデータを特定の染色体のみに絞り込む (デフォルト: chr1-22, chrX, chrY)
        
        Parameters:
        -----------
        bed_df : pandas.DataFrame
            フィルタリングするBEDデータ
        chromosomes : list of str, optional
            保持する染色体のリスト
            
        Returns:
        --------
        tuple
            (filtered_df, original_df, removed_df) - フィルタリング後、元のデータ、除外されたデータ
        """
        if bed_df is None:
            self.logger.warning("BED data is not loaded. Nothing to filter.")
            return None, None, None

        # 指定されていない場合は標準のヒト染色体をデフォルトとする
        if chromosomes is None:
            chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']

        self.logger.info(f"Filtering BED data for chromosomes: {', '.join(chromosomes)}")

        # 元のデータを保存
        original_df = bed_df.copy()

        # データをフィルタリング
        filtered_df = bed_df.loc[bed_df['se_chr'].isin(chromosomes)].copy()
        removed_df = bed_df.loc[~bed_df['se_chr'].isin(chromosomes)].copy()

        original_size = len(original_df)
        filtered_size = len(filtered_df)
        removed_size = len(removed_df)

        self.logger.info(f"BED data filtering completed:")
        self.logger.info(f"  Original size: {original_size} rows")
        self.logger.info(f"  Removed size: {removed_size} rows")
        self.logger.info(f"  Filtered size: {filtered_size} rows")

        return filtered_df, original_df, removed_df

    def count_cell_id_occurrences(self, bed_df, sample_info):
        """
        bed_dfにおける'cell_id'の出現回数をカウント
        
        Parameters:
        -----------
        bed_df : pandas.DataFrame
            BEDデータ
        sample_info : pandas.DataFrame
            サンプル情報
            
        Returns:
        --------
        pandas.DataFrame
            'cell_id'とその出現'count'を含むDataFrame
        """
        if sample_info is None:
            self.logger.warning("Sample information is not loaded. Cannot count cell_id occurrences.")
            return None
        if bed_df is None:
            self.logger.warning("BED data is not loaded. Cannot count cell_id occurrences.")
            return None

        # bed_dfにおける'cell_id'の出現回数を取得
        bed_counts = bed_df['cell_id'].value_counts()

        result = []
        # sample_infoのcell_idに基づいてカウント
        for cell_id in sample_info['cell_id']:
            count = bed_counts.get(cell_id, 0)
            result.append({'cell_id': cell_id, 'count': count})
            if count == 0:
                message = f"Warning: cell_id '{cell_id}' not found in BED data."
                self.logger.warning(message)
                warnings.warn(message)

        # 結果をDataFrameに変換
        cell_id_counts = pd.DataFrame(result)
        self.logger.info("cell_id occurrence counting completed.")
        return cell_id_counts