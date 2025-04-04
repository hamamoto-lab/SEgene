import pandas as pd
import json
import os
import logging
from pathlib import Path


from sedb_common.logging import get_module_logger

# モジュールレベルのロガー
logger = get_module_logger(__name__)

class SEdbDataPackager:
    """Saves processed SEDB data as a standardized package."""
    
    def __init__(self, logger=None):
        self.logger = logger or get_module_logger(f"{__name__}.loader")
    
    def _setup_default_logger(self):
        """Setup default logger."""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        return logger
    
    def save_processed_package(self, output_dir, bed_df, sample_info, 
                              cell_id_counts=None, statistics=None, 
                              original_columns=None, partition_by_chr=False,
                              source_files=None):
        """
        Save processed data as a standardized package.
        
        Parameters:
        -----------
        output_dir : str
            Output directory path.
        bed_df : pandas.DataFrame
            Processed BED data.
        sample_info : pandas.DataFrame
            Processed sample information.
        cell_id_counts : pandas.DataFrame, optional
            Counts of cell IDs.
        statistics : dict, optional
            Calculated statistics.
        original_columns : list, optional
            List of original BED file column names.
        partition_by_chr : bool, optional
            Flag to partition by chromosome.
        source_files : dict, optional
            Dictionary containing source file paths and MD5 hashes.
        
        Returns:
        --------
        dict : 
            Information about saved files.
        """
        # ディレクトリ作成
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # サブディレクトリ作成
        data_dir = output_path / "data"
        stats_dir = output_path / "statistics"
        data_dir.mkdir(exist_ok=True)
        stats_dir.mkdir(exist_ok=True)
        
        saved_files = {}
        
        # BEDデータの保存 (染色体ごとのパーティショニングあり/なし)
        if partition_by_chr and bed_df is not None:
            self.logger.info("Saving BED data with chromosome partitioning")
            chr_dir = data_dir / "bed_by_chr"
            chr_dir.mkdir(exist_ok=True)
            
            # 染色体ごとにグループ化して保存
            for chrom, group in bed_df.groupby('se_chr'):
                chr_file = chr_dir / f"{chrom}.parquet"
                group.to_parquet(chr_file, index=False)
                saved_files[f"bed_data_{chrom}"] = str(chr_file)
                self.logger.info(f"Saved {len(group)} rows to {chr_file}")
            
            # インデックスファイルの保存（高速検索用）
            index_df = bed_df[['se_chr', 'se_start', 'se_end', 'cell_id']]
            index_file = data_dir / "bed_index.parquet"
            index_df.to_parquet(index_file, index=False)
            saved_files["bed_index"] = str(index_file)
            
        # 通常の単一ファイル保存
        if bed_df is not None:
            bed_file = data_dir / "bed_data.parquet"
            bed_df.to_parquet(bed_file, index=False)
            saved_files["bed_data"] = str(bed_file)
            self.logger.info(f"Saved {len(bed_df)} rows to {bed_file}")
        
        # サンプル情報の保存
        if sample_info is not None:
            sample_file = data_dir / "sample_info.parquet"
            sample_info.to_parquet(sample_file, index=False)
            saved_files["sample_info"] = str(sample_file)
            self.logger.info(f"Saved {len(sample_info)} sample records to {sample_file}")
        
        # cell_id出現回数の保存
        if cell_id_counts is not None:
            cell_id_file = data_dir / "cell_id_counts.parquet"
            cell_id_counts.to_parquet(cell_id_file, index=False)
            saved_files["cell_id_counts"] = str(cell_id_file)
            self.logger.info(f"Saved cell_id counts to {cell_id_file}")
        
        # 統計情報の保存
        if statistics:
            stats_file = stats_dir / "statistics.json"
            with open(stats_file, 'w') as f:
                json.dump(statistics, f, indent=2)
            saved_files["statistics"] = str(stats_file)
            self.logger.info(f"Saved statistics to {stats_file}")
            
            # 主要な統計は個別ファイルとしても保存
            for stat_name, stat_data in statistics.items():
                if isinstance(stat_data, dict) and len(stat_data) > 0:
                    stat_file = stats_dir / f"{stat_name}.json"
                    with open(stat_file, 'w') as f:
                        json.dump(stat_data, f, indent=2)
                    saved_files[f"stat_{stat_name}"] = str(stat_file)
        
        # メタデータの保存
        metadata = {
            "version": "1.0",
            "created_at": pd.Timestamp.now().isoformat(),
            "total_samples": len(sample_info) if sample_info is not None else 0,
            "total_regions": len(bed_df) if bed_df is not None else 0,
            "original_bed_columns": original_columns,
            "has_chr_partitioning": partition_by_chr,
            "files": list(saved_files.keys()),
        }
        
        # ソースファイル情報を追加（存在する場合）
        if source_files:
            metadata["source_files"] = source_files
            self.logger.info("Added source files information with MD5 hashes to metadata")
        
        metadata_file = output_path / "metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        saved_files["metadata"] = str(metadata_file)
        
        self.logger.info(f"Completed saving processed package to {output_dir}")
        return {
            "output_dir": str(output_path),
            "files": saved_files,
            "metadata": metadata
        }
