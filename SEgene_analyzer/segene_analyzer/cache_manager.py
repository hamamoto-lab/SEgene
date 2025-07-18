"""
Cache management for SEdb database statistics
"""

import pickle
import json
import os
from datetime import datetime
import hashlib
import logging
from typing import Dict, Any, Optional
import pandas as pd


class SEdbCacheManager:
    """SEdbデータベース統計のキャッシュ管理"""
    
    def __init__(self, logger=None):
        self.logger = logger or self._setup_default_logger()
        
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
    
    def calculate_file_hash(self, filepath: str) -> str:
        """ファイルのハッシュ値を計算"""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def save_statistics(self, analyzer, cache_file: str, 
                       bed_file: str, sample_info_file: str) -> Dict[str, Any]:
        """
        データベース統計をキャッシュファイルに保存
        
        Parameters:
        -----------
        analyzer : SEdbRegionAnalyzer
            統計が計算済みのアナライザーインスタンス
        cache_file : str
            保存先キャッシュファイル
        bed_file : str
            元のBEDファイルパス（検証用）
        sample_info_file : str
            元のサンプル情報ファイルパス（検証用）
            
        Returns:
        --------
        dict
            保存した統計情報のメタデータ
        """
        # 統計情報の収集
        stats = {
            'metadata': {
                'created_at': datetime.now().isoformat(),
                'bed_file': bed_file,
                'sample_info_file': sample_info_file,
                'bed_file_hash': self.calculate_file_hash(bed_file),
                'sample_info_hash': self.calculate_file_hash(sample_info_file),
                'total_samples': len(analyzer.sample_info) if analyzer.sample_info is not None else 0,
                'total_se_records': len(analyzer.bed_df) if analyzer.bed_df is not None else 0,
            },
            'statistics': {
                'tissue_distribution': analyzer._tissue_distribution if hasattr(analyzer, '_tissue_distribution') else None,
                'cell_id_counts': analyzer._cell_id_counts if hasattr(analyzer, '_cell_id_counts') else None,
                'filtered_chromosomes': hasattr(analyzer, '_original_bed_df') and analyzer._original_bed_df is not None,
            }
        }
        
        # 追加の統計情報があれば保存
        if hasattr(analyzer, '_biosample_distribution'):
            stats['statistics']['biosample_distribution'] = analyzer._biosample_distribution
        
        # ディレクトリが存在しない場合は作成
        os.makedirs(os.path.dirname(cache_file) if os.path.dirname(cache_file) else '.', exist_ok=True)
        
        # Pickleで保存
        with open(cache_file, 'wb') as f:
            pickle.dump(stats, f)
        
        # メタデータをJSONでも保存（人が読める形式）
        metadata_file = cache_file.replace('.pkl', '_metadata.json')
        with open(metadata_file, 'w') as f:
            json.dump(stats['metadata'], f, indent=2)
        
        self.logger.info(f"Statistics cached to: {cache_file}")
        self.logger.info(f"Metadata saved to: {metadata_file}")
        
        return stats['metadata']
    
    def load_statistics(self, analyzer, cache_file: str, 
                       verify_files: bool = True,
                       bed_file: Optional[str] = None,
                       sample_info_file: Optional[str] = None) -> bool:
        """
        キャッシュファイルから統計を読み込み
        
        Parameters:
        -----------
        analyzer : SEdbRegionAnalyzer
            統計を設定するアナライザーインスタンス
        cache_file : str
            キャッシュファイルパス
        verify_files : bool
            元ファイルのハッシュ値を検証するか
        bed_file : str, optional
            検証用のBEDファイルパス
        sample_info_file : str, optional
            検証用のサンプル情報ファイルパス
            
        Returns:
        --------
        bool
            読み込みが成功したかどうか
        """
        if not os.path.exists(cache_file):
            self.logger.error(f"Cache file not found: {cache_file}")
            return False
        
        try:
            # キャッシュを読み込み
            with open(cache_file, 'rb') as f:
                stats = pickle.load(f)
            
            metadata = stats['metadata']
            self.logger.info(f"Loading cached statistics from {metadata['created_at']}")
            
            # ファイルの検証（オプション）
            if verify_files and bed_file and sample_info_file:
                bed_hash = self.calculate_file_hash(bed_file)
                sample_hash = self.calculate_file_hash(sample_info_file)
                
                if bed_hash != metadata['bed_file_hash']:
                    self.logger.warning("BED file has changed since cache was created")
                    return False
                
                if sample_hash != metadata['sample_info_hash']:
                    self.logger.warning("Sample info file has changed since cache was created")
                    return False
            
            # 統計情報をアナライザーに設定
            statistics = stats['statistics']
            
            if statistics['tissue_distribution'] is not None:
                analyzer._tissue_distribution = statistics['tissue_distribution']
                self.logger.info(f"Loaded tissue distribution: {len(analyzer._tissue_distribution)} types")
            
            if statistics['cell_id_counts'] is not None:
                analyzer._cell_id_counts = statistics['cell_id_counts']
                self.logger.info(f"Loaded cell_id counts: {len(analyzer._cell_id_counts)} samples")
            
            if 'biosample_distribution' in statistics and statistics['biosample_distribution'] is not None:
                analyzer._biosample_distribution = statistics['biosample_distribution']
                self.logger.info(f"Loaded biosample distribution: {len(analyzer._biosample_distribution)} types")
            
            self.logger.info("Successfully loaded cached statistics")
            return True
            
        except Exception as e:
            self.logger.error(f"Error loading cache: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False
    
    def is_cache_valid(self, cache_file: str, max_age_days: int = 30) -> bool:
        """
        キャッシュが有効かチェック
        
        Parameters:
        -----------
        cache_file : str
            キャッシュファイルパス
        max_age_days : int
            キャッシュの最大有効日数
            
        Returns:
        --------
        bool
            キャッシュが有効かどうか
        """
        if not os.path.exists(cache_file):
            return False
        
        try:
            # メタデータファイルを確認
            metadata_file = cache_file.replace('.pkl', '_metadata.json')
            if os.path.exists(metadata_file):
                with open(metadata_file, 'r') as f:
                    metadata = json.load(f)
                
                created_at = datetime.fromisoformat(metadata['created_at'])
                age = datetime.now() - created_at
                
                if age.days > max_age_days:
                    self.logger.warning(f"Cache is {age.days} days old (max: {max_age_days})")
                    return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error checking cache validity: {e}")
            return False