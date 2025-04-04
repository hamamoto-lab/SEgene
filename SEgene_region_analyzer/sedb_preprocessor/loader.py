# sedb_tools_package/data_loader/loaders.py

import pandas as pd
import logging
import hashlib
from pathlib import Path

# 修正後:
from sedb_common.logging import get_module_logger

# モジュールレベルのロガー
logger = get_module_logger(__name__)

class SEdbDataLoader:
    """SEdb BED file and metadata loader"""
    
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
    
    def calculate_file_md5(self, file_path):
        """Calculate the MD5 hash of the file"""
        self.logger.info(f"Calculating MD5 hash for: {file_path}")
        hash_md5 = hashlib.md5()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        md5_hash = hash_md5.hexdigest()
        self.logger.debug(f"MD5 hash: {md5_hash}")
        return md5_hash
    
    def load_bed_file(self, file_path):
        """Load the BED file and return a standardized DataFrame"""
        self.logger.info(f"Reading BED file: {file_path}")
        
        # MD5ハッシュを計算
        md5_hash = self.calculate_file_md5(file_path)
        
        try:
            # Read file with header, preserving original column names
            se_df = pd.read_csv(file_path, sep='\s+', header=0)
            
            # Store original column names for reference
            original_columns = list(se_df.columns)
            self.logger.debug(f"BED file original columns: {', '.join(original_columns)}")
            
            # Log basic information about the data
            self.logger.info(f"BED file loaded with {len(se_df)} rows and {len(se_df.columns)} columns")
            
            # Validate essential columns
            required_cols = ['se_chr', 'se_start', 'se_end', 'cell_id']
            missing_cols = [col for col in required_cols if col not in se_df.columns]
            if missing_cols:
                self.logger.warning(f"Missing essential columns: {', '.join(missing_cols)}")
            
            return se_df, original_columns, md5_hash  # MD5ハッシュも返す
        except Exception as e:
            self.logger.exception(f"Error reading BED file: {e}")
            raise
    
    def load_sample_info(self, file_path):
        """Load the sample information file (supports multiple encodings)"""
        self.logger.info(f"Reading sample information file: {file_path}")
        
        # MD5ハッシュを計算
        md5_hash = self.calculate_file_md5(file_path)
        
        # 複数のエンコーディングを試す
        encodings = ['utf-16', 'utf-16-le', 'utf-8', 'utf-8-sig']
        
        for encoding in encodings:
            try:
                df = pd.read_csv(file_path, sep='\t', encoding=encoding)
                self.logger.info(f"Successfully read sample info with {encoding} encoding")
                
                # 必須カラムの検証
                required_cols = ['Sample ID', 'Species']
                missing_cols = [col for col in required_cols if col not in df.columns]
                if missing_cols:
                    self.logger.warning(f"Missing expected columns: {', '.join(missing_cols)}")
                
                return df, md5_hash  # MD5ハッシュも返す
            except Exception as e:
                self.logger.debug(f"{encoding} encoding failed: {str(e)}")
                continue
        
        # すべてのエンコーディングが失敗した場合
        self.logger.error("Failed to read sample information with any encoding")
        raise ValueError(f"Could not read file {file_path} with any supported encoding")
