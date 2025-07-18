# sedb_tools/data_loader/loaders.py
import pandas as pd
import logging
import time


class SEdbDataLoader:
    """Super-enhancer (SE) データローダー"""
    
    def __init__(self, logger=None):
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbDataLoader")
    
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


    # def read_bed_file(self, file_path):
    #     """
    #     Read BED file with original column names - internal method

    #     Parameters:
    #     file_path (str): Path to the BED file

    #     Returns:
    #     pandas.DataFrame: Loaded BED data with original column names
    #     """
    #     self.logger.info(f"Reading BED file: {file_path}")
    #     try:
    #         # Read file with header, preserving original column names
    #         se_df = pd.read_csv(file_path, sep='\s+', header=0)

    #         # Store original column names for reference
    #         self._original_bed_columns = list(se_df.columns)
    #         self.logger.debug(f"BED file original columns: {', '.join(self._original_bed_columns)}")

    #         # Log basic information about the data
    #         self.logger.info(f"BED file loaded with {len(se_df)} rows and {len(se_df.columns)} columns")
    #         self.logger.debug(f"First few sample IDs: {', '.join(se_df['cell_id'].head().tolist())}")

    #         # Log column data types
    #         column_types = {col: str(se_df[col].dtype) for col in se_df.columns}
    #         self.logger.debug(f"Column data types: {column_types}")

    #         return se_df
    #     except Exception as e:
    #         self.logger.exception(f"Error reading BED file: {e}")
    #         raise


    def read_bed_file(self, file_path):
        """
        Read BED file with original column names
        
        Parameters:
        file_path (str): Path to the BED file
        
        Returns:
        tuple: (DataFrame, list of column names)
            - pandas.DataFrame: Loaded BED data
            - list: Original column names
        """
        self.logger.info(f"Reading BED file: {file_path}")
        try:
            # Read file with header, preserving original column names
            se_df = pd.read_csv(file_path, sep='\s+', header=0)
            
            # Store original column names for reference
            original_columns = list(se_df.columns)
            self.logger.debug(f"BED file original columns: {', '.join(original_columns)}")
            
            # Log basic information about the data
            self.logger.info(f"BED file loaded with {len(se_df)} rows and {len(se_df.columns)} columns")
            self.logger.debug(f"First few sample IDs: {', '.join(se_df['cell_id'].head().tolist())}")
            
            # Log column data types
            column_types = {col: str(se_df[col].dtype) for col in se_df.columns}
            self.logger.debug(f"Column data types: {column_types}")
            
            # 明示的に2つの値をタプルとして返す
            return se_df, original_columns
        except Exception as e:
            self.logger.exception(f"Error reading BED file: {e}")
            raise


    def read_sample_info(self, file_path):
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