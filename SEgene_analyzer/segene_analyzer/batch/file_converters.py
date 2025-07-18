# SEgene_analyzer/batch/file_converters.py
import pandas as pd
import os
import logging
from datetime import datetime
from typing import Optional, Tuple

class FileConverter:
    """
    ゲノム領域ファイルの形式変換を担当するクラス
    """
    
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
    
    def detect_file_format(self, file_path: str) -> str:
        """
        ファイル形式を判別（BED or TSV）
        
        Returns:
        --------
        str: 'bed' または 'tsv'
        """
        # 拡張子チェック
        if file_path.lower().endswith('.bed'):
            return 'bed'
        elif file_path.lower().endswith('.tsv'):
            return 'tsv'
        
        # ファイル内容チェック
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                
            # ヘッダー行の存在でTSV判別
            if 'se_data' in first_line.lower():
                return 'tsv'
            else:
                return 'bed'
                
        except Exception as e:
            self.logger.warning(f"Could not detect file format: {e}. Assuming TSV.")
            return 'tsv'
    
    def convert_bed_to_tsv(self, bed_file: str, output_dir: str, 
                          output_filename: Optional[str] = None) -> str:
        """
        BEDファイルを専用TSV形式に変換
        
        Parameters:
        -----------
        bed_file : str
            入力BEDファイルのパス
        output_dir : str
            出力ディレクトリ
        output_filename : str, optional
            出力ファイル名（Noneの場合は自動生成）
            
        Returns:
        --------
        str: 変換されたTSVファイルのパス
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # BEDファイル読み込み
        try:
            bed_df = pd.read_csv(bed_file, sep='\t', header=None)
            
            # 最低限の列数チェック
            if len(bed_df.columns) < 3:
                raise ValueError("BED file must have at least 3 columns (chrom, start, end)")
            
            # 列名を設定
            bed_df.columns = ['chrom', 'start', 'end'] + [f'col{i}' for i in range(3, len(bed_df.columns))]
            
            self.logger.info(f"Loaded BED file with {len(bed_df)} regions")
            
        except Exception as e:
            self.logger.error(f"Error reading BED file {bed_file}: {e}")
            raise
        
        # se_data列を作成（chr_start_end形式）
        bed_df['se_data'] = (bed_df['chrom'].astype(str) + '_' + 
                            bed_df['start'].astype(str) + '_' + 
                            bed_df['end'].astype(str))
        
        # 専用TSV形式のDataFrameを作成
        tsv_df = pd.DataFrame({
            'se_data': bed_df['se_data'],
            'se_count': 0,           # デフォルト値
            'sample_count': 0,       # デフォルト値  
            'gene_list': '[]',       # 空リスト
            'gene_count': 0          # デフォルト値
        })
        
        # 出力ファイル名の設定
        if output_filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            base_name = os.path.splitext(os.path.basename(bed_file))[0]
            output_filename = f"converted_{base_name}_{timestamp}.tsv"
        
        # TSVファイルとして保存
        output_path = os.path.join(output_dir, output_filename)
        tsv_df.to_csv(output_path, sep='\t', index=False)
        
        self.logger.info(f"Converted BED to TSV: {bed_file} → {output_path}")
        self.logger.info(f"Generated {len(tsv_df)} region entries")
        
        return output_path