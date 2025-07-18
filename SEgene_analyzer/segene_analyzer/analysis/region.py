# sedb_tools/analysis/region.py
import pandas as pd
import numpy as np
import logging
import os
import warnings
from collections import Counter
from typing import Optional, Dict, List, Tuple, Any
from pybedtools import BedTool


class SEdbRegionOverlapAnalyzer:
    """
    Super-enhancer (SE) 領域分析クラス
    特定のゲノム領域に重複するスーパーエンハンサーを抽出・分析する
    """
    
    def __init__(self, logger=None):
        """
        Initialize the region overlap analyzer
        
        Parameters:
        logger (Logger): Optional logger instance to use. If None, creates a basic console logger.
        """
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbRegionOverlapAnalyzer")
    
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
    
    def extract_overlapping_se(self, bed_df, sample_info, region_chr, region_start, region_end, 
                               region_name=None, extraction_method='pybedtools'):
        """
        Extract super enhancers that overlap with the specified genomic region
        
        Parameters:
        ----------
        bed_df : pandas.DataFrame
            BEDデータ
        sample_info : pandas.DataFrame
            サンプル情報
        region_chr : str
            染色体名 (例: 'chr1')
        region_start : int
            領域の開始位置
        region_end : int
            領域の終了位置
        region_name : str, optional
            領域の名前（プロットやレポートで使用）
        extraction_method : str
            抽出に使用するメソッド ('pybedtools' or 'pandas')
        
        Returns:
        -------
        dict:
            抽出・分析結果を含む辞書
            {
                'region': {'chr': str, 'start': int, 'end': int, 'name': str},
                'extraction_method': str,
                'overlapping_se': pandas.DataFrame,
                'unique_se_data': pandas.DataFrame,
                'final_data': pandas.DataFrame,
                'unique_sample_count': int,
                'region_tissue_counts': pandas.Series,
                'region_biosample_counts': pandas.Series,
                'region_gene_counts': pandas.Series
            }
        """
        if bed_df is None:
            self.logger.warning("BED data not loaded. Cannot extract overlapping SEs.")
            return None
        
        self.logger.info(f"Extracting SEs overlapping with region: {region_chr}:{region_start}-{region_end}")
        
        # 領域情報を保存
        region = {
            'chr': region_chr,
            'start': region_start,
            'end': region_end,
            'name': region_name or f"{region_chr}:{region_start}-{region_end}"
        }
        
        # 結果を格納する辞書
        results = {
            'region': region,
            'extraction_method': extraction_method,
            'overlapping_se': None,
            'unique_se_data': None,
            'final_data': None,
            'unique_sample_count': 0,
            'region_tissue_counts': None,
            'region_biosample_counts': None,
            'region_gene_counts': None
        }
        
        try:
            # 抽出メソッドに基づいて処理
            if extraction_method == 'pybedtools':
                overlapping_se = self._extract_with_pybedtools(bed_df, region_chr, region_start, region_end)
            elif extraction_method == 'pandas':
                overlapping_se = self._extract_with_pandas(bed_df, region_chr, region_start, region_end)
            else:
                self.logger.warning(f"Unknown extraction method: {extraction_method}. Using pybedtools.")
                overlapping_se = self._extract_with_pybedtools(bed_df, region_chr, region_start, region_end)
            
            # 抽出結果を保存
            results['overlapping_se'] = overlapping_se
            
            # 結果の分析（overlapping_seが正常に作成された場合）
            if overlapping_se is not None and not overlapping_se.empty:
                analysis_results = self._analyze_overlapping_se(overlapping_se, sample_info, region)
                # 分析結果を辞書にマージ
                results.update(analysis_results)
                
                # ユニークSEデータを取得
                unique_se_data = self.get_unique_se_data(overlapping_se)
                results['unique_se_data'] = unique_se_data
            
            return results
            
        except Exception as e:
            self.logger.exception(f"Error extracting overlapping SEs: {e}")
            return results

    def _extract_with_pybedtools(self, bed_df, region_chr, region_start, region_end):
        """
        Extract overlapping SEs using pybedtools
        
        Parameters:
        ----------
        bed_df : pandas.DataFrame
            BEDデータ
        region_chr : str
            染色体名
        region_start : int
            領域の開始位置
        region_end : int
            領域の終了位置
            
        Returns:
        -------
        pandas.DataFrame
            重複するSEのデータフレーム
        """
        try:
            import pybedtools
            from pybedtools import BedTool
            
            # 必要なカラムを確認
            required_columns = ['se_chr', 'se_start', 'se_end']
            
            # 必要なカラムが存在するか確認
            missing_cols = [col for col in required_columns if col not in bed_df.columns]
            if missing_cols:
                self.logger.warning(f"Missing columns in BED data: {missing_cols}")
                self.logger.warning("Trying to find alternative column names...")
                
                # 代替カラム名を検索
                alt_cols = {'se_chr': ['chr'], 'se_start': ['start'], 'se_end': ['end']}
                col_mapping = {}
                
                for req_col in missing_cols:
                    for alt_col in alt_cols.get(req_col, []):
                        if alt_col in bed_df.columns:
                            col_mapping[alt_col] = req_col
                            self.logger.info(f"Using '{alt_col}' for '{req_col}'")
                            break
                
                # 代替カラムが見つからない場合は処理を中止
                still_missing = [col for col in missing_cols if col not in col_mapping.values()]
                if still_missing:
                    self.logger.error(f"Could not find alternatives for columns: {still_missing}. Cannot proceed.")
                    return None
                
                # 一時的なコピーを作成し、カラム名を変更
                temp_df = bed_df.copy()
                temp_df.rename(columns=col_mapping, inplace=True)
            else:
                temp_df = bed_df
            
            # BedToolオブジェクトに変換
            try:
                # BEDフォーマット用に必要なカラムのみ抽出
                bed_cols = ['se_chr', 'se_start', 'se_end', 'cell_id']
                se_df_bed = temp_df[bed_cols].copy()
                # 標準BEDフォーマットにカラム名を変更
                se_df_bed.columns = ['chrom', 'start', 'end', 'name']
                
                # BedToolに変換
                se_bedtool = BedTool.from_dataframe(se_df_bed)
                
                # 対象領域を定義
                interval_str = f"{region_chr}\t{region_start}\t{region_end}"
                interval_bedtool = BedTool(interval_str, from_string=True)
                
                # 重複するSEを検索
                intersect_result = se_bedtool.intersect(interval_bedtool, wa=True)
                
                # 結果をDataFrameに変換
                result_df = pd.read_table(
                    intersect_result.fn, 
                    names=['chrom', 'start', 'end', 'name']
                )
                
                # 重複するcell_idを抽出
                overlapping_cell_ids = set(result_df['name'])
                
                # 元のデータから重複するSEに対応する行を取得
                overlapping_se = temp_df[temp_df['cell_id'].isin(overlapping_cell_ids)].copy()
                
                self.logger.info(f"Found {len(overlapping_se)} SEs overlapping with the region")
                
                return overlapping_se
                
            except Exception as e:
                self.logger.exception(f"Error using pybedtools: {e}")
                return None
                
        except ImportError:
            self.logger.error("pybedtools module not available. Please install it or use 'pandas' method.")
            return self._extract_with_pandas(bed_df, region_chr, region_start, region_end)
    
    def _extract_with_pandas(self, bed_df, region_chr, region_start, region_end):
        """
        Extract overlapping SEs using pandas (no external dependencies)
        
        Parameters:
        ----------
        bed_df : pandas.DataFrame
            BEDデータ
        region_chr : str
            染色体名
        region_start : int
            領域の開始位置
        region_end : int
            領域の終了位置
            
        Returns:
        -------
        pandas.DataFrame
            重複するSEのデータフレーム
        """
        # 必要なカラムを確認
        required_columns = ['se_chr', 'se_start', 'se_end']
        alt_columns = {'se_chr': 'chr', 'se_start': 'start', 'se_end': 'end'}
        
        col_mapping = {}
        for req_col in required_columns:
            if req_col in bed_df.columns:
                col_mapping[req_col] = req_col
            elif alt_columns[req_col] in bed_df.columns:
                col_mapping[alt_columns[req_col]] = req_col
            else:
                self.logger.error(f"Could not find column for {req_col}. Cannot proceed.")
                return None
        
        # カラム名を取得
        chr_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_chr')
        start_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_start')
        end_col = next(col for col, mapped in col_mapping.items() if mapped == 'se_end')
        
        # 領域と重複するSEをフィルタリング
        overlapping_se = bed_df[
            (bed_df[chr_col] == region_chr) & 
            (bed_df[start_col] <= region_end) & 
            (bed_df[end_col] >= region_start)
        ].copy()
        
        self.logger.info(f"Found {len(overlapping_se)} SEs overlapping with the region")
        return overlapping_se
    
    def _analyze_overlapping_se(self, overlapping_se, sample_info, region):
        """
        Analyze the overlapping SEs to calculate tissues, biosamples, and gene counts
        
        Parameters:
        ----------
        overlapping_se : pandas.DataFrame
            重複するSEのデータフレーム
        sample_info : pandas.DataFrame
            サンプル情報のデータフレーム
        region : dict
            領域情報の辞書
            
        Returns:
        -------
        dict
            分析結果を含む辞書
        """
        if overlapping_se is None or overlapping_se.empty or sample_info is None:
            self.logger.warning("No overlapping SEs or sample info available.")
            return {
                'final_data': None,
                'unique_sample_count': 0,
                'region_tissue_counts': None,
                'region_biosample_counts': None,
                'region_gene_counts': None
            }
        
        self.logger.info("Analyzing overlapping SEs based on unique samples")
        
        # 重複するSEから一意のcell_idを取得
        overlapping_cell_ids = set(overlapping_se['cell_id'])
        unique_sample_count = len(overlapping_cell_ids)
        self.logger.info(f"Found {unique_sample_count} unique samples with overlapping SEs")
        
        # sample_infoから一意の重複するセルIDに対応する行だけをフィルタリング
        # 各サンプルはこのフィルタリングされたデータフレームで正確に1回だけ表現される
        overlapping_samples = sample_info[sample_info['cell_id'].isin(overlapping_cell_ids)]
        
        # tissue_countsを計算 - 各サンプルが正確に1回カウントされる
        region_tissue_counts = None
        if 'Tissue type' in overlapping_samples.columns:
            region_tissue_counts = overlapping_samples['Tissue type'].value_counts()
            self.logger.info(f"Found {len(region_tissue_counts)} tissue types across {unique_sample_count} unique samples")
        else:
            self.logger.warning("'Tissue type' column not found in sample information")
        
        # biosample_countsを計算 - 各サンプルが正確に1回カウントされる
        region_biosample_counts = None
        if 'Biosample type' in overlapping_samples.columns:
            region_biosample_counts = overlapping_samples['Biosample type'].value_counts()
            self.logger.info(f"Found {len(region_biosample_counts)} biosample types across {unique_sample_count} unique samples")
        else:
            self.logger.warning("'Biosample type' column not found in sample information")
        
        # gene_countsを計算
        # 遺伝子カウントについては、同じサンプルからの複数のSEが異なる遺伝子に関連付けられる可能性があるため、
        # すべての重複するSEを使用します
        region_gene_counts = None
        gene_column = next((col for col in overlapping_se.columns if 'gene' in col.lower()), None)
        if gene_column:
            all_genes = []
            for entry in overlapping_se[gene_column].dropna():
                if isinstance(entry, str):
                    genes = entry.split(',')
                    all_genes.extend([g.strip() for g in genes if g.strip()])
            
            region_gene_counts = pd.Series(Counter(all_genes))
            self.logger.info(f"Found {len(region_gene_counts)} unique genes in overlapping SEs")
        else:
            self.logger.warning("No gene column found in BED data")
            
        # 統合データフレームを作成
        try:
            # overlapping_seとsample_infoをマージして組織タイプなどを取得
            final_data = pd.merge(
                overlapping_se, 
                sample_info,
                on='cell_id',
                how='left'
            )
            
            # 領域情報を追加
            for key, value in region.items():
                final_data[f'region_{key}'] = value
                
            self.logger.info(f"Created consolidated data with {len(final_data)} rows from {unique_sample_count} unique samples")
            
            # サンプルごとの概要データフレームも作成（各サンプルが1回だけ表示される）
            sample_summary = overlapping_samples.copy()
            
            # サンプルごとのSEカウントを追加
            se_counts = overlapping_se['cell_id'].value_counts().reset_index()
            se_counts.columns = ['cell_id', 'se_count']
            sample_summary = pd.merge(
                sample_summary,
                se_counts,
                on='cell_id',
                how='left'
            )
            
            self.logger.info(f"Created sample summary data with {len(sample_summary)} rows (one per unique sample)")
            
        except Exception as e:
            self.logger.warning(f"Could not create consolidated data: {e}")
            final_data = None
            sample_summary = None
            
        return {
            'final_data': final_data,
            'sample_summary': sample_summary,
            'unique_sample_count': unique_sample_count,
            'region_tissue_counts': region_tissue_counts,
            'region_biosample_counts': region_biosample_counts,
            'region_gene_counts': region_gene_counts
        }
    
    def get_unique_se_data(self, overlapping_se):
        """
        Get unique SE data for the region (one row per unique SE)
        
        Parameters:
        ----------
        overlapping_se : pandas.DataFrame
            重複するSEのデータフレーム
            
        Returns:
        -------
        pandas.DataFrame
            ユニークなSEデータ
        """
        if overlapping_se is None or overlapping_se.empty:
            self.logger.warning("No overlapping SEs. Cannot extract unique SE data.")
            return None
        
        # SEの座標でグループ化してユニークなSEを取得
        try:
            if 'se_id' in overlapping_se.columns:
                # SE IDがある場合はそれを使用
                unique_se = overlapping_se.drop_duplicates('se_id')
            else:
                # そうでなければ座標を使用
                coords_cols = [col for col in ['se_chr', 'se_start', 'se_end'] 
                            if col in overlapping_se.columns]
                if len(coords_cols) < 3:
                    self.logger.warning("Missing coordinate columns, using fallback.")
                    coords_cols = [col for col in ['chr', 'start', 'end'] 
                                if col in overlapping_se.columns]
                
                if len(coords_cols) < 3:
                    self.logger.error("Cannot identify coordinate columns. Aborting.")
                    return None
                    
                unique_se = overlapping_se.drop_duplicates(coords_cols)
            
            self.logger.info(f"Found {len(unique_se)} unique SEs in the region")
            
            return unique_se
            
        except Exception as e:
            self.logger.exception(f"Error getting unique SE data: {e}")
            return None
        

    # def get_overlapping_regions(self, bed_df, region_chr, region_start, region_end):
    #     """
    #     pybedtoolsを使用して指定された領域と重複するゲノム領域を抽出
    #     現状は４カラム設計
    #     Parameters:
    #     -----------
    #     bed_df : pandas.DataFrame
    #         BEDデータフレーム
    #     region_chr : str
    #         染色体名 (例: 'chr1')
    #     region_start : int
    #         領域の開始位置
    #     region_end : int
    #         領域の終了位置
        
    #     Returns:
    #     --------
    #     pandas.DataFrame
    #         重複する領域のデータフレーム
    #     """
        
    #     self.logger.info(f"Extracting regions overlapping with {region_chr}:{region_start}-{region_end}")
        
    #     try:
    #         # BEDフォーマット用にデータ準備
    #         bed_cols = ['se_chr', 'se_start', 'se_end', 'cell_id']
    #         se_df_bed = bed_df[bed_cols].copy()
    #         se_df_bed.columns = ['chrom', 'start', 'end', 'name']
            
    #         # BedToolオブジェクトに変換
    #         se_bedtool = BedTool.from_dataframe(se_df_bed)
            
    #         # 対象領域を定義
    #         interval_str = f"{region_chr}\t{region_start}\t{region_end}"
    #         interval_bedtool = BedTool(interval_str, from_string=True)
            
    #         # 重複検出
    #         intersect_result = se_bedtool.intersect(interval_bedtool, wa=True)
            
    #         # 結果をDataFrameに変換
    #         result_df = pd.read_table(
    #             intersect_result.fn, 
    #             names=['chrom', 'start', 'end', 'name']
    #         )
            
    #         # 重複するcell_idを抽出
    #         overlapping_cell_ids = set(result_df['name'])
            
    #         # 元のデータから重複するSEに対応する行を取得
    #         overlapping_df = bed_df[bed_df['cell_id'].isin(overlapping_cell_ids)].copy()
            
    #         self.logger.info(f"Found {len(overlapping_df)} regions from {len(overlapping_cell_ids)} unique samples")
            
    #         return overlapping_df
            
    #     except Exception as e:
    #         self.logger.exception(f"Error in get_overlapping_regions: {e}")
    #         return pd.DataFrame()  # 空のデータフレームを返す
        

    # def calculate_sample_statistics(self, overlapping_df, sample_info=None):
    #     """
    #     重複領域の基本サンプル統計を計算します。
        
    #     Parameters:
    #     -----------
    #     overlapping_df : pandas.DataFrame
    #         get_overlapping_regions()によって返された重複領域
    #     sample_info : pandas.DataFrame, optional
    #         サンプル情報データ
            
    #     Returns:
    #     --------
    #     dict
    #         {
    #             'unique_sample_ids': set,  # ユニークなサンプルIDのセット
    #             'sample_count': int,       # ユニークなサンプル数
    #             'region_count': int,       # 領域の総数
    #             'samples_df': DataFrame    # サンプル情報がある場合のみ
    #         }
    #     """
    #     if overlapping_df is None or overlapping_df.empty:
    #         self.logger.warning("Empty overlapping regions data")
    #         return {
    #             'unique_sample_ids': set(),
    #             'sample_count': 0,
    #             'region_count': 0,
    #             'samples_df': None
    #         }
        
    #     # ユニークなサンプルIDと数
    #     unique_sample_ids = set(overlapping_df['cell_id'])
        
    #     stats = {
    #         'unique_sample_ids': unique_sample_ids,
    #         'sample_count': len(unique_sample_ids),
    #         'region_count': len(overlapping_df),
    #         'samples_df': None
    #     }
        
    #     # サンプル情報がある場合、サンプルデータフレームを作成
    #     if sample_info is not None:
    #         samples_df = sample_info[sample_info['cell_id'].isin(unique_sample_ids)].copy()
    #         stats['samples_df'] = samples_df
        
    #     return stats