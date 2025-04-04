"""
SEdb Region Extraction and Analysis Module
Provides functions for extracting and analyzing super-enhancers overlapping a specified genomic region.
"""

import pandas as pd
import numpy as np
import logging
from collections import Counter
from typing import Dict, Any, Optional, List, Set
from pathlib import Path

from sedb_common.logging import get_module_logger

import pybedtools
from pybedtools import BedTool

import os

class RegionOverlapAnalyzer:
    """Class for analyzing super-enhancers overlapping a specified genomic region."""
    
    def __init__(self, logger=None):
        """Initialize the RegionOverlapAnalyzer instance."""
        self.logger = logger or get_module_logger(__name__)
    

    # def extract_region(self, region_chr: str, region_start: int, region_end: int, 
    #                   region_name: Optional[str] = None, bed_df: Optional[pd.DataFrame] = None, 
    #                   sample_info: Optional[pd.DataFrame] = None, 
    #                   extraction_method: str = 'pybedtools',
    #                   analysis_level: str = 'basic') -> Dict[str, Any]:

    def extract_region(self, region_chr: str, region_start: int, region_end: int,
                      region_name: Optional[str] = None, bed_df: Optional[pd.DataFrame] = None,
                      sample_info: Optional[pd.DataFrame] = None,
                      extraction_method: str = 'pybedtools',
                      analysis_level: str = 'basic',
                      output_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Extract super-enhancers overlapping the specified region.
        
        Parameters:
        -----------
        region_chr : str
            Chromosome name (e.g., 'chr1').
        region_start : int
            Start position of the region.
        region_end : int
            End position of the region.
        region_name : str, optional
            Name of the region.
        bed_df : pandas.DataFrame, optional
            BED data (required).
        sample_info : pandas.DataFrame, optional
            Sample information (required).
        extraction_method : str, optional
            Extraction method ('pybedtools' or 'pandas').
        analysis_level : str, optional
            Analysis level ('none' or 'basic').
        output_dir : str, optional
            Output directory for saving debug files.
            
        Returns:
        --------
        dict
            Dictionary containing extraction results.
        """
        if bed_df is None or sample_info is None:
            self.logger.error("BED data and sample information are required")
            raise ValueError("BED data and sample information are required")
        
        # 領域情報を準備
        region = {
            'chr': region_chr,
            'start': region_start,
            'end': region_end,
            'name': region_name or f"{region_chr}:{region_start}-{region_end}",
            'size': region_end - region_start
        }
        
        self.logger.info(f"Extracting overlapping super-enhancers for region {region['name']} ({region_chr}:{region_start}-{region_end})...")
        
        # 抽出メソッドの選択
        intersect_result_bedtool_debug = None # デバッグ用ファイルが作られたかの確認用（必須ではない）
        if extraction_method == 'pybedtools':

            # overlapping_se = self._extract_with_pybedtools(bed_df, region_chr, region_start, region_end)
            overlapping_se = self._extract_with_pybedtools(bed_df, region_chr, region_start, region_end, output_dir=output_dir)


        else:  # pandas
            overlapping_se = self._extract_with_pandas(bed_df, region_chr, region_start, region_end)
        
        if overlapping_se is None or len(overlapping_se) == 0:
            self.logger.warning(f"No overlapping super-enhancers found in region {region['name']}")
            return {
                'region': region,
                'extraction_method': extraction_method,
                'overlapping_se': pd.DataFrame(),
                'sample_summary': pd.DataFrame(),
                'stats': {
                    'total_regions': 0,
                    'unique_samples': 0,
                    'tissue_types': 0,
                    'biosample_types': 0
                }
            }
        
        self.logger.info(f"Found {len(overlapping_se)} overlapping super-enhancers in the region.")
        
        # 結果をまとめる基本情報
        results = {
            'region': region,
            'extraction_method': extraction_method,
            'overlapping_se': overlapping_se,
            'stats': {
                'total_regions': len(overlapping_se)
            }
        }
        
        # 分析レベルに応じた処理
        if analysis_level == 'basic':
            # 重複するSEのサンプル情報を分析
            self.logger.info("Performing basic analysis...")
            analysis_results = self._analyze_overlapping_samples(overlapping_se, sample_info, region)
            # ファイル保存処理
            overlapping_ids_set = analysis_results.get('overlapping_cell_id_set')
            if output_dir and overlapping_ids_set: # output_dir が渡されていて、セットが存在する場合
                ids_filename = "debug_overlapping_sample_ids.txt" # 保存ファイル名
                ids_filepath = os.path.join(output_dir, ids_filename)
                try:
                    # セットをリストに変換してソートし、ファイルに書き出す
                    sorted_ids = sorted(list(overlapping_ids_set))
                    with open(ids_filepath, 'w') as f:
                        for sample_id in sorted_ids:
                            f.write(f"{sample_id}\n")
                    self.logger.info(f"Debug: Saved {len(sorted_ids)} overlapping sample IDs to {ids_filepath}")
                except Exception as e:
                    self.logger.error(f"Debug: Failed to save overlapping sample IDs to {ids_filepath}: {e}")





            # 分析結果を追加
            results['sample_summary'] = analysis_results['sample_summary']
            results['stats'].update({
                'unique_samples': analysis_results['unique_sample_count'],
                'tissue_types': len(analysis_results['tissue_counts']) if analysis_results['tissue_counts'] is not None else 0,
                'biosample_types': len(analysis_results['biosample_counts']) if analysis_results['biosample_counts'] is not None else 0,
                'gene_count': len(analysis_results['gene_counts']) if analysis_results['gene_counts'] is not None else 0
            })
            
            # 分布データを追加
            for dist_key in ['tissue_counts', 'biosample_counts', 'gene_counts']:
                if analysis_results[dist_key] is not None:
                    results[dist_key] = analysis_results[dist_key]
        else:
            # basic分析なしの場合は最小限の統計情報のみ設定
            self.logger.info("Running in no-analysis mode...")
            
            # サンプルのユニーク数だけは基本情報として取得
            cell_id_col = 'cell_id' if 'cell_id' in overlapping_se.columns else 'name'
            unique_samples = set(overlapping_se[cell_id_col])
            results['stats']['unique_samples'] = len(unique_samples)
        
        return results
    
    # def _extract_with_pybedtools(self, bed_df, region_chr, region_start, region_end):
    #     """pybedtoolsを使用して重複するSEを抽出"""
    #     try:
            
    #         self.logger.info("pybedtoolsを使用して重複領域を抽出")
            
    #         # BEDフォーマット用にデータを準備
    #         bed_cols = ['se_chr', 'se_start', 'se_end', 'cell_id']
    #         # カラム名を確認
    #         actual_columns = {}
    #         for needed_col in ['se_chr', 'se_start', 'se_end', 'cell_id']:
    #             if needed_col in bed_df.columns:
    #                 actual_columns[needed_col] = needed_col
    #             else:
    #                 alt_names = {'se_chr': ['chr'], 'se_start': ['start'], 
    #                             'se_end': ['end'], 'cell_id': ['name', 'sample_id']}
    #                 found = False
    #                 for alt in alt_names.get(needed_col, []):
    #                     if alt in bed_df.columns:
    #                         actual_columns[needed_col] = alt
    #                         found = True
    #                         break
    #                 if not found:
    #                     raise ValueError(f"カラム '{needed_col}' またはその代替が見つかりません")
            
    #         # データフレームを準備
    #         se_df_bed = bed_df[[actual_columns[col] for col in bed_cols]].copy()
    #         se_df_bed.columns = ['chrom', 'start', 'end', 'name']
            
    #         # BedToolオブジェクトに変換
    #         se_bedtool = BedTool.from_dataframe(se_df_bed)
            
    #         # 対象領域を定義
    #         interval_str = f"{region_chr}\t{region_start}\t{region_end}"
    #         interval_bedtool = BedTool(interval_str, from_string=True)
    #         
    #         # 重複するSEを検索
    #         intersect_result = se_bedtool.intersect(interval_bedtool, wa=True)
    #         
    #         # 結果をDataFrameに変換
    #         result_df = pd.read_table(
    #             intersect_result.fn, 
    #             names=['chrom', 'start', 'end', 'name']
    #         )
    #         
    #         # 重複するcell_idを抽出
    #         overlapping_cell_ids = set(result_df['name'])
    #         
    #         # 元のデータから重複するSEに対応する行を取得
    #         overlapping_se = bed_df[bed_df[actual_columns['cell_id']].isin(overlapping_cell_ids)].copy()
    #         
    #         return overlapping_se
    #         
    #     except ImportError:
    #         self.logger.warning("pybedtoolsモジュールが利用できません。pandasメソッドに切り替えます")
    #         return self._extract_with_pandas(bed_df, region_chr, region_start, region_end)
    #     except Exception as e:
    #         self.logger.exception(f"pybedtoolsによる抽出中にエラーが発生しました: {e}")
    #         return None


    # def _extract_with_pybedtools(self, bed_df, region_chr, region_start, region_end):
    #     """pybedtoolsを使用して重複するSEを抽出 (修正案)"""
    def _extract_with_pybedtools(self, bed_df, region_chr, region_start, region_end, output_dir=None):
        """Extracts overlapping super-enhancers using pybedtools (with debug file output)."""

        try:
            self.logger.info("Extracting overlapping regions using pybedtools")

            # 対象領域を BedTool オブジェクトとして定義
            interval_str = f"{region_chr}\t{region_start}\t{region_end}"
            interval_bedtool = BedTool(interval_str, from_string=True)
            self.logger.debug(f"Target interval: {interval_str}")

            # --- ここから修正 ---
            # 元の bed_df から BedTool を作成する際、必要な列のみにする
            # pybedtools が必要とするのは chrom, start, end (必須), name (任意だがcell_idを使う)
            # 実際のカラム名を取得 (以前と同様のロジック)
            actual_columns = {}
            required_for_bedtool = ['se_chr', 'se_start', 'se_end', 'cell_id']
            missing_cols = []
            for col in required_for_bedtool:
                 potential_names = [col]
                 if col == 'se_chr': potential_names.append('chr')
                 if col == 'se_start': potential_names.append('start')
                 if col == 'se_end': potential_names.append('end')
                 if col == 'cell_id': potential_names.extend(['name', 'sample_id'])
                 found_name = next((name for name in potential_names if name in bed_df.columns), None)
                 if found_name:
                     actual_columns[col] = found_name
                 else:
                     missing_cols.append(col)
            if missing_cols:
                 raise ValueError(f"BED DataFrame missing columns needed for BedTool: {', '.join(missing_cols)}")

            # BedTool 用の一時的な DataFrame を作成 (必須列のみ)
            temp_bed_df = bed_df[[
                actual_columns['se_chr'],
                actual_columns['se_start'],
                actual_columns['se_end'],
                actual_columns['cell_id']
            ]].copy()
            # pybedtools 用に列名を変更
            temp_bed_df.columns = ['chrom', 'start', 'end', 'name']

            # DataFrame から BedTool オブジェクトを作成
            se_bedtool = BedTool.from_dataframe(temp_bed_df)
            self.logger.debug(f"Created BedTool from DataFrame with {len(se_bedtool)} features.")

            # intersect を実行 (-wa オプションで重複した元の特徴量を取得)
            # intersect_result は指定領域にオーバーラップする特徴量のみを含む BedTool オブジェクト
            intersect_result = se_bedtool.intersect(interval_bedtool, wa=True)
            self.logger.debug(f"Intersection found {len(intersect_result)} overlapping features.")

            # --- ★デバッグ用ファイル出力処理を追加 ---
            if output_dir and intersect_result and len(intersect_result) > 0:
                debug_filename = "debug_intersect_result.bed" # デバッグ用ファイル名
                debug_filepath = os.path.join(output_dir, debug_filename)
                try:
                    saved_path = intersect_result.saveas(debug_filepath)
                    self.logger.info(f"Debug: Direct intersection result saved to: {saved_path.fn}")
                except Exception as e:
                    self.logger.error(f"Debug: Failed to save direct intersection result to {debug_filepath}: {e}")
            elif output_dir:
                 self.logger.info("Debug: No intersection results found, skipping debug file save.")
            # --- ★デバッグ用出力ここまで ---

            if len(intersect_result) > 0:
                # intersect の結果 (BedTool オブジェクト) を直接 DataFrame に変換する
                # この DataFrame は指定領域内の SE の基本情報 (chrom, start, end, name=cell_id) のみを含む
                overlapping_coords_df = intersect_result.to_dataframe()
                self.logger.debug(f"Overlapping coordinates DataFrame head:\n{overlapping_coords_df.head()}")

                # 元の bed_df と intersect 結果の座標 DataFrame をマージして、
                # 指定領域内の SE に対応する元の全カラム情報を取得する
                # マージキーとして使うために、intersect 結果のカラム名を元の名前に戻す
                overlapping_coords_df = overlapping_coords_df.rename(columns={
                    'chrom': actual_columns['se_chr'],
                    'start': actual_columns['se_start'],
                    'end': actual_columns['se_end'],
                    'name': actual_columns['cell_id']
                })

                # マージ実行 (inner join)
                # これにより、座標と cell_id が一致する行のみが抽出される
                overlapping_se = pd.merge(
                    bed_df,
                    overlapping_coords_df,
                    on=[actual_columns['se_chr'], actual_columns['se_start'], actual_columns['se_end'], actual_columns['cell_id']],
                    how='inner'
                )
                self.logger.info(f"Successfully merged to get {len(overlapping_se)} full SE records for the target region.")

            else:
                 # 重複がなければ空の DataFrame を返す (カラムは元の bed_df に合わせる)
                 overlapping_se = pd.DataFrame(columns=bed_df.columns)
                 self.logger.info("No overlapping SE records found in the target region.")

            # --- 修正ここまで ---

            # 正しくフィルタリングされた DataFrame を返す
            return overlapping_se

        except ImportError:
            self.logger.warning("pybedtools module is not available. Switching to pandas method.")
            return self._extract_with_pandas(bed_df, region_chr, region_start, region_end)
        except Exception as e:
            self.logger.exception(f"An error occurred during extraction with pybedtools: {e}")
            return pd.DataFrame(columns=bed_df.columns) # エラー時も空のDataFrame


    
    def _extract_with_pandas(self, bed_df, region_chr, region_start, region_end):
        """Extracts overlapping super-enhancers using pandas."""
        try:
            self.logger.info("Extracting overlapping regions using pandas")
            
            # カラム名を確認
            chr_col = 'se_chr' if 'se_chr' in bed_df.columns else 'chr'
            start_col = 'se_start' if 'se_start' in bed_df.columns else 'start'
            end_col = 'se_end' if 'se_end' in bed_df.columns else 'end'
            
            # 重複領域をフィルタリング
            overlapping_se = bed_df[
                (bed_df[chr_col] == region_chr) & 
                (bed_df[start_col] <= region_end) & 
                (bed_df[end_col] >= region_start)
            ].copy()
            
            return overlapping_se
            
        except Exception as e:
            self.logger.exception(f"An error occurred during extraction with pandas: {e}")
            return None
    
    def _analyze_overlapping_samples(self, overlapping_se, sample_info, region):
        """Analyzes sample information for overlapping super-enhancers."""
        # 重複するSEからcell_idを取得
        cell_id_col = 'cell_id' if 'cell_id' in overlapping_se.columns else 'name'
        overlapping_cell_ids = set(overlapping_se[cell_id_col])
        unique_sample_count = len(overlapping_cell_ids)
        
        self.logger.info(f"Found {unique_sample_count} unique samples overlapping the region.")
        
        # サンプル情報からcell_idに対応する行をフィルタリング
        sample_id_col = 'cell_id' if 'cell_id' in sample_info.columns else 'Sample ID'
        sample_summary = sample_info[sample_info[sample_id_col].isin(overlapping_cell_ids)].copy()
        
        # サンプルごとのSE数をカウント
        se_counts = overlapping_se[cell_id_col].value_counts().reset_index()
        se_counts.columns = [sample_id_col, 'se_count']
        
        # サンプル情報にSE数を結合
        sample_summary = pd.merge(
            sample_summary,
            se_counts,
            on=sample_id_col,
            how='left'
        )
        

        self.logger.debug(f"Sample Summary DataFrame shape after filtering: {sample_summary.shape}") # 474行あるか？
        self.logger.debug(f"Sample Summary columns: {list(sample_summary.columns)}") # カラム名を確認
        if not sample_summary.empty:
            self.logger.debug(f"Sample Summary head:\n{sample_summary.head()}") # 最初の数行を表示
            # Tissue type 列の確認
            if 'Tissue type' in sample_summary.columns:
                # ↓ NaN でない値の数をカウント
                self.logger.debug(f"Sample Summary 'Tissue type' non-null count: {sample_summary['Tissue type'].count()}")
                # ↓ ユニークな値（NaNも含む可能性あり）の最初のいくつかを表示
                self.logger.debug(f"Sample Summary 'Tissue type' unique values (first 10): {sample_summary['Tissue type'].unique()[:10]}")
            else:
                self.logger.warning("Column 'Tissue type' not found in sample_summary!")
            # Biosample type 列も同様に確認
            if 'Biosample type' in sample_summary.columns:
                self.logger.debug(f"Sample Summary 'Biosample type' non-null count: {sample_summary['Biosample type'].count()}")
                self.logger.debug(f"Sample Summary 'Biosample type' unique values (first 10): {sample_summary['Biosample type'].unique()[:10]}")
            else:
                self.logger.warning("Column 'Biosample type' not found in sample_summary!")



        
        # 組織タイプのカウント
        tissue_counts = None
        if 'Tissue type' in sample_summary.columns:
            tissue_counts = sample_summary['Tissue type'].value_counts()
            self.logger.info(f"Found {len(tissue_counts)} tissue types.")
        
        # バイオサンプルタイプのカウント
        biosample_counts = None
        if 'Biosample type' in sample_summary.columns:
            biosample_counts = sample_summary['Biosample type'].value_counts()
            self.logger.info(f"Found {len(biosample_counts)} biosample types.")
        
        # 遺伝子カウント
        gene_counts = None
        gene_column = next((col for col in overlapping_se.columns if 'gene' in col.lower()), None)
        if gene_column:
            all_genes = []
            for entry in overlapping_se[gene_column].dropna():
                if isinstance(entry, str):
                    genes = entry.split(',')
                    all_genes.extend([g.strip() for g in genes if g.strip()])
            
            gene_counts = pd.Series(Counter(all_genes))
            self.logger.info(f"Found {len(gene_counts)} unique genes.")
        
        return {
            'sample_summary': sample_summary,
            'unique_sample_count': unique_sample_count,
            'tissue_counts': tissue_counts,
            'biosample_counts': biosample_counts,
            'gene_counts': gene_counts,
            'overlapping_cell_id_set': overlapping_cell_ids # <<<  (ID のセット)

        }
    
    def export_to_bed(self, overlapping_se, output_path, format='classic'):
        """
        Export overlapping super-enhancers as a BED file.
        
        Parameters:
        -----------
        overlapping_se : pandas.DataFrame
            Overlapping super-enhancer data.
        output_path : str
            Output file path.
        format : str
            Output format ('classic' or 'extended').
        
        Returns:
        --------
        str
            File path of the exported file.
        """
        if overlapping_se is None or len(overlapping_se) == 0:
            self.logger.warning("No super-enhancer data available for export")
            return None
        
        self.logger.info(f"Exporting BED file in {format} format: {output_path}")
        
        try:
            # カラム名を確認
            chr_col = 'se_chr' if 'se_chr' in overlapping_se.columns else 'chr'
            start_col = 'se_start' if 'se_start' in overlapping_se.columns else 'start'
            end_col = 'se_end' if 'se_end' in overlapping_se.columns else 'end'
            name_col = 'cell_id' if 'cell_id' in overlapping_se.columns else 'name'
            
            if format == 'classic':
                # 基本的なBEDフォーマット (chrom, start, end, name)
                bed_data = overlapping_se[[chr_col, start_col, end_col, name_col]]
                # BEDファイルとして出力
                bed_data.to_csv(output_path, sep='\t', header=False, index=False)
            else:  # 'extended'
                # 拡張BEDフォーマット (より多くの情報を含む)
                # 基本列を含める
                bed_cols = [chr_col, start_col, end_col, name_col]
                
                # 利用可能な追加の列を探す
                extra_cols = []
                gene_col = next((col for col in overlapping_se.columns if 'gene' in col.lower()), None)
                if gene_col:
                    extra_cols.append(gene_col)
                
                # スコア列を検索
                score_col = next((col for col in overlapping_se.columns if 'score' in col.lower()), None)
                if score_col:
                    extra_cols.append(score_col)
                
                # 必要に応じてその他の列も追加
                all_cols = bed_cols + extra_cols
                
                # BEDファイルとして出力
                overlapping_se[all_cols].to_csv(output_path, sep='\t', header=False, index=False)
            
            self.logger.info(f"BED file exported: {output_path}")
            return output_path
        
        except Exception as e:
            self.logger.exception(f"An error occurred during BED file export: {e}")
            return None
