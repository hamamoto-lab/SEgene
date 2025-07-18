# visualization/distribution_plots.py
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import os
from typing import Optional, Dict, List, Tuple, Any
import logging
import seaborn as sns

class SEdbDistributionPlotter:
    """
    分布の可視化に特化したクラス。
    様々な種類の分布プロットを生成するための機能を提供します。
    """
    
    def __init__(self, logger=None):
        """
        分布プロッターを初期化
        
        Parameters:
        logger (Logger): オプションのロガーインスタンス。Noneの場合、基本的なコンソールロガーを作成。
        """
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbDistributionPlotter")

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
    


    def visualize_tissue_distribution(self, sample_info, top_n=20, horizontal=True, store_all=True, figsize=(10, 6),
                                    save_dir=None, save_filename=None, 
                                    save_png=False, save_svg=False, save_pdf=False, save_eps=False,
                                    save_dpi=300, save_transparent=False):
        """
        Calculate and visualize tissue type distribution.
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            サンプル情報を含むデータフレーム
        top_n : int, optional
            表示する上位組織タイプの数（デフォルト: 20）
        horizontal : bool, optional
            水平バーチャートを使用するか（デフォルト: True）
        store_all : bool, optional
            全ての分布データを保存するか（デフォルト: True）
        figsize : tuple, optional
            図のサイズ（デフォルト: (10, 6)）
        save_dir : str, optional
            図を保存するディレクトリ（指定しない場合は保存しない）
        save_filename : str, optional
            保存するファイル名（拡張子なし）
        save_png, save_svg, save_pdf, save_eps : bool, optional
            各フォーマットで保存するかどうか（デフォルト: False）
        save_dpi : int, optional
            保存時の解像度（デフォルト: 300）
        save_transparent : bool, optional
            透明な背景で保存するかどうか（デフォルト: False）
            
        Returns:
        --------
        pandas.DataFrame
            組織タイプの分布データ
        """
        if sample_info is None:
            self.logger.warning("Sample information is not provided. Cannot calculate tissue distribution.")
            return None
        
#        print("組織型分布を計算...")
        
        # 組織タイプのカラムを確認
        tissue_col = 'Tissue type' if 'Tissue type' in sample_info.columns else None
        
        if tissue_col is None:
            self.logger.warning("Tissue type column not found in sample information.")
            return None
        
        # 組織タイプの分布を計算
        tissue_counts = sample_info[tissue_col].value_counts()
        total_samples = len(sample_info)
        
        # 計算結果をDataFrameに変換
        tissue_distribution = pd.DataFrame({
            'Tissue Type': tissue_counts.index,
            'Sample Count': tissue_counts.values,
            'Percentage (%)': (tissue_counts / total_samples * 100).round(2)
        })
        
        # トップNを表示用に選択
        if top_n > 0 and len(tissue_distribution) > top_n:
            display_data = tissue_distribution.head(top_n)
        else:
            display_data = tissue_distribution
        
        # 可視化用のフィギュアを作成
        plt.figure(figsize=figsize)
        
        if horizontal:
            # 水平バーチャート（長い名前のリストに適している）
            plt.barh(display_data['Tissue Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='skyblue')
            plt.xlabel('Sample Count')
            plt.ylabel('Tissue Type')
        else:
            # 垂直バーチャート
            plt.bar(display_data['Tissue Type'], display_data['Sample Count'], color='skyblue')
            plt.xlabel('Tissue Type')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
        
        plt.title(f'Tissue Type Distribution (Top {len(display_data)} of {len(tissue_distribution)})')
        plt.tight_layout()
        
        # 図を保存（オプション）
        if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
            import os
            from datetime import datetime
            
            # ディレクトリが存在しなければ作成
            os.makedirs(save_dir, exist_ok=True)
            
            # タイムスタンプを生成（ファイル名がない場合）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # ファイル名がなければデフォルト名を生成
            if save_filename is None:
                save_filename = f'tissue_distribution_{timestamp}'
                
            # 指定された各フォーマットで保存
            saved_files = []
            
            if save_png:
                png_path = os.path.join(save_dir, f'{save_filename}.png')
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_svg:
                svg_path = os.path.join(save_dir, f'{save_filename}.svg')
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_pdf:
                pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = os.path.join(save_dir, f'{save_filename}.eps')
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # 保存されたファイルをログに記録
            if saved_files:
                self.logger.info(f"Tissue distribution visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return tissue_distribution
    


    def visualize_biosample_distribution(self, sample_info, top_n=20, horizontal=True, figsize=(10, 6),
                                        save_dir=None, save_filename=None, 
                                        save_png=False, save_svg=False, save_pdf=False, save_eps=False,
                                        save_dpi=300, save_transparent=False):
        """
        Calculate and visualize biosample type distribution.
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            サンプル情報を含むデータフレーム
        top_n : int, optional
            表示する上位バイオサンプルタイプの数（デフォルト: 20）
        horizontal : bool, optional
            水平バーチャートを使用するか（デフォルト: True）
        figsize : tuple, optional
            図のサイズ（デフォルト: (10, 6)）
        save_dir : str, optional
            図を保存するディレクトリ（指定しない場合は保存しない）
        save_filename : str, optional
            保存するファイル名（拡張子なし）
        save_png, save_svg, save_pdf, save_eps : bool, optional
            各フォーマットで保存するかどうか（デフォルト: False）
        save_dpi : int, optional
            保存時の解像度（デフォルト: 300）
        save_transparent : bool, optional
            透明な背景で保存するかどうか（デフォルト: False）
            
        Returns:
        --------
        pandas.DataFrame
            バイオサンプルタイプの分布データ
        """
        if sample_info is None:
            self.logger.warning("Sample information is not provided. Cannot calculate biosample distribution.")
            return None
        
        # バイオサンプルタイプのカラムを確認
        biosample_col = 'Biosample type' if 'Biosample type' in sample_info.columns else None
        
        if biosample_col is None:
            self.logger.warning("Biosample type column not found in sample information.")
            return None
        
        # バイオサンプルタイプの分布を計算
        biosample_counts = sample_info[biosample_col].value_counts()
        total_samples = len(sample_info)
        
        # 計算結果をDataFrameに変換
        biosample_distribution = pd.DataFrame({
            'Biosample Type': biosample_counts.index,
            'Sample Count': biosample_counts.values,
            'Percentage (%)': (biosample_counts / total_samples * 100).round(2)
        })
        
        # トップNを表示用に選択
        if top_n > 0 and len(biosample_distribution) > top_n:
            display_data = biosample_distribution.head(top_n)
        else:
            display_data = biosample_distribution
        
        # 可視化用のフィギュアを作成
        plt.figure(figsize=figsize)
        
        if horizontal:
            # 水平バーチャート（長い名前のリストに適している）
            plt.barh(display_data['Biosample Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='lightgreen')
            plt.xlabel('Sample Count')
            plt.ylabel('Biosample Type')
        else:
            # 垂直バーチャート
            plt.bar(display_data['Biosample Type'], display_data['Sample Count'], color='lightgreen')
            plt.xlabel('Biosample Type')
            plt.ylabel('Sample Count')
            plt.xticks(rotation=45, ha='right')
        
        plt.title(f'Biosample Type Distribution (Top {len(display_data)} of {len(biosample_distribution)})')
        plt.tight_layout()
        
        # 図を保存（オプション）
        if save_dir is not None and (save_png or save_svg or save_pdf or save_eps):
            import os
            from datetime import datetime
            
            # ディレクトリが存在しなければ作成
            os.makedirs(save_dir, exist_ok=True)
            
            # タイムスタンプを生成（ファイル名がない場合）
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # ファイル名がなければデフォルト名を生成
            if save_filename is None:
                save_filename = f'biosample_distribution_{timestamp}'
                
            # 指定された各フォーマットで保存
            saved_files = []
            
            if save_png:
                png_path = os.path.join(save_dir, f'{save_filename}.png')
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_svg:
                svg_path = os.path.join(save_dir, f'{save_filename}.svg')
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_pdf:
                pdf_path = os.path.join(save_dir, f'{save_filename}.pdf')
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = os.path.join(save_dir, f'{save_filename}.eps')
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight', transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # 保存されたファイルをログに記録
            if saved_files:
                self.logger.info(f"Biosample distribution visualization saved in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return biosample_distribution
    

    def visualize_tissue_biosample_relationship(self, sample_info, figsize=(15, 12), min_count=5):
        """
        組織タイプとバイオサンプルタイプの関係を可視化します。
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            サンプル情報を含むデータフレーム
        figsize : tuple, optional
            図のサイズ
        min_count : int, optional
            ヒートマップに含める最小カウント数
        
        Returns:
        --------
        pandas.DataFrame
            組織タイプとバイオサンプルタイプのクロスタブ
        """
        if sample_info is None:
            self.logger.warning("Sample information is not provided. Cannot analyze relationship.")
            return None
        
        if 'Tissue type' not in sample_info.columns or 'Biosample type' not in sample_info.columns:
            self.logger.warning("'Tissue type' or 'Biosample type' column not found in sample information.")
            return None
        
        # クロスタブを作成
        cross_tab = pd.crosstab(sample_info['Tissue type'], sample_info['Biosample type'])
        
        # 行と列の合計を追加
        cross_tab['Total'] = cross_tab.sum(axis=1)
        cross_tab.loc['Total'] = cross_tab.sum()
        
        # ヒートマップ用データの準備（合計行と列を除外）
        heatmap_data = cross_tab.iloc[:-1, :-1].copy()
        
        # 行と列を合計値で降順に並べ替え
        row_order = heatmap_data.sum(axis=1).sort_values(ascending=False).index
        col_order = heatmap_data.sum(axis=0).sort_values(ascending=False).index
        
        heatmap_data = heatmap_data.loc[row_order, col_order]
        
        # 少数のサンプルを持つ組織タイプとバイオサンプルタイプをフィルタリング
        if min_count > 0:
            row_mask = heatmap_data.sum(axis=1) >= min_count
            col_mask = heatmap_data.sum(axis=0) >= min_count
            filtered_data = heatmap_data.loc[row_mask, col_mask]
            
            if filtered_data.empty:
                self.logger.warning(f"No tissue type and biosample type combinations meet the min_count={min_count} criteria.")
                filtered_data = heatmap_data
        else:
            filtered_data = heatmap_data
        
        # 可視化
        plt.figure(figsize=figsize)
        
        # カラーマップ（0の値は灰色で表示）
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        mask = filtered_data == 0
        
        # ヒートマップ
        sns.heatmap(filtered_data, annot=True, fmt='d', cmap=cmap, mask=mask,
                    linewidths=0.5, cbar_kws={'label': 'Sample Count'})
        
        plt.title('Relationship Between Tissue Type and Biosample Type')
        plt.xlabel('Biosample Type')
        plt.ylabel('Tissue Type')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        return cross_tab
    


    def visualize_data_sources(self, sample_info, figsize=(12, 8)):
        """
        データソースの分布を可視化します。
        
        Parameters:
        -----------
        sample_info : pandas.DataFrame
            サンプル情報を含むデータフレーム
        figsize : tuple, optional
            図のサイズ
        
        Returns:
        --------
        tuple
            (matplotlib.figure.Figure, データソース集計情報)
        """
        if sample_info is None:
            self.logger.warning("Sample information is not provided. Cannot analyze data sources.")
            return None, None
        
        if 'Data source' not in sample_info.columns:
            self.logger.warning("'Data source' column not found in sample information.")
            return None, None
        
        # データソースごとのサンプル数を計算
        source_counts = sample_info['Data source'].value_counts()
        
        # 組織タイプとデータソースのクロスタブ
        tissue_source = None
        if 'Tissue type' in sample_info.columns:
            tissue_source = pd.crosstab(sample_info['Tissue type'], sample_info['Data source'])
            tissue_source['Total'] = tissue_source.sum(axis=1)
            tissue_source = tissue_source.sort_values('Total', ascending=False)
        
        # バイオサンプルタイプとデータソースのクロスタブ
        biosample_source = None
        if 'Biosample type' in sample_info.columns:
            biosample_source = pd.crosstab(sample_info['Biosample type'], sample_info['Data source'])
            biosample_source['Total'] = biosample_source.sum(axis=1)
            biosample_source = biosample_source.sort_values('Total', ascending=False)
        
        # 可視化
        fig, axs = plt.subplots(2, 2, figsize=figsize)
        
        # 1. データソースごとのサンプル数
        bars = axs[0, 0].bar(source_counts.index, source_counts.values, 
                    color=plt.cm.Set3(np.arange(len(source_counts))))
        axs[0, 0].set_title('Sample Count by Data Source')
        axs[0, 0].set_ylabel('Sample Count')
        axs[0, 0].tick_params(axis='x', rotation=45, labelsize=8)
        
        # バーの上に数値を表示
        for bar in bars:
            height = bar.get_height()
            axs[0, 0].text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{height}', ha='center', va='bottom')
        
        # 2. データソースの割合（円グラフ）
        axs[0, 1].pie(source_counts.values, labels=source_counts.index, autopct='%1.1f%%', 
                    colors=plt.cm.Set3(np.arange(len(source_counts))))
        axs[0, 1].set_title('Percentage by Data Source')
        axs[0, 1].axis('equal')
        
        # 3. 組織タイプ別のデータソース内訳
        if tissue_source is not None:
            # 上位10組織タイプのみ表示
            top_tissues = tissue_source.iloc[:10, :-1]  # 'Total'列を除外
            top_tissues.plot(kind='bar', stacked=True, ax=axs[1, 0])
            axs[1, 0].set_title('Data Source by Major Tissue Types')
            axs[1, 0].set_ylabel('Sample Count')
            axs[1, 0].set_xlabel('Tissue Type')
            axs[1, 0].tick_params(axis='x', rotation=45, labelsize=8)
            axs[1, 0].legend(loc='upper right', fontsize='small')
        
        # 4. バイオサンプルタイプ別のデータソース内訳
        if biosample_source is not None:
            # 上位10バイオサンプルタイプのみ表示
            top_biosamples = biosample_source.iloc[:10, :-1]  # 'Total'列を除外
            top_biosamples.plot(kind='bar', stacked=True, ax=axs[1, 1])
            axs[1, 1].set_title('Data Source by Major Biosample Types')
            axs[1, 1].set_ylabel('Sample Count')
            axs[1, 1].set_xlabel('Biosample Type')
            axs[1, 1].tick_params(axis='x', rotation=45, labelsize=8)
            axs[1, 1].legend(loc='upper right', fontsize='small')
        
        plt.tight_layout()
        
        # 集計情報を構築
        summary = {
            'source_counts': source_counts,
            'tissue_source': tissue_source,
            'biosample_source': biosample_source
        }
        
        return fig, summary
    


    def visualize_cell_id_distribution(self, cell_id_counts, top_n=10, bottom_n=10, 
                                    show_histogram=True, bins=30, figsize=(15, 10)):
        """
        cell_idの分布を可視化します。
        
        Parameters:
        -----------
        cell_id_counts : pandas.DataFrame
            cell_idとそのカウント情報を含むデータフレーム
        top_n : int, optional
            表示する最頻出cell_idの数
        bottom_n : int, optional
            表示する最少出現cell_idの数
        show_histogram : bool, optional
            ヒストグラムを表示するかどうか
        bins : int, optional
            ヒストグラムのビン数
        figsize : tuple, optional
            図のサイズ
        
        Returns:
        --------
        tuple
            (matplotlib.figure.Figure, 集計統計情報)
        """
        if cell_id_counts is None or len(cell_id_counts) == 0:
            self.logger.warning("cell_id counts not available.")
            return None, None
        
        counts_df = cell_id_counts
        
        # 基本統計情報を計算
        stats = counts_df['count'].describe().to_dict()
        
        # ゼロカウントのcell_id数
        zero_counts = (counts_df['count'] == 0).sum()
        
        # 最頻出cell_ids
        top_counts = counts_df.nlargest(top_n, 'count')
        
        # 最低頻出cell_ids（ゼロを除く）
        non_zero_df = counts_df[counts_df['count'] > 0]
        if len(non_zero_df) > 0:
            bottom_counts = non_zero_df.nsmallest(bottom_n, 'count')
        else:
            bottom_counts = pd.DataFrame()
        
        # 可視化
        if show_histogram:
            fig, axs = plt.subplots(2, 2, figsize=figsize)
            
            # 1. ヒストグラム
            axs[0, 0].hist(counts_df['count'], bins=bins, edgecolor='black', color='skyblue')
            axs[0, 0].set_title('Histogram of cell_id Occurrences')
            axs[0, 0].set_xlabel('Occurrence Count')
            axs[0, 0].set_ylabel('Number of cell_ids')
            
            # 2. 対数スケールのヒストグラム
            axs[0, 1].hist(counts_df['count'], bins=bins, edgecolor='black', log=True, color='skyblue')
            axs[0, 1].set_title('Histogram of cell_id Occurrences (Log Scale)')
            axs[0, 1].set_xlabel('Occurrence Count')
            axs[0, 1].set_ylabel('Number of cell_ids (log)')
            
            # 3. ボックスプロット
            axs[1, 0].boxplot(counts_df['count'])
            axs[1, 0].set_title('Box Plot of cell_id Occurrences')
            axs[1, 0].set_ylabel('Occurrence Count')
            
            # 4. 累積分布
            counts_sorted = np.sort(counts_df['count'])
            p = 1. * np.arange(len(counts_df)) / (len(counts_df) - 1)
            axs[1, 1].plot(counts_sorted, p)
            axs[1, 1].set_title('Cumulative Distribution of cell_id Occurrences')
            axs[1, 1].set_xlabel('Occurrence Count')
            axs[1, 1].set_ylabel('Cumulative Probability')
            
            plt.tight_layout()
        else:
            fig = None
        
        # 集計情報をまとめる
        summary = {
            'stats': stats,
            'zero_counts': zero_counts,
            'zero_percent': zero_counts/len(counts_df)*100 if len(counts_df) > 0 else 0,
            'top_counts': top_counts,
            'bottom_counts': bottom_counts
        }
        
        return fig, summary