# sedb_tools/reporting/report_generator.py

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import logging
import warnings
import traceback

class SEdbReportGenerator:
    """
    SEdbデータのHTMLレポート生成を担当するクラス
    """
    
    def __init__(self, analyzer):
        """
        初期化
        
        Parameters:
        -----------
        analyzer : SEdbRegionAnalyzer
            元となるアナライザークラスのインスタンス
        """
        self.analyzer = analyzer
        self.logger = analyzer.logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbReportGenerator")
    
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
    # ▲▲▲ 新規追加部分：初期化処理 ▲▲▲
    
    def generate_sedb_report(self, output_dir='results/sedb_report',
                          top_n=20, figure_formats=None, save_tables=True,
                          dpi=300, fig_width=12, fig_height=8):
        """
        SEdbメタデータの総合HTMLレポートを生成します
        
        Parameters:
        -----------
        output_dir : str, optional
            出力ディレクトリ（デフォルト: 'results/sedb_report'）
        top_n : int, optional
            可視化で表示するトップカテゴリの数（デフォルト: 20）
        figure_formats : list, optional
            保存する図のフォーマット（デフォルト: ['png', 'svg', 'pdf']）
        save_tables : bool, optional
            分布テーブルをCSVファイルとして保存するかどうか（デフォルト: True）
        dpi : int, optional
            ラスター画像のDPI（デフォルト: 300）
        fig_width : int, optional
            図の幅（インチ単位）（デフォルト: 12）
        fig_height : int, optional
            図の高さ（インチ単位）（デフォルト: 8）
            
        Returns:
        --------
        str
            生成されたHTMLレポートのパス
        """
        # 必要なデータの存在チェック
        # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
        if self.analyzer.bed_df is None or self.analyzer.sample_info is None:
            self.logger.warning("SE BED data or sample information not loaded. Call load_databases() first.")
            return None
        # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
        
        # Default figure formats if not provided
        if figure_formats is None:
            figure_formats = ['png', 'svg', 'pdf']
        
        # Create timestamp for filenames
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create directory structure with relative paths
        figures_dir = 'figures'  # Use relative paths for portable HTML
        tables_dir = 'tables'
        
        # Create full directory paths
        figures_full_path = os.path.join(output_dir, figures_dir)
        tables_full_path = os.path.join(output_dir, tables_dir)
        
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(figures_full_path, exist_ok=True)
        os.makedirs(tables_full_path, exist_ok=True)
        
        self.logger.info(f"Creating SEdb metadata report in directory: {output_dir}")
        
        # Dictionary to track all generated files
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Collect basic database stats
        # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
        db_stats = {
            'total_samples': len(self.analyzer.sample_info) if self.analyzer.sample_info is not None else 0,
            'total_se_records': len(self.analyzer.bed_df) if self.analyzer.bed_df is not None else 0,
            'bed_columns': list(self.analyzer.bed_df.columns) if self.analyzer.bed_df is not None else [],
            'sample_info_columns': list(self.analyzer.sample_info.columns) if self.analyzer.sample_info is not None else [],
            'timestamp': timestamp
        }
        
        # Calculate missing values in sample_info
        if self.analyzer.sample_info is not None:
            missing_values = self.analyzer.sample_info.isnull().sum()
            missing_values = missing_values[missing_values > 0]
            db_stats['missing_values'] = missing_values.to_dict()
        # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
        
        # Set Matplotlib params for better quality figures
        plt.rcParams['figure.figsize'] = (fig_width, fig_height)
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.dpi'] = dpi
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # Step 1: Tissue Distribution Analysis
        self.logger.info("Calculating and visualizing tissue type distribution...")
        
        tissue_dist_figures = []
        tissue_csv_path = None
        
        try:
            # 1. 直接保存オプションを使用してtissue distributionを生成
            tissue_viz_basename = f'tissue_distribution_direct_{timestamp}'
            
            # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
            # 重要な修正：直接ファイル保存オプションを使用
            self.analyzer.visualize_tissue_distribution(
                top_n=top_n, 
                horizontal=True, 
                store_all=True,
                figsize=(fig_width, fig_height),
                save_dir=figures_full_path,
                save_filename=tissue_viz_basename,
                save_png=True,
                save_svg=('svg' in figure_formats),
                save_pdf=('pdf' in figure_formats),
                save_eps=('eps' in figure_formats),
                save_dpi=dpi
            )
            # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
            
            # 保存されたファイルを登録
            for fmt in figure_formats:
                if fmt == 'png' or f'save_{fmt}' in locals():
                    figure_filename = f'{tissue_viz_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Tissue Distribution (Direct) ({fmt.upper()})', figure_relative_path))
                        
                        if fmt == 'png':
                            tissue_dist_figures.append(('Tissue Distribution (Direct)', figure_relative_path))
            
            # Now create supplementary visualizations using the stored distribution data
            # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
            if hasattr(self.analyzer, '_tissue_distribution') and self.analyzer._tissue_distribution is not None:
                # Save tissue data as CSV
                if save_tables:
                    tissue_table_filename = f'tissue_distribution_{timestamp}.csv'
                    tissue_table_path = os.path.join(tables_full_path, tissue_table_filename)
                    self.analyzer._tissue_distribution.to_csv(tissue_table_path, index=False)
                    generated_files['tables'].append(('Tissue Distribution', os.path.join(tables_dir, tissue_table_filename)))
                    tissue_csv_path = os.path.join(tables_dir, tissue_table_filename)
                
                # 1. Bar chart with sorted_bar helper (directly save)
                tissue_bar_basename = f'tissue_distribution_bar_{timestamp}'
                self.analyzer.plot_sorted_bar(
                    data=self.analyzer._tissue_distribution,
                    category='Tissue Type',
                    count_col='Sample Count',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='skyblue',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=tissue_bar_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Register the figure
                for fmt in figure_formats:
                    if fmt == 'png' or f'save_{fmt}' in locals():
                        figure_filename = f'{tissue_bar_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_dist_figures.append(('Tissue Distribution Bar Chart', figure_relative_path))
                
                # 2. Pie chart with simple_pie helper (directly save)
                tissue_pie_basename = f'tissue_distribution_pie_{timestamp}'
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                self.analyzer.plot_simple_pie(
                    data=self.analyzer._tissue_distribution,
                    category='Tissue Type',
                    count_col='Sample Count',
                    title=f'Tissue Type Distribution (Top {top_n})',
                    limit=top_n,
                    figsize=(fig_width, fig_height),
                    colormap='tab20',
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Register pie chart
                for fmt in figure_formats:
                    if fmt == 'png' or f'save_{fmt}' in locals():
                        figure_filename = f'{tissue_pie_basename}.{fmt}'
                        figure_relative_path = os.path.join(figures_dir, figure_filename)
                        figure_full_path = os.path.join(figures_full_path, figure_filename)
                        
                        if os.path.exists(figure_full_path):
                            generated_files['figures'].append((f'Tissue Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                            if fmt == 'png':
                                tissue_dist_figures.append(('Tissue Distribution Pie Chart', figure_relative_path))
                
                print(f"Successfully created tissue distribution visualizations")
            else:
                self.logger.warning("No tissue distribution data stored after visualization.")
                
        except Exception as e:
            self.logger.error(f"Error in tissue distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

# Step 2: Biosample Distribution Analysis
        self.logger.info("Calculating and visualizing biosample type distribution...")
        
        biosample_dist_figures = []
        biosample_csv_path = None
        
        try:
            # visualize_biosample_distributionメソッドが存在しないため、直接実装
            # バイオサンプルタイプの分布を計算
            # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
            if 'Biosample type' in self.analyzer.sample_info.columns:
                biosample_counts = self.analyzer.sample_info['Biosample type'].value_counts()
                total_samples = len(self.analyzer.sample_info)
            # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # 計算結果をDataFrameに変換
                biosample_data = pd.DataFrame({
                    'Biosample Type': biosample_counts.index,
                    'Sample Count': biosample_counts.values,
                    'Percentage (%)': (biosample_counts / total_samples * 100).round(2)
                })
                
                # CSV出力
                if save_tables:
                    biosample_table_filename = f'biosample_distribution_{timestamp}.csv'
                    biosample_table_path = os.path.join(tables_full_path, biosample_table_filename)
                    biosample_data.to_csv(biosample_table_path, index=False)
                    generated_files['tables'].append(('Biosample Distribution', os.path.join(tables_dir, biosample_table_filename)))
                    biosample_csv_path = os.path.join(tables_dir, biosample_table_filename)
                
                # バイオサンプルタイプの分布を直接可視化
                plt.figure(figsize=(fig_width, fig_height))
                
                # トップNを表示用に選択
                if top_n > 0 and len(biosample_data) > top_n:
                    display_data = biosample_data.head(top_n)
                else:
                    display_data = biosample_data
                
                # 水平バーチャート
                plt.barh(display_data['Biosample Type'].iloc[::-1], display_data['Sample Count'].iloc[::-1], color='lightgreen')
                plt.xlabel('Sample Count')
                plt.ylabel('Biosample Type')
                plt.title(f'Biosample Type Distribution (Top {len(display_data)} of {len(biosample_data)})')
                plt.tight_layout()
                
                # 分布図を保存
                biosample_viz_basename = f'biosample_distribution_direct_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{biosample_viz_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # 相対パスで登録
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Biosample Distribution (Direct) ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        biosample_dist_figures.append(('Biosample Distribution (Direct)', figure_relative_path))
                
                plt.close()
                
                # 追加の可視化 - バーチャート
                biosample_bar_basename = f'biosample_distribution_bar_{timestamp}'
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                self.analyzer.plot_sorted_bar(
                    data=biosample_data,
                    category='Biosample Type',
                    count_col='Sample Count',
                    title=f'Biosample Type Distribution (Top {top_n})',
                    limit=top_n,
                    horizontal=True,
                    color='lightgreen',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=biosample_bar_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Register the figure
                for fmt in figure_formats:
                    figure_filename = f'{biosample_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Biosample Distribution Bar Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            biosample_dist_figures.append(('Biosample Distribution Bar Chart', figure_relative_path))
                
                # 追加の可視化 - 円グラフ
                biosample_pie_basename = f'biosample_distribution_pie_{timestamp}'
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                self.analyzer.plot_simple_pie(
                    data=biosample_data,
                    category='Biosample Type',
                    count_col='Sample Count',
                    title=f'Biosample Type Distribution (Top {top_n})',
                    limit=top_n,
                    figsize=(fig_width, fig_height),
                    colormap='Set3',
                    save_dir=figures_full_path,
                    save_filename=biosample_pie_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Register pie chart
                for fmt in figure_formats:
                    figure_filename = f'{biosample_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Biosample Distribution Pie Chart ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            biosample_dist_figures.append(('Biosample Distribution Pie Chart', figure_relative_path))
                
                print(f"Successfully created biosample distribution visualizations")
            else:
                self.logger.warning("'Biosample type' column not found in sample information")
        except Exception as e:
            self.logger.error(f"Error in biosample distribution analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Step 3: Tissue-Biosample Relationship Analysis
        self.logger.info("Analyzing relationship between tissue types and biosample types...")
        
        relationship_figures = []
        relationship_csv_path = None
        
        try:
            # 手動でクロスタブとヒートマップを作成
            # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
            if 'Tissue type' in self.analyzer.sample_info.columns and 'Biosample type' in self.analyzer.sample_info.columns:
                # クロスタブを作成
                cross_tab = pd.crosstab(
                    self.analyzer.sample_info['Tissue type'], 
                    self.analyzer.sample_info['Biosample type']
                )
            # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # クロスタブデータを保存
                if save_tables:
                    cross_tab_filename = f'tissue_biosample_crosstab_{timestamp}.csv'
                    cross_tab_path = os.path.join(tables_full_path, cross_tab_filename)
                    cross_tab.to_csv(cross_tab_path)
                    generated_files['tables'].append(('Tissue-Biosample Cross-tabulation', os.path.join(tables_dir, cross_tab_filename)))
                    relationship_csv_path = os.path.join(tables_dir, cross_tab_filename)
                
                # 上位の組織とバイオサンプルタイプを取得
                row_sums = cross_tab.sum(axis=1)
                col_sums = cross_tab.sum(axis=0)
                
                row_order = row_sums.sort_values(ascending=False).index[:min(top_n, len(row_sums))]
                col_order = col_sums.sort_values(ascending=False).index[:min(top_n, len(col_sums))]
                
                # 表示用にデータをフィルタリング
                filtered_data = cross_tab.loc[row_order, col_order]
                
                # ヒートマップ図作成のためseabornをインポート
                import seaborn as sns
                
                # 新しい図を作成
                plt.figure(figsize=(fig_width, fig_height * 1.2))
                
                # ヒートマップ作成 - mask パラメータを使用して0の値をマスク
                cmap = sns.diverging_palette(220, 10, as_cmap=True)
                mask = filtered_data == 0  # 値が0のセルをマスク
                
                # ヒートマップを描画
                heatmap = sns.heatmap(
                    filtered_data, 
                    annot=True,  # セルに値を表示
                    fmt='d',     # 整数フォーマット
                    cmap=cmap,   # カラーマップ
                    linewidths=0.5, 
                    mask=mask,   # マスク適用
                    cbar_kws={'label': 'Sample Count'}
                )
                
                plt.title('Relationship Between Tissue Type and Biosample Type')
                plt.xlabel('Biosample Type')
                plt.ylabel('Tissue Type')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()  # レイアウトを調整
                
                # ヒートマップを保存
                heatmap_basename = f'tissue_biosample_heatmap_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{heatmap_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # 相対パスで登録
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Tissue-Biosample Heatmap ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        relationship_figures.append(('Tissue-Biosample Heatmap', figure_relative_path))
                
                plt.close()
                print(f"Successfully created tissue-biosample heatmap visualization")
            else:
                self.logger.warning("'Tissue type' or 'Biosample type' column not found in sample information")
        except Exception as e:
            self.logger.error(f"Error in tissue-biosample relationship analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())

# Step 4: Cell ID Occurrence Analysis (specific to SEdb)
        self.logger.info("Analyzing cell_id occurrences...")
        
        cell_id_figures = []
        cell_id_csv_path = None
        
        try:
            # Count cell_id occurrences if not already done
            # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
            if not hasattr(self.analyzer, '_cell_id_counts') or self.analyzer._cell_id_counts is None:
                self.analyzer.count_cell_id_occurrences()
            
            if hasattr(self.analyzer, '_cell_id_counts') and self.analyzer._cell_id_counts is not None:
                # Save cell_id counts data
                if save_tables:
                    cellid_table_filename = f'cell_id_counts_{timestamp}.csv'
                    cellid_table_path = os.path.join(tables_full_path, cellid_table_filename)
                    self.analyzer._cell_id_counts.to_csv(cellid_table_path, index=False)
                    generated_files['tables'].append(('Cell ID Counts', os.path.join(tables_dir, cellid_table_filename)))
                    cell_id_csv_path = os.path.join(tables_dir, cellid_table_filename)
                
                # Generate histogram of cell_id occurrences
                plt.figure(figsize=(fig_width, fig_height))
                
                # Create histogram
                plt.hist(self.analyzer._cell_id_counts['count'], bins=30, color='skyblue', edgecolor='black')
                plt.title('Distribution of cell_id Occurrences')
                plt.xlabel('Super Enhancers per Sample')
                plt.ylabel('Number of Samples')
                plt.grid(True, linestyle='--', alpha=0.7)
                plt.tight_layout()
            # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Save histogram
                cellid_hist_basename = f'cell_id_histogram_{timestamp}'
                
                for fmt in figure_formats:
                    figure_filename = f'{cellid_hist_basename}.{fmt}'
                    figure_path = os.path.join(figures_full_path, figure_filename)
                    plt.savefig(figure_path, format=fmt, dpi=dpi, bbox_inches='tight')
                    
                    # Register with relative path
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    generated_files['figures'].append((f'Cell ID Occurrence Histogram ({fmt.upper()})', figure_relative_path))
                    
                    if fmt == 'png':
                        cell_id_figures.append(('Cell ID Occurrence Histogram', figure_relative_path))
                
                plt.close()
                
                # Generate summary statistics for display in report
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                count_stats = self.analyzer._cell_id_counts['count'].describe().to_dict()
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                db_stats['cell_id_stats'] = count_stats
                
                # Display counts per sample using helper method
                cellid_bars_basename = f'cell_id_top_counts_{timestamp}'
                
                # Create a DataFrame of top samples by SE count
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                top_samples = self.analyzer._cell_id_counts.sort_values('count', ascending=False).head(top_n)
                
                # Create bar chart of top samples
                self.analyzer.plot_sorted_bar(
                    data=top_samples,
                    category='cell_id',
                    count_col='count',
                    title=f'Top {top_n} Samples by SE Count',
                    limit=top_n,
                    horizontal=True,
                    color='lightcoral',
                    figsize=(fig_width, fig_height),
                    save_dir=figures_full_path,
                    save_filename=cellid_bars_basename,
                    save_png=True,
                    save_svg=('svg' in figure_formats),
                    save_pdf=('pdf' in figure_formats),
                    save_eps=('eps' in figure_formats),
                    save_dpi=dpi
                )
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                
                # Register the figure
                for fmt in figure_formats:
                    figure_filename = f'{cellid_bars_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        generated_files['figures'].append((f'Top Samples by SE Count ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            cell_id_figures.append(('Top Samples by SE Count', figure_relative_path))
                
                print(f"Successfully created cell_id occurrence visualizations")
        except Exception as e:
            self.logger.error(f"Error in cell_id occurrence analysis: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
        
        # Generate HTML Report
        self.logger.info("Generating HTML report...")
        
        html_filename = f"sedb_metadata_report_{timestamp}.html"
        html_path = os.path.join(output_dir, html_filename)
        
        try:
            with open(html_path, 'w', encoding='utf-8') as f:
                # Start HTML document
                f.write("<!DOCTYPE html>\n")
                f.write("<html>\n")
                f.write("<head>\n")
                f.write("  <title>SEdb Metadata Report</title>\n")
                f.write("  <style>\n")
                f.write("    body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }\n")
                f.write("    h1 { color: #2c3e50; }\n")
                f.write("    h2 { color: #3498db; margin-top: 30px; border-bottom: 2px solid #3498db; padding-bottom: 5px; }\n")
                f.write("    h3 { color: #2980b9; margin-top: 20px; }\n")
                f.write("    h4 { color: #27ae60; margin-top: 15px; }\n")
                f.write("    table { border-collapse: collapse; width: 100%; margin: 20px 0; }\n")
                f.write("    th { background-color: #3498db; color: white; text-align: left; padding: 8px; }\n")
                f.write("    td { border: 1px solid #ddd; padding: 8px; }\n")
                f.write("    tr:nth-child(even) { background-color: #f2f2f2; }\n")
                f.write("    img { max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; box-shadow: 0 0 5px rgba(0,0,0,0.1); }\n")
                f.write("    .file-formats { margin: 10px 0; }\n")
                f.write("    .file-formats a { margin-right: 10px; color: #3498db; text-decoration: none; }\n")
                f.write("    .file-formats a:hover { text-decoration: underline; }\n")
                f.write("    .chart-section { display: flex; flex-wrap: wrap; justify-content: space-around; }\n")
                f.write("    .chart-container { margin: 15px; max-width: 48%; }\n")
                f.write("    .highlight { background-color: #f8f9fa; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }\n")
                f.write("    .explanation { background-color: #e8f4f8; padding: 15px; border-radius: 5px; margin: 10px 0; }\n")
                f.write("    .data-link { background-color: #f0f7ff; padding: 8px; margin: 5px 0; border-radius: 3px; }\n")
                f.write("    .data-link a { color: #0066cc; text-decoration: none; }\n")
                f.write("    .data-link a:hover { text-decoration: underline; }\n")
                f.write("  </style>\n")
                f.write("</head>\n")
                f.write("<body>\n")
                
                # Title and timestamp
                f.write(f"  <h1>SEdb Metadata Analysis Report</h1>\n")
                f.write(f"  <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>\n")
                
                # Database overview
                f.write("  <h2>1. Database Overview</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This section provides an overview of the SEdb database.</p>\n")
                f.write("  </div>\n")
                
                f.write("  <table>\n")
                f.write("    <tr><th>Metric</th><th>Value</th></tr>\n")
                f.write(f"    <tr><td>Total samples</td><td>{db_stats['total_samples']:,}</td></tr>\n")
                f.write(f"    <tr><td>Total SE records</td><td>{db_stats['total_se_records']:,}</td></tr>\n")
                
                # Add cell_id stats if available
                if 'cell_id_stats' in db_stats:
                    stats = db_stats['cell_id_stats']
                    f.write(f"    <tr><td>Average SE per sample</td><td>{stats['mean']:.2f}</td></tr>\n")
                    f.write(f"    <tr><td>Maximum SE per sample</td><td>{int(stats['max'])}</td></tr>\n")
                    f.write(f"    <tr><td>Minimum SE per sample</td><td>{int(stats['min'])}</td></tr>\n")
                
                f.write("  </table>\n")
                
                # Data schema section
                f.write("  <h3>Data Schema</h3>\n")
                
                # BED file columns
                f.write("  <h4>BED File Columns</h4>\n")
                f.write("  <ul>\n")
                for col in db_stats['bed_columns']:
                    f.write(f"    <li>{col}</li>\n")
                f.write("  </ul>\n")
                
                # Sample info columns
                f.write("  <h4>Sample Info Columns</h4>\n")
                f.write("  <ul>\n")
                for col in db_stats['sample_info_columns']:
                    f.write(f"    <li>{col}</li>\n")
                f.write("  </ul>\n")
                
                # Missing values section
                if 'missing_values' in db_stats and db_stats['missing_values']:
                    f.write("  <h3>Missing Values</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Column</th><th>Missing Count</th><th>Percentage</th></tr>\n")
                    
                    for col, count in db_stats['missing_values'].items():
                        percentage = (count / db_stats['total_samples']) * 100
                        f.write(f"    <tr><td>{col}</td><td>{count}</td><td>{percentage:.2f}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                
                # Tissue Distribution section
                f.write("  <h2>2. Tissue Type Distribution</h2>\n")
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This section shows the distribution of tissue types among samples in the SEdb database.</p>\n")
                f.write("    <p>The visualizations below help identify which tissues are most represented in the database.</p>\n")
                f.write("  </div>\n")
                
                # Direct visualization from visualize_tissue_distribution
                tissue_direct_png = [path for name, path in tissue_dist_figures if 'Direct' in name]
                if tissue_direct_png:
                    f.write("  <h3>Tissue Type Distribution (Direct Visualization)</h3>\n")
                    f.write(f"  <img src='{tissue_direct_png[0]}' alt='Tissue Distribution Direct Visualization'>\n")
                    
                    # Add data link for CSV if available
                    if tissue_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{tissue_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    tissue_direct_formats = [tup for tup in generated_files['figures'] if 'Tissue Distribution (Direct)' in tup[0]]
                    if len(tissue_direct_formats) > 1:
                        f.write("  <div class='file-formats'>Download in other formats: ")
                        for name, path in tissue_direct_formats:
                            if 'PNG' not in name:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[-1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Tissue distribution table
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                if hasattr(self.analyzer, '_tissue_distribution') and self.analyzer._tissue_distribution is not None:
                    # Table of top tissues
                    f.write(f"  <h3>Top {min(top_n, len(self.analyzer._tissue_distribution))} Tissue Types</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Tissue Type</th><th>Sample Count</th><th>Percentage</th></tr>\n")
                    
                    top_tissues = self.analyzer._tissue_distribution.head(top_n)
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                    for _, row in top_tissues.iterrows():
                        f.write(f"    <tr><td>{row['Tissue Type']}</td><td>{row['Sample Count']}</td><td>{row['Percentage (%)']}%</td></tr>\n")
                    
                    f.write("  </table>\n")
                    
                    # Additional visualizations
                    f.write("  <h3>Alternative Visualizations</h3>\n")
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    tissue_bar_png = [path for name, path in tissue_dist_figures if 'Bar Chart' in name]
                    if tissue_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{tissue_bar_png[0]}' alt='Tissue Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    tissue_pie_png = [path for name, path in tissue_dist_figures if 'Pie Chart' in name]
                    if tissue_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{tissue_pie_png[0]}' alt='Tissue Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    tissue_figures_all = [tup for tup in generated_files['figures'] if 'Tissue Distribution' in tup[0] and 'Direct' not in tup[0]]
                    if len(tissue_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        # Group by format
                        format_groups = {}
                        for name, path in tissue_figures_all:
                            format_name = name.split('(')[1].split(')')[0]
                            if format_name not in format_groups:
                                format_groups[format_name] = []
                            format_groups[format_name].append((name, path))
                        
                        # Add links by format
                        for format_name, items in sorted(format_groups.items()):
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Biosample Distribution section
                f.write("  <h2>3. Biosample Type Distribution</h2>\n")
                f.write("  <div class='explanation'>\n")
                f.write("    <p>This section shows the distribution of biosample types (cell lines, primary cells, tissues, etc.) among samples in the SEdb database.</p>\n")
                f.write("    <p>These visualizations help understand the types of biological samples represented in the database.</p>\n")
                f.write("  </div>\n")
                
                # Direct visualization from biosample distribution
                biosample_direct_png = [path for name, path in biosample_dist_figures if 'Direct' in name]
                if biosample_direct_png:
                    f.write("  <h3>Biosample Type Distribution (Direct Visualization)</h3>\n")
                    f.write(f"  <img src='{biosample_direct_png[0]}' alt='Biosample Distribution Direct Visualization'>\n")
                    
                    # Add data link for CSV if available
                    if biosample_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{biosample_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    biosample_direct_formats = [tup for tup in generated_files['figures'] if 'Biosample Distribution (Direct)' in tup[0]]
                    if len(biosample_direct_formats) > 1:
                        f.write("  <div class='file-formats'>Download in other formats: ")
                        for name, path in biosample_direct_formats:
                            if 'PNG' not in name:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[-1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Additional biosample visualizations
                other_biosample_png = [path for name, path in biosample_dist_figures if 'Direct' not in name]
                if other_biosample_png:
                    f.write("  <h3>Alternative Visualizations</h3>\n")
                    f.write("  <div class='chart-section'>\n")
                    
                    # Bar chart
                    biosample_bar_png = [path for name, path in biosample_dist_figures if 'Bar Chart' in name]
                    if biosample_bar_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Bar Chart</h4>\n")
                        f.write(f"      <img src='{biosample_bar_png[0]}' alt='Biosample Distribution Bar Chart'>\n")
                        f.write("    </div>\n")
                    
                    # Pie chart
                    biosample_pie_png = [path for name, path in biosample_dist_figures if 'Pie Chart' in name]
                    if biosample_pie_png:
                        f.write("    <div class='chart-container'>\n")
                        f.write("      <h4>Pie Chart</h4>\n")
                        f.write(f"      <img src='{biosample_pie_png[0]}' alt='Biosample Distribution Pie Chart'>\n")
                        f.write("    </div>\n")
                    
                    f.write("  </div>\n")
                    
                    # Add links to other formats
                    biosample_figures_all = [tup for tup in generated_files['figures'] if 'Biosample Distribution' in tup[0] and 'Direct' not in tup[0]]
                    if len(biosample_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        # Group by format
                        format_groups = {}
                        for name, path in biosample_figures_all:
                            format_name = name.split('(')[1].split(')')[0]
                            if format_name not in format_groups:
                                format_groups[format_name] = []
                            format_groups[format_name].append((name, path))
                        
                        # Add links by format
                        for format_name, items in sorted(format_groups.items()):
                            for name, path in items:
                                chart_type = "Bar Chart" if "Bar Chart" in name else "Pie Chart"
                                f.write(f"<a href='{path}'>{format_name} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Tissue-Biosample Relationship section
                if relationship_figures:
                    f.write("  <h2>4. Tissue-Biosample Relationship</h2>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>This section shows the relationship between tissue types and biosample types in the SEdb database.</p>\n")
                    f.write("    <p>The heatmap indicates the number of samples with both the corresponding tissue type and biosample type.</p>\n")
                    f.write("    <p>Areas with darker colors represent more samples with that specific tissue and biosample combination.</p>\n")
                    f.write("  </div>\n")
                    
                    # Find PNG relationship image
                    relationship_png = relationship_figures[0][1] if relationship_figures else None
                    if relationship_png:
                        f.write(f"  <img src='{relationship_png}' alt='Tissue-Biosample Relationship'>\n")
                        
                        # Add data link for CSV if available
                        if relationship_csv_path:
                            f.write(f"  <div class='data-link'>📊 <a href='{relationship_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Add links to other formats
                    relationship_figures_all = [tup for tup in generated_files['figures'] if ('Relationship' in tup[0] or 'Heatmap' in tup[0])]
                    if len(relationship_figures_all) > 1:
                        f.write("  <div class='file-formats'>Download figure as: ")
                        for name, path in relationship_figures_all:
                            if '.png' not in path:  # Skip PNG as it's already displayed
                                format_name = name.split('(')[1].split(')')[0]
                                f.write(f"<a href='{path}'>{format_name}</a> ")
                        f.write("  </div>\n")
                
                # Cell ID Occurrence Analysis section (SEdb specific)
                if cell_id_figures and 'cell_id_stats' in db_stats:
                    f.write("  <h2>5. Cell ID Occurrence Analysis</h2>\n")
                    f.write("  <div class='explanation'>\n")
                    f.write("    <p>This section analyzes the distribution of super enhancer occurrences per sample (cell_id).</p>\n")
                    f.write("    <p>It shows how many super enhancers (SEs) are found in each sample, helping to identify")
                    f.write("       samples with unusual SE counts that might be of special interest.</p>\n")
                    f.write("  </div>\n")
                    
                    # Display statistics
                    stats = db_stats['cell_id_stats']
                    f.write("  <h3>Statistics</h3>\n")
                    f.write("  <table>\n")
                    f.write("    <tr><th>Metric</th><th>Value</th></tr>\n")
                    f.write(f"    <tr><td>Mean</td><td>{stats['mean']:.2f} SEs per sample</td></tr>\n")
                    f.write(f"    <tr><td>Standard Deviation</td><td>{stats['std']:.2f}</td></tr>\n")
                    f.write(f"    <tr><td>Minimum</td><td>{int(stats['min'])} SEs</td></tr>\n")
                    f.write(f"    <tr><td>25th Percentile</td><td>{stats['25%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>Median (50th Percentile)</td><td>{stats['50%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>75th Percentile</td><td>{stats['75%']:.2f} SEs</td></tr>\n")
                    f.write(f"    <tr><td>Maximum</td><td>{int(stats['max'])} SEs</td></tr>\n")
                    f.write("  </table>\n")
                    
                    # Add data link for CSV if available
                    if cell_id_csv_path:
                        f.write(f"  <div class='data-link'>📊 <a href='{cell_id_csv_path}'>Download raw data (CSV)</a></div>\n")
                    
                    # Find PNG histogram image
                    histogram_png = [path for name, path in cell_id_figures if 'Histogram' in name]
                    if histogram_png:
                        f.write("  <h3>Distribution</h3>\n")
                        f.write(f"  <img src='{histogram_png[0]}' alt='Cell ID Occurrence Histogram'>\n")
                    
                    # Top samples bar chart
                    top_samples_png = [path for name, path in cell_id_figures if 'Top Samples' in name]
                    if top_samples_png:
                        f.write("  <h3>Top Samples by Super Enhancer Count</h3>\n")
                        f.write(f"  <img src='{top_samples_png[0]}' alt='Top Samples by SE Count'>\n")
                    
                    # Add links to other formats
                    if len(cell_id_figures) > 1:
                        f.write("  <div class='file-formats'>Download figures as: ")
                        for fmt in ['SVG', 'PDF', 'EPS']:  # Using common formats
                            format_files = [tup for tup in generated_files['figures'] if 'Cell ID' in tup[0] and fmt in tup[0]]
                            for name, path in format_files:
                                chart_type = "Histogram" if "Histogram" in name else "Bar Chart"
                                f.write(f"<a href='{path}'>{fmt} ({chart_type})</a> ")
                        f.write("  </div>\n")
                
                # Generated Files section
                f.write("  <h2>6. Generated Files</h2>\n")
                
                # List of all data files (tables)
                if generated_files['tables']:
                    f.write("  <h3>Data Files</h3>\n")
                    f.write("  <ul>\n")
                    for name, path in generated_files['tables']:
                        f.write(f"    <li><a href='{path}'>{name} (CSV)</a></li>\n")
                    f.write("  </ul>\n")
                
                # List of figures by category
                if generated_files['figures']:
                    f.write("  <h3>Figures</h3>\n")
                    
                    # Group figures by category
                    figure_groups = {}
                    for name, path in generated_files['figures']:
                        if 'Tissue' in name and 'Biosample' not in name:
                            key = 'Tissue Type Distribution'
                        elif 'Biosample' in name:
                            key = 'Biosample Type Distribution'
                        elif 'Relationship' in name or 'Heatmap' in name:
                            key = 'Tissue-Biosample Relationship'
                        elif 'Cell ID' in name or 'Samples by SE Count' in name:
                            key = 'Cell ID Analysis'
                        else:
                            key = 'Other Figures'
                        
                        if key not in figure_groups:
                            figure_groups[key] = []
                        figure_groups[key].append((name, path))
                    
                    # Display figures by group
                    for group_name, figures in figure_groups.items():
                        f.write(f"  <h4>{group_name}</h4>\n")
                        f.write("  <ul>\n")
                        for name, path in figures:
                            f.write(f"    <li><a href='{path}'>{name}</a></li>\n")
                        f.write("  </ul>\n")
                
                # Summary
                f.write("  <h2>7. Summary</h2>\n")
                f.write("  <div class='highlight'>\n")
                f.write("    <p>This report provides a comprehensive overview of the SEdb database metadata.</p>\n")
                f.write("    <p>Key findings include:</p>\n")
                f.write("    <ul>\n")
                f.write(f"      <li>The database contains {db_stats['total_samples']:,} samples with a total of {db_stats['total_se_records']:,} super enhancer records.</li>\n")
                
                # ▼▼▼ 改変：self -> self.analyzer ▼▼▼
                if hasattr(self.analyzer, '_tissue_distribution') and self.analyzer._tissue_distribution is not None and len(self.analyzer._tissue_distribution) > 0:
                    top_tissue = self.analyzer._tissue_distribution.iloc[0]
                # ▲▲▲ 改変：self -> self.analyzer ▲▲▲
                    f.write(f"      <li>The most common tissue type is '{top_tissue['Tissue Type']}' with {top_tissue['Sample Count']} samples ({top_tissue['Percentage (%)']}%).</li>\n")
                
                if 'cell_id_stats' in db_stats:
                    f.write(f"      <li>On average, each sample has {db_stats['cell_id_stats']['mean']:.2f} super enhancers.</li>\n")
                
                f.write("    </ul>\n")
                f.write("    <p>The visualizations show the distribution of samples across different biological categories using direct SEdb methods and additional representations.</p>\n")
                f.write("  </div>\n")
                
                # Footer
                f.write("  <hr>\n")
                f.write(f"  <p style='text-align:center; font-size:0.8em; color:#777;'>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by SEdbRegionAnalyzer</p>\n")
                
                # End HTML document
                f.write("</body>\n")
                f.write("</html>\n")
            
            self.logger.info(f"SEdb metadata report generated: {html_path}")
            return html_path
            
        except Exception as e:
            self.logger.error(f"Error generating HTML report: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None
        

    # def generate_region_report(self, analysis_results=None, chrom=None, start=None, end=None, 
    #                         region_name=None, output_dir='results/region_report',
    #                         figure_formats=None, save_tables=True, dpi=300, 
    #                         fig_width=12, fig_height=8, pvalue_threshold=0.05, 
    #                         fc_threshold=1.5):
    #     """
    #     Generate an HTML report from comprehensive region analysis results.
        
    #     [元のdocstringをコピー]
    #     """
    #     # 元のコードをほぼそのままコピー
    #     # self → self.analyzer に変更
        
    #     # 例:
    #     if analysis_results is None:
    #         if chrom is None or start is None or end is None:
    #             self.logger.error("Must provide either analysis_results or region coordinates (chrom, start, end)")
    #             return None
            
    #         self.logger.info(f"No analysis results provided. Running comprehensive analysis for {chrom}:{start}-{end}")
    #         analysis_results = self.analyzer.analyze_region_comprehensive(
    #             chrom=chrom, 
    #             start=start, 
    #             end=end, 
    #             region_name=region_name,
    #             output_dir=output_dir
    #         )
        
    #     # 以下、同様にコードを移植...