"""
Extended report generation features - Development code

This module contains development code for enhanced report generation.
Currently preserved for future development but not integrated into main application.

Note: This code contains experimental features and is not production-ready.
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime
import logging


class ExtendedReportGenerator:
    """
    Development version of enhanced report generator
    
    This code is preserved for future development but is not currently
    integrated into the main application.
    """
    
    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.logger = analyzer.logger if hasattr(analyzer, 'logger') else logging.getLogger(__name__)
    
    def generate_enhanced_report_dev(self, analysis_results=None, chrom=None, start=None, end=None,
                                   region_name=None, output_dir='results/region_report',
                                   figure_formats=None, save_tables=True, dpi=300, 
                                   fig_width=12, fig_height=8, pvalue_threshold=0.05, 
                                   fc_threshold=1.5):
        """
        DEVELOPMENT VERSION: Enhanced report with multiple visualizations
        
        This is development code preserved from the original extended report.
        Currently not integrated but preserved for future development.
        """
        

    # ver3
    def generate_enhanced_report_dev(self, analysis_results=None, chrom=None, start=None, end=None,
                                region_name=None, output_dir='results/region_report',
                                figure_formats=None, save_tables=True, dpi=300, 
                                fig_width=12, fig_height=8, pvalue_threshold=0.05, 
                                fc_threshold=1.5):
        """
        DEVELOPMENT VERSION: Enhanced report with multiple visualizations
        
        This is development code preserved from the original extended report.
        Currently not integrated but preserved for future development.
        """
        # 必要に応じて解析を実行
        if analysis_results is None:
            self.logger.info(f"レポート生成のための解析を実行: {chrom}:{start}-{end}")
            analysis_results = self.analyze_region_comprehensive(
                chrom=chrom, 
                start=start, 
                end=end, 
                region_name=region_name,
                output_dir=output_dir
            )
        
        # まず標準のレポートを生成
        report_path = self.generate_region_report(
            analysis_results=analysis_results,
            chrom=chrom,
            start=start,
            end=end,
            region_name=region_name,
            output_dir=output_dir,
            figure_formats=figure_formats,
            save_tables=save_tables,
            dpi=dpi,
            fig_width=fig_width,
            fig_height=fig_height,
            pvalue_threshold=pvalue_threshold,
            fc_threshold=fc_threshold
        )
        
        if report_path is None:
            self.logger.error("基本レポートの生成に失敗しました")
            # 結果辞書を作成して返す（エラー時）
            return {
                'report_path': None,
                'has_overlaps': False,
                'overlap_sample_count': 0,
                'overlap_read_count': 0,
                'significant_tissues': []
            }
        
        # 必要なパラメータを設定
        if figure_formats is None:
            figure_formats = ['png', 'svg', 'pdf']
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        figures_dir = 'figures'
        tables_dir = 'tables'
        figures_full_path = os.path.join(output_dir, figures_dir)
        tables_full_path = os.path.join(output_dir, tables_dir)
        os.makedirs(tables_full_path, exist_ok=True)
        
        # HTMLファイルを編集する前に読み込み
        with open(report_path, 'r', encoding='utf-8') as f:
            html_content = f.read()
        
        # 図の保存オプション
        save_opts = {}
        for fmt in figure_formats:
            save_opts[f'save_{fmt}'] = True
        
        # region_infoの取得
        if analysis_results:
            region_info = analysis_results.get('region_info', {})
        else:
            region_info = {}
        
        if not region_info and (chrom is None or start is None or end is None):
            if hasattr(self, 'region') and self.region:
                chrom = self.region.get('chr')
                start = self.region.get('start')
                end = self.region.get('end')
                region_name = self.region.get('name', f"{chrom}:{start}-{end}")
        else:
            chrom = region_info.get('chrom', chrom)
            start = region_info.get('start', start)
            end = region_info.get('end', end)
            region_name = region_info.get('name', region_name or f"{chrom}:{start}-{end}")
        
        # 新しいセクション用の画像と表を生成
        bar_figure_path = None
        pie_figure_path = None
        generated_files = {
            'figures': [],
            'tables': []
        }
        
        # Tissue Data取得
        tissue_data = None
        if hasattr(self, 'region_tissue_counts') and self.region_tissue_counts is not None:
            if isinstance(self.region_tissue_counts, pd.Series):
                # シリーズをDataFrameに変換
                tissue_data = pd.DataFrame({
                    'Tissue Type': self.region_tissue_counts.index,
                    'Region Count': self.region_tissue_counts.values
                })
                
                # Region Percentageを計算
                if hasattr(self, 'unique_sample_count'):
                    total_samples = self.unique_sample_count
                else:
                    total_samples = tissue_data['Region Count'].sum()
                    
                tissue_data['Region Percentage'] = (tissue_data['Region Count'] / total_samples * 100).round(2)
            else:
                tissue_data = self.region_tissue_counts
                
                # Percentageが存在しなければ追加
                if 'Region Percentage' not in tissue_data.columns and 'Region Count' in tissue_data.columns:
                    total_samples = tissue_data['Region Count'].sum()
                    tissue_data['Region Percentage'] = (tissue_data['Region Count'] / total_samples * 100).round(2)
        
        if tissue_data is not None:
            # 1. Bar chart
            tissue_bar_basename = f'region_{chrom}_{start}_{end}_tissue_bar_test_{timestamp}'
            try:
                self.logger.info(f"テスト用バーチャートを生成: {tissue_bar_basename}")
                tissue_bar_fig = self.plot_generator.plot_sorted_bar(
                    data=tissue_data,
                    category='Tissue Type',
                    count_col='Region Count',
                    title=f'Tissue Type Distribution in {region_name} (Test Bar)',
                    limit=15,
                    horizontal=True,
                    color='skyblue',
                    save_dir=figures_full_path,
                    save_filename=tissue_bar_basename,
                    **save_opts
                )
                
                # 生成された画像パスを保存
                for fmt in figure_formats:
                    figure_filename = f'{tissue_bar_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        self.logger.info(f"バーチャート ({fmt}) 生成確認: {figure_full_path}")
                        generated_files['figures'].append((f'Tissue Distribution Test Bar ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            bar_figure_path = figure_relative_path
            except Exception as e:
                self.logger.error(f"Test Tissue Bar Chart生成中にエラー発生: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
            
            # 2. Pie chart
            tissue_pie_basename = f'region_{chrom}_{start}_{end}_tissue_pie_test_{timestamp}'
            try:
                self.logger.info(f"テスト用パイチャートを生成: {tissue_pie_basename}")
                tissue_pie_fig = self.plot_generator.plot_simple_pie(
                    data=tissue_data,
                    category='Tissue Type',
                    count_col='Region Count',  # countカラムを明示的に指定
                    title=f'Tissue Type Distribution in {region_name} (Test Pie)',
                    limit=10,
                    save_dir=figures_full_path,
                    save_filename=tissue_pie_basename,
                    **save_opts
                )
                
                # 生成された画像パスを保存
                for fmt in figure_formats:
                    figure_filename = f'{tissue_pie_basename}.{fmt}'
                    figure_relative_path = os.path.join(figures_dir, figure_filename)
                    figure_full_path = os.path.join(figures_full_path, figure_filename)
                    
                    if os.path.exists(figure_full_path):
                        self.logger.info(f"パイチャート ({fmt}) 生成確認: {figure_full_path}")
                        generated_files['figures'].append((f'Tissue Distribution Test Pie ({fmt.upper()})', figure_relative_path))
                        if fmt == 'png':
                            pie_figure_path = figure_relative_path
            except Exception as e:
                self.logger.error(f"Test Tissue Pie Chart生成中にエラー発生: {e}")
                import traceback
                self.logger.error(traceback.format_exc())
                
            # 3. データテーブルをCSVとして保存
            if save_tables:
                tissue_table_test_filename = f'region_{chrom}_{start}_{end}_tissue_test_dist_{timestamp}.csv'
                tissue_table_test_path = os.path.join(tables_full_path, tissue_table_test_filename)
                tissue_data.to_csv(tissue_table_test_path, index=False)
                relative_table_path = os.path.join(tables_dir, tissue_table_test_filename)
                generated_files['tables'].append(('Tissue Test Distribution', relative_table_path))
        
        # 3.1セクションの終わりを見つけて新しいセクションを挿入
        section_3_1_marker = """  <h3>3.1 Tissue Type Distribution</h3>"""
        section_3_2_marker = """  <h3>3.2 Biosample Type Distribution</h3>"""
        
        section_3_1b_html = """
        <h3>3.1B(test)</h3>
        <div class='explanation'>
            <p>このセクションはテスト用に追加されたもので、eRNA分析レポートと同様のグラフ表示を実装しています。</p>
            <p>以下のグラフはSEdbデータに基づく組織分布を表示しています。</p>
        </div>
        """
        
        # 組織分布テーブルをHTMLに直接表示
        if tissue_data is not None:
            # 上位10の組織をテーブルで表示
            top_tissues = tissue_data.sort_values('Region Count', ascending=False).head(10)
            
            section_3_1b_html += """  <h4>Top Tissue Types in Overlapping Samples (Test)</h4>\n"""
            section_3_1b_html += """  <table>\n"""
            section_3_1b_html += """    <tr><th>Tissue Type</th><th>Count</th><th>Percentage</th></tr>\n"""
            
            for _, row in top_tissues.iterrows():
                tissue = row['Tissue Type']
                count = row['Region Count']
                percentage = row['Region Percentage'] if 'Region Percentage' in row else (count / tissue_data['Region Count'].sum() * 100)
                section_3_1b_html += f"""    <tr><td>{tissue}</td><td>{count}</td><td>{percentage:.2f}%</td></tr>\n"""
                
            section_3_1b_html += """  </table>\n"""
            section_3_1b_html += """  <p class='notes' style='margin-bottom:20px;'>このテストテーブルはSEdbデータから直接生成されています。より詳細なデータはCSVファイルでダウンロードできます。</p>\n"""
        
        # チャートセクションを明示的に作成（パイチャートとバーチャート両方）
        section_3_1b_html += """  <div class='chart-section'>\n"""
        
        # バーチャートを追加
        if bar_figure_path:
            section_3_1b_html += """    <div class='chart-container'>\n"""
            section_3_1b_html += """      <h4>Tissue Bar Chart (Test)</h4>\n"""
            section_3_1b_html += f"""      <img src='{bar_figure_path}' alt='Tissue Bar Chart (Test)'>\n"""
            section_3_1b_html += """    </div>\n"""
        
        # パイチャートを追加
        if pie_figure_path:
            section_3_1b_html += """    <div class='chart-container'>\n"""
            section_3_1b_html += """      <h4>Tissue Pie Chart (Test)</h4>\n"""
            section_3_1b_html += f"""      <img src='{pie_figure_path}' alt='Tissue Pie Chart (Test)'>\n"""
            section_3_1b_html += """    </div>\n"""
        
        section_3_1b_html += """  </div>\n"""
        

        # ダウンロードリンクを明示的に追加（図とCSVファイルを同じ行に）
        section_3_1b_html += """  <div class='file-formats'>Download test figures as: """

        # 各形式のリンクを追加
        for fmt in figure_formats:
            # バーチャートリンク
            bar_filename = f'{tissue_bar_basename}.{fmt}'
            bar_path = os.path.join(figures_dir, bar_filename)
            bar_full_path = os.path.join(figures_full_path, bar_filename)
            if os.path.exists(bar_full_path):
                section_3_1b_html += f"""<a href='{bar_path}'>Bar ({fmt.upper()})</a> """
            
            # パイチャートリンク
            pie_filename = f'{tissue_pie_basename}.{fmt}'
            pie_path = os.path.join(figures_dir, pie_filename)
            pie_full_path = os.path.join(figures_full_path, pie_filename)
            if os.path.exists(pie_full_path):
                section_3_1b_html += f"""<a href='{pie_path}'>Pie ({fmt.upper()})</a> """

        # CSVファイルへのリンクを同じ行に追加
        if save_tables and 'tables' in generated_files and generated_files['tables']:
            # 区切りを追加
            section_3_1b_html += """ | Data: """
            
            # CSVファイルへのリンク
            for table_name, table_path in generated_files['tables']:
                section_3_1b_html += f"""<a href='{table_path}'>CSV</a> """

        section_3_1b_html += """  </div>\n"""


        # 3.1と3.2セクションの間に新しいセクションを挿入
        if section_3_1_marker in html_content and section_3_2_marker in html_content:
            split_pos = html_content.find(section_3_2_marker)
            
            if split_pos > 0:
                new_html_content = (
                    html_content[:split_pos] + 
                    section_3_1b_html + 
                    html_content[split_pos:]
                )
                
                # 変更されたHTMLを保存
                with open(report_path, 'w', encoding='utf-8') as f:
                    f.write(new_html_content)
                
                self.logger.info(f"拡張レポートが生成されました: {report_path} (テストセクション3.1B追加済み)")
                
                # 生成確認
                figure_counts = {}
                for fmt in figure_formats:
                    bar_path = os.path.join(figures_full_path, f'{tissue_bar_basename}.{fmt}')
                    pie_path = os.path.join(figures_full_path, f'{tissue_pie_basename}.{fmt}')
                    figure_counts[f'bar_{fmt}'] = os.path.exists(bar_path)
                    figure_counts[f'pie_{fmt}'] = os.path.exists(pie_path)
                
                self.logger.info(f"図ファイル生成確認: {figure_counts}")
            else:
                self.logger.warning("次のセクションが見つからないため、テストセクションの追加に失敗しました")
        else:
            self.logger.warning("必要なセクションマーカーが見つからないため、テストセクションの追加に失敗しました")
        
        # 戻り値を作成
        result_dict = {
            'report_path': report_path,
            'analysis_results': analysis_results,  # analysis_resultsをそのまま含める
            'has_overlaps': False,
            'overlap_sample_count': 0,
            'overlap_read_count': 0,
            'significant_tissues': []
        }

        # analyze_and_report_regionと同様に必要な情報を抽出
        # 重複する領域があるかどうか
        has_overlaps = False
        overlap_sample_count = 0
        overlap_read_count = 0
        significant_tissues = []
        
        # overlapping_samples_countから抽出
        if analysis_results and 'overlapping_samples_count' in analysis_results:
            overlap_sample_count = analysis_results['overlapping_samples_count']
            has_overlaps = overlap_sample_count > 0
            result_dict['has_overlaps'] = has_overlaps
            result_dict['overlap_sample_count'] = overlap_sample_count
        
        # overlapping_resultsからリード数を抽出
        if analysis_results and 'overlapping_results' in analysis_results and analysis_results['overlapping_results'] is not None:
            overlap_read_count = len(analysis_results['overlapping_results'])
            result_dict['overlap_read_count'] = overlap_read_count
        
        # 有意な組織のリスト
        if analysis_results and 'tissue_enrichment' in analysis_results and analysis_results['tissue_enrichment'] is not None:
            tissue_df = analysis_results['tissue_enrichment']
            # pまたはFDR値列を特定
            fdr_col = 'adjusted p-value' if 'adjusted p-value' in tissue_df.columns else 'p-value'
            
            # 有意な組織のみをフィルター
            sig_tissues = tissue_df[(tissue_df[fdr_col] < pvalue_threshold) & 
                                (tissue_df['Fold Change'] > fc_threshold)]
            
            # 有意性でソート
            sig_tissues = sig_tissues.sort_values(fdr_col)
            
            # 結果を辞書のリストに変換
            for _, row in sig_tissues.iterrows():
                tissue_info = {
                    'tissue': row['Tissue Type'],
                    'fold_change': row['Fold Change'],
                    'pvalue': row['p-value'],
                    'fdr': row[fdr_col] if fdr_col in row else row['p-value']
                }
                significant_tissues.append(tissue_info)
        
        result_dict['significant_tissues'] = significant_tissues
        
        # 代替方法: インスタンス変数から情報を直接取得
        if not has_overlaps and hasattr(self, 'unique_sample_count'):
            result_dict['overlap_sample_count'] = self.unique_sample_count
            result_dict['has_overlaps'] = self.unique_sample_count > 0
        
        if overlap_read_count == 0 and hasattr(self, 'overlapping_se') and self.overlapping_se is not None:
            result_dict['overlap_read_count'] = len(self.overlapping_se)
        
        # 最終ログ出力
        self.logger.info(f"戻り値の確認: has_overlaps={result_dict['has_overlaps']}, "
                        f"sample_count={result_dict['overlap_sample_count']}, "
                        f"read_count={result_dict['overlap_read_count']}, "
                        f"significant_tissues_count={len(result_dict['significant_tissues'])}")
        
        print(report_path)
        return result_dict
