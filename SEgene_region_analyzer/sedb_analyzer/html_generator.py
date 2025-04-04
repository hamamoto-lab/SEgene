"""
sedb_analyzer.html_generator - HTML Report Generation Module

This module provides functionality for generating HTML reports of SEdb analysis results.
"""

import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple

def get_safe_relative_path(file_path, base_dir=None):
    """
    Safely converts a file path to a relative path.

    Parameters
    ----------
    file_path : str
        The file path to be converted
    base_dir : str or Path, optional
        Base directory for the relative path. If None, returns the path as is.

    Returns
    -------
    str
        Relative path or the original path string
    """
    if base_dir is None:
        return str(file_path)
    
    try:
        # パスオブジェクトに変換
        p_file = Path(file_path)
        p_base = Path(base_dir)
        
        # 両方を絶対パスに変換
        if not p_file.is_absolute():
            p_file = Path.cwd() / p_file
        if not p_base.is_absolute():
            p_base = Path.cwd() / p_base
        
        # 可能であれば相対パスを計算
        try:
            return str(p_file.relative_to(p_base))
        except ValueError:
            # 相対化できない場合はファイル名のみを使用
            return p_file.name
    except Exception:
        # エラー発生時は元のパスを返す
        return str(file_path)

def generate_html_report(region_info, output_base_dir=None, enrichment_results=None, sample_analysis=None, 
                      gene_analysis=None, figures=None, saved_files=None):
    """
    Generates an HTML report of analysis results.
    
    Parameters
    ----------
    region_info : dict
        Region information
    output_base_dir : str or Path, optional
        Base path for the output directory (reference for relative path calculation)
    enrichment_results : pd.DataFrame, optional
        Enrichment analysis results
    sample_analysis : dict, optional
        Sample distribution analysis results
    gene_analysis : dict, optional
        Gene distribution analysis results
    figures : dict, optional
        Dictionary of figure objects
    saved_files : dict, optional
        Dictionary of saved files
        
    Returns
    -------
    str
        Generated HTML string
    """
    try:
        from jinja2 import Template
        has_jinja = True
    except ImportError:
        has_jinja = False
    
    # 地域情報の整理
    region_name = region_info.get('name', 'Unknown')
    chrom = region_info.get('chr', region_info.get('chrom', 'Unknown'))
    start = region_info.get('start', 0)
    end = region_info.get('end', 0)
    size = region_info.get('size', end - start)
    
    # 日時
    now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    # エンリッチメント情報の整理
    enrichment_data = {}
    if enrichment_results is not None and not enrichment_results.empty:
        if 'Significant' in enrichment_results.columns:
            significant = enrichment_results[enrichment_results['Significant'] == True]
            enrichment_data['total'] = len(enrichment_results)
            enrichment_data['significant'] = len(significant)
            
            # 上位5件
            if len(significant) > 0:
                top_enriched = []
                for i, (_, row) in enumerate(significant.sort_values('p-value').head(5).iterrows()):
                    top_enriched.append({
                        'tissue': row.get('Tissue Type', 'Unknown'),
                        'fc': row.get('Fold Change', 0),
                        'pvalue': row.get('p-value', 1),
                        'adj_pvalue': row.get('adjusted p-value', row.get('p-value', 1))
                    })
                enrichment_data['top_enriched'] = top_enriched
    
    # サンプル情報の整理
    sample_data = {}
    if sample_analysis:
        sample_data['count'] = sample_analysis.get('sample_count', 0)
        
        # SE数の統計
        se_stats = sample_analysis.get('se_count_stats', {})
        if se_stats:
            sample_data['se_stats'] = {
                'total': se_stats.get('total', 0),
                'mean': se_stats.get('mean', 0),
                'median': se_stats.get('median', 0),
                'std': se_stats.get('std', 0),
                'min': se_stats.get('min', 0),
                'max': se_stats.get('max', 0)
            }
        
        # 組織タイプの分布
        tissue_counts = sample_analysis.get('tissue_counts')
        if tissue_counts is not None and len(tissue_counts) > 0:
            top_tissues = []
            total_samples = tissue_counts.sum()
            for tissue, count in tissue_counts.nlargest(5).items():
                percent = (count / total_samples * 100) if total_samples > 0 else 0
                top_tissues.append({
                    'name': tissue,
                    'count': count,
                    'percent': percent
                })
            sample_data['top_tissues'] = top_tissues
            sample_data['tissue_count'] = len(tissue_counts)
    
    # 遺伝子情報の整理
    gene_data = {}
    if gene_analysis:
        gene_data['count'] = gene_analysis.get('gene_count', 0)
        
        # 出現回数の統計
        occ_stats = gene_analysis.get('occurrence_stats', {})
        if occ_stats:
            gene_data['occ_stats'] = {
                'total': occ_stats.get('total', 0),
                'mean': occ_stats.get('mean', 0),
                'median': occ_stats.get('median', 0),
                'std': occ_stats.get('std', 0),
                'min': occ_stats.get('min', 0),
                'max': occ_stats.get('max', 0)
            }
        
        # 上位遺伝子
        top_genes = gene_analysis.get('top_genes')
        if top_genes is not None and len(top_genes) > 0:
            top_gene_list = []
            for gene, count in top_genes.nlargest(10).items():
                top_gene_list.append({
                    'name': gene,
                    'count': count
                })
            gene_data['top_genes'] = top_gene_list
    
    # 図ファイル情報の整理
    figure_files = []
    if saved_files and 'human_figures' in saved_files:
        for file_key, file_path in saved_files['human_figures'].items():
            if file_key.endswith('_png'):
                # 安全に相対パスを計算
                rel_path = get_safe_relative_path(file_path, output_base_dir)
                fig_name = file_key.replace('_png', '')
                figure_files.append({
                    'name': fig_name,
                    'path': rel_path
                })
    
    # ファイルリンク情報の整理
    file_links = {
        'tables': [],
        'machine_data': [],
        'summary': None
    }
    
    if saved_files:
        # 人間可読テーブル
        if 'human_tables' in saved_files:
            for file_key, file_path in saved_files['human_tables'].items():
                rel_path = get_safe_relative_path(file_path, output_base_dir)
                file_links['tables'].append({
                    'name': file_key,
                    'path': rel_path
                })
        
        # 機械可読データ
        if 'machine_data' in saved_files:
            for file_key, file_path in saved_files['machine_data'].items():
                rel_path = get_safe_relative_path(file_path, output_base_dir)
                file_links['machine_data'].append({
                    'name': file_key,
                    'path': rel_path
                })
        
        # サマリーテキスト
        if 'summary_txt' in saved_files:
            rel_path = get_safe_relative_path(saved_files['summary_txt'], output_base_dir)
            file_links['summary'] = rel_path
    
    # Jinja2を使用してテンプレートからHTMLを生成
    if has_jinja:
        # Jinja2テンプレート
        template_str = """<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SEdb Analysis Report - {{ region_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; line-height: 1.6; margin: 0; padding: 20px; }
        h1, h2, h3 { color: #333; }
        .container { max-width: 1200px; margin: 0 auto; }
        .section { margin-bottom: 30px; border: 1px solid #eee; padding: 20px; border-radius: 5px; }
        table { border-collapse: collapse; width: 100%; margin-bottom: 15px; }
        table, th, td { border: 1px solid #ddd; }
        th, td { padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .file-link { margin-bottom: 5px; }
        img { max-width: 100%; height: auto; margin-bottom: 20px; box-shadow: 0 0 10px rgba(0,0,0,0.1); }
        .summary-box { background-color: #f8f9fa; padding: 15px; border-radius: 5px; margin-bottom: 20px; }
        .stat-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin-bottom: 15px; }
        .stat-item { background-color: #f2f2f2; padding: 10px; border-radius: 3px; }
        .stat-value { font-weight: bold; font-size: 1.2em; }
        .stat-label { font-size: 0.9em; color: #666; }
    </style>
</head>
<body>
    <div class="container">
        <div class="section">
            <h1>SEdb Analysis Report</h1>
            <p>Generated on: {{ now }}</p>
            
            <div class="summary-box">
                <h3>Target Region</h3>
                <p><strong>Region name:</strong> {{ region_name }}</p>
                <p><strong>Coordinates:</strong> {{ chrom }}:{{ '{:,}'.format(start) }}-{{ '{:,}'.format(end) }} ({{ '{:,}'.format(size) }} bp)</p>
            </div>
        </div>
        
        {% if enrichment_data %}
        <div class="section">
            <h2>Tissue Enrichment Analysis</h2>
            
            <div class="stat-grid">
                <div class="stat-item">
                    <div class="stat-value">{{ enrichment_data.total }}</div>
                    <div class="stat-label">Total tissues</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">{{ enrichment_data.significant }}</div>
                    <div class="stat-label">Significantly enriched tissues</div>
                </div>
            </div>
            
            {% if enrichment_data.top_enriched %}
            <h3>Top Enriched Tissues</h3>
            <table>
                <thead>
                    <tr>
                        <th>Tissue Type</th>
                        <th>Fold Change</th>
                        <th>p-value</th>
                        <th>adjusted p-value</th>
                    </tr>
                </thead>
                <tbody>
                    {% for item in enrichment_data.top_enriched %}
                    <tr>
                        <td>{{ item.tissue }}</td>
                        <td>{{ '{:.2f}'.format(item.fc) }}</td>
                        <td>{{ '{:.2e}'.format(item.pvalue) }}</td>
                        <td>{{ '{:.2e}'.format(item.adj_pvalue) }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
            {% endif %}
        </div>
        {% endif %}
        
        {% if sample_data %}
        <div class="section">
            <h2>Sample Distribution</h2>
            
            <div class="stat-grid">
                <div class="stat-item">
                    <div class="stat-value">{{ '{:,}'.format(sample_data.count) }}</div>
                    <div class="stat-label">Total samples</div>
                </div>
                {% if sample_data.tissue_count %}
                <div class="stat-item">
                    <div class="stat-value">{{ sample_data.tissue_count }}</div>
                    <div class="stat-label">Total tissue types</div>
                </div>
                {% endif %}
                {% if sample_data.se_stats %}
                <div class="stat-item">
                    <div class="stat-value">{{ '{:.2f}'.format(sample_data.se_stats.mean) }}</div>
                    <div class="stat-label">Average SEs per sample</div>
                </div>
                {% endif %}
            </div>
            
            {% if sample_data.se_stats %}
            <h3>SE Count Statistics</h3>
            <table>
                <tr>
                    <th>Total</th>
                    <th>Mean</th>
                    <th>Median</th>
                    <th>Std Dev</th>
                    <th>Min</th>
                    <th>Max</th>
                </tr>
                <tr>
                    <td>{{ '{:,}'.format(sample_data.se_stats.total) }}</td>
                    <td>{{ '{:.2f}'.format(sample_data.se_stats.mean) }}</td>
                    <td>{{ '{:.2f}'.format(sample_data.se_stats.median) }}</td>
                    <td>{{ '{:.2f}'.format(sample_data.se_stats.std) }}</td>
                    <td>{{ sample_data.se_stats.min }}</td>
                    <td>{{ sample_data.se_stats.max }}</td>
                </tr>
            </table>
            {% endif %}
            
            {% if sample_data.top_tissues %}
            <h3>Major Tissue Types</h3>
            <table>
                <thead>
                    <tr>
                        <th>Tissue Type</th>
                        <th>Sample Count</th>
                        <th>Percentage (%)</th>
                    </tr>
                </thead>
                <tbody>
                    {% for tissue in sample_data.top_tissues %}
                    <tr>
                        <td>{{ tissue.name }}</td>
                        <td>{{ tissue.count }}</td>
                        <td>{{ '{:.2f}'.format(tissue.percent) }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
            {% endif %}
        </div>
        {% endif %}
        
        {% if gene_data %}
        <div class="section">
            <h2>Gene Distribution</h2>
            
            <div class="stat-grid">
                <div class="stat-item">
                    <div class="stat-value">{{ '{:,}'.format(gene_data.count) }}</div>
                    <div class="stat-label">Detected genes</div>
                </div>
                {% if gene_data.occ_stats %}
                <div class="stat-item">
                    <div class="stat-value">{{ '{:.2f}'.format(gene_data.occ_stats.mean) }}</div>
                    <div class="stat-label">Average occurrences per gene</div>
                </div>
                <div class="stat-item">
                    <div class="stat-value">{{ gene_data.occ_stats.max }}</div>
                    <div class="stat-label">Maximum occurrences</div>
                </div>
                {% endif %}
            </div>
            
            {% if gene_data.top_genes %}
            <h3>Most Frequent Genes</h3>
            <table>
                <thead>
                    <tr>
                        <th>Gene Name</th>
                        <th>Occurrence Count</th>
                    </tr>
                </thead>
                <tbody>
                    {% for gene in gene_data.top_genes %}
                    <tr>
                        <td>{{ gene.name }}</td>
                        <td>{{ gene.count }}</td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
            {% endif %}
        </div>
        {% endif %}
        
        {% if figure_files %}
        <div class="section">
            <h2>Analysis Figures</h2>
            
            {% for figure in figure_files %}
            <div>
                <h3>{{ figure.name }}</h3>
                <img src="{{ figure.path }}" alt="{{ figure.name }}">
            </div>
            {% endfor %}
        </div>
        {% endif %}
        
        <div class="section">
            <h2>Analysis Result Files</h2>
            
            {% if file_links.summary %}
            <h3>Overall Summary</h3>
            <p><a href="{{ file_links.summary }}" target="_blank">Summary Text</a></p>
            {% endif %}
            
            {% if file_links.tables %}
            <h3>Tables and Summaries</h3>
            <ul>
                {% for file in file_links.tables %}
                <li class="file-link"><a href="{{ file.path }}" target="_blank">{{ file.name }}</a></li>
                {% endfor %}
            </ul>
            {% endif %}
            
            {% if file_links.machine_data %}
            <h3>Machine Readable Data</h3>
            <ul>
                {% for file in file_links.machine_data %}
                <li class="file-link"><a href="{{ file.path }}" target="_blank">{{ file.name }}</a></li>
                {% endfor %}
            </ul>
            {% endif %}
        </div>
    </div>
</body>
</html>
"""
        
        # テンプレートレンダリング
        template = Template(template_str)
        html = template.render(
            region_name=region_name,
            chrom=chrom,
            start=start,
            end=end,
            size=size,
            now=now,
            enrichment_data=enrichment_data,
            sample_data=sample_data,
            gene_data=gene_data,
            figure_files=figure_files,
            file_links=file_links
        )
        
        return html
    else:
        # Jinja2が利用できない場合は例外を発生させる
        raise ImportError("Jinja2 is not installed. Jinja2 is required to generate HTML reports.")