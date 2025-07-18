# sedb_tools/analysis/statistics.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import warnings
from collections import Counter
from typing import Optional, Dict, List, Tuple, Any
import networkx as nx
from scipy import stats
from statsmodels.stats.multitest import multipletests

class SEdbStatisticsAnalyzer:
    """
    Super-enhancer (SE) 統計分析クラス
    エンリッチメント分析やサンプル分布分析などの統計関連機能を提供
    """
    
    def __init__(self, logger=None):
        """
        Initialize the statistics analyzer
        
        Parameters:
        logger (Logger): Optional logger instance to use. If None, creates a basic console logger.
        """
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbStatisticsAnalyzer")
    
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
    
    def test_tissue_enrichment(self, region_tissue_counts, background_distribution, method='fisher', correction='fdr_bh'):
        """
        Test for tissue type enrichment in a specific region.
        
        Parameters:
        -----------
        region_tissue_counts : dict or pandas.Series
            Tissue type counts in the specific region
        background_distribution : pandas.DataFrame or pandas.Series
            Background distribution of tissue types - if DataFrame, must have 'Tissue Type' and 'Sample Count' columns
        method : str, optional
            Statistical test method ('fisher' or 'chi2')
        correction : str, optional
            Multiple testing correction method. Same as the method parameter in statsmodels.stats.multitest.multipletests.
            
        Returns:
        --------
        pandas.DataFrame
            Test results for each tissue type (p-values, corrected p-values, odds ratios, etc.)
        """
        self.logger.info("Testing tissue enrichment...")
        
        # Convert background distribution to Series if it's a DataFrame
        if isinstance(background_distribution, pd.DataFrame):
            if 'Tissue Type' in background_distribution.columns and 'Sample Count' in background_distribution.columns:
                background = background_distribution.set_index('Tissue Type')['Sample Count']
            else:
                self.logger.error("Background distribution DataFrame must have 'Tissue Type' and 'Sample Count' columns")
                return None
        else:
            background = background_distribution
            
        total_background = background.sum()
        
        # Convert region counts to Series (if it's a dictionary)
        if isinstance(region_tissue_counts, dict):
            region_counts = pd.Series(region_tissue_counts)
        else:
            region_counts = region_tissue_counts
        
        total_region = region_counts.sum()
        
        # List to store results
        results = []
        
        # Perform test for each tissue type
        for tissue in background.index:
            if tissue in region_counts.index:
                a = region_counts[tissue]  # Count of tissue in region
            else:
                a = 0
            
            b = total_region - a  # Count of other tissues in region
            c = background[tissue]  # Count of tissue in overall data
            d = total_background - c  # Count of other tissues in overall data
            
            # 2x2 contingency table
            contingency_table = np.array([[a, b], [c, d]])
            
            if method == 'fisher':
                # Fisher's exact test
                odds_ratio, p_value = stats.fisher_exact(contingency_table)
            elif method == 'chi2':
                # Chi-squared test
                chi2, p_value, dof, expected = stats.chi2_contingency(contingency_table)
                odds_ratio = (a * d) / (b * c) if (b * c != 0) else np.nan
            else:
                raise ValueError(f"Unsupported method: {method}")
            
            # Observed and expected values
            observed = a
            expected = (c / total_background) * total_region
            
            # Fold change
            fold_change = (a / total_region) / (c / total_background) if (c != 0) else np.nan
            
            results.append({
                'Tissue Type': tissue,
                'Region Count': a,
                'Region Percentage': (a / total_region * 100) if total_region > 0 else 0,
                'Total Count': c,
                'Total Percentage': (c / total_background * 100),
                'Expected Value': expected,
                'Odds Ratio': odds_ratio,
                'Fold Change': fold_change,
                'p-value': p_value
            })
        
        # Convert results to DataFrame
        result_df = pd.DataFrame(results)
        
        # Sort by p-value
        result_df = result_df.sort_values('p-value')
        
        # Multiple testing correction
        if len(result_df) > 1:
            reject, p_corr, _, _ = multipletests(result_df['p-value'].values, method=correction)
            result_df['adjusted p-value'] = p_corr
            result_df['Significant'] = reject
        
        return result_df
    
    def analyze_region_tissue_enrichment(self, region_tissue_counts, background_distribution, region_name=None, 
                                         method='fisher', correction='fdr_bh', **kwargs):
        """
        Analyze tissue enrichment in the overlapping SEs compared to the background distribution.
        Counts are based on unique samples (each sample is counted only once).
        
        Parameters:
        -----------
        region_tissue_counts : dict or pandas.Series
            Region tissue counts
        background_distribution : pandas.DataFrame or pandas.Series
            Background distribution of tissue types
        region_name : str, optional
            Region name for the plot title
        method : str, optional
            Statistical test method ('fisher' or 'chi2')
        correction : str, optional
            Multiple testing correction method
        **kwargs : dict
            Additional parameters for visualization
        
        Returns:
        --------
        pandas.DataFrame
            Results of the enrichment analysis
        """
        if region_tissue_counts is None:
            self.logger.warning("No region tissue counts provided.")
            return None
            
        region_name_str = region_name or "Unknown Region"
        self.logger.info(f"Analyzing tissue enrichment for region {region_name_str} based on unique samples")
        
        # Perform enrichment test
        enrichment_results = self.test_tissue_enrichment(
            region_tissue_counts,
            background_distribution,
            method=method,
            correction=correction
        )
        
        self.logger.info("Tissue enrichment analysis completed based on unique samples")
        
        return enrichment_results
    
    def analyze_samples_distribution(self, sample_summary, figsize=(12, 8)):
        """
        Analyze and visualize the distribution of samples and their SEs in the region.
        
        Parameters:
        -----------
        sample_summary : pandas.DataFrame
            DataFrame with sample information and SE counts
        figsize : tuple, optional
            Figure size for the plots
        
        Returns:
        --------
        dict
            Dictionary with analysis results
        """
        if sample_summary is None or len(sample_summary) == 0:
            self.logger.warning("Sample summary data not available. Cannot analyze distribution.")
            return None
        
        self.logger.info("Analyzing sample distribution...")
        
        results = {
            'sample_count': len(sample_summary),
            'se_count_stats': None,
            'tissue_counts': None,
            'biosample_counts': None,
            'top_samples': None
        }
        
        # SE count per sample statistics
        se_count_stats = sample_summary['se_count'].describe().to_dict()
        results['se_count_stats'] = se_count_stats
        
        # Tissue distribution
        if 'Tissue type' in sample_summary.columns:
            tissue_counts = sample_summary['Tissue type'].value_counts()
            tissue_percent = (tissue_counts / len(sample_summary) * 100).round(2)
            tissue_df = pd.DataFrame({
                'Sample Count': tissue_counts,
                'Percentage': tissue_percent
            })
            results['tissue_counts'] = tissue_df
        
        # Biosample distribution
        if 'Biosample type' in sample_summary.columns:
            biosample_counts = sample_summary['Biosample type'].value_counts()
            biosample_percent = (biosample_counts / len(sample_summary) * 100).round(2)
            biosample_df = pd.DataFrame({
                'Sample Count': biosample_counts,
                'Percentage': biosample_percent
            })
            results['biosample_counts'] = biosample_df
        
        # Top samples by SE count
        top_samples = sample_summary.sort_values('se_count', ascending=False).head(10)
        results['top_samples'] = top_samples
        
        # Generate plots
        if hasattr(plt, 'figure'):  # Check if matplotlib is available
            try:
                fig, axes = plt.subplots(2, 2, figsize=figsize)
                
                # Plot 1: SE count distribution
                sns.histplot(data=sample_summary, x='se_count', bins=20, ax=axes[0, 0])
                axes[0, 0].set_title('Distribution of SE Counts per Sample')
                axes[0, 0].set_xlabel('Number of SEs per Sample')
                axes[0, 0].set_ylabel('Number of Samples')
                
                # Plot 2: Tissue type pie chart
                if 'Tissue type' in sample_summary.columns:
                    tissue_counts = sample_summary['Tissue type'].value_counts()
                    top_tissues = tissue_counts.head(8)
                    
                    # Group others if needed
                    if len(tissue_counts) > 8:
                        other_count = tissue_counts[8:].sum()
                        top_tissues['Other'] = other_count
                    
                    top_tissues.plot.pie(autopct='%1.1f%%', ax=axes[0, 1])
                    axes[0, 1].set_title('Distribution of Samples by Tissue Type')
                    axes[0, 1].set_ylabel('')
                
                # Plot 3: Top samples by SE count
                if len(top_samples) > 0:
                    sns.barplot(data=top_samples, y='cell_id', x='se_count', ax=axes[1, 0])
                    axes[1, 0].set_title('Top 10 Samples by SE Count')
                    axes[1, 0].set_xlabel('Number of SEs')
                    axes[1, 0].set_ylabel('Sample ID')
                
                # Plot 4: SE count vs Tissue Type
                if 'Tissue type' in sample_summary.columns:
                    sns.boxplot(data=sample_summary, y='Tissue type', x='se_count', ax=axes[1, 1])
                    axes[1, 1].set_title('SE Count by Tissue Type')
                    axes[1, 1].set_xlabel('Number of SEs')
                    axes[1, 1].set_ylabel('Tissue Type')
                
                plt.tight_layout()
                results['figure'] = fig
                
            except Exception as e:
                self.logger.warning(f"Error generating sample distribution plots: {e}")
        
        return results
    
    def refine_tissue_enrichment_results(self, enrichment_results, fdr_threshold=0.05, fc_threshold=1.5, 
                                         top_n=None, rename_columns=True):
        """
        Refine tissue enrichment analysis results and add significance flags.
        
        Parameters:
        -----------
        enrichment_results : pandas.DataFrame
            Analysis results DataFrame
        fdr_threshold : float, optional
            FDR threshold for significance
        fc_threshold : float, optional
            Fold change threshold for significance
        top_n : int, optional
            Number of top tissues to display. If None, display all
        rename_columns : bool, optional
            Whether to convert column names to English
            
        Returns:
        --------
        pandas.DataFrame
            Refined analysis results
        """
        if enrichment_results is None or len(enrichment_results) == 0:
            self.logger.warning("No enrichment analysis results to refine.")
            return None
        
        # Copy the results for manipulation
        results = enrichment_results.copy()
        
        # Check and rename FDR column
        if '補正p値' in results.columns:
            fdr_col = '補正p値'
        elif 'adjusted p-value' in results.columns:
            fdr_col = 'adjusted p-value'
        else:
            # Use regular p-value if FDR column is not available
            fdr_col = 'p値' if 'p値' in results.columns else 'p-value'
            self.logger.warning(f"FDR-corrected p-value not found. Using {fdr_col} instead.")
        
        # Check p-value column
        if 'p値' in results.columns:
            p_col = 'p値'
        else:
            p_col = 'p-value'
        
        # Check fold change column
        if 'フォールドチェンジ' in results.columns:
            fc_col = 'フォールドチェンジ'
        else:
            fc_col = 'Fold Change'
        
        # Check tissue type column
        if '組織型' in results.columns:
            tissue_col = '組織型'
        else:
            tissue_col = 'Tissue Type'
        
        # Check count columns
        if 'リージョンカウント' in results.columns:
            region_count_col = 'リージョンカウント'
        else:
            region_count_col = 'Region Count'
            
        if '全体カウント' in results.columns:
            total_count_col = '全体カウント'
        else:
            total_count_col = 'Total Count'
        
        # Add significance column
        results['Significant'] = (results[fdr_col] < fdr_threshold) & (results[fc_col] > fc_threshold)
        
        # Check for negative enrichment (depletion)
        results['Depleted'] = (results[fdr_col] < fdr_threshold) & (results[fc_col] < (1/fc_threshold))
        
        # Sort for display (by FDR)
        results = results.sort_values(by=fdr_col)
        
        # Limit to top N entries
        if top_n is not None and top_n > 0 and len(results) > top_n:
            results = results.head(top_n)
        
        # Rename columns to English
        if rename_columns:
            column_mapping = {
                '組織型': 'Tissue Type',
                'リージョンカウント': 'Region Count',
                '全体カウント': 'Total Count',
                'フォールドチェンジ': 'Fold Change',
                '補正p値': 'FDR',
                'p値': 'p-value',
                '有意': 'Significant',
                'デプリート': 'Depleted'
            }
            
            # Only rename existing columns
            rename_cols = {col: column_mapping[col] for col in results.columns if col in column_mapping}
            results = results.rename(columns=rename_cols)
        
        # Select and order main columns for display
        main_columns = [
            tissue_col if not rename_columns else 'Tissue Type',
            region_count_col if not rename_columns else 'Region Count',
            total_count_col if not rename_columns else 'Total Count',
            fc_col if not rename_columns else 'Fold Change',
            fdr_col if not rename_columns else 'FDR',
            'Significant' if rename_columns else '有意',
            'Depleted' if rename_columns else 'デプリート'
        ]
        
        # Select only columns that exist
        available_columns = [col for col in main_columns if col in results.columns]
        results = results[available_columns]
        
        return results
    
    def compare_multiple_regions(self, region_results_dict, plot_type='heatmap', 
                                 figsize=(14, 10), p_threshold=0.05):
        """
        Compare tissue enrichment results across multiple regions.
        
        Parameters:
        -----------
        region_results_dict : dict
            Dictionary with region names as keys and results from test_tissue_enrichment() as values
        plot_type : str, optional
            Plot type ('heatmap', 'barplot', 'network')
        figsize : tuple, optional
            Figure size
        p_threshold : float, optional
            P-value threshold for significance
        
        Returns:
        --------
        tuple
            (figure, dataframe) - Created figure object and cross-comparison dataframe
        """
        if not region_results_dict or len(region_results_dict) == 0:
            self.logger.warning("Region results dictionary is empty.")
            return None, None
        
        # List of region names
        region_names = list(region_results_dict.keys())
        
        # Find all unique tissue types
        all_tissues = set()
        for results in region_results_dict.values():
            all_tissues.update(results['Tissue Type'].values)
        
        # Collect significant tissue types for each region
        significant_tissues = {}
        for region, results in region_results_dict.items():
            p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
            sig_tissues = results[results[p_col] <= p_threshold]['Tissue Type'].values
            significant_tissues[region] = set(sig_tissues)
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # 比較結果を保存するためのデータ構造
        comparison_data = None
        
        if plot_type == 'heatmap':
            # Create dataframe for heatmap
            heatmap_data = []
            
            for tissue in sorted(all_tissues):
                row = {'Tissue Type': tissue}
                
                for region in region_names:
                    results = region_results_dict[region]
                    tissue_row = results[results['Tissue Type'] == tissue]
                    
                    if len(tissue_row) > 0:
                        fc = tissue_row['Fold Change'].values[0]
                        p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                        p_val = tissue_row[p_col].values[0]
                        
                        # Set fold change to 1 (no change) for non-significant results
                        if p_val > p_threshold:
                            row[region] = 1.0
                        else:
                            row[region] = fc
                    else:
                        row[region] = 1.0  # No data means fold change = 1 (no change)
                
                heatmap_data.append(row)
            
            # Convert to DataFrame
            heatmap_df = pd.DataFrame(heatmap_data).set_index('Tissue Type')
            comparison_data = heatmap_df
            
            # Filter by significance (significant in at least one region)
            sig_in_any = set()
            for sig_set in significant_tissues.values():
                sig_in_any.update(sig_set)
            
            if len(sig_in_any) > 0:
                filtered_df = heatmap_df.loc[list(sig_in_any)]
            else:
                filtered_df = heatmap_df
            
            # Sort by maximum fold change
            max_fc = filtered_df.max(axis=1)
            filtered_df = filtered_df.iloc[max_fc.argsort(ascending=False)]
            
            # Create heatmap
            cmap = sns.diverging_palette(220, 10, as_cmap=True)
            
            # Create mask for values = 1 (to display in gray)
            mask = filtered_df == 1.0
            
            # Draw heatmap
            sns.heatmap(filtered_df, cmap=cmap, center=1, annot=True, fmt='.2f', linewidths=0.5,
                    mask=mask, cbar_kws={'label': 'Fold Change'})
            
            plt.title(f'Tissue Enrichment Comparison Across Multiple Regions (p < {p_threshold})')
            
        elif plot_type == 'barplot':
            # Collect all significant tissue types
            all_sig_tissues = set()
            for sig_set in significant_tissues.values():
                all_sig_tissues.update(sig_set)
            
            # Limit to top 15 tissues
            if len(all_sig_tissues) > 15:
                # Get maximum fold change for each tissue across all regions
                max_fc_by_tissue = {}
                for tissue in all_sig_tissues:
                    max_fc = 0
                    for region, results in region_results_dict.items():
                        tissue_row = results[results['Tissue Type'] == tissue]
                        if len(tissue_row) > 0:
                            fc = abs(tissue_row['Fold Change'].values[0] - 1.0)  # Absolute change from 1.0
                            max_fc = max(max_fc, fc)
                    max_fc_by_tissue[tissue] = max_fc
                
                # Sort by maximum fold change and get top 15 tissues
                sorted_tissues = sorted(max_fc_by_tissue.items(), key=lambda x: x[1], reverse=True)
                all_sig_tissues = {t for t, _ in sorted_tissues[:15]}
            
            # Plot fold changes for significant tissue types in each region
            n_regions = len(region_names)
            region_colors = plt.cm.tab10(np.arange(n_regions) % 10)
            
            for i, (region, color) in enumerate(zip(region_names, region_colors)):
                results = region_results_dict[region]
                p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                
                # Filter to get significant tissue types only
                sig_results = results[
                    (results[p_col] <= p_threshold) & 
                    (results['Tissue Type'].isin(all_sig_tissues))
                ]
                
                if len(sig_results) > 0:
                    plt.subplot(n_regions, 1, i+1)
                    bars = plt.barh(
                        y=sig_results['Tissue Type'],
                        width=sig_results['Fold Change'],
                        color=color,
                        alpha=0.7
                    )
                    
                    plt.axvline(x=1, linestyle='--', color='gray')
                    plt.xlim(0, max(5, sig_results['Fold Change'].max() * 1.1))
                    plt.title(f'{region} - Significant Tissue Enrichment')
                    
                    if i == n_regions - 1:
                        plt.xlabel('Fold Change')
                    
                    # Display p-values
                    for bar, p_val in zip(bars, sig_results[p_col]):
                        plt.text(
                            x=max(bar.get_width() + 0.1, 1.1),
                            y=bar.get_y() + bar.get_height() * 0.5,
                            s=f"p={p_val:.1e}",
                            va='center',
                            ha='left',
                            fontsize=8
                        )
            
            # Collect data for comparison
            comparison_data = pd.DataFrame()
            for region, results in region_results_dict.items():
                p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                sig_results = results[results[p_col] <= p_threshold].copy()
                if len(sig_results) > 0:
                    sig_results['Region'] = region
                    comparison_data = pd.concat([comparison_data, sig_results])
            
        elif plot_type == 'network':
            # Prepare nodes and edges for network plot
            
            G = nx.Graph()
            
            # Add region nodes
            for region in region_names:
                G.add_node(region, type='region')
            
            # Add significant tissue type nodes
            all_sig_tissues = set()
            for sig_set in significant_tissues.values():
                all_sig_tissues.update(sig_set)
            
            # Add tissue nodes
            for tissue in all_sig_tissues:
                G.add_node(tissue, type='tissue')
            
            # Add edges (between regions and their significant tissue types)
            for region, results in region_results_dict.items():
                p_col = 'adjusted p-value' if 'adjusted p-value' in results.columns else 'p-value'
                
                for _, row in results[(results[p_col] <= p_threshold) & 
                                    (results['Tissue Type'].isin(all_sig_tissues))].iterrows():
                    tissue = row['Tissue Type']
                    fc = row['Fold Change']
                    G.add_edge(region, tissue, weight=fc)
            
            # Calculate node positions
            pos = nx.spring_layout(G, k=0.3, seed=42)
            
            # Separate region nodes and tissue nodes
            region_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'region']
            tissue_nodes = [n for n, attr in G.nodes(data=True) if attr.get('type') == 'tissue']
            
            # Edge colors based on fold change
            edge_colors = []
            for u, v in G.edges():
                fc = G[u][v]['weight']
                if fc > 1:
                    edge_colors.append('red')
                else:
                    edge_colors.append('blue')
            
            # Edge width proportional to fold change
            edge_widths = [abs(G[u][v]['weight'] - 1) * 2 + 0.5 for u, v in G.edges()]
            
            # Draw network
            plt.figure(figsize=figsize)
            
            # Draw region nodes
            nx.draw_networkx_nodes(G, pos, nodelist=region_nodes, node_size=800,
                                node_color='lightgreen', alpha=0.8, 
                                node_shape='s', label='Region')
            
            # Draw tissue nodes
            nx.draw_networkx_nodes(G, pos, nodelist=tissue_nodes, node_size=600,
                                node_color='lightblue', alpha=0.8,
                                node_shape='o', label='Tissue Type')
            
            # Draw edges
            nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.7,
                                edge_color=edge_colors)
            
            # Draw labels
            nx.draw_networkx_labels(G, pos, font_size=8)
            
            plt.title('Network of Regions and Significant Tissue Types')
            plt.legend()
            plt.axis('off')
            
            # Create comparison data
            edge_data = []
            for u, v, data in G.edges(data=True):
                if u in region_nodes and v in tissue_nodes:
                    edge_data.append({
                        'Region': u,
                        'Tissue Type': v,
                        'Fold Change': data['weight']
                    })
            comparison_data = pd.DataFrame(edge_data)
            
        else:
            self.logger.warning(f"Unknown plot type: {plot_type}")
            plt.close(fig)
            return None, None
        
        plt.tight_layout()
        
        return fig, comparison_data