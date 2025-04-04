"""
sedb_analyzer.statistics - SEdb Statistical Analysis Module

This module provides functionalities for statistical and enrichment analysis on SEdb data.
The computational logic and visualization logic are separated, and each function has clearly defined responsibilities.
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Any, Optional, Union, Tuple
from scipy import stats
from statsmodels.stats.multitest import multipletests
from collections import Counter, defaultdict
import warnings

class StatisticsAnalyzer:
    """
    A class that provides statistical analysis of SEdb data.
    """
    
    def __init__(self, logger=None):
        """
        Initialization
        
        Parameters
        ----------
        logger : Logger, optional
            The logger instance to be used.
        """
        # ロガーの設定
        self.logger = logger or self._setup_default_logger()
    
    def _setup_default_logger(self):
        """Set up the default logger."""
        logger = logging.getLogger(f"{__name__}.{id(self)}")
        
        # 既存のハンドラを確認
        if not logger.handlers:
            logger.setLevel(logging.INFO)
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.propagate = False
        
        return logger
    
    def analyze_samples_distribution(self, sample_df: pd.DataFrame, 
                                   se_count_col: str = 'se_count') -> Dict[str, Any]:
        """
        Analyze the distribution of samples and SE counts.
        
        Parameters
        ----------
        sample_df : pd.DataFrame
            DataFrame containing sample information.
        se_count_col : str, optional
            The column name for SE counts.
            
        Returns
        -------
        Dict[str, Any]
            A dictionary containing the analysis results (excluding figures).
        """
        self.logger.info("Starting sample distribution analysis")
        
        if sample_df is None or len(sample_df) == 0:
            self.logger.warning("Sample data is empty. Skipping analysis.")
            return {
                'sample_count': 0,
                'se_count_stats': {},
                'tissue_counts': None,
                'biosample_counts': None
            }
        
        # SEカウントがない場合は0で埋める
        if se_count_col not in sample_df.columns:
            self.logger.warning(f"Column '{se_count_col}' not found. Filling with 0.")
            sample_df[se_count_col] = 0
        
        # 基本統計の計算
        sample_count = len(sample_df)
        
        # SEカウントの統計
        se_count_stats = {
            'total': int(sample_df[se_count_col].sum()),
            'mean': float(sample_df[se_count_col].mean()),
            'median': float(sample_df[se_count_col].median()),
            'std': float(sample_df[se_count_col].std()),
            'min': int(sample_df[se_count_col].min()),
            'max': int(sample_df[se_count_col].max())
        }
        
        # 組織タイプのカウント
        tissue_counts = None
        for col_name in ['Tissue type', 'Tissue']:
            if col_name in sample_df.columns:
                tissue_counts = sample_df[col_name].value_counts()
                self.logger.debug(f"Tissue type counts: {len(tissue_counts)} categories")
                break
        
        # バイオサンプルタイプのカウント
        biosample_counts = None
        for col_name in ['Biosample type', 'Biosample']:
            if col_name in sample_df.columns:
                biosample_counts = sample_df[col_name].value_counts()
                self.logger.debug(f"Biosample type counts: {len(biosample_counts)} categories")
                break
        
        # 結果を辞書としてまとめる
        result = {
            'sample_count': sample_count,
            'se_count_stats': se_count_stats,
            'tissue_counts': tissue_counts,
            'biosample_counts': biosample_counts
        }
        
        self.logger.info("Sample distribution analysis completed")
        return result

    
    def plot_samples_distribution(self, analysis_result: Dict[str, Any],
                                # figsize=(12, 8) # 元のサイズ
                                figsize: Tuple[int, int] = (12, 5)) -> Optional[Any]:  # 1x2 用に調整
        """
        Plot the analysis results of sample distribution (displaying only tissue and biosample types).
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns

            # プロットに必要なデータを取得
            tissue_counts = analysis_result.get('tissue_counts')
            biosample_counts = analysis_result.get('biosample_counts')

            # 実際にプロットできるデータがあるかチェック
            can_plot_tissue = tissue_counts is not None and not tissue_counts.empty
            can_plot_biosample = biosample_counts is not None and not biosample_counts.empty

            if not can_plot_tissue and not can_plot_biosample:
                self.logger.warning("Neither tissue type nor biosample type data is available. Unable to generate distribution plots.")
                return None

            # Figure を作成
            fig = plt.figure(figsize=figsize)

            plot_index = 1  # subplot のインデックス用

            # --- 組織タイプのプロット ---
            ax1 = fig.add_subplot(1, 2, plot_index)  # 1行2列の plot_index 番目
            if can_plot_tissue:
                top_tissues = tissue_counts.nlargest(10)
                sns.barplot(ax=ax1, x=top_tissues.values, y=top_tissues.index)
                ax1.set_title('Top 10 Tissue Types')
                ax1.set_xlabel('Sample Count')
                ax1.set_ylabel('Tissue type')
            else:
                # データがない場合の表示
                ax1.text(0.5, 0.5, 'No Tissue Type data available', horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes)
                ax1.set_title('Top 10 Tissue Types')
                ax1.set_xticks([])  # 軸を消すなど
                ax1.set_yticks([])
            plot_index += 1

            # --- バイオサンプルタイプのプロット ---
            ax2 = fig.add_subplot(1, 2, plot_index)  # 1行2列の plot_index 番目
            if can_plot_biosample:
                top_biosamples = biosample_counts.nlargest(10)
                sns.barplot(ax=ax2, x=top_biosamples.values, y=top_biosamples.index)
                ax2.set_title('Top 10 Biosample Types')
                ax2.set_xlabel('Sample Count')
                ax2.set_ylabel('Biosample type')
            else:
                # データがない場合の表示
                ax2.text(0.5, 0.5, 'No Biosample Type data available', horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)
                ax2.set_title('Top 10 Biosample Types')
                ax2.set_xticks([])
                ax2.set_yticks([])
            plot_index += 1

            plt.tight_layout()  # プロット間のスペースを調整
            return fig

        except ImportError:
            self.logger.debug("matplotlib/seaborn is not installed, plot cannot be generated")
            return None
        except Exception as e:
            self.logger.error(f"An error occurred during the generation of sample distribution plots: {e}", exc_info=True)
            # エラー発生時に figure をクローズする試み
            try:
                if 'fig' in locals() and fig:
                    plt.close(fig)
            except Exception:
                pass
            return None








    def test_tissue_enrichment(self, region_tissue_counts, 
                            background_tissue_counts,
                            method: str = 'fisher',
                            correction: str = 'fdr_bh',
                            min_region_count: int = 1,
                            pvalue_threshold: float = 0.05,
                            fc_threshold: float = 1.5) -> pd.DataFrame:
        """
        Perform enrichment analysis for tissue types.
        
        Parameters
        ----------
        region_tissue_counts : pd.Series or pd.DataFrame
            Tissue type counts within the region.
            Can be a Series, or a DataFrame containing columns 'tissue_type' and 'count' or 'Tissue Type' and 'Count'.
        background_tissue_counts : pd.Series or pd.DataFrame
            Overall tissue type counts.
            Can be a Series, or a DataFrame containing columns 'tissue_type' and 'count' or 'Tissue Type' and 'Count'.
        method : str, optional
            Statistical test method ('fisher' or 'chi2').
        correction : str, optional
            Method for multiple testing correction.
        min_region_count : int, optional
            Minimum region count to include in the analysis.
        pvalue_threshold : float, optional
            Significance threshold.
        fc_threshold : float, optional
            Fold change threshold to consider enrichment significant.
            
        Returns
        -------
        pd.DataFrame
            Enrichment analysis results.
        """
        self.logger.info("Starting tissue type enrichment analysis")
        
        # region_tissue_countsの型チェックと変換
        if isinstance(region_tissue_counts, pd.DataFrame):
            self.logger.info(f"Detected region_tissue_counts as DataFrame: column names = {list(region_tissue_counts.columns)}")
            
            # カラム名のパターンを確認
            if 'tissue_type' in region_tissue_counts.columns and 'count' in region_tissue_counts.columns:
                self.logger.info("Detected pattern with 'tissue_type' and 'count' columns")
                # この形式なら数値型に変換
                try:
                    original_dtype = region_tissue_counts['count'].dtype
                    region_tissue_counts_converted = pd.Series(
                        pd.to_numeric(region_tissue_counts['count'], errors='coerce'), 
                        index=region_tissue_counts['tissue_type']
                    )
                    self.logger.info(f"Converted region_tissue_counts from DataFrame to Series: original dtype = {original_dtype} → converted dtype = {region_tissue_counts_converted.dtype}")
                    region_tissue_counts = region_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert region_tissue_counts: {e}")
            
            elif 'Tissue Type' in region_tissue_counts.columns and 'Count' in region_tissue_counts.columns:
                self.logger.info("Detected pattern with 'Tissue Type' and 'Count' columns")
                # この形式なら数値型に変換
                try:
                    original_dtype = region_tissue_counts['Count'].dtype
                    region_tissue_counts_converted = pd.Series(
                        pd.to_numeric(region_tissue_counts['Count'], errors='coerce'), 
                        index=region_tissue_counts['Tissue Type']
                    )
                    self.logger.info(f"Converted region_tissue_counts from DataFrame to Series: original dtype = {original_dtype} → converted dtype = {region_tissue_counts_converted.dtype}")
                    region_tissue_counts = region_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert region_tissue_counts: {e}")
            
            else:
                self.logger.warning(f"Unrecognized column pattern: {list(region_tissue_counts.columns)}")
        
        # region_tissue_countsがSeriesだが値が数値でない場合
        elif isinstance(region_tissue_counts, pd.Series):
            if not pd.api.types.is_numeric_dtype(region_tissue_counts):
                self.logger.info(f"Detected region_tissue_counts as Series with non-numeric dtype ({region_tissue_counts.dtype})")
                try:
                    original_dtype = region_tissue_counts.dtype
                    region_tissue_counts_converted = pd.to_numeric(region_tissue_counts, errors='coerce')
                    self.logger.info(f"Converted region_tissue_counts to numeric type: original dtype = {original_dtype} → converted dtype = {region_tissue_counts_converted.dtype}")
                    region_tissue_counts = region_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert region_tissue_counts to numeric: {e}")
            else:
                self.logger.debug(f"region_tissue_counts is already a numeric Series ({region_tissue_counts.dtype})")
        
        # background_tissue_countsの型チェックと変換
        if isinstance(background_tissue_counts, pd.DataFrame):
            self.logger.info(f"Detected background_tissue_counts as DataFrame: column names = {list(background_tissue_counts.columns)}")
            
            # カラム名のパターンを確認
            if 'tissue_type' in background_tissue_counts.columns and 'count' in background_tissue_counts.columns:
                self.logger.info("Detected pattern with 'tissue_type' and 'count' columns")
                # この形式なら数値型に変換
                try:
                    original_dtype = background_tissue_counts['count'].dtype
                    background_tissue_counts_converted = pd.Series(
                        pd.to_numeric(background_tissue_counts['count'], errors='coerce'), 
                        index=background_tissue_counts['tissue_type']
                    )
                    self.logger.info(f"Converted background_tissue_counts from DataFrame to Series: original dtype = {original_dtype} → converted dtype = {background_tissue_counts_converted.dtype}")
                    background_tissue_counts = background_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert background_tissue_counts: {e}")
            
            elif 'Tissue Type' in background_tissue_counts.columns and 'Count' in background_tissue_counts.columns:
                self.logger.info("Detected pattern with 'Tissue Type' and 'Count' columns")
                # この形式なら数値型に変換
                try:
                    original_dtype = background_tissue_counts['Count'].dtype
                    background_tissue_counts_converted = pd.Series(
                        pd.to_numeric(background_tissue_counts['Count'], errors='coerce'), 
                        index=background_tissue_counts['Tissue Type']
                    )
                    self.logger.info(f"Converted background_tissue_counts from DataFrame to Series: original dtype = {original_dtype} → converted dtype = {background_tissue_counts_converted.dtype}")
                    background_tissue_counts = background_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert background_tissue_counts: {e}")
            
            else:
                self.logger.warning(f"Unrecognized column pattern: {list(background_tissue_counts.columns)}")
        
        # background_tissue_countsがSeriesだが値が数値でない場合
        elif isinstance(background_tissue_counts, pd.Series):
            if not pd.api.types.is_numeric_dtype(background_tissue_counts):
                self.logger.info(f"Detected background_tissue_counts as Series with non-numeric dtype ({background_tissue_counts.dtype})")
                try:
                    original_dtype = background_tissue_counts.dtype
                    background_tissue_counts_converted = pd.to_numeric(background_tissue_counts, errors='coerce')
                    self.logger.info(f"Converted background_tissue_counts to numeric type: original dtype = {original_dtype} → converted dtype = {background_tissue_counts_converted.dtype}")
                    background_tissue_counts = background_tissue_counts_converted
                except Exception as e:
                    self.logger.error(f"Failed to convert background_tissue_counts to numeric: {e}")
            else:
                self.logger.debug(f"background_tissue_counts is already a numeric Series ({background_tissue_counts.dtype})")
        
        
        if region_tissue_counts is None or len(region_tissue_counts) == 0:
            self.logger.warning("Region tissue data is empty. Skipping analysis.")
            return pd.DataFrame()
        
        if background_tissue_counts is None or len(background_tissue_counts) == 0:
            self.logger.warning("Background tissue data is empty. Skipping analysis.")
            return pd.DataFrame()
        
        # 解析対象の組織タイプ（最小カウント以上のもの）
        target_tissues = region_tissue_counts[region_tissue_counts >= min_region_count].index
        self.logger.debug(f"Target tissue types for analysis: {len(target_tissues)} categories")
        
        # デバッグ情報
        self.logger.info(f"Type of region_tissue_counts: {type(region_tissue_counts)}")
        self.logger.info(f"Data type of region_tissue_counts: {region_tissue_counts.dtype}")
        self.logger.info(f"First 5 entries of region_tissue_counts: {region_tissue_counts.head() if hasattr(region_tissue_counts, 'head') else list(region_tissue_counts)[:5]}")
        self.logger.info(f"Non-NaN count in region_tissue_counts: {region_tissue_counts.count() if hasattr(region_tissue_counts, 'count') else 'Uncountable'}")
        self.logger.info(f"Count for region_tissue_counts >= {min_region_count}: {(region_tissue_counts >= min_region_count).sum() if hasattr(region_tissue_counts, 'sum') else 'Uncomputable'}")

        if len(target_tissues) == 0:
            self.logger.warning(f"No tissue types with at least {min_region_count} counts.")
            return pd.DataFrame()
        
        # 解析結果の保存先
        results = []
        
        # 全サンプル数
        region_total = region_tissue_counts.sum()
        bg_total = background_tissue_counts.sum()
        
        # *** 修正: 全組織タイプの検定結果を格納するリスト ***
        all_tissue_results = []
        
        # *** 修正：すべての組織タイプに対して検定を実行 ***
        # background_tissue_counts.index を使用して全組織タイプをループする
        all_tissues = background_tissue_counts.index
        self.logger.info(f"Total number of tissue types: {len(all_tissues)}")
        
        # 各組織タイプごとに検定を実行
        for tissue in all_tissues:
            # 2×2分割表の作成
            region_count = region_tissue_counts.get(tissue, 0)
            bg_count = background_tissue_counts.get(tissue, 0)
            
            region_other = region_total - region_count
            bg_other = bg_total - bg_count
            
            # フォールドチェンジの計算
            region_freq = region_count / region_total if region_total > 0 else 0
            bg_freq = bg_count / bg_total if bg_total > 0 else 0
            fold_change = region_freq / bg_freq if bg_freq > 0 else float('inf')
            
            # 検定の実行
            if method == 'fisher':
                # Fisherの正確確率検定
                try:
                    contingency_table = [[region_count, region_other], [bg_count, bg_other]]
                    _, p_value = stats.fisher_exact(contingency_table)
                except Exception as e:
                    self.logger.warning(f"Fisher test failed (tissue: {tissue}): {e}")
                    p_value = 1.0
                    
            elif method == 'chi2':
                # カイ二乗検定
                try:
                    contingency_table = [[region_count, bg_count], [region_other, bg_other]]
                    chi2, p_value, _, _ = stats.chi2_contingency(contingency_table)
                except Exception as e:
                    self.logger.warning(f"Chi-squared test failed (tissue: {tissue}): {e}")
                    p_value = 1.0
                    
            else:
                raise ValueError(f"Unsupported test method: {method}")
            
            # 結果を保存（全組織）
            result = {
                'Tissue Type': tissue,
                'Region Count': region_count,
                'Total Count': bg_count,
                'Fold Change': fold_change,
                'p-value': p_value,
                'is_target': tissue in target_tissues  # ターゲット組織かどうかのフラグ
            }
            
            all_tissue_results.append(result)
            
            # ターゲット組織の場合は結果リストにも追加
            if tissue in target_tissues:
                results.append({
                    'Tissue Type': tissue,
                    'Region Count': region_count,
                    'Total Count': bg_count,
                    'Fold Change': fold_change,
                    'p-value': p_value
                })
        
        # 結果をDataFrameに変換
        if not results:
            self.logger.warning("No results for target tissue types. Returning an empty DataFrame.")
            return pd.DataFrame()
            
        # 全組織の結果をDataFrameに変換（多重検定補正用）
        all_tissue_df = pd.DataFrame(all_tissue_results)
        
        # ターゲット組織の結果をDataFrameに変換（最終的な返り値用）
        result_df = pd.DataFrame(results)
        
        # p値の多重検定補正（全組織のp値を使用）
        if len(all_tissue_df) > 0:
            all_p_values = all_tissue_df['p-value'].tolist()
            
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', category=RuntimeWarning)
                    # 全組織のp値を使って多重検定補正を実行
                    reject, all_adjusted_p_values, _, _ = multipletests(all_p_values, method=correction)
                
                # 補正後のp値を全組織のDataFrameに追加
                all_tissue_df['adjusted p-value'] = all_adjusted_p_values
                
                # ターゲット組織だけの結果DataFrame（result_df）に、補正後のp値をマージ
                # まず、組織タイプと補正後p値のマッピングを作成
                adjusted_p_dict = dict(zip(all_tissue_df['Tissue Type'], all_tissue_df['adjusted p-value']))
                
                # 結果DFに補正後p値を追加
                result_df['adjusted p-value'] = result_df['Tissue Type'].map(adjusted_p_dict)
                
                # 有意性フラグの設定
                result_df['Significant'] = (
                    (result_df['adjusted p-value'] <= pvalue_threshold) & 
                    (result_df['Fold Change'] >= fc_threshold)
                )
                
                # デプリートフラグの設定
                result_df['Depleted'] = (
                    (result_df['adjusted p-value'] <= pvalue_threshold) & 
                    (result_df['Fold Change'] <= 1/fc_threshold)
                )
                
                # 有意な結果の数をカウント
                significant_count = result_df['Significant'].sum()
                depleted_count = result_df['Depleted'].sum()
                
            except Exception as e:
                self.logger.warning(f"Multiple testing correction failed: {e}")
                # 多重検定補正に失敗した場合は、元のp値を使用する（フォールバック）
                result_df['adjusted p-value'] = result_df['p-value']
                result_df['Significant'] = (
                    (result_df['p-value'] <= pvalue_threshold) & 
                    (result_df['Fold Change'] >= fc_threshold)
                )
                result_df['Depleted'] = (
                    (result_df['p-value'] <= pvalue_threshold) & 
                    (result_df['Fold Change'] <= 1/fc_threshold)
                )
                
                # 有意な結果の数をカウント (フォールバック時)
                significant_count = result_df['Significant'].sum()
                depleted_count = result_df['Depleted'].sum()
        else:
            # 結果が空の場合は、カウントを0に設定
            significant_count = 0
            depleted_count = 0
        
        # 有意な結果の数を報告
        self.logger.info(f"Enrichment analysis completed: Total {len(result_df)} tissue types, {significant_count} enriched, {depleted_count} depleted")
        
        # p値でソートして返す
        return result_df.sort_values('p-value')










    def plot_tissue_enrichment(self, enrichment_result: pd.DataFrame, 
                             figsize: Tuple[int, int] = (10, 8),
                             top_n: int = 20,
                             show_depleted: bool = True) -> Optional[Any]:
        """
        Plot the tissue type enrichment analysis results.
        
        Parameters
        ----------
        enrichment_result : pd.DataFrame
            The result from test_tissue_enrichment().
        figsize : Tuple[int, int], optional
            The size of the generated figure.
        top_n : int, optional
            The number of top results to display.
        show_depleted : bool, optional
            Whether to display depleted tissue types as well.
            
        Returns
        -------
        Optional[Figure]
            The Matplotlib figure object generated, or None if the required libraries are not available.
        """
        if enrichment_result is None or len(enrichment_result) == 0:
            self.logger.warning("Enrichment result is empty. Skipping plot generation.")
            return None
            
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            fig = plt.figure(figsize=figsize)
            
            # エンリッチされた組織
            enriched = enrichment_result[enrichment_result['Significant'] == True].sort_values('p-value').head(top_n)
            
            # 1. 上位エンリッチ組織のプロット
            plt.subplot(2, 1, 1)
            if len(enriched) > 0:
                sns.barplot(x='Fold Change', y='Tissue Type', data=enriched)
                plt.title(f'Top {min(len(enriched), top_n)} Enriched Tissue Types')
                plt.xlabel('Fold Change')
            else:
                plt.title('No Significantly Enriched Tissue Types Found')
                plt.text(0.5, 0.5, 'No enriched tissue types', 
                        horizontalalignment='center', verticalalignment='center',
                        transform=plt.gca().transAxes)
            
            # 2. デプリートされた組織のプロット（オプション）
            if show_depleted:
                plt.subplot(2, 1, 2)
                depleted = enrichment_result[enrichment_result['Depleted'] == True].sort_values('Fold Change').head(top_n)
                
                if len(depleted) > 0:
                    sns.barplot(x='Fold Change', y='Tissue Type', data=depleted)
                    plt.title(f'Top {min(len(depleted), top_n)} Depleted Tissue Types')
                    plt.xlabel('Fold Change')
                else:
                    plt.title('No Significantly Depleted Tissue Types Found')
                    plt.text(0.5, 0.5, 'No depleted tissue types', 
                            horizontalalignment='center', verticalalignment='center',
                            transform=plt.gca().transAxes)
            
            plt.tight_layout()
            return fig
            
        except ImportError:
            self.logger.debug("matplotlib/seaborn is not installed, plot cannot be generated")
            return None

    def analyze_gene_distribution(self, gene_counts: pd.Series, 
                               top_n: int = 20) -> Dict[str, Any]:
        """
        Analyze gene distribution.
        
        Parameters
        ----------
        gene_counts : pd.Series
            A Series of genes and their occurrence counts.
        top_n : int, optional
            The number of top genes to return.
            
        Returns
        -------
        Dict[str, Any]
            A dictionary containing the analysis results (excluding figures).
        """
        self.logger.info("Starting gene distribution analysis")
        
        if gene_counts is None or len(gene_counts) == 0:
            self.logger.warning("Gene data is empty. Skipping analysis.")
            return {
                'gene_count': 0,
                'occurrence_stats': {},
                'top_genes': pd.Series()
            }
        
        # 基本統計の計算
        gene_count = len(gene_counts)
        total_occurrences = gene_counts.sum()
        
        # 出現回数の統計
        occurrence_stats = {
            'total': int(total_occurrences),
            'mean': float(gene_counts.mean()),
            'median': float(gene_counts.median()),
            'std': float(gene_counts.std()),
            'min': int(gene_counts.min()),
            'max': int(gene_counts.max())
        }
        
        # 上位の遺伝子
        top_genes = gene_counts.nlargest(top_n)
        
        # 結果を辞書としてまとめる
        result = {
            'gene_count': gene_count,
            'occurrence_stats': occurrence_stats,
            'top_genes': top_genes
        }
        
        self.logger.info("Gene distribution analysis completed")
        return result
    
    def plot_gene_distribution(self, analysis_result: Dict[str, Any],
                             figsize: Tuple[int, int] = (10, 8)) -> Optional[Any]:
        """
        Plot the gene distribution analysis results.
        
        Parameters
        ----------
        analysis_result : Dict[str, Any]
            The result from analyze_gene_distribution().
        figsize : Tuple[int, int], optional
            The size of the generated figure.
            
        Returns
        -------
        Optional[Figure]
            The Matplotlib figure object generated, or None if the required libraries are not available.
        """
        # プロット作成（matplotlibがある場合のみ）
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            fig = plt.figure(figsize=figsize)
            
            top_genes = analysis_result.get('top_genes')
            occurrence_stats = analysis_result.get('occurrence_stats', {})
            
            # 1. 上位遺伝子の棒グラフ
            plt.subplot(2, 1, 1)
            if top_genes is not None and len(top_genes) > 0:
                sns.barplot(x=top_genes.values, y=top_genes.index)
                plt.title(f'Top {len(top_genes)} Genes by Occurrence')
                plt.xlabel('Occurrence Count')
            else:
                plt.title('No Gene Data Available')
                plt.text(0.5, 0.5, 'No gene data', 
                        horizontalalignment='center', verticalalignment='center',
                        transform=plt.gca().transAxes)
            
            # 2. 出現回数の分布
            plt.subplot(2, 1, 2)
            plt.title('Gene Occurrence Distribution')
            plt.xlabel('Occurrence Count')
            plt.ylabel('Frequency')
            
            if 'mean' in occurrence_stats and 'std' in occurrence_stats:
                # ヒストグラムを直接描画するためのデータはないので、情報テキストを表示
                info_text = (
                    f"Total genes: {analysis_result.get('gene_count', 'N/A')}\n"
                    f"Total occurrences: {occurrence_stats.get('total', 'N/A')}\n"
                    f"Mean: {occurrence_stats.get('mean', 'N/A'):.2f}\n"
                    f"Median: {occurrence_stats.get('median', 'N/A'):.2f}\n"
                    f"Std Dev: {occurrence_stats.get('std', 'N/A'):.2f}\n"
                    f"Min: {occurrence_stats.get('min', 'N/A')}\n"
                    f"Max: {occurrence_stats.get('max', 'N/A')}"
                )
                plt.text(0.5, 0.5, info_text, 
                        horizontalalignment='center', verticalalignment='center',
                        transform=plt.gca().transAxes)
            else:
                plt.text(0.5, 0.5, 'Distribution data not available', 
                        horizontalalignment='center', verticalalignment='center',
                        transform=plt.gca().transAxes)
            
            plt.tight_layout()
            return fig
            
        except ImportError:
            self.logger.debug("matplotlib/seaborn is not installed, plot cannot be generated")
            return None
