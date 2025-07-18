# sedb_tools/visualization/plots.py
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import os
from typing import Optional, Dict, List, Tuple, Any
import logging

class SEdbPlotGenerator:
    """
    基本的なプロット生成機能を提供するクラス。
    様々な種類のプロットを生成するためのユーティリティ関数を含む。
    """
    
    def __init__(self, logger=None):
        """
        プロット生成器を初期化
        
        Parameters:
        logger (Logger): オプションのロガーインスタンス。Noneの場合、基本的なコンソールロガーを作成。
        """
        self.logger = logger or self._setup_default_logger()
        self.logger.info("Initializing SEdbPlotGenerator")
    
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
    
    def plot_simple_pie(self, data: pd.DataFrame, category: str, count_col: str = 'Count',
                    title: Optional[str] = None, limit: int = 10,
                    figsize: tuple = (12, 8),
                    colormap: str = 'tab20',
                    title_y: float = 1.1,  # タイトル位置を上に調整
                    save_dir: Optional[str] = None, 
                    save_filename: Optional[str] = None,
                    save_svg: bool = False, 
                    save_png: bool = False, 
                    save_pdf: bool = False, 
                    save_eps: bool = False,
                    save_dpi: int = 600, 
                    save_transparent: bool = False) -> plt.Figure:
        """
        Create a simple, clear pie chart
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame containing category and count columns
        category : str
            Category column name to display in pie chart
        count_col : str, default='Count'
            Column name for the count values
        title : str, optional
            Chart title
        limit : int, default=10
            Number of top categories to display
        figsize : tuple, default=(12, 8)
            Figure size
        colormap : str, default='tab20'
            Matplotlib colormap to use ('tab20', 'tab20b', 'tab20c', 'Set3', etc.)
        title_y : float, default=1.1
            Y position of the title (higher value moves title up)
        save_dir : str, optional
            Directory to save output files. If None, current directory is used
        save_filename : str, optional
            Base filename for saved files (without extension). If None, auto-generated
        save_svg : bool, optional
            Whether to save as SVG format (default: False)
        save_png : bool, optional
            Whether to save as PNG format (default: False)
        save_pdf : bool, optional
            Whether to save as PDF format (default: False)
        save_eps : bool, optional
            Whether to save as EPS format (default: False)
        save_dpi : int, optional
            DPI (resolution) for raster formats like PNG (default: 600)
        save_transparent : bool, optional
            Whether to save with transparent background (default: False)
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        # Check if required columns exist
        if category not in data.columns:
            self.logger.warning(f"Column '{category}' does not exist in the dataframe")
            return None
        
        if count_col not in data.columns:
            self.logger.warning(f"Column '{count_col}' does not exist in the dataframe")
            return None
        
        # Sort data by count (descending)
        sorted_data = data.sort_values(count_col, ascending=False)
        total_categories = len(sorted_data)
        
        # Set title
        if title is None:
            title = f'{category} Distribution'
            
        # Divide into top N and others
        if total_categories > limit:
            top_data = sorted_data.head(limit)
            others_count = sorted_data.iloc[limit:][count_col].sum()
            
            # Create a DataFrame with 'Others' row
            others_row = pd.DataFrame({category: ['Others'], count_col: [others_count]})
            plot_data = pd.concat([top_data, others_row], ignore_index=True)
            title_suffix = f" (Top {limit}, other {total_categories-limit} types grouped as 'Others')"
        else:
            plot_data = sorted_data
            title_suffix = ""
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # Get labels and values
        labels = plot_data[category]
        values = plot_data[count_col]
        
        # Generate distinct colors for pie segments
        cmap = plt.cm.get_cmap(colormap)
        num_colors = len(plot_data)
        
        # 重要な修正: Othersカテゴリと他のカテゴリの色が被らないようにする
        colors = []
        others_idx = None
        
        # まずOthersのインデックスを見つける（存在する場合）
        for i, label in enumerate(labels):
            if label == 'Others':
                others_idx = i
                break
        
        # カラーマップから色を生成
        for i in range(num_colors):
            if i == others_idx:
                # Othersには特別な色（灰色）を使用
                colors.append('grey')
            else:
                # カラーマップからインデックスを計算（Othersがある場合はずらす）
                if others_idx is not None and i > others_idx:
                    # Othersの後のインデックスは-1して詰める
                    colors.append(cmap((i-1) % cmap.N))
                else:
                    colors.append(cmap(i % cmap.N))
        
        # Create pie chart (maintaining original design)
        plt.pie(
            values, 
            labels=labels,
            autopct='%1.1f%%',
            startangle=90,
            shadow=False,
            wedgeprops={'edgecolor': 'none'},  # Remove edge lines
            textprops={'fontsize': 10},
            colors=colors  # 修正した色のリストを使用
        )
        
        # タイトル位置を上に調整し、余白を確保
        plt.title(f'{title}{title_suffix}', y=title_y)
        plt.axis('equal')
        
        # タイトル用に十分なスペースを確保
        plt.subplots_adjust(top=0.85)  # これによりタイトル用のスペースを確保
        plt.tight_layout(pad=1.5)  # パディングを少し増やす
        
        # Save figure in requested formats
        if save_svg or save_png or save_pdf or save_eps:
            import os
            
            # Create directory if needed
            if save_dir is None:
                save_dir = os.getcwd()  # Current directory
            elif not os.path.exists(save_dir):
                os.makedirs(save_dir, exist_ok=True)
                self.logger.info(f"Created output directory: {save_dir}")
            
            # Generate base filename if not provided
            if save_filename is None:
                save_filename = f"pie_chart_{category.lower().replace(' ', '_')}"
            
            # Create full path without extension
            base_path = os.path.join(save_dir, save_filename)
            
            # Save in each requested format
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved pie chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        # Display detailed information
        total_count = values.sum()
        print(f"\n{category} distribution:")
        print(f"Total {total_categories} types, showing top {min(limit, total_categories)}")
        
        for i, (label, count) in enumerate(zip(labels, values)):
            if label != 'Others':
                print(f"  {label}: {int(count)} samples ({count/total_count*100:.1f}%)")
        
        if 'Others' in labels.values:
            others_idx = labels.tolist().index('Others')
            others_count = values.iloc[others_idx]
            print(f"  Others ({total_categories - limit} types): {int(others_count)} samples ({others_count/total_count*100:.1f}%)")
        
        return fig

    def plot_sorted_bar(self, data: pd.DataFrame, category: str, count_col: str = 'Count', 
                    title: Optional[str] = None, limit: int = 15,
                    horizontal: bool = True, color: str = 'skyblue', 
                    figsize: tuple = (12, 8),
                    save_dir: Optional[str] = None, 
                    save_filename: Optional[str] = None,
                    save_svg: bool = False, 
                    save_png: bool = False, 
                    save_pdf: bool = False, 
                    save_eps: bool = False,
                    save_dpi: int = 600, 
                    save_transparent: bool = False) -> plt.Figure:
        """
        Helper method to plot sorted bar chart for distribution visualization.
        
        Parameters:
        -----------
        data : pandas.DataFrame
            DataFrame containing category and count columns
        category : str
            Column name for the category
        count_col : str, optional
            Column name for the count values (default: 'Count')
        title : str, optional
            Plot title
        limit : int, optional
            Number of top categories to display
        horizontal : bool, optional
            Whether to use horizontal bar chart (default: True)
        color : str, optional
            Bar color
        figsize : tuple, optional
            Figure size
        save_dir : str, optional
            Directory to save the figure
        save_filename : str, optional
            Base filename for saved figure (without extension)
        save_svg, save_png, save_pdf, save_eps : bool, optional
            Whether to save in respective formats
        save_dpi : int, optional
            DPI for saved images
        save_transparent : bool, optional
            Whether to save with transparent background
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure
        """
        # Sort data and limit to top N categories
        if data is None or len(data) == 0:
            self.logger.warning(f"No data provided for plotting {category} distribution")
            return None
        
        sorted_data = data.sort_values(count_col, ascending=False)
        if limit > 0 and len(sorted_data) > limit:
            plot_data = sorted_data.head(limit)
        else:
            plot_data = sorted_data
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        if horizontal:
            # Horizontal bar chart (better for many categories with long names)
            # Reverse order for better display (highest at top)
            plot_data = plot_data.iloc[::-1].reset_index(drop=True)
            bars = ax.barh(plot_data[category], plot_data[count_col], color=color, alpha=0.8)
            
            # Add count values to the right of bars
            for i, bar in enumerate(bars):
                ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                    f"{int(bar.get_width())}", va='center')
            
            ax.set_xlabel('Count')
            ax.set_ylabel(category)
        else:
            # Vertical bar chart
            bars = ax.bar(plot_data[category], plot_data[count_col], color=color, alpha=0.8)
            
            # Add count values above bars
            for i, bar in enumerate(bars):
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f"{int(bar.get_height())}", ha='center', va='bottom')
            
            ax.set_xlabel(category)
            ax.set_ylabel('Count')
            plt.xticks(rotation=45, ha='right')
        
        # Set title
        if title:
            ax.set_title(title)
        else:
            ax.set_title(f'{category} Distribution')
        
        plt.tight_layout()
        
        # Save figure if requested
        if save_dir and save_filename:
            import os
            os.makedirs(save_dir, exist_ok=True)
            base_path = os.path.join(save_dir, save_filename)
            
            saved_files = []
            
            if save_svg:
                svg_path = f"{base_path}.svg"
                plt.savefig(svg_path, format='svg', dpi=save_dpi, bbox_inches='tight', 
                        transparent=save_transparent)
                saved_files.append(f"SVG: {svg_path}")
            
            if save_png:
                png_path = f"{base_path}.png"
                plt.savefig(png_path, format='png', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PNG: {png_path}")
            
            if save_pdf:
                pdf_path = f"{base_path}.pdf"
                plt.savefig(pdf_path, format='pdf', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"PDF: {pdf_path}")
            
            if save_eps:
                eps_path = f"{base_path}.eps"
                plt.savefig(eps_path, format='eps', dpi=save_dpi, bbox_inches='tight',
                        transparent=save_transparent)
                saved_files.append(f"EPS: {eps_path}")
            
            # Log saved files
            if saved_files:
                self.logger.info(f"Saved bar chart in the following formats:")
                for file_info in saved_files:
                    self.logger.info(f"  - {file_info}")
        
        return fig