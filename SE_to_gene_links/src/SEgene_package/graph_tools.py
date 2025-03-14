

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import japanize_matplotlib
import pathlib
import os
import collections

from itertools import islice
from pybedtools import BedTool

from typing import Optional, Tuple, List

import matplotlib.ticker as ticker

from SEgene_package.tools import return_bed_obj_tupple_from_rose_all_table
from SEgene_package.tools import sort_se_merge_SAMPLE_count_df
from SEgene_package.tools import sort_se_merge_SE_count_df






def make_sort_se_merge_SAMPLE_count_graph(se_count_df,graph_max_number):
    x_list=se_count_df["se_data"].to_list()
    y_list=se_count_df["sample_count"].to_list()
    
    max_number=graph_max_number
    title_str="merge-se sample count"

    x_list_graph=x_list[0:max_number]
    y_list_graph=y_list[0:max_number]
    plt.bar(x_list_graph,y_list_graph)
    plt.tight_layout()    

    plt.xticks(rotation=90)

    plt.gca().get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))

    plt.xlabel('merge SE')
    plt.ylabel('sample count')
    plt.title(title_str)   
    plt.show()    
    plt.clf()
    plt.close()


def make_sort_se_merge_SE_count_graph(
    se_count_df: pd.DataFrame,
    graph_max_number: int,
    sort_column: str = "se_count",
    title: Optional[str] = "merge-se SE count",
    xlabel: Optional[str] = "merge SE",
    ylabel: Optional[str] = "SE count",
    rotation: int = 90,
    save_svg: Optional[str] = None,
    save_tsv: Optional[str] = None,
    save_and_show: bool = False
) -> None:

    # Validate sort_column
    if sort_column not in se_count_df.columns:
        raise ValueError(f"Invalid sort_column: '{sort_column}'. Available columns: {list(se_count_df.columns)}")
    
    # Extract data
    x_list = se_count_df["se_data"].to_list()
    y_list = se_count_df[sort_column].to_list()

    # Limit the number of data points for display
    max_number = graph_max_number
    x_list_graph = x_list[:max_number]
    y_list_graph = y_list[:max_number]

    # Create a DataFrame for the displayed data
    display_data = pd.DataFrame({"se_data": x_list_graph, sort_column: y_list_graph})

    # Plot the bar chart
    plt.bar(x_list_graph, y_list_graph)

    # Adjust layout
    plt.tight_layout()

    # Rotate x-axis labels
    plt.xticks(rotation=rotation)

    # Ensure y-axis has integer ticks only if y is integer
    if pd.api.types.is_integer_dtype(se_count_df[sort_column]):
        plt.gca().get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))

    # Add labels and title
    plt.xlabel(xlabel)
    plt.ylabel(ylabel if ylabel else sort_column)
    plt.title(title)

    # Save as SVG if path is provided
    if save_svg:
        plt.savefig(save_svg, format="svg")
        print(f"Graph saved as SVG at {save_svg}.")

    # Show the graph if save_and_show is True or save_svg is None
    if save_and_show or not save_svg:
        plt.show()

    # Clear the plot to avoid overlap
    plt.clf()
    plt.close()

    # Save displayed data as TSV if path is provided
    if save_tsv:
        display_data.to_csv(save_tsv, sep="\t", index=False)
        print(f"Data saved as TSV at {save_tsv}.")



def se_merge_SAMPLE_count_graph(se_count_df,graph_max_number):

    sort_df=sort_se_merge_SAMPLE_count_df(se_count_df)

    make_sort_se_merge_SAMPLE_count_graph(sort_df, graph_max_number)

def se_merge_SE_count_graph(se_count_df,graph_max_number):

    sort_df=sort_se_merge_SE_count_df(se_count_df)

    make_sort_se_merge_SE_count_graph(sort_df, graph_max_number)




def create_ROSE_summary(
    bed_obj_bed_filter: BedTool, 
    path_enhancer_table: str,
    show_plot: bool = True,
    save_path: Optional[str] = None,
    save_text: Optional[str] = None,
    save_and_show: bool = False,
    include_pie_chart: bool = True,
    title: str = "ROSE SE Summary",
    color_all: str = "blue",
    color_super: str = "red",
    line_color: str = "green",
    line_style: str = "-",
    dpi: int = 300,
    figsize: Tuple[float, float] = (10, 8),
    output_format: str = "svg",
    # Font size parameters with defaults
    title_fontsize: int = 16,
    axis_label_fontsize: int = 14,
    legend_fontsize: int = 12,
    pie_percent_fontsize: int = 16,
    pie_legend_fontsize: int = 14,
    pie_title_fontsize: int = 18
) -> None:
    """
    Create a summary visualization of ROSE Super Enhancers analysis.
    
    Parameters:
    -----------
    bed_obj_bed_filter : BedTool
        BedTool object containing filter information
    path_enhancer_table : str
        Path to the enhancer table
    show_plot : bool, default=True
        Whether to display the plot
    save_path : Optional[str], default=None
        Path to save the figure (without extension)
    save_text : Optional[str], default=None
        Path to save the summary text
    save_and_show : bool, default=False
        Whether to save and show the plot
    include_pie_chart : bool, default=True
        Whether to include a pie chart
    title : str, default="ROSE SE Summary"
        Plot title
    color_all : str, default="blue"
        Color for all enhancers
    color_super : str, default="red"
        Color for super enhancers
    line_color : str, default="green"
        Color for threshold line
    line_style : str, default="-"
        Style for threshold line
    dpi : int, default=300
        Resolution in dots per inch for saved figure
    figsize : Tuple[float, float], default=(10, 8)
        Figure size (width, height) in inches
    output_format : str, default="svg"
        Output format for the figure: 'svg', 'png', 'jpg', 'pdf', 'eps'
    title_fontsize : int, default=16
        Font size for the main plot title
    axis_label_fontsize : int, default=14
        Font size for axis labels
    legend_fontsize : int, default=12
        Font size for the legend
    pie_percent_fontsize : int, default=16
        Font size for pie chart percentage labels
    pie_legend_fontsize : int, default=14
        Font size for pie chart legend
    pie_title_fontsize : int, default=18
        Font size for pie chart title
    """
    
    # Load enhancer table
    bed_table_all, bed_table_super, _ = return_bed_obj_tupple_from_rose_all_table(path_enhancer_table)
    len_super = len(bed_table_super)

    # Check coverage for all enhancers
    bed_all_covercheck = bed_table_all.intersect(bed_obj_bed_filter, c=True)
    df_all_covercheck = bed_all_covercheck.to_dataframe()

    # Extract covered super enhancers
    bed_super_covercheck = bed_table_super.intersect(bed_obj_bed_filter, c=True)
    df_super_covercheck = bed_super_covercheck.to_dataframe()

    # Add flag and rename columns
    df_super_covercheck.rename(columns={"thickStart": "count"}, inplace=True)
    df_super_covercheck["Filtered"] = df_super_covercheck["count"] > 0

    if "strand" in df_super_covercheck.columns:
        df_super_covercheck.drop(columns=["strand"], inplace=True)

    # Extract filtered entries
    df_super_covered = df_super_covercheck[df_super_covercheck["Filtered"]]
    len_super_covered = len(df_super_covered)

    # Print summary messages
    print("Output Graph: Super Table - Coverage Ratio of SEs by Peak to Gene")
    print(f"Total SEs: {len_super}")
    print(f"SEs Not Covered: {len_super - len_super_covered}")
    print(f"SEs Covered: {len_super_covered}")

    # Save summary text
    if save_text:
        with open(save_text, "w") as f:
            f.write(f"Summary for {path_enhancer_table}\n")
            f.write("=" * 60 + "\n")
            f.write(f"Total SEs: {len_super}\n")
            f.write(f"SEs Not Covered: {len_super - len_super_covered}\n")
            f.write(f"SEs Covered: {len_super_covered}\n")
            f.write("=" * 60 + "\n\n")
            f.write("Enhanced Table with Filter Flag:\n")
            f.write(df_super_covercheck.to_string(index=True, header=True))
            f.write("\n")
        print(f"Saved summary text with filter flags to {save_text}.")

    # Set figure size
    plt.figure(figsize=figsize)
    
    # Main graph plotting
    plt.plot(
        df_all_covercheck.index[::-1], df_all_covercheck["score"][::-1], 
        '.', c=color_all, label='All Enhancers'
    )
    plt.plot(
        df_super_covered.index, df_super_covered["score"], 
        ".", c=color_super, label='Filtered Super Enhancers'
    )
    plt.axvline(
        x=len_super - 0.5, color=line_color, linestyle=line_style, label='Threshold'
    )
    plt.gca().invert_xaxis()
    
    # Set font sizes for main plot
    plt.xlabel('SE rank', fontsize=axis_label_fontsize)
    plt.ylabel('Enhancer signal', fontsize=axis_label_fontsize)
    plt.title(title, fontsize=title_fontsize)
    plt.legend(fancybox=False, frameon=True, fontsize=legend_fontsize)  

    # Save figure if path is provided
    if save_path:
        # Standardize output format
        output_format = output_format.lower().strip('.')
        
        # Create full file path with extension
        figure_path = f"{save_path}.{output_format}"
        
        # Save the figure with the specified dpi
        plt.savefig(figure_path, format=output_format, dpi=dpi, bbox_inches='tight')
        print(f"Saved figure to {figure_path} with DPI={dpi}.")
    
    if show_plot or save_and_show:
        plt.show()

    plt.clf()
    plt.close()

    # Pie chart generation and saving
    if include_pie_chart:
        # Create new figure for pie chart with the specified figure size
        plt.figure(figsize=figsize)
        
        pie_label_covered_se = ["SE_COVERED", "SE_NOT_COVERED"]
        pie_super_and_covered_se = np.array([len_super_covered, len_super - len_super_covered])

        # Create pie chart with custom font sizes
        wedges, texts, autotexts = plt.pie(
            pie_super_and_covered_se, 
            labels=None,  # Disable default labels to use custom legend
            counterclock=False, 
            startangle=90, 
            autopct="%1.1f%%",
            textprops={'fontsize': pie_percent_fontsize, 'fontweight': 'bold'},  # Set percentage font size
            wedgeprops={'edgecolor': 'w', 'linewidth': 1.5},  # Add wedge borders
            pctdistance=0.65  # Adjust percentage label position
        )
        
        # Add custom legend with count values
        plt.legend(
            wedges, 
            [f"{pie_label_covered_se[0]} ({len_super_covered})", 
             f"{pie_label_covered_se[1]} ({len_super - len_super_covered})"],
            loc="upper right",
            bbox_to_anchor=(1.1, 1),
            fontsize=pie_legend_fontsize,
            frameon=True,
            framealpha=0.7,
            edgecolor='black'
        )

        # Set pie chart title with custom font size
        pie_chart_title = f"{title} - Pie Chart"
        plt.title(pie_chart_title, fontsize=pie_title_fontsize)
        
        # Keep the pie chart circular
        plt.axis('equal')

        if save_path:
            # Create a pie chart filename
            pie_path = f"{save_path}_pie.{output_format}"
            plt.savefig(pie_path, format=output_format, dpi=dpi, bbox_inches='tight')
            print(f"Saved pie chart to {pie_path} with DPI={dpi}.")

        if show_plot or save_and_show:
            plt.show()

        plt.clf()
        plt.close()
