import os
import sys
import numpy as np
import pandas as pd
import pathlib
import pickle 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker 
import json
import random

import datetime


import seaborn as sns

from scipy import stats

from pathlib import Path

from IPython.display import display


from logging import getLogger, Formatter, StreamHandler, FileHandler, DEBUG, INFO, WARNING, Logger
from logging.handlers import RotatingFileHandler

import logging
import io
import sys
import contextlib
from io import StringIO



from typing import Optional, List, Dict, Tuple, Any

from pybedtools import BedTool
import networkx as nx
from networkx.readwrite import json_graph

from networkx.drawing.nx_agraph import to_agraph
from SEgene_package.tools import (
    bed_filter_from_p2g_output,
    return_all_se_concat_full_df,
    return_merge_se_count_df_full_edit,
    count_genes_from_bed,
    p2gl_path_to_filter_df,
    return_merge_se_count_df_full_with_info,
    find_gene_linked_super_enhancers
)

from SEgene_package.graph_tools import create_ROSE_summary
from SEgene_package.graph_tools import make_sort_se_merge_SE_count_graph


from SEgene_package.networkx_tools import create_3_layer_DiGraph_mergeSE_to_gene_custom

from SEgene_package.networkx_tools import csv_path_to_df_RNA
from SEgene_package.networkx_tools import rna_data_to_dic_metadata

from SEgene_package.region_visualize_tools import plot_stacked_reads_bed

import glob
from SEgene_package.tools import find_gene_linked_super_enhancers_in_directory


def sort_se_merge_count_by_column(
    merge_se_count_df: pd.DataFrame, 
    column: str = "se_count", 
    ascending: bool = False
) -> pd.DataFrame:

    if column not in merge_se_count_df.columns:
        raise ValueError(f"Invalid column name: {column}")

    sorted_df = merge_se_count_df.copy().sort_values(by=column, ascending=ascending)
    return sorted_df






class LoggerManager:

    def __init__(
        self, 
        logger_name: str, 
        enable_console_logging: bool = False, 
        console_verbose: bool = False,
        enable_file_logging: bool = False,
        log_dir: str = './logs',
        log_file_name: Optional[str] = None
    ):

        self.logger = getLogger(logger_name)
        self.logger.setLevel(DEBUG)  


        formatter = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')


        if enable_console_logging:
            console_handler = StreamHandler(sys.stdout)
            console_level = INFO if console_verbose else WARNING
            console_handler.setLevel(console_level)
            console_handler.setFormatter(formatter)

            if not any(isinstance(handler, StreamHandler) for handler in self.logger.handlers):
                self.logger.addHandler(console_handler)


        if enable_file_logging:

            if not os.path.exists(log_dir):
                os.makedirs(log_dir, exist_ok=True)

                if enable_console_logging:
                    self.logger.info(f"Log directory '{log_dir}' was created.")


            if log_file_name is None:
                log_file_name = f"{logger_name}.log"
            log_file_path = os.path.join(log_dir, log_file_name)


            file_handler = RotatingFileHandler(
                log_file_path, 
                maxBytes=10**6,       
                backupCount=5        
            )
            file_handler.setLevel(DEBUG)  
            file_handler.setFormatter(formatter)

            if not any(isinstance(handler, RotatingFileHandler) and handler.baseFilename == os.path.abspath(log_file_path) for handler in self.logger.handlers):
                self.logger.addHandler(file_handler)


        self.logger.propagate = False

    def get_logger(self) -> Logger:

        return self.logger




class SEgeneJupyter:
    """
    SEGENE Analysis Class for SE Merging and Counting
    """

    def __init__(
        self, 
        rose_file: Optional[str] = None,
        rose_files_dir: Optional[str] = None,
        p2g_file: Optional[str] = None,
        rna_info_file: Optional[str] = None,
        FDR: float = 0.05,  
        r: float = 0.5,        
        auto_filter_p2g: bool = True,  
        verbose: bool = False,         
        enable_console_logging: bool = False,  
        console_verbose: bool = False,         
        enable_file_logging: bool = False,     
        log_dir: str = './logs',               
        log_file_name: Optional[str] = None,     
        log_relative_paths: bool = True  
    ) -> None:
        
        self.verbose = verbose 

        
        self.logger_manager = LoggerManager(
            logger_name=self.__class__.__name__,
            enable_console_logging=enable_console_logging,
            console_verbose=console_verbose,
            enable_file_logging=enable_file_logging,
            log_dir=log_dir,
            log_file_name=log_file_name
        )
        self.logger = self.logger_manager.get_logger() 

        self.logger.info("Initializing SEGENE_jupyter instance.")


        
        self.rose_file = rose_file
        self.rose_files_dir = rose_files_dir
        self.p2g_file = p2g_file
        self.rna_info_file = rna_info_file

        
        self._bed_filter: Optional[BedTool] = None
        self._temp_concat_full: Optional[pd.DataFrame] = None
        self._temp_full_df_edit: Optional[pd.DataFrame] = None
        self._sorted_df: Optional[pd.DataFrame] = None

        
        self._FDR: float = FDR  
        self._r: float = r        

        
        self._network_FDR: float = FDR  
        self._network_r: float = r        

        
        self._bed_after_threshold_network: Optional[BedTool] = None  
        self._bed_before_threshold_network: Optional[BedTool] = None  
        self._bed_df_concat_sort_network: Optional[BedTool] = None  

        
        self._df_bed_mergeSE_to_gene_network: Optional[pd.DataFrame] = None  

        
        self._DG_network: Optional[nx.DiGraph] = None  

        
        self._selected_se_index: Optional[int] = None
        self._selected_se_region: Optional[str] = None


        self.log_relative_paths = log_relative_paths  


        
        self._input_se_file_list = self._generate_se_file_list()

        self._selected_subgraph = None
        
        
        
        
        
        
        
        
        
        



        
        self.default_settings = {
            'layout': 'spring',
            'node_size': 300,
            'with_labels': True,
            'font_size': 10,
            'figsize': (12, 9),
            'verbose': False,
            'color_map': {
                'gene': 'red',
                'mergeSE': 'blue',
                'enhancer': 'green',
                'default': 'gray'
            }
        }
        
        
        self.settings_3layer = {
            
            'seed': self._generate_seed('3layer')
        }
        
        self.settings_2layer = {
            'edge_scale': 1.0,
            'node_color_map': {'mergeSE': 'blue', 'gene': 'red'},
            'seed': self._generate_seed('2layer')
            
        }



        
        
        
        
        


        
        if self.rose_file is None and self.rose_files_dir is not None and len(self._input_se_file_list) > 0:
            self.rose_file = str(self._input_se_file_list[0])
            if self.verbose:
                relative_rose_file = self._get_log_path(self.rose_file)
                self.logger.info(f"rose_file not specified. Automatically selected the first file from rose_files_dir: {relative_rose_file}")


        
        if self.verbose:
            self.datafile_exists()  


        
        if auto_filter_p2g:
            self.filter_p2g_file(self._FDR, self._r)  


        
        self._current_sort_key: Optional[str] = None  


        # Variables for directory search results
        self._se_directory_results = {}  # Dictionary to store directory search results by gene symbol
        self._last_searched_gene_directory = None  # Last searched gene symbol in directory

        self.logger.info("SEGENE_jupyter instance initialized successfully.")

        #for search gene 
        self._gene_enhancer_results = {}  # Dictionary to store search results by gene symbol
        self._last_searched_gene = None   # Last searched gene symbol
        self._gene_bed_filter = None      # BedTool object for the last searched gene

    

    def __enter__(self):
        self.logger.info("Entering SEGENE_jupyter context.")  
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:
            self.logger.error(f"An exception occurred: {exc_val}", exc_info=True)  
        self.logger.info("Exiting SEGENE_jupyter context.")  

    def __del__(self):
        self.logger.info("SEGENE_jupyter instance is being deleted. Cleaning up resources.")  
        try:
            handlers = self.logger.handlers[:]
            for handler in handlers:
                handler.close()
                self.logger.removeHandler(handler)
            self.logger.info("Cleanup completed successfully.")  
        except Exception as e:
            self.logger.exception(f"Error during cleanup in __del__: {e}")










    
    
    
    
    


    def _generate_seed(self, graph_type: str) -> int:
        self.logger.info(f"Generating random seed for {graph_type} graph.")  
        try:
            seed = random.randint(0, 1000000)
            self.logger.info(f"Generated random seed for {graph_type}: {seed}")  
            return seed
        except Exception as e:
            self.logger.exception(f"Error in _generate_seed: {e}")
            raise
        self.logger.info(f"Random seed for {graph_type} graph generated successfully.")  




    def _set_seed(self, graph_type: str, seed: Optional[int] = None) -> int:

        if graph_type == '3layer':
            if seed is not None:
                self.settings_3layer['seed'] = seed
                self.logger.info(f"Set random seed for 3layer graph to: {self.settings_3layer['seed']}")  
            else:
                self.logger.info(f"Using existing seed for 3layer graph: {self.settings_3layer['seed']}")  
            return self.settings_3layer['seed']
        elif graph_type == '2layer':
            if seed is not None:
                self.settings_2layer['seed'] = seed
                self.logger.info(f"Set random seed for 2layer graph to: {self.settings_2layer['seed']}")  
            else:
                self.logger.info(f"Using existing seed for 2layer graph: {self.settings_2layer['seed']}")  
            return self.settings_2layer['seed']
        else:
            self.logger.error(f"Unsupported graph type for seed setting: {graph_type}")  
            raise ValueError(f"Unsupported graph type: {graph_type}")



    def _prepare_graph(self, subgraph_id: Optional[int] = None, layers: int = 3) -> Optional[nx.DiGraph]:

        self.logger.info(f"Preparing graph with subgraph_id={subgraph_id} and layers={layers}.")  
        try:
            
            stats = self.get_subgraph_statistics(subgraph_id=subgraph_id)
            subgraph = self._selected_subgraph if subgraph_id is None else self._DG_network.subgraph([
                n for n, d in self._DG_network.nodes(data=True) if d.get('subgraph_id') == subgraph_id
            ]).copy()
            
            if subgraph is None or subgraph.number_of_nodes() == 0:
                self.logger.warning(f"Subgraph ID {subgraph_id} has no nodes to display.")
                return None
            
            
            filtered_graph = nx.DiGraph()
            
            
            if layers == 3:
                
                for node, attr in subgraph.nodes(data=True):
                    filtered_graph.add_node(node, **attr)
            elif layers == 2:
                
                for node, attr in subgraph.nodes(data=True):
                    node_type = attr.get('node_type', 'default')
                    if node_type in ['mergeSE', 'gene']:
                        filtered_graph.add_node(node, **attr)
            else:
                self.logger.error(f"Unsupported number of layers: {layers}. Supported: 2, 3.")
                return None
            
            
            if layers == 3:
                
                filtered_graph.add_edges_from(subgraph.edges(data=True))
            elif layers == 2:
                
                edge_weights = {}
                for u, v, attr in subgraph.edges(data=True):
                    u_type = subgraph.nodes[u].get('node_type', 'default')
                    v_type = subgraph.nodes[v].get('node_type', 'default')
                    
                    if u_type == 'mergeSE' and v_type == 'enhancer':
                        for _, target_gene, _ in subgraph.out_edges(v, data=True):
                            target_type = subgraph.nodes[target_gene].get('node_type', 'default')
                            if target_type == 'gene':
                                key = (u, target_gene)
                                edge_weights[key] = edge_weights.get(key, 0) + 1
                    elif u_type == 'mergeSE' and v_type == 'gene':
                        
                        key = (u, v)
                        edge_weights[key] = edge_weights.get(key, 0) + 1
                
                for (u, v), weight in edge_weights.items():
                    filtered_graph.add_edge(u, v, weight=weight)
                    self.logger.debug(f"Added edge from {u} to {v} with weight {weight}")
            
        except Exception as e:
            self.logger.exception(f"Error in _prepare_graph: {e}")
            raise
        self.logger.info("Graph preparation completed successfully.")  
        return filtered_graph




    def _draw_graph(self, 
                graph: nx.DiGraph,
                layout: str,
                node_size: int,
                with_labels: bool,
                font_size: int,
                figsize: Tuple[int, int],
                node_colors: list,
                edge_colors: list,
                edge_widths: list,
                seed: Optional[int],
                verbose: bool) -> None:

        self.logger.info(f"Drawing graph with layout={layout}, node_size={node_size}, with_labels={with_labels}, font_size={font_size}, figsize={figsize}, seed={seed}.")  
        try:
            
            if layout == 'spring':
                pos = nx.spring_layout(graph, seed=seed)
            elif layout == 'circular':
                pos = nx.circular_layout(graph)
            elif layout == 'shell':
                pos = nx.shell_layout(graph)
            elif layout == 'kamada_kawai':
                pos = nx.kamada_kawai_layout(graph, seed=seed)
            else:
                self.logger.warning(f"Unsupported layout '{layout}'. Falling back to 'spring'.")
                pos = nx.spring_layout(graph, seed=seed)
            
            if verbose:
                self.logger.debug(f"Node positions: {pos}")
            
            
            plt.figure(figsize=figsize)
            
            
            nx.draw_networkx_nodes(graph, pos, node_size=node_size, node_color=node_colors, alpha=0.8)
            
            
            if edge_widths:
                nx.draw_networkx_edges(graph, pos, edge_color=edge_colors, width=edge_widths, arrowstyle='->', arrowsize=15, alpha=0.5)
            else:
                self.logger.warning("No edges to draw.")
            
            
            if with_labels:
                nx.draw_networkx_labels(graph, pos, font_size=font_size, font_color='black')
            
            
            plt.axis('off')
            plt.tight_layout()
            plt.show()
        except Exception as e:
            self.logger.exception(f"Error in _draw_graph: {e}")
            raise
        self.logger.info("Graph drawn successfully.")  



    def display_subgraph(
        self,
        subgraph_id: Optional[int] = None,
        layout: Optional[str] = None,
        node_size: Optional[int] = None,
        with_labels: Optional[bool] = None,
        font_size: Optional[int] = None,
        figsize: Optional[Tuple[int, int]] = None,
        verbose: Optional[bool] = None,
        seed: Optional[int] = None
    ) -> None:

        self.logger.info(f"Displaying subgraph with subgraph_id={subgraph_id}.")  

        try:
            
            current_seed = self._set_seed('3layer', seed)
            
            
            settings = self.default_settings.copy()
            if layout is not None:
                settings['layout'] = layout
            if node_size is not None:
                settings['node_size'] = node_size
            if with_labels is not None:
                settings['with_labels'] = with_labels
            if font_size is not None:
                settings['font_size'] = font_size
            if figsize is not None:
                settings['figsize'] = figsize
            if verbose is not None:
                settings['verbose'] = verbose
            
            
            graph = self._prepare_graph(subgraph_id=subgraph_id, layers=3)
            if graph is None:
                self.logger.warning("No graph to display.")
                return
            
            
            node_colors = []
            color_map = self.default_settings.get('color_map', {
                'gene': 'red',
                'mergeSE': 'blue',
                'enhancer': 'green',
                'default': 'gray'
            })
            for node, attr in graph.nodes(data=True):
                node_type = attr.get('node_type', 'default')
                node_colors.append(color_map.get(node_type, color_map['default']))
            
            
            edge_colors = ['gray' for _ in graph.edges()]
            edge_widths = [1 for _ in graph.edges()]
            
            if settings['verbose']:
                self.logger.info(f"Filtered Graph has {graph.number_of_edges()} edges.")
                for u, v, attr in graph.edges(data=True):
                    self.logger.debug(f"Edge from {u} to {v} with attributes {attr}")
            
            
            self._draw_graph(
                graph=graph,
                layout=settings['layout'],
                node_size=settings['node_size'],
                with_labels=settings['with_labels'],
                font_size=settings['font_size'],
                figsize=settings['figsize'],
                node_colors=node_colors,
                edge_colors=edge_colors,
                edge_widths=edge_widths,
                seed=current_seed,
                verbose=settings['verbose']
            )
            
            if settings['verbose']:
                stats = self.get_subgraph_statistics(subgraph_id=subgraph_id)
                self.logger.info(f"Displayed Subgraph ID {stats['subgraph_id']}:")
                self.logger.info(f"  Number of Nodes: {stats['number_of_nodes']}")
                self.logger.info(f"  Number of Edges: {stats['number_of_edges']}")
                self.logger.info("  Node Types:")
                for nt, count in stats['node_types'].items():
                    self.logger.info(f"    {nt}: {count}")
                self.logger.info("  Edge Types:")
                for et, count in stats['edge_types'].items():
                    self.logger.info(f"    {et}: {count}")
        except Exception as e:
            self.logger.exception(f"Error in display_subgraph: {e}")
            raise
        self.logger.info("Subgraph displayed successfully.")  





    def display_two_layer_subgraph(
        self,
        subgraph_id: Optional[int] = None,
        layout: Optional[str] = None,
        node_size: Optional[int] = None,
        with_labels: Optional[bool] = None,
        font_size: Optional[int] = None,
        figsize: Optional[Tuple[int, int]] = None,
        edge_scale: Optional[float] = None,
        node_color_map: Optional[Dict[str, str]] = None,
        verbose: Optional[bool] = None,
        seed: Optional[int] = None
    ) -> None:

        self.logger.info(f"Displaying two-layer subgraph with subgraph_id={subgraph_id}.")  
        try:
            
            current_seed = self._set_seed('2layer', seed)
            
            
            settings = self.default_settings.copy()
            
            settings.update({
                'edge_scale': self.settings_2layer.get('edge_scale', 1.0),
                'node_color_map': self.settings_2layer.get('node_color_map', {'mergeSE': 'blue', 'gene': 'red'}),
                'verbose': False  
            })
            if layout is not None:
                settings['layout'] = layout
            if node_size is not None:
                settings['node_size'] = node_size
            if with_labels is not None:
                settings['with_labels'] = with_labels
            if font_size is not None:
                settings['font_size'] = font_size
            if figsize is not None:
                settings['figsize'] = figsize
            if verbose is not None:
                settings['verbose'] = verbose
            if edge_scale is not None:
                settings['edge_scale'] = edge_scale
            if node_color_map is not None:
                settings['node_color_map'] = node_color_map
            
            
            graph = self._prepare_graph(subgraph_id=subgraph_id, layers=2)
            if graph is None:
                self.logger.warning("No graph to display.")
                return
            
            
            node_colors = []
            color_map = settings['node_color_map']
            for node, attr in graph.nodes(data=True):
                node_type = attr.get('node_type', 'default')
                node_colors.append(color_map.get(node_type, 'gray'))
            
            
            edge_colors = ['black' for _ in graph.edges()]
            edge_widths = [max(attr.get('weight', 1) * settings['edge_scale'], 1) for _, _, attr in graph.edges(data=True)]
            
            if settings['verbose']:
                self.logger.info(f"Filtered Graph has {graph.number_of_edges()} edges.")
                for u, v, attr in graph.edges(data=True):
                    self.logger.debug(f"Edge from {u} to {v} with weight {attr.get('weight', 1)}")
            
            
            self._draw_graph(
                graph=graph,
                layout=settings['layout'],
                node_size=settings['node_size'],
                with_labels=settings['with_labels'],
                font_size=settings['font_size'],
                figsize=settings['figsize'],
                node_colors=node_colors,
                edge_colors=edge_colors,
                edge_widths=edge_widths,
                seed=current_seed,
                verbose=settings['verbose']
            )
            
            if settings['verbose']:
                stats = self.get_subgraph_statistics(subgraph_id=subgraph_id)
                self.logger.info(f"Displayed Subgraph ID {stats['subgraph_id']}:")
                self.logger.info(f"  Number of Nodes: {stats['number_of_nodes']}")
                self.logger.info(f"  Number of Edges: {stats['number_of_edges']}")
                self.logger.info("  Node Types:")
                for nt, count in stats['node_types'].items():
                    if nt in ['mergeSE', 'gene']:
                        self.logger.info(f"    {nt}: {count}")
                self.logger.info("  Edge Types:")
                for et, count in stats['edge_types'].items():
                    self.logger.info(f"    {et}: {count}")
        except Exception as e:
            self.logger.exception(f"Error in display_two_layer_subgraph: {e}")
            raise
        self.logger.info("Two-layer subgraph displayed successfully.")  





    def update_default_settings(self, new_settings: Dict[str, Any]) -> None:

        self.logger.info(f"Updating default settings with: {new_settings}")  
        try:
            self.default_settings.update(new_settings)
            self.logger.info(f"Default settings updated successfully.")
        except Exception as e:
            self.logger.exception(f"Error in update_default_settings: {e}")
            raise


    def update_settings(self, 
                        graph_type: str, 
                        new_settings: Dict[str, Any]) -> None:

        self.logger.info(f"Updating settings for graph_type='{graph_type}' with: {new_settings}")  
        try:
            if graph_type == '3layer':
                self.settings_3layer.update(new_settings)
                self.logger.info(f"Settings for 3layer graph updated successfully.")
            elif graph_type == '2layer':
                self.settings_2layer.update(new_settings)
                self.logger.info(f"Settings for 2layer graph updated successfully.")
            else:
                self.logger.error(f"Unsupported graph type for settings update: {graph_type}")
                raise ValueError(f"Unsupported graph type: {graph_type}")
        except Exception as e:
            self.logger.exception(f"Error in update_settings: {e}")
            raise





    def get_current_settings(self, graph_type: Optional[str] = None) -> Dict[str, Any]:

        if graph_type == '3layer':
            return {
                'default_settings': self.default_settings,
                '3layer_settings': self.settings_3layer
            }
        elif graph_type == '2layer':
            return {
                'default_settings': self.default_settings,
                '2layer_settings': self.settings_2layer
            }
        elif graph_type is None:
            return {
                'default_settings': self.default_settings,
                '3layer_settings': self.settings_3layer,
                '2layer_settings': self.settings_2layer
            }
        else:
            self.logger.error(f"Unsupported graph type for settings retrieval: {graph_type}")
        raise ValueError(f"Unsupported graph type: {graph_type}")



    def display_current_settings(self, graph_type: Optional[str] = None) -> None:

        try:
            settings = self.get_current_settings(graph_type)
        except ValueError as ve:
            self.logger.error(f"Error retrieving settings: {ve}")
            return
        
        self.logger.info("----- Current Settings -----")
        for key, value in settings.items():
            self.logger.info(f"{key}: {value}")
        self.logger.info("----------------------------")




















    def _generate_se_file_list(self) -> list:
        self.logger.info("Generating SE file list from directory.")  
        try:
            if not self.rose_files_dir:
                raise ValueError("ROSE files directory is not specified.")
            
            se_path = pathlib.Path(self.rose_files_dir).resolve()
            if not se_path.exists() or not se_path.is_dir():
                raise FileNotFoundError(f"Specified ROSE directory does not exist: {se_path}")

            se_file_list = sorted(list(se_path.iterdir()))
            if not se_file_list:
                raise FileNotFoundError(f"No SE files found in the directory: {se_path}")
            
            self.logger.info(f"Generated SE file list with {len(se_file_list)} files.")  
            return se_file_list
        except Exception as e:
            self.logger.exception(f"Error in _generate_se_file_list: {e}")
            raise

    




    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    


    def datafile_exists(self) -> None:
        """Check if specified data files exist."""
        self.logger.info("Checking existence of data files.")  
        try:
            self.logger.info("Data Files Information:")
            files_to_check = {
                "P2GL File": self.p2g_file,
                "RNA Info File": self.rna_info_file,
                "ROSE Files Directory": self.rose_files_dir,
                "ROSE File": self.rose_file,
            }

            for label, file in files_to_check.items():
                if file:
                    if os.path.exists(file):
                        
                        log_path = self._get_log_path(file)
                        self.logger.info(f"{label}: {log_path}")
                    else:
                        self.logger.warning(f"{label} does not exist: {file}")
                else:
                    self.logger.warning(f"{label} is not provided.")
        except Exception as e:
            self.logger.exception(f"Error in datafile_exists: {e}")
            raise
        self.logger.info("Data file existence check completed.")  





    def count_genes(
        self, 
        output_file: Optional[str] = None, 
        show_item_number: bool = True, 
        show_genes: bool = True, 
        show_gene_count: bool = False, 
        sort_genes: bool = False, 
        save_and_print: bool = False
    ) -> None:

        self.logger.info("Starting gene count process.")  
        try:
            if not self._bed_filter:
                self.logger.error("Bed filter is not set. Please run filter_p2g_file first.")
                raise ValueError("Bed filter is not set. Please run filter_p2g_file first.")
            
            count_genes_from_bed(
                bed=self._bed_filter, 
                output_file=output_file, 
                show_item_number=show_item_number, 
                show_genes=show_genes, 
                show_gene_count=show_gene_count, 
                sort_genes=sort_genes, 
                save_and_print=save_and_print
            )
        except Exception as e:
            self.logger.exception(f"Error in count_genes: {e}")
            raise
        self.logger.info("Gene count process completed successfully.")  




    def create_ROSE_summary(
        self,
        show_plot: bool = True,
        save_svg: Optional[str] = None,
        save_text: Optional[str] = None,
        save_and_show: bool = False,
        include_pie_chart: bool = True,
        title: str = "ROSE SE Summary",
        color_all: str = "blue",
        color_super: str = "red",
        line_color: str = "green",
        line_style: str = "-"
    ) -> None:

        self.logger.info("Starting ROSE summary creation.")  
        try:
            if not self._bed_filter:
                raise ValueError("Bed filter is not set. Please run filter_p2g_file first.")
            
            if not self.rose_file or not os.path.exists(self.rose_file):
                raise FileNotFoundError(f"Enhancer table not found: {self.rose_file}")

            
            create_ROSE_summary(
                bed_obj_bed_filter=self._bed_filter,
                path_enhancer_table=self.rose_file,
                show_plot=show_plot,
                save_svg=save_svg,
                save_text=save_text,
                save_and_show=save_and_show,
                include_pie_chart=include_pie_chart,
                title=title,
                color_all=color_all,
                color_super=color_super,
                line_color=line_color,
                line_style=line_style
            )
        except Exception as e:
            self.logger.exception(f"Error in create_ROSE_summary: {e}")
            raise
        self.logger.info("ROSE summary created successfully.")  







    def filter_p2g_file(self, FDR: float = 0.05, r: float = 0.5) -> None:

        self.logger.info(f"Filtering P2G file with FDR={FDR} and r={r}.")  
        try:
            if not self.p2g_file or not os.path.exists(self.p2g_file):
                raise FileNotFoundError(f"P2G file not found: {self.p2g_file}")

            
            self._FDR = FDR
            self._r = r

            
            self._bed_filter = bed_filter_from_p2g_output(
                path_p2g_output=self.p2g_file, fdr_filter=FDR, r_filter=r
            )

            self.logger.info(f"Filtered Bedfile created with {len(self._bed_filter.to_dataframe())} entries.")
            self.logger.info(f"Current Filter Parameters: FDR={self._FDR}, r={self._r}")
        except Exception as e:
            self.logger.exception(f"Error in filter_p2g_file: {e}")
            raise
        self.logger.info("P2G file filtered successfully.")  









    def create_network(
        self,
        threshold: int = 10,
        network_FDR: Optional[float] = None,
        network_r: Optional[float] = None,
        verbose: bool = False,  
        verbose_gene: bool = False,  
        save_subgraph_info: bool = False,  
        subgraph_csv_path: Optional[str] = None  
    ) -> None:

        self.logger.info("Creating 3-layer DiGraph network...")  
        
        try:
            
            self._threshold_network = threshold
            self.logger.info(f"Using threshold: {self._threshold_network}")  

            
            self._network_FDR = network_FDR if network_FDR is not None else self._FDR
            self._network_r = network_r if network_r is not None else self._r
            self.logger.info(f"Using network filter parameters: FDR={self._network_FDR}, r={self._network_r}")

            
            df_p2gl_filtered_network = p2gl_path_to_filter_df(
                self.p2g_file, 
                fdr_filter=self._network_FDR, 
                r_filter=self._network_r
            )
            self.logger.info(f"Network-specific filtered DataFrame created with {len(df_p2gl_filtered_network)} entries.")
            
            
            if self.rna_info_file:
                self.logger.info(f"Processing RNA information from {self.rna_info_file}...")
                df_rna = csv_path_to_df_RNA(self.rna_info_file)
                dic_rna = rna_data_to_dic_metadata(df_rna, gtf=True)
                self.logger.info(f"RNA metadata dictionary created with {len(dic_rna)} entries.")
            else:
                self.logger.warning("No RNA info file provided. Proceeding without RNA metadata.")
                dic_rna = {}  
            
            
            list_se_files = [str(path) for path in self._input_se_file_list]
            
            
            DG, bed_after_threshold, bed_before_threshold, bed_df_concat_sort = create_3_layer_DiGraph_mergeSE_to_gene_custom(
                list_path=list_se_files,
                df_p2gl_filtered=df_p2gl_filtered_network,
                dic_rna=dic_rna,
                threshold=threshold
            )
            
            
            self._DG_network = DG
            self._bed_after_threshold_network = bed_after_threshold
            self._bed_before_threshold_network = bed_before_threshold
            self._bed_df_concat_sort_network = bed_df_concat_sort
            
            
            edge_data = list(DG.edges(data=True))
            if edge_data:
                
                df_edges_network = pd.DataFrame(edge_data, columns=['source', 'target', 'attributes'])
            else:
                df_edges_network = pd.DataFrame(columns=['source', 'target', 'attributes'])
            
            self._df_bed_mergeSE_to_gene_network = df_edges_network
            self.logger.info(f"df_bed_mergeSE_to_gene_network has been set with {len(df_edges_network)} rows.")
            
            
            self.logger.info("Network created successfully!")
            self.logger.info(f"Number of nodes: {DG.number_of_nodes()}")
            self.logger.info(f"Number of edges: {DG.number_of_edges()}")
            
            
            self.assign_subgraph_ids(verbose=verbose, verbose_gene=verbose_gene)
            
            
            if save_subgraph_info:
                self.save_subgraph_info(output_path=subgraph_csv_path)
        except Exception as e:
            self.logger.exception(f"Error in create_network: {e}")
            raise
        self.logger.info("Network creation process completed successfully.")  




    




    def assign_subgraph_ids(self, verbose: bool = True, verbose_gene: bool = True) -> None:

        self.logger.info("Assigning subgraph IDs to weakly connected components.")  
        try:
            if self._DG_network is None:
                raise ValueError("Network graph is not set. Please run create_network first.")
            
            
            weakly_connected_components = list(nx.weakly_connected_components(self._DG_network))
            if verbose:
                self.logger.info(f"Total weakly connected components to process: {len(weakly_connected_components)}")
            
            
            subgraphs_info = []

            for component in weakly_connected_components:
                subgraph = self._DG_network.subgraph(component)
                
                mergeSE_in_subgraph = [node for node, data in subgraph.nodes(data=True) if data.get('node_type') == 'mergeSE']
                num_mergeSE_nodes = len(mergeSE_in_subgraph)
                if mergeSE_in_subgraph:
                    
                    sample_counts = [subgraph.nodes[node].get('sample_count', 0) for node in mergeSE_in_subgraph]
                    SE_counts = [subgraph.nodes[node].get('SE_count', 0) for node in mergeSE_in_subgraph]
                    max_sample_count = max(sample_counts)
                    min_sample_count = min(sample_counts)
                    max_SE_count = max(SE_counts)
                    min_SE_count = min(SE_counts)
                    
                    max_sample_count_of_mergeSE_name = mergeSE_in_subgraph[sample_counts.index(max_sample_count)]
                    min_sample_count_of_mergeSE_name = mergeSE_in_subgraph[sample_counts.index(min_sample_count)]
                    max_SE_count_of_mergeSE_name = mergeSE_in_subgraph[SE_counts.index(max_SE_count)]
                    min_SE_count_of_mergeSE_name = mergeSE_in_subgraph[SE_counts.index(min_SE_count)]
                    
                    total_SE_count = sum(SE_counts)
                else:
                    max_sample_count = 0
                    min_sample_count = 0
                    max_SE_count = 0
                    min_SE_count = 0
                    max_sample_count_of_mergeSE_name = 'None'
                    min_sample_count_of_mergeSE_name = 'None'
                    max_SE_count_of_mergeSE_name = 'None'
                    min_SE_count_of_mergeSE_name = 'None'
                    total_SE_count = 0
                
                
                num_nodes = subgraph.number_of_nodes()
                
                
                first_node_name = sorted(component)[0] if component else ""
                
                
                mergeSE_names = mergeSE_in_subgraph
                mergeSE_data = []
                for node in mergeSE_in_subgraph:
                    sample_count = subgraph.nodes[node].get('sample_count', 0)
                    SE_count = subgraph.nodes[node].get('SE_count', 0)
                    mergeSE_data.append({
                        'name': node,
                        'sample_count': sample_count,
                        'SE_count': SE_count
                    })
                
                
                subgraphs_info.append({
                    'subgraph': subgraph,
                    'max_sample_count': max_sample_count,
                    'min_sample_count': min_sample_count,
                    'max_SE_count': max_SE_count,
                    'min_SE_count': min_SE_count,
                    'total_SE_count': total_SE_count,  
                    'max_sample_count_of_mergeSE_name': max_sample_count_of_mergeSE_name,
                    'min_sample_count_of_mergeSE_name': min_sample_count_of_mergeSE_name,
                    'max_SE_count_of_mergeSE_name': max_SE_count_of_mergeSE_name,
                    'min_SE_count_of_mergeSE_name': min_SE_count_of_mergeSE_name,
                    'num_mergeSE_nodes': num_mergeSE_nodes,
                    'num_nodes': num_nodes,
                    'first_node_name': first_node_name,
                    'mergeSE_names': mergeSE_names,  
                    'mergeSE_data': mergeSE_data     
                })
            
            
            
            
            
            sorted_subgraphs = sorted(
                subgraphs_info,
                key=lambda x: (-x['max_sample_count'], -x['num_nodes'], x['first_node_name'])
            )
            
            if verbose:
                self.logger.info("Subgraphs sorted based on sample count, number of nodes, and node names.")
            
            
            self.sorted_subgraphs_info = sorted_subgraphs

            
            for idx, subgraph_info in enumerate(sorted_subgraphs, start=1):
                subgraph = subgraph_info['subgraph']
                subgraph_id = idx
                
                for node in subgraph.nodes():
                    self._DG_network.nodes[node]['subgraph_id'] = subgraph_id
                
                for u, v, data in subgraph.edges(data=True):
                    self._DG_network.edges[u, v].update({'subgraph_id': subgraph_id})
                
                if verbose:
                    
                    num_mergeSE_nodes = subgraph_info['num_mergeSE_nodes']
                    enhancer_nodes = len([n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'enhancer'])
                    gene_nodes = len([n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'gene'])
                    
                    self.logger.info(f"Assigned subgraph ID{subgraph_id}: {subgraph.number_of_nodes()} nodes "
                        f"({num_mergeSE_nodes} mergeSE, {enhancer_nodes} enhancer, {gene_nodes} gene)")
                    
                    
                    self.logger.info(f"    - Max sample_count used for sorting: {subgraph_info['max_sample_count']}")
                    
                    
                    self.logger.info(f"    - Total SE_count in subgraph: {subgraph_info['total_SE_count']}")
                    
                    
                    self.logger.info(f"    - Max sample_count of mergeSE: {subgraph_info['max_sample_count']} (Node: {subgraph_info['max_sample_count_of_mergeSE_name']})")
                    self.logger.info(f"    - Min sample_count of mergeSE: {subgraph_info['min_sample_count']} (Node: {subgraph_info['min_sample_count_of_mergeSE_name']})")
                    self.logger.info(f"    - Max SE_count of mergeSE: {subgraph_info['max_SE_count']} (Node: {subgraph_info['max_SE_count_of_mergeSE_name']})")
                    self.logger.info(f"    - Min SE_count of mergeSE: {subgraph_info['min_SE_count']} (Node: {subgraph_info['min_SE_count_of_mergeSE_name']})")
                    
                    
                    is_have_mergeSE = num_mergeSE_nodes > 0
                    self.logger.info(f"    - Is have mergeSE: {is_have_mergeSE}")

                    
                    if verbose_gene:
                        genes = [n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'gene']
                        self.logger.info(f"    - Genes in subgraph ID{subgraph_id}: {', '.join(genes) if genes else 'None'}")
                        
                        self.logger.debug(f"        Debug: Number of genes found: {len(genes)}")
                        self.logger.debug(f"        Genes list: {genes}")

            if verbose:
                self.logger.info(f"Assigned subgraph IDs from 1 to {len(sorted_subgraphs)}.")
            
            
            self._subgraph_info_df = self._generate_subgraph_info_df()
        except Exception as e:
            self.logger.exception(f"Error in assign_subgraph_ids: {e}")
            raise
        self.logger.info("Subgraph ID assignment completed successfully.")  






    def _generate_subgraph_info_df(self) -> pd.DataFrame:

        self.logger.info("Generating subgraph information DataFrame.")  
        try:
            subgraph_data = []
            for subgraph_info in self.sorted_subgraphs_info:
                subgraph = subgraph_info['subgraph']
                
                subgraph_id = subgraph.nodes[next(iter(subgraph.nodes()))].get('subgraph_id', None)
                if subgraph_id is None:
                    continue  
                
                total_nodes = subgraph_info['num_nodes']
                mergeSE_nodes = subgraph_info['num_mergeSE_nodes']  
                enhancer_nodes = len([n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'enhancer'])
                gene_nodes = len([n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'gene'])
                genes = [n for n, d in subgraph.nodes(data=True) if d.get('node_type') == 'gene']
                genes_str = ', '.join(genes)
                max_sample_count = subgraph_info['max_sample_count']
                min_sample_count = subgraph_info['min_sample_count']
                max_SE_count = subgraph_info['max_SE_count']
                min_SE_count = subgraph_info['min_SE_count']
                max_sample_count_of_mergeSE_name = subgraph_info['max_sample_count_of_mergeSE_name']
                min_sample_count_of_mergeSE_name = subgraph_info['min_sample_count_of_mergeSE_name']
                max_SE_count_of_mergeSE_name = subgraph_info['max_SE_count_of_mergeSE_name']
                min_SE_count_of_mergeSE_name = subgraph_info['min_SE_count_of_mergeSE_name']
                is_have_mergeSE = subgraph_info['num_mergeSE_nodes'] > 0
                
                
                mergeSE_names = subgraph_info['mergeSE_names']
                mergeSE_names_str = ', '.join(mergeSE_names) if mergeSE_names else 'None'
                
                
                
                mergeSE_data_entries = [
                    f"{entry['name']}: {entry['sample_count']}, {entry['SE_count']}" 
                    for entry in subgraph_info['mergeSE_data']
                ]
                mergeSE_data_str = '; '.join(mergeSE_data_entries) if mergeSE_data_entries else 'None'
                
                subgraph_data.append({
                    'subgraph_id': subgraph_id,
                    'total_nodes': total_nodes,
                    'mergeSE_nodes': mergeSE_nodes,
                    'enhancer_nodes': enhancer_nodes,
                    'gene_nodes': gene_nodes,
                    'genes': genes_str,
                    'max_sample_count': max_sample_count,
                    'min_sample_count': min_sample_count,
                    'max_SE_count': max_SE_count,
                    'min_SE_count': min_SE_count,
                    'max_SAMPLE_count_of_mergeSE_name': max_sample_count_of_mergeSE_name,
                    'max_SAMPLE_count_of_mergeSE': max_sample_count,
                    'min_SAMPLE_count_of_mergeSE_name': min_sample_count_of_mergeSE_name,
                    'min_SAMPLE_count_of_mergeSE': min_sample_count,
                    'max_SE_count_of_mergeSE_name': max_SE_count_of_mergeSE_name,
                    'max_SE_count_of_mergeSE': max_SE_count,
                    'min_SE_count_of_mergeSE_name': min_SE_count_of_mergeSE_name,
                    'min_SE_count_of_mergeSE': min_SE_count,
                    'is_have_mergeSE': is_have_mergeSE,
                    'mergeSE_name': mergeSE_names_str,
                    'mergeSE_data': mergeSE_data_str
                })
            
            
            df = pd.DataFrame(subgraph_data)
            self.logger.info("Subgraph information DataFrame generated successfully.")  
            return df
        except Exception as e:
            self.logger.exception(f"Error in _generate_subgraph_info_df: {e}")
            raise



        

                

    def save_subgraph_info(self, output_path: Optional[str] = None) -> None:
        self.logger.info("Saving subgraph information to CSV.")  
        try:
            if not hasattr(self, '_subgraph_info_df') or self._subgraph_info_df is None:
                raise ValueError("Subgraph information is not available. Please run assign_subgraph_ids first.")

            if output_path:
                
                output_dir = os.path.dirname(output_path)
                if output_dir and not os.path.exists(output_dir):
                    raise FileNotFoundError(f"Specified output directory does not exist: {output_dir}")
            else:
                
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)

                
                p2g_filename = os.path.splitext(os.path.basename(self.p2g_file))[0]

                
                rose_dir_name = os.path.basename(os.path.normpath(self.rose_files_dir))

                
                sort_by = self.current_sort_key or "se_count"  
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  
                filename = f"subgraph_{p2g_filename}_{rose_dir_name}_{self._network_r}_{self._threshold_network}_{timestamp}.csv"  
                output_path = os.path.join(output_dir, filename)

            
            self._subgraph_info_df.to_csv(output_path, index=False)
            self.logger.info(f"Subgraph information saved successfully to {output_path} with {len(self._subgraph_info_df)} records.")
        except Exception as e:
            self.logger.exception(f"Error in save_subgraph_info: {e}")
            raise
        self.logger.info("Subgraph information saved successfully.")  




    def draw_network(
        self,
        output_path: Optional[str] = None,
        prog: str = 'fdp',
        format: str = 'svg',
        italic_genes: bool = False  
    ) -> None:

        self.logger.info("Drawing network graph.")  

        try:
            if self._DG_network is None:
                raise ValueError("Network graph is not set. Please run create_network first.")
            
            
            agraph_DG = to_agraph(self._DG_network)
            
            
            if italic_genes:
                for node in agraph_DG.nodes():
                    node_type = self._DG_network.nodes[node].get('node_type', 'default')
                    if node_type == 'gene':
                        agraph_DG.get_node(node).attr['fontstyle'] = 'italic'
            
            
            if output_path is None:
                
                try:
                    default_filename = f"network_FDR{self.network_FDR}_r{self.network_r}_threshold{self.network_threshold}.{format}"
                except AttributeError as e:
                    self.logger.error(f"Failed to generate default filename: {e}")
                    raise ValueError("Network parameters are not properly set.")
                
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, default_filename)
            else:
                
                output_dir = os.path.dirname(output_path)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
            
            
            self.logger.info(f"Saving network graph to {output_path} with prog={prog} and format={format}.")
            
            
            
            
            
            agraph_DG.draw(output_path, prog=prog, format=format)

            self.logger.info(f"Network graph saved successfully to {output_path} in {format} format.")
        except Exception as e:
            self.logger.exception(f"Error in draw_network: {e}")
            raise
        self.logger.info("Network graph drawing process completed successfully.")  



    def save_network(
        self,
        output_path: Optional[str] = None,
        format: str = "pickle"
    ) -> None:
        """
        Save the networkx graph object to a file in the specified format.
        """
        self.logger.info(f"Saving network graph in format '{format}' to '{output_path}'." if output_path else "Saving network graph with default path and format.")  
        try:
            if self._DG_network is None:
                raise ValueError("Network graph is not set. Please run create_network first.")

            supported_formats = ['pickle', 'gml', 'json', 'edgelist']
            if format not in supported_formats:
                self.logger.error(f"Unsupported format '{format}'. Supported formats are: {supported_formats}")
                raise ValueError(f"Unsupported format '{format}'. Supported formats are: {supported_formats}")

            if output_path:
                
                output_dir = os.path.dirname(output_path)
                if output_dir and not os.path.exists(output_dir):
                    os.makedirs(output_dir, exist_ok=True)
            else:
                
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)

                
                p2g_filename = os.path.splitext(os.path.basename(self.p2g_file))[0]

                
                rose_dir_name = os.path.basename(os.path.normpath(self.rose_files_dir))

                
                extension_map = {
                    'pickle': 'pickle',
                    'gml': 'gml',
                    'json': 'json',
                    'edgelist': 'edgelist'
                }
                extension = extension_map.get(format, 'pickle')  

                
                filename = f"network_{p2g_filename}_{rose_dir_name}_{self.network_r}_threshold{self.network_threshold}.{extension}"
                output_path = os.path.join(output_dir, filename)

            
            if format == 'pickle':
                with open(output_path, 'wb') as f:
                    pickle.dump(self._DG_network, f, protocol=pickle.HIGHEST_PROTOCOL)
            elif format == 'gml':
                nx.write_gml(self._DG_network, output_path)
            elif format == 'json':
                data = json_graph.node_link_data(self._DG_network)
                with open(output_path, 'w') as f:
                    json.dump(data, f, indent=4)  
            elif format == 'edgelist':
                nx.write_edgelist(self._DG_network, output_path, data=True)
            
            self.logger.info(f"Network graph saved successfully to {output_path} in {format} format.")
        except Exception as e:
            self.logger.exception(f"Error in save_network: {e}")
            raise
        self.logger.info("Network saving process completed successfully.")  




    def search_subnetwork(
        self, 
        name: str, 
        node_type: str,
        set_subgraph: bool = True,  
        verbose: bool = True        
    ) -> Optional[nx.DiGraph]:

        self.logger.info(f"Searching for subnetwork with name='{name}' and node_type='{node_type}'.")  
        try:
            
            valid_node_types = ['mergeSE', 'gene', 'enhancer']
            if node_type not in valid_node_types:
                raise ValueError(f"Invalid node type '{node_type}' specified. Supported types: {valid_node_types}")

            
            if self._DG_network is None:
                raise ValueError("Network graph is not set. Please run create_network first.")

            
            node_type_normalized = node_type.strip().lower()
            name_normalized = name.strip().lower()

            
            self.logger.debug(f"Normalized search criteria - Name: {name_normalized}, Node Type: {node_type_normalized}")

            
            matching_nodes = [
                n for n, attr in self._DG_network.nodes(data=True)
                if attr.get('node_type', '').strip().lower() == node_type_normalized and str(n).strip().lower() == name_normalized
            ]

            
            self.logger.debug(f"Matching nodes: {matching_nodes}")

            if not matching_nodes:
                self.logger.info(f"No nodes found matching name '{name}' and node type '{node_type}'.")
                return None  

            
            subgraph_ids = set()
            for node in matching_nodes:
                subgraph_id = self._DG_network.nodes[node].get('subgraph_id')
                if subgraph_id is not None:
                    subgraph_ids.add(subgraph_id)
                    self.logger.debug(f"Node '{node}' belongs to Subgraph ID {subgraph_id}.")

            if len(subgraph_ids) > 1:
                self.logger.error(f"Multiple subgraphs matched the specified criteria. Subgraph IDs: {subgraph_ids}")
                raise ValueError(f"Multiple subgraphs matched the specified criteria. Subgraph IDs: {subgraph_ids}")

            subgraph_id = subgraph_ids.pop()
            self.logger.info(f"Target Subgraph ID: {subgraph_id}")

            
            nodes_in_subgraph = [
                n for n, attr in self._DG_network.nodes(data=True) if attr.get('subgraph_id') == subgraph_id
            ]

            subgraph = self._DG_network.subgraph(nodes_in_subgraph).copy()

            
            self._print_subgraph_summary(subgraph, verbose=verbose)  

            
            if set_subgraph:
                self._selected_subgraph_id = subgraph_id
                self._selected_subgraph = subgraph
                self.logger.info(f"Subgraph ID {subgraph_id} and the subgraph have been saved internally.")

            self.logger.info(f"Subgraph matching name '{name}' and node type '{node_type}' has been found.")
        except Exception as e:
            self.logger.exception(f"Error in search_subnetwork: {e}")
            raise
        return subgraph





    def _print_subgraph_summary(self, subgraph: nx.DiGraph, verbose: bool = False) -> None:

        num_nodes = subgraph.number_of_nodes()
        num_edges = subgraph.number_of_edges()
        node_types = {}
        mergeSE_sample_counts = []
        mergeSE_SE_counts = []
        mergeSE_details = []
        genes = []

        
        for node, attr in subgraph.nodes(data=True):
            nt = attr.get('node_type', 'unknown')
            node_types[nt] = node_types.get(nt, 0) + 1
            
            if nt == 'mergeSE':
                sample_count = attr.get('sample_count', 0)
                SE_count = attr.get('SE_count', 0)
                mergeSE_sample_counts.append(sample_count)
                mergeSE_SE_counts.append(SE_count)
                mergeSE_details.append({
                    'name': node,
                    'sample_count': sample_count,
                    'SE_count': SE_count
                })
            elif nt == 'gene':
                genes.append(node)

        
        self.logger.info("----- Subgraph Summary -----")
        self.logger.info(f"Number of Nodes: {num_nodes}")
        self.logger.info(f"Number of Edges: {num_edges}")
        self.logger.info("Node Types Breakdown:")
        for nt, count in node_types.items():
            self.logger.info(f"  {nt}: {count}")

        
        if mergeSE_sample_counts and mergeSE_SE_counts:
            max_sample = max(mergeSE_sample_counts)
            min_sample = min(mergeSE_sample_counts)
            max_SE = max(mergeSE_SE_counts)
            min_SE = min(mergeSE_SE_counts)
            
            self.logger.info(f"MergeSE Sample Counts: Max = {max_sample}, Min = {min_sample}")
            self.logger.info(f"MergeSE SE Counts: Max = {max_SE}, Min = {min_SE}")
        else:
            self.logger.info("No 'mergeSE' nodes present in the subgraph.")

        
        if verbose:
            self.logger.info("----- Detailed MergeSE Information -----")
            if mergeSE_details:
                for ms in mergeSE_details:
                    self.logger.info(f"MergeSE Name: {ms['name']}, Sample Count: {ms['sample_count']}, SE Count: {ms['SE_count']}")
            else:
                self.logger.info("No 'mergeSE' nodes to display.")
            self.logger.info("----------------------------------------")

            self.logger.info("----- Genes in Subgraph -----")
            if genes:
                self.logger.info(", ".join(genes))
            else:
                self.logger.info("No genes present in the subgraph.")
            self.logger.info("-----------------------------")



    def get_selected_subgraph_info(self) -> Dict[str, any]:
        """
        Returns detailed information about the selected subgraph as a dictionary.
        
        Returns:
            Dict[str, any]: Detailed information about the subgraph.
        
        Raises:
            ValueError: If no subgraph has been selected.
        """
        if self._selected_subgraph is None:
            raise ValueError("No subgraph has been selected. Please run search_subnetwork first.")
        
        subgraph = self._selected_subgraph
        info = {
            "subgraph_id": self._selected_subgraph_id,
            "number_of_nodes": subgraph.number_of_nodes(),
            "number_of_edges": subgraph.number_of_edges(),
            "node_types": {}
        }
        
        for _, attr in subgraph.nodes(data=True):
            nt = attr.get('node_type', 'unknown')
            info["node_types"][nt] = info["node_types"].get(nt, 0) + 1
        
        return info




    def get_subgraph_statistics(self, subgraph_id: Optional[int] = None) -> Dict[str, any]:

        self.logger.info(f"Retrieving statistics for subgraph_id={subgraph_id}.")  
        try:
            if subgraph_id is None:
                subgraph_id = self.selected_subgraph_id
                if subgraph_id is None:
                    raise ValueError("No subgraph is currently selected. Please set a subgraph first.")
            
            
            if not hasattr(self, '_subgraph_info_df') or self._subgraph_info_df is None:
                raise ValueError("Subgraph information is not available. Please run assign_subgraph_ids first.")
            
            if subgraph_id not in self._subgraph_info_df['subgraph_id'].values:
                raise ValueError(f"Subgraph ID {subgraph_id} does not exist.")
            
            
            nodes_in_subgraph = [
                n for n, d in self._DG_network.nodes(data=True) if d.get('subgraph_id') == subgraph_id
            ]
            
            if not nodes_in_subgraph:
                raise ValueError(f"No nodes found for Subgraph ID {subgraph_id}.")
            
            
            subgraph = self._DG_network.subgraph(nodes_in_subgraph).copy()
            
            
            stats = {
                "subgraph_id": subgraph_id,
                "number_of_nodes": subgraph.number_of_nodes(),
                "number_of_edges": subgraph.number_of_edges(),
                "node_types": {},
                "edge_types": {},
                "degree_centrality": {},
                "betweenness_centrality": {}
            }
            
            
            for _, attr in subgraph.nodes(data=True):
                nt = attr.get('node_type', 'unknown')
                stats["node_types"][nt] = stats["node_types"].get(nt, 0) + 1
            
            
            for _, _, attr in subgraph.edges(data=True):
                et = attr.get('edge_type', 'unknown')
                stats["edge_types"][et] = stats["edge_types"].get(et, 0) + 1
            
            
            stats["degree_centrality"] = nx.degree_centrality(subgraph)
            stats["betweenness_centrality"] = nx.betweenness_centrality(subgraph)
            
            self.logger.info(f"Statistics for subgraph_id={subgraph_id} retrieved successfully.")  
            return stats
        except Exception as e:
            self.logger.exception(f"Error in get_subgraph_statistics: {e}")
            raise





    def print_subgraph_statistics(self, subgraph_id: Optional[int] = None, verbose: bool = False) -> None:

        self.logger.info(f"Printing statistics for subgraph_id={subgraph_id}.")  
        try:
            stats = self.get_subgraph_statistics(subgraph_id=subgraph_id)
        except ValueError as ve:
            self.logger.error(f"Error: {ve}")
            return
        
        self.logger.info(f"----- Subgraph Statistics (ID: {stats['subgraph_id']}) -----")
        self.logger.info(f"Number of Nodes: {stats['number_of_nodes']}")
        self.logger.info(f"Number of Edges: {stats['number_of_edges']}")
        
        self.logger.info("Node Types Breakdown:")
        for nt, count in stats["node_types"].items():
            self.logger.info(f"  {nt}: {count}")
        
        self.logger.info("Edge Types Breakdown:")
        for et, count in stats["edge_types"].items():
            self.logger.info(f"  {et}: {count}")
        
        if verbose:
            self.logger.info("\nDegree Centrality:")
            for node, centrality in stats["degree_centrality"].items():
                self.logger.info(f"  {node}: {centrality:.4f}")
            
            self.logger.info("\nBetweenness Centrality:")
            for node, centrality in stats["betweenness_centrality"].items():
                self.logger.info(f"  {node}: {centrality:.4f}")
        
        self.logger.info("--------------------------------------------------")
        self.logger.info("Subgraph statistics printed successfully.")  




    def set_subgraph_id(self, subgraph_id: int, verbose: bool = True) -> None:

        self.logger.info(f"Setting subgraph ID to {subgraph_id}.")  
        try:
            if not hasattr(self, '_subgraph_info_df') or self._subgraph_info_df is None:
                raise ValueError("Subgraph information is not available. Please run assign_subgraph_ids first.")

            
            if subgraph_id not in self._subgraph_info_df['subgraph_id'].values:
                self.logger.error(f"Subgraph ID {subgraph_id} does not exist.")
                raise ValueError(f"Subgraph ID {subgraph_id} does not exist.")

            
            nodes_in_subgraph = [
                n for n, d in self._DG_network.nodes(data=True) if d.get('subgraph_id') == subgraph_id
            ]

            
            if not nodes_in_subgraph:
                self.logger.error(f"No nodes found for Subgraph ID {subgraph_id}.")
                raise ValueError(f"No nodes found for Subgraph ID {subgraph_id}.")

            
            subgraph = self._DG_network.subgraph(nodes_in_subgraph).copy()

            
            self._selected_subgraph_id = subgraph_id
            self._selected_subgraph = subgraph

            if verbose:
                self.logger.info(f"Subgraph ID {subgraph_id} has been set as the current selected subgraph.")
                self._print_subgraph_summary(subgraph, verbose=True)
        except Exception as e:
            self.logger.exception(f"Error in set_subgraph_id: {e}")
            raise
        self.logger.info(f"Subgraph ID {subgraph_id} has been set successfully.")  






    def analyze_merge_SE(
        self, 
        save_concat_tsv: Optional[str] = None, 
        save_merge_tsv: Optional[str] = None, 
        save_sorted_tsv: Optional[str] = None, 
        save_svg: Optional[str] = None, 
        save_and_show: bool = True, 
        graph_max_number: int = 30, 
        sort_by: str = "se_count"
    ) -> None:

        self.logger.info("Starting analyze_merge_SE process.")  

        
        self.logger.info("Step 1: Concatenating SE data.")  
        try:
            self._temp_concat_full, _ = return_all_se_concat_full_df(
                self._input_se_file_list, self.bed_filter
            )
            self.logger.info(f"Concatenated SE data created with {len(self._temp_concat_full)} entries.")
            if save_concat_tsv:
                self._temp_concat_full.to_csv(save_concat_tsv, sep="\t", index=False)
                self.logger.info(f"Concatenated SE data saved to {save_concat_tsv}.")
        except Exception as e:
            self.logger.exception(f"Error in analyze_merge_SE during concatenating SE data: {e}")
            raise

        
        self.logger.info("Step 2: Merging SE counts.")  
        try:
            self._temp_full_df_edit = return_merge_se_count_df_full_edit(
                self._temp_concat_full, self.bed_filter
            )
            self.logger.info(f"Merged SE data created with {len(self._temp_full_df_edit)} entries.")
            if save_merge_tsv:
                self._temp_full_df_edit.to_csv(save_merge_tsv, sep="\t", index=False)
                self.logger.info(f"Merged SE data saved to {save_merge_tsv}.")
        except Exception as e:
            self.logger.exception(f"Error in analyze_merge_SE during merging SE counts: {e}")
            raise

        
        self.logger.info(f"Step 3: Sorting SE data by column '{sort_by}'.")  

        if sort_by not in self._temp_full_df_edit.columns:
            self.logger.error(f"Invalid sort_by column: '{sort_by}'. Available columns: {list(self._temp_full_df_edit.columns)}")
            raise ValueError(f"Invalid sort_by column: '{sort_by}'.")

        try:
            self._sorted_df = sort_se_merge_count_by_column(
                self._temp_full_df_edit, column=sort_by
            )
            self.logger.info(f"Sorted SE data created by '{sort_by}' column.")
            
            
            self._sorted_df = self._sorted_df.reset_index(drop=True)
            self.logger.debug("DataFrame index has been reset to start from 0.")
            
            if save_sorted_tsv:
                self._sorted_df.to_csv(save_sorted_tsv, sep="\t", index=True)  
                self.logger.info(f"Sorted SE data saved to {save_sorted_tsv}.")
        
            
            self.current_sort_key = sort_by  
            self.logger.debug(f"Current sort key set to '{self.current_sort_key}'.")  
        except Exception as e:
            self.logger.exception(f"Error in analyze_merge_SE during sorting SE data: {e}")
            raise

        
        self.logger.info("Step 4: Generating and saving SE count graph.")  
        try:
            self._make_sort_se_merge_count_graph(
                se_count_df=self._sorted_df, 
                graph_max_number=graph_max_number, 
                count_type=sort_by,
                title=f"Merged SE Count (Sorted by {sort_by})",
                save_svg=save_svg,
                save_and_show=save_and_show
            )
            self.logger.info("SE count graph generated successfully.")
            if save_svg:
                self.logger.info(f"SE count graph saved to {save_svg}.")
        except Exception as e:
            self.logger.exception(f"Error in analyze_merge_SE during graph generation: {e}")
            raise

        self.logger.info("Completed analyze_merge_SE process.")  





    def select_se_region(self, index: int) -> None:

        self.logger.info(f"Selecting SE region at index {index}.")


        if self._sorted_df is None:
            self.logger.error("SE regions are not available. Run analyze_merge_SE first.")
            raise ValueError("SE regions are not available. Run analyze_merge_SE first.")
        if index < 1 or index > len(self._sorted_df):
            self.logger.error(f"Index {index} out of range. Valid range: 1 to {len(self._sorted_df)}.")
            raise IndexError(f"Index {index} out of range.")
        try:
            self._selected_se_index = index
            self._selected_se_region = self._sorted_df.iloc[index - 1]["se_data"]
            self.logger.info(f"Selected SE Region: {self._selected_se_region} (Index: {index}).")
        except Exception as e:
            self.logger.error(f"Error selecting SE region at index {index}: {e}")
            raise



    def visualize_selected_se_region(
        self, 
        title: Optional[str] = None, 
        save_svg: Optional[str] = None, 
        save_region_bed: Optional[str] = None, 
        save_full_bed: Optional[str] = None, 
        show_plot: bool = True
    ) -> None:

        self.logger.info("Visualizing selected SE region.")  

        if self._selected_se_region is None:
            self.logger.error("No SE region selected. Use select_se_region first.")
            raise ValueError("No SE region selected. Use select_se_region first.")
        
        
        self.logger.debug(f"Parsing selected SE region: {self._selected_se_region}.")  
        try:
            chrom, start, end = self._selected_se_region.split("_")
            start = int(start)
            end = int(end)
            self.logger.debug(f"Parsed SE region - Chromosome: {chrom}, Start: {start}, End: {end}.")
        except ValueError as e:
            self.logger.error(f"Selected SE region format is incorrect: {self._selected_se_region}. Error: {e}")
            raise ValueError(f"Selected SE region format is incorrect: {self._selected_se_region}")
        
        
        self.logger.debug("Converting SE data to BedTool object.")
        try:
            bed_data = BedTool.from_dataframe(self._temp_concat_full)
            self.logger.info("Converted SE data to BedTool object successfully.")
        except Exception as e:
            self.logger.exception(f"Error converting SE data to BedTool: {e}")
            raise
        
        
        self.logger.debug("Generating SE region visualization.")
        try:
            plot_stacked_reads_bed(
                bed=bed_data, 
                chrom=chrom, 
                start=start, 
                end=end, 
                title=title or f"SE Region: {self._selected_se_region}", 
                save_path=save_svg,
                save_region_bed=save_region_bed,
                save_full_bed=save_full_bed,
                show_plot=show_plot
            )
            self.logger.info("SE region visualization completed successfully.")
            if save_svg:
                self.logger.info(f"SE region plot saved to {save_svg}.")
            if save_region_bed:
                self.logger.info(f"Selected SE region BED saved to {save_region_bed}.")
            if save_full_bed:
                self.logger.info(f"Full SE BED data saved to {save_full_bed}.")
        except Exception as e:
            self.logger.exception(f"Error during SE region visualization: {e}")
            raise







    def display_temp_concat_full(self, top_n: int = 30) -> None:

        self.logger.info(f"Displaying top {top_n} rows of concatenated SE data.")
        try:
            df = self.temp_concat_full.head(top_n)
            display(df)
            self.logger.info("Displayed concatenated SE data successfully.")
        except Exception as e:
            self.logger.exception(f"Error in display_temp_concat_full: {e}")
            raise

    def display_temp_full_df_edit(self, top_n: int = 30) -> None:

        self.logger.info(f"Displaying top {top_n} rows of merged SE data.")
        try:
            df = self.temp_full_df_edit.head(top_n)
            display(df)
            self.logger.info("Displayed merged SE data successfully.")
        except Exception as e:
            self.logger.exception(f"Error in display_temp_full_df_edit: {e}")
            raise

    def _apply_styling(self, styler: 'pd.Styler', df: pd.DataFrame) -> 'pd.Styler':

        
        styler = styler.set_table_styles([
            {'selector': 'th', 'props': [('text-align', 'right')]},
            {'selector': 'td', 'props': [('text-align', 'right')]}
        ])
        
        
        if 'peak_list' in df.columns:
            styler = styler.set_properties(
                subset=['peak_list'],
                **{
                    'white-space': 'nowrap',
                    'overflow': 'hidden',
                    'text-overflow': 'ellipsis',
                    'max-width': '150px'
                }
            )
        
        return styler

    def display_sorted_df(self, top_n: int = 30) -> None:

        self.logger.info(f"Displaying top {top_n} rows of sorted SE data.")
        try:
            
            df = self.sorted_df.head(top_n).copy()
            
            
            df.index += 1
            df.index.name = "rank"
            
            
            for col in ['gene_list', 'peak_list']:
                if col in df.columns:
                    df[col] = df[col].apply(
                        lambda x: ', '.join(x) if isinstance(x, list) else str(x)
                    )
                    if col == 'peak_list':  
                        df[col] = df[col].str.strip('[]')
            
            
            styled_df = df.style.pipe(self._apply_styling, df)
            
            
            display(styled_df)
            
            self.logger.info("Displayed sorted SE data successfully.")
        except Exception as e:
            self.logger.exception(f"Error in display_sorted_df: {e}")
            raise




    def search_sorted_df(self, search_term: str, search_column: str = 'gene_list', case_sensitive: bool = False) -> None:

        self.logger.info(f"Searching for term '{search_term}' in column '{search_column}'.")
        try:
            if self._sorted_df is None:
                self.logger.error("Sorted data is not available. Please run analyze_merge_SE first.")
                raise ValueError("Sorted data is not available. Please run analyze_merge_SE first.")
            
            if search_column not in self._sorted_df.columns:
                self.logger.error(f"Column '{search_column}' not found in sorted_df.")
                raise ValueError(f"Column '{search_column}' not found in sorted_df.")
            
            
            if case_sensitive:
                mask = self._sorted_df[search_column].apply(
                    lambda x: search_term in (x if isinstance(x, list) else str(x).split(', '))
                )
            else:
                search_term_lower = search_term.lower()
                mask = self._sorted_df[search_column].apply(
                    lambda x: search_term_lower in (
                        [g.lower() for g in x] if isinstance(x, list) 
                        else [g.lower() for g in str(x).split(', ')]
                    )
                )
            
            
            matching_rows = self._sorted_df[mask].copy()
            
            if matching_rows.empty:
                self.logger.info(f"No matches found for term '{search_term}'.")
                print(f"No matches found for term '{search_term}'.")
                return
            
            
            original_indices = np.where(mask)[0] + 1
            matching_rows['rank'] = original_indices
            
            
            cols = matching_rows.columns.tolist()
            cols.remove('rank')
            cols = ['rank'] + cols
            matching_rows = matching_rows[cols]
            
            
            for col in ['gene_list', 'peak_list']:
                if col in matching_rows.columns:
                    matching_rows[col] = matching_rows[col].apply(
                        lambda x: ', '.join(x) if isinstance(x, list) else str(x)
                    )
                    if col == 'peak_list':  
                        matching_rows[col] = matching_rows[col].str.strip('[]')
            
            
            matching_rows = matching_rows.set_index('rank')
            
            
            styled_df = matching_rows.style.pipe(self._apply_styling, matching_rows)
            
            
            self.logger.info(f"Found {len(matching_rows)} matches for term '{search_term}'.")
            print(f"Found {len(matching_rows)} matches for term '{search_term}'.")
            display(styled_df)
            
        except Exception as e:
            self.logger.exception(f"Error in search_sorted_df: {e}")
            raise





    def save_temp_concat_full(self, output_path: Optional[str] = None) -> None:

        self.logger.info("Saving concatenated SE data to TSV.")
        try:
            if output_path is None:
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                p2g_filename = os.path.splitext(os.path.basename(self.p2g_file))[0]
                filename = f"concat_SE_data_{p2g_filename}_{timestamp}.tsv"
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, filename)
            
            self.temp_concat_full.to_csv(output_path, sep="\t", index=False)
            self.logger.info(f"Concatenated SE data saved to {output_path}.")
        except Exception as e:
            self.logger.exception(f"Error in save_temp_concat_full: {e}")
            raise

    def save_temp_full_df_edit(self, output_path: Optional[str] = None) -> None:
        """
        Save the merged SE count DataFrame to a TSV file.

        Args:
            output_path (Optional[str]): Path to save the TSV file. If None, uses a default path.
        """
        self.logger.info("Saving merged SE data to TSV.")
        try:
            if output_path is None:
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                p2g_filename = os.path.splitext(os.path.basename(self.p2g_file))[0]
                filename = f"merged_SE_data_{p2g_filename}_{timestamp}.tsv"
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, filename)
            
            self.temp_full_df_edit.to_csv(output_path, sep="\t", index=False)
            self.logger.info(f"Merged SE data saved to {output_path}.")
        except Exception as e:
            self.logger.exception(f"Error in save_temp_full_df_edit: {e}")
            raise

    def save_sorted_df(self, output_path: Optional[str] = None) -> None:
        """
        Save the sorted SE count DataFrame to a TSV file.

        Args:
            output_path (Optional[str]): Path to save the TSV file. If None, uses a default path.
        """
        self.logger.info("Saving sorted SE data to TSV.")
        try:
            if output_path is None:
                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                p2g_filename = os.path.splitext(os.path.basename(self.p2g_file))[0]
                
                sort_by = self.current_sort_key or "se_count"  
                filename = f"sorted_SE_data_{p2g_filename}_{sort_by}_{timestamp}.tsv"
                output_dir = './output'
                os.makedirs(output_dir, exist_ok=True)
                output_path = os.path.join(output_dir, filename)
            
            self.sorted_df.to_csv(output_path, sep="\t", index=False)
            self.logger.info(f"Sorted SE data saved to {output_path}.")
        except Exception as e:
            self.logger.exception(f"Error in save_sorted_df: {e}")
            raise



    def resort_sorted_df(self, new_sort_key: str, ascending: bool = False) -> None:

        self.logger.info(f"Resorting sorted_df by new_sort_key='{new_sort_key}' with ascending={ascending}.")
        try:
            if self._sorted_df is None:
                self.logger.error("Sorted SE data is not available. Please run analyze_merge_SE first.")
                raise ValueError("Sorted SE data is not available. Please run analyze_merge_SE first.")
            
            if new_sort_key not in self._sorted_df.columns:
                self.logger.error(f"Invalid new_sort_key: '{new_sort_key}'. Available columns: {list(self._sorted_df.columns)}")
                raise ValueError(f"Invalid new_sort_key: '{new_sort_key}'.")
            
            
            self._sorted_df = self._sorted_df.sort_values(by=new_sort_key, ascending=ascending).reset_index(drop=True)
            self.logger.info(f"Sorted sorted_df by '{new_sort_key}' successfully.")
            
            
            self.current_sort_key = new_sort_key
            self.logger.debug(f"Current sort key set to '{self.current_sort_key}'.")
            
        except Exception as e:
            self.logger.exception(f"Error in resort_sorted_df: {e}")
            raise

    
        


    def plot_sorted_se_count_graph(
        self, 
        graph_max_number: int = 30, 
        y_column: Optional[str] = None,  
        save_svg: Optional[str] = None, 
        save_and_show: bool = True,
        ylabel: Optional[str] = None,  
        title: Optional[str] = None    
    ) -> None:

        self.logger.info("Generating SE count graph based on the current sorted_df.")
        try:
            if self._sorted_df is None:
                self.logger.error("Sorted SE data is not available. Please run analyze_merge_SE first.")
                raise ValueError("Sorted SE data is not available. Please run analyze_merge_SE first.")
            
            
            sort_by = self.current_sort_key
            if sort_by is None:
                self.logger.warning("Current sort key is not set. Using the first available sort key.")
                sort_by = self._sorted_df.columns[0]
                self.current_sort_key = sort_by
                self.logger.info(f"Defaulting to sort key '{sort_by}'.")
            
            
            if y_column is None:
                y_column = sort_by
                self.logger.debug(f"y_column is not specified. Using current_sort_key '{sort_by}' as y_column.")
            
            
            if title is None:
                title = f"Merged SE Count (Sorted by {sort_by})"
            
            
            if ylabel is None:
                ylabel = y_column.replace("_", " ").capitalize()
            
            self.logger.info(f"Generating SE count graph sorted by '{sort_by}' with y_column='{y_column}'.")
            
            
            self._make_sort_se_merge_count_graph(
                se_count_df=self._sorted_df, 
                graph_max_number=graph_max_number, 
                count_type=y_column,  
                title=title,
                xlabel="SE Data",      
                ylabel=ylabel,
                save_svg=save_svg,
                save_and_show=save_and_show
            )
            self.logger.info("SE count graph generated successfully.")
            if save_svg:
                self.logger.info(f"SE count graph saved to {save_svg}.")
        except Exception as e:
            self.logger.exception(f"Error in plot_sorted_se_count_graph: {e}")
            raise










    def _make_sort_se_merge_count_graph(
        self,
        se_count_df: pd.DataFrame,
        graph_max_number: int,
        count_type: str = "se_count",
        title: Optional[str] = None,
        xlabel: Optional[str] = "merge SE",
        ylabel: Optional[str] = None,
        rotation: int = 90,
        save_svg: Optional[str] = None,
        save_tsv: Optional[str] = None,
        save_and_show: bool = False
    ) -> None:

        self.logger.info("Starting _make_sort_se_merge_count_graph.")
        
        
        if count_type not in se_count_df.columns:
            self.logger.error(f"Invalid count_type: '{count_type}'. Available columns: {list(se_count_df.columns)}")
            raise ValueError(f"Invalid count_type: '{count_type}'. Available columns: {list(se_count_df.columns)}")
        
        
        x_list = se_count_df["se_data"].tolist()
        y_list = se_count_df[count_type].tolist()

        
        x_list_graph = x_list[:graph_max_number]
        y_list_graph = y_list[:graph_max_number]

        
        display_data = pd.DataFrame({
            "se_data": x_list_graph,
            count_type: y_list_graph
        })

        

        plt.bar(x_list_graph, y_list_graph)  
        plt.tight_layout()

        
        plt.xticks(rotation=rotation)

        
        if pd.api.types.is_integer_dtype(se_count_df[count_type]):
            plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        
        plt.xlabel(xlabel)
        plt.ylabel(ylabel if ylabel else count_type.replace("_", " ").capitalize())
        plt.title(title if title else f"Merge SE {count_type.replace('_', ' ').capitalize()}")

        
        if save_svg:
            plt.savefig(save_svg, format="svg")
            self.logger.info(f"Graph saved as SVG at {save_svg}.")

        
        if save_and_show or not save_svg:
            plt.show()

        
        plt.clf()
        plt.close()

        
        if save_tsv:
            display_data.to_csv(save_tsv, sep="\t", index=False)
            self.logger.info(f"Data saved as TSV at {save_tsv}.")
        
        self.logger.info("Completed _make_sort_se_merge_count_graph.")






    def _get_log_path(self, path: str, base_dir: Optional[str] = None) -> str:

        p = Path(path)
        if self.log_relative_paths:
            try:
                if base_dir:
                    relative_path = os.path.relpath(p, base_dir)
                    return relative_path
                else:
                    relative_path = os.path.relpath(p, Path.cwd())
                    return relative_path
            except Exception as e:
                
                self.logger.warning(f"Failed to compute relative path for {path}: {e}")
                return str(p.resolve())
        else:
            return str(p.resolve())



    def analyze_merge_SE_info(
        self,
        save_info_tsv: Optional[str] = None,
        save_info_pkl: Optional[str] = None
    ) -> None:
        """
        Create and store additional merged SE information DataFrame.
        
        This method uses the internally stored `temp_concat_full` DataFrame and `bed_filter`
        (the p2g filtered BedTool object) as inputs to call the function
        `return_merge_se_count_df_full_with_info` from tools.py.
        
        If `temp_concat_full` is not available, a warning is logged and the method
        `analyze_merge_SE()` is automatically executed to generate it before continuing.
        
        Optionally, the resulting DataFrame can be saved as a TSV file and/or as a pickle file.
        """
        self.logger.info("Starting analyze_merge_SE_info process.")

        # Check if _temp_concat_full is available; if not, generate it.
        if self._temp_concat_full is None:
            self.logger.warning("temp_concat_full is not available. Running analyze_merge_SE() to generate temp_concat_full.")
            self.analyze_merge_SE()  # Assumes analyze_merge_SE() populates self._temp_concat_full.
            if self._temp_concat_full is None:
                self.logger.error("Failed to generate temp_concat_full using analyze_merge_SE().")
                raise ValueError("temp_concat_full is not available even after running analyze_merge_SE().")

        try:
            self._temp_full_df_info = return_merge_se_count_df_full_with_info(self._temp_concat_full, self.bed_filter)
            self.logger.info("Temporary full DataFrame info computed successfully.")
            
            if save_info_tsv:
                self._temp_full_df_info.to_csv(save_info_tsv, sep="\t", index=True)
                self.logger.info(f"Temporary full DataFrame info saved to {save_info_tsv}.")
            
            if save_info_pkl:
                self._temp_full_df_info.to_pickle(save_info_pkl)
                self.logger.info(f"Temporary full DataFrame info saved as pickle to {save_info_pkl}.")
        except Exception as e:
            self.logger.exception(f"Error in analyze_merge_SE_info: {e}")
            raise


    def display_se_files(self, use_absolute_paths: bool = False) -> None:
        """
        Display the list of SE files being analyzed.
        
        Args:
            use_absolute_paths: Whether to display absolute paths (True) or
                            relative paths (False, default).
        """
        if not self._input_se_file_list:
            self.logger.warning("No SE files found.")
            print("No SE files found.")
            return
        
        path_type = "absolute" if use_absolute_paths else "relative"
        self.logger.info(f"Displaying list of {len(self._input_se_file_list)} SE files using {path_type} paths.")
        
        original_setting = self.log_relative_paths
        
        if use_absolute_paths:
            self.log_relative_paths = False
        
        try:
            for i, file_path in enumerate(self._input_se_file_list, 1):
                print(f"{i}. {self._get_log_path(str(file_path))}")
        finally:
            self.log_relative_paths = original_setting

    def search_gene_enhancer_links(
        self, 
        gene_symbol: str, 
        display_full_info: bool = True,
        save_tsv: bool = False,
        output_dir: str = './output',
        output_prefix: Optional[str] = None,
        store_bed_filter: bool = True
    ) -> None:
        """
        Search for enhancer links associated with a specific gene symbol.
        
        Args:
            gene_symbol (str): Gene symbol to search for.
            display_full_info (bool): Whether to display all columns or just the key information.
            save_tsv (bool): Whether to save results as TSV files.
            output_dir (str): Directory to save TSV files (default: './output').
            output_prefix (Optional[str]): Prefix for output filenames. If None, uses gene_symbol.
            store_bed_filter (bool): Whether to store the filtered BED object internally for later use.
        """
        self.logger.info(f"Searching for enhancer links for gene '{gene_symbol}'.")
        
        try:
            if not self._bed_filter:
                self.logger.error("P2G data is not filtered. Please run filter_p2g_file first.")
                print("Error: P2G data is not filtered. Please run filter_p2g_file first.")
                return
            
            # Get gene location information
            gene_info = self._get_gene_info(gene_symbol)
            
            bed_df = self._bed_filter.to_dataframe()
            if "score" in bed_df.columns and "symbol" not in bed_df.columns:
                bed_df = bed_df.rename(columns={
                    "name": "PeakID", 
                    "score": "symbol",
                    "strand": "strand"
                })
            
            gene_enhancers = bed_df[bed_df["symbol"] == gene_symbol]
            
            if gene_enhancers.empty:
                self.logger.info(f"No enhancer links found for gene '{gene_symbol}' with current filter settings (FDR≤{self._FDR}, r≥{self._r}).")
                
                # Store search information even if no results found
                self._gene_enhancer_results[gene_symbol] = {
                    "gene_symbol": gene_symbol,
                    "gene_enhancers": None,
                    "enriched_df": None,
                    "search_params": {
                        "FDR": self._FDR,
                        "r": self._r
                    },
                    "timestamp": datetime.datetime.now().isoformat(),
                    "found": False,
                    "gene_info": gene_info
                }
                self._last_searched_gene = gene_symbol
                
                # Display gene information even if no enhancers found
                if gene_info:
                    print(f"Search gene: {gene_symbol} - Location: {gene_info.get('chr', 'Unknown')}:{gene_info.get('start', 'Unknown')}-{gene_info.get('end', 'Unknown')}")
                
                print(f"No enhancer links found for gene '{gene_symbol}' with current filter settings (FDR≤{self._FDR}, r≥{self._r}).")
                return


            enriched_df = None
            if display_full_info and hasattr(self, 'p2g_file') and self.p2g_file:
                try:
                    original_p2gl = pd.read_table(self.p2g_file)
                    
                    peak_ids = gene_enhancers["PeakID"].tolist()
                    
                    detailed_data = original_p2gl[
                        (original_p2gl["symbol"] == gene_symbol) & 
                        (original_p2gl["PeakID"].isin(peak_ids))
                    ].copy()
                    
                    if not detailed_data.empty:
                        self.logger.debug(f"Found {len(detailed_data)} entries in P2GL file matching gene '{gene_symbol}'")
                        self.logger.debug(f"Using filter settings from filter_p2g_file: FDR≤{self._FDR}, r≥{self._r}")
                        
                        display_columns = ["chr", "Start", "End", "PeakID", "symbol", "FDR", "r"]
                        
                        self._add_signed_distance(detailed_data, gene_info)
                        
                        if 'signed_distance' in detailed_data.columns:
                            display_columns.append("signed_distance")
                            
                        min_fdr = detailed_data["FDR"].min() if "FDR" in detailed_data.columns else "N/A"
                        max_fdr = detailed_data["FDR"].max() if "FDR" in detailed_data.columns else "N/A"
                        min_r = detailed_data["r"].min() if "r" in detailed_data.columns else "N/A"
                        max_r = detailed_data["r"].max() if "r" in detailed_data.columns else "N/A"
                        self.logger.debug(f"FDR range in matched data: {min_fdr} to {max_fdr}")
                        self.logger.debug(f"r value range in matched data: {min_r} to {max_r}")
                        
                        enriched_df = detailed_data[display_columns].sort_values(by=["chr", "Start"])
                except Exception as e:
                    self.logger.warning(f"Could not load additional information from original P2GL file: {e}")
                    self.logger.debug(f"Exception details: {str(e)}")

            # enriched_df = None
            # if display_full_info and hasattr(self, 'p2g_file') and self.p2g_file:
            #     try:
            #         original_p2gl = pd.read_table(self.p2g_file)
                    
            #         original_gene_data = original_p2gl[original_p2gl["symbol"] == gene_symbol]
                    
            #         original_gene_data = original_gene_data[
            #             (original_gene_data["FDR"] <= self._FDR) & 
            #             (original_gene_data["r"] >= self._r)
            #         ]
                    
            #         if "PeakID" in original_gene_data.columns:
            #             peak_ids = gene_enhancers["PeakID"].tolist()
            #             original_gene_data = original_gene_data[original_gene_data["PeakID"].isin(peak_ids)]
                        
            #             if not original_gene_data.empty:
            #                 display_columns = ["chr", "Start", "End", "PeakID", "symbol", "FDR", "r"]
                            
            #                 # Always calculate and add signed distance
            #                 self._add_signed_distance(original_gene_data, gene_info)
                            
            #                 # Add signed_distance column to display if it was successfully added
            #                 if 'signed_distance' in original_gene_data.columns:
            #                     display_columns.append("signed_distance")
                                
            #                 enriched_df = original_gene_data[display_columns].sort_values(by=["chr", "Start"])
            #     except Exception as e:
            #         self.logger.warning(f"Could not load additional information from original P2GL file: {e}")

            # Store search results
            self._gene_enhancer_results[gene_symbol] = {
                "gene_symbol": gene_symbol,
                "gene_enhancers": gene_enhancers.copy() if not gene_enhancers.empty else None,
                "enriched_df": enriched_df.copy() if enriched_df is not None and not enriched_df.empty else None,
                "search_params": {
                    "FDR": self._FDR,
                    "r": self._r
                },
                "timestamp": datetime.datetime.now().isoformat(),
                "found": True,
                "gene_info": gene_info
            }
            self._last_searched_gene = gene_symbol

            # Display results with gene information first
            self._display_gene_enhancer_results(
                gene_symbol=gene_symbol,
                gene_enhancers=gene_enhancers,
                enriched_df=enriched_df,
                display_full_info=display_full_info,
                gene_info=gene_info
            )
            
            # Save results as TSV if requested
            if save_tsv:
                self._save_gene_enhancer_results_tsv(
                    gene_symbol=gene_symbol,
                    gene_enhancers=gene_enhancers,
                    enriched_df=enriched_df,
                    output_dir=output_dir,
                    output_prefix=output_prefix
                )
            
            # Store filtered BED object if requested
            if store_bed_filter and not gene_enhancers.empty:
                # Create a BED filter containing only this gene's enhancers
                try:
                    self._gene_bed_filter = BedTool.from_dataframe(gene_enhancers)
                    self.logger.info(f"Stored filtered BED object for gene '{gene_symbol}' with {len(gene_enhancers)} enhancers.")
                except Exception as e:
                    self.logger.error(f"Error creating gene-specific BED filter: {e}")
                    print(f"Error creating gene-specific BED filter: {e}")
                        
        except Exception as e:
            self.logger.exception(f"Error in search_gene_enhancer_links: {e}")
            print(f"Error searching for enhancer links: {e}")

    def _get_gene_info(self, gene_symbol: str) -> dict:
        """
        Get gene location information from RNA info file.
        
        Args:
            gene_symbol (str): Gene symbol to look up
            
        Returns:
            dict: Gene location information
            
        Raises:
            ValueError: If RNA info file is not available or gene not found
        """
        # Check if RNA info file is available
        if not hasattr(self, 'rna_info_file') or not self.rna_info_file:
            self.logger.error("RNA info file is not specified. Cannot get gene location information.")
            raise ValueError("RNA info file is not specified. Cannot get gene location information.")
        
        try:
            
            df_rna = csv_path_to_df_RNA(self.rna_info_file)
            dic_rna = rna_data_to_dic_metadata(df_rna, gtf=True)
            
            if gene_symbol not in dic_rna:
                self.logger.error(f"Gene '{gene_symbol}' not found in RNA info file.")
                raise ValueError(f"Gene '{gene_symbol}' not found in RNA info file.")
            
            gene_data = dic_rna[gene_symbol]
            gene_info = {
                'chr': gene_data.get('chr'),
                'start': gene_data.get('start'),
                'end': gene_data.get('end'),
                'strand': gene_data.get('strand', '.'),
                'source': 'RNA info file'
            }
            self.logger.debug(f"Gene info for {gene_symbol} found in RNA info file.")
            return gene_info
        except Exception as e:
            self.logger.exception(f"Error loading gene info from RNA info file: {e}")
            raise ValueError(f"Error loading gene info from RNA info file: {e}")

    def _add_signed_distance(self, df: pd.DataFrame, gene_info: dict) -> None:
        """
        Add signed distance information to enhancer DataFrame based on gene location.
        Always calculates distance dynamically from genomic coordinates.
        
        Args:
            df (pd.DataFrame): DataFrame containing enhancer information
            gene_info (dict): Gene location information
        """
        if 'Start' not in df.columns or 'End' not in df.columns:
            self.logger.warning("Cannot add signed distance: enhancer position columns not found in DataFrame.")
            return
        
        if not gene_info or 'start' not in gene_info or 'end' not in gene_info or 'strand' not in gene_info:
            self.logger.warning("Cannot add signed distance: gene location information incomplete.")
            return
        
        # Calculate enhancer midpoint
        df['peak_mid'] = (df['Start'] + df['End']) // 2
        
        # Get gene TSS based on strand
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        gene_strand = gene_info['strand']
        
        # TSS is at start for '+' strand, at end for '-' strand
        tss = gene_start if gene_strand == '+' else gene_end
        
        # Calculate signed distance for each enhancer
        signed_distance = []
        
        for _, row in df.iterrows():
            peak_mid = row['peak_mid']
            
            # Always calculate absolute distance from TSS to enhancer midpoint
            dist_value = abs(peak_mid - tss)
            
            # Determine sign based on strand and relative position
            if gene_strand == '+':
                # For + strand: positive if enhancer is downstream of TSS, negative if upstream
                sign = 1 if peak_mid > tss else -1
            else:  # strand == '-'
                # For - strand: negative if enhancer is downstream of TSS, positive if upstream
                sign = -1 if peak_mid > tss else 1
            
            signed_distance.append(sign * dist_value)
        
        # Add the calculated distances to the DataFrame
        df['signed_distance'] = signed_distance
        
        # Remove temporary column
        if 'peak_mid' in df.columns:
            df.drop('peak_mid', axis=1, inplace=True)



    def _display_gene_enhancer_results(
        self, 
        gene_symbol: str, 
        gene_enhancers: pd.DataFrame,
        enriched_df: Optional[pd.DataFrame] = None,
        display_full_info: bool = True,
        gene_info: Optional[dict] = None
    ) -> None:
        """
        Helper method to display gene enhancer results.
        This separates display logic from data processing logic.
        
        Args:
            gene_symbol (str): Gene symbol
            gene_enhancers (pd.DataFrame): DataFrame containing enhancer information
            enriched_df (Optional[pd.DataFrame]): DataFrame with detailed enhancer information
            display_full_info (bool): Whether to display all columns
            gene_info (Optional[dict]): Gene location information
        """
        # Display gene information on separate lines for better readability
        print(f"Search gene: {gene_symbol}")
        
        if gene_info:
            location = f"{gene_info.get('chr', 'Unknown')}:{gene_info.get('start', 'Unknown')}-{gene_info.get('end', 'Unknown')}"
            strand = gene_info.get('strand', 'Unknown')
            strand_symbol = "+" if strand == "+" else "-" if strand == "-" else strand
            
            # Determine TSS based on strand
            if strand == "+":
                tss = gene_info.get('start', 'Unknown')
            elif strand == "-":
                tss = gene_info.get('end', 'Unknown')
            else:
                tss = 'Unknown'
                
            print(f"Location: {location}")
            print(f"TSS: {tss}")
            print(f"Strand: {strand_symbol}")
        else:
            print("Location information not available")
        
        if gene_enhancers is None or gene_enhancers.empty:
            print(f"No enhancer links found for gene '{gene_symbol}'.")
            return
            
        print(f"\nFOUND {len(gene_enhancers)} ENHANCERS LINKED TO GENE '{gene_symbol}':")
        print("-" * 50)
        
        # Display detailed information if available and requested
        if display_full_info and enriched_df is not None and not enriched_df.empty:
            # Rename columns for better readability
            renamed_df = enriched_df.copy()
            if 'signed_distance' in renamed_df.columns:
                renamed_df.rename(columns={'signed_distance': 'distance(bp)'}, inplace=True)
            elif 'distance' in renamed_df.columns:
                renamed_df.rename(columns={'distance': 'distance(bp)'}, inplace=True)
            
            # Reset index and add 1 to make it 1-based
            display_df = renamed_df.reset_index(drop=True)
            display_df.index += 1
            
            display(display_df)
        else:
            # Display basic enhancer information if detailed info not available or not requested
            gene_enhancers_sorted = gene_enhancers.sort_values(by=['chrom', 'start'])
            
            display_columns = ['chrom', 'start', 'end', 'PeakID']
            if display_full_info:
                display_df = gene_enhancers_sorted
            else:
                display_df = gene_enhancers_sorted[display_columns]
            
            # Reset index and add 1 to make it 1-based
            display_df = display_df.reset_index(drop=True)
            display_df.index += 1
            
            display(display_df)
        
        # Always show genomic coordinates for copy/paste
        print("\nGENOMIC COORDINATES FOR COPY/PASTE:")
        print("-" * 50)
        
        # Use sorted data frame for consistent order with displayed table
        sorted_enhancers = gene_enhancers.sort_values(by=['chrom', 'start']).reset_index(drop=True)
        
        # Add index (1-based) to each coordinate for reference
        for i, row in enumerate(sorted_enhancers.itertuples(), 1):
            # Format index with leading zeros for better readability when copying
            formatted_index = f"{i:02d}"
            print(f"{formatted_index}: {row.chrom}:{row.start}-{row.end}")





    def _save_gene_enhancer_results_tsv(
        self,
        gene_symbol: str,
        gene_enhancers: pd.DataFrame,
        enriched_df: Optional[pd.DataFrame] = None,
        output_dir: str = './output',
        output_prefix: Optional[str] = None
    ) -> None:
        """
        Helper method to save gene enhancer results as TSV files.
        
        Args:
            gene_symbol (str): Gene symbol.
            gene_enhancers (pd.DataFrame): Basic enhancer information.
            enriched_df (Optional[pd.DataFrame]): Detailed enhancer information with FDR and r values.
            output_dir (str): Directory to save TSV files.
            output_prefix (Optional[str]): Prefix for output filenames.
        """
        try:
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate timestamp for unique filenames
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Use provided prefix or default to gene_symbol
            prefix = output_prefix if output_prefix else gene_symbol
            
            # Save basic enhancer information
            if gene_enhancers is not None and not gene_enhancers.empty:
                basic_file = os.path.join(output_dir, f"{prefix}_enhancers_{timestamp}.tsv")
                gene_enhancers.to_csv(basic_file, sep='\t', index=False)
                self.logger.info(f"Basic enhancer data saved to {basic_file}")
                print(f"Basic enhancer data saved to {basic_file}")
            
            # Save detailed enhancer information if available
            if enriched_df is not None and not enriched_df.empty:
                # Rename distance column for clarity
                save_df = enriched_df.copy()
                if 'signed_distance' in save_df.columns:
                    save_df.rename(columns={'signed_distance': 'distance(bp)'}, inplace=True)
                    
                detailed_file = os.path.join(output_dir, f"{prefix}_enhancers_detailed_{timestamp}.tsv")
                save_df.to_csv(detailed_file, sep='\t', index=False)
                self.logger.info(f"Detailed enhancer data saved to {detailed_file}")
                print(f"Detailed enhancer data saved to {detailed_file}")
            
            # Save genomic coordinates as a text file for convenient copy-paste
            if gene_enhancers is not None and not gene_enhancers.empty:
                coords_file = os.path.join(output_dir, f"{prefix}_coordinates_{timestamp}.txt")
                with open(coords_file, 'w') as f:
                    for _, row in gene_enhancers.iterrows():
                        f.write(f"{row['chrom']}:{row['start']}-{row['end']}\n")
                self.logger.info(f"Genomic coordinates saved to {coords_file}")
                print(f"Genomic coordinates saved to {coords_file}")
        
        except Exception as e:
            self.logger.error(f"Error saving enhancer results to TSV: {e}")
            print(f"Error saving enhancer results to TSV: {e}")

    def save_gene_search_results(
        self,
        gene_symbol: Optional[str] = None,
        output_dir: str = './output',
        output_prefix: Optional[str] = None
    ) -> None:
        """
        Save previously searched gene enhancer results to TSV files.
        
        Args:
            gene_symbol (Optional[str]): Gene symbol to save results for. 
                                        If None, uses the last searched gene.
            output_dir (str): Directory to save TSV files.
            output_prefix (Optional[str]): Prefix for output filenames.
        """
        self.logger.info(f"Saving gene search results for gene_symbol={gene_symbol if gene_symbol else 'last searched gene'}.")
        
        try:
            if gene_symbol is None:
                if self._last_searched_gene is None:
                    self.logger.error("No gene has been searched. Please specify a gene symbol.")
                    print("Error: No gene has been searched. Please specify a gene symbol.")
                    return
                gene_symbol = self._last_searched_gene
            
            if gene_symbol not in self._gene_enhancer_results:
                self.logger.error(f"No results found for gene '{gene_symbol}'. Please search for this gene first.")
                print(f"Error: No results found for gene '{gene_symbol}'. Please search for this gene first.")
                return
            
            result_data = self._gene_enhancer_results[gene_symbol]
            
            if not result_data.get("found", False):
                self.logger.warning(f"No enhancer links were found for gene '{gene_symbol}' during the search.")
                print(f"Warning: No enhancer links were found for gene '{gene_symbol}' during the search.")
                return
            
            self._save_gene_enhancer_results_tsv(
                gene_symbol=result_data["gene_symbol"],
                gene_enhancers=result_data["gene_enhancers"],
                enriched_df=result_data["enriched_df"],
                output_dir=output_dir,
                output_prefix=output_prefix
            )
            
            # Save search parameters
            params_file = os.path.join(output_dir, f"{output_prefix or gene_symbol}_search_params_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
            with open(params_file, 'w') as f:
                f.write(f"Gene: {gene_symbol}\n")
                f.write(f"Search timestamp: {result_data['timestamp']}\n")
                f.write(f"FDR threshold: {result_data['search_params']['FDR']}\n")
                f.write(f"r threshold: {result_data['search_params']['r']}\n")
                if result_data["gene_enhancers"] is not None:
                    f.write(f"Number of enhancers found: {len(result_data['gene_enhancers'])}\n")
                else:
                    f.write("Number of enhancers found: 0\n")
            
            self.logger.info(f"Search parameters saved to {params_file}")
            print(f"Search parameters saved to {params_file}")
            
        except Exception as e:
            self.logger.exception(f"Error in save_gene_search_results: {e}")
            print(f"Error saving gene search results: {e}")





    def find_gene_name_linked_super_enhancers(self, gene_symbol):
        """
        Find super enhancers linked to a specific gene and display the results.
        
        Parameters:
        -----------
        gene_symbol : str
            Gene symbol to search for
        
        Returns:
        --------
        pd.DataFrame
            DataFrame containing information about super enhancers linked to the specified gene
        """
        self.logger.info(f"Finding super enhancers linked to gene '{gene_symbol}'")
        
        try:
            # 1. Create gene-specific BED filter using the existing method
            self.search_gene_enhancer_links(
                gene_symbol=gene_symbol,
                display_full_info=False,  # Suppress extra information display
                save_tsv=False,
                store_bed_filter=True     # Important: This sets _gene_bed_filter
            )
            
            # Check if _gene_bed_filter was properly set
            if not hasattr(self, '_gene_bed_filter') or self._gene_bed_filter is None:
                self.logger.error(f"Could not create BED filter for gene '{gene_symbol}'")
                print(f"Error: No enhancers found for gene '{gene_symbol}'")
                return None
            
            # 2. Call find_gene_linked_super_enhancers function
            
            if not self.rose_file:
                self.logger.error("No ROSE file specified")
                print("Error: No ROSE file specified. Please set a ROSE file first.")
                return None
            
            self.logger.info(f"Filtering ROSE file '{self.rose_file}' for gene-linked super enhancers")
            result_df = find_gene_linked_super_enhancers(self.rose_file, self._gene_bed_filter)
            
            # 3. Display results
            if result_df.empty:
                self.logger.info(f"No super enhancers linked to gene '{gene_symbol}' found")
                print(f"No super enhancers linked to gene '{gene_symbol}' found in the current ROSE file.")
                return result_df
            
            # Create display dataframe with only necessary columns
            display_cols = ['enhancerRank', 'total_super_enhancers', 'REGION_ID', 
                            'CHROM', 'START', 'STOP', 'linked_genes']
            
            if 'linked_peaks' in result_df.columns:
                display_cols.append('linked_peaks')
            
            display_df = result_df[display_cols].copy()
            
            # Display summary information
            total_ses = result_df['total_super_enhancers'].iloc[0]
            found_ses = len(result_df)
            
            print(f"Search results for super enhancers linked to gene '{gene_symbol}':")
            print(f"Total super enhancers in ROSE file: {total_ses}")
            print(f"Super enhancers linked to gene '{gene_symbol}': {found_ses} ({found_ses/total_ses*100:.1f}%)")
            print("")
            
            # Display the dataframe
            display(display_df)
            
            return result_df
            
        except Exception as e:
            self.logger.exception(f"Error in find_gene_name_linked_super_enhancers: {e}")
            print(f"Error: {e}")
            return None













    ##### find_gene_linked_super_enhancers #####

    def find_gene_linked_super_enhancers_in_directory(
        self,
        gene_symbol: str,
        output_dir: Optional[str] = None,
        save_results: bool = False,
        output_prefix: Optional[str] = None,
        file_format: str = "csv",
        verbose: bool = True
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
        """
        Finds super enhancers linked to the specified gene across all samples in the ROSE directory.
        
        Parameters:
        -----------
        gene_symbol : str
            Gene symbol to search for.
        output_dir : Optional[str]
            Directory to save results (default: './output').
        save_results : bool
            Whether to save results to files.
        output_prefix : Optional[str]
            Prefix for output filenames (default: gene_symbol).
        file_format : str
            File format for results: "csv" or "tsv" (default: "csv").
        verbose : bool
            Whether to print progress and summary information.
            
        Returns:
        --------
        Tuple[pd.DataFrame, pd.DataFrame, Dict]
            - Combined results from all samples
            - Summary information for all samples
            - Statistics summary with sample counts and lists
        """
        self.logger.info(f"Finding super enhancers linked to gene '{gene_symbol}' across all samples in ROSE directory.")
        
        try:
            # Create gene-specific BED filter
            self.search_gene_enhancer_links(
                gene_symbol=gene_symbol,
                display_full_info=False,  # Suppress additional information display
                save_tsv=False,
                store_bed_filter=True     # This sets _gene_bed_filter
            )
            
            # Check if gene_bed_filter was created
            if not hasattr(self, '_gene_bed_filter') or self._gene_bed_filter is None:
                self.logger.error(f"Could not create BED filter for gene '{gene_symbol}'")
                print(f"Error: No enhancers found for gene '{gene_symbol}'")
                return pd.DataFrame(), pd.DataFrame(), {}
            
            # Verify ROSE files directory is valid
            if not self.rose_files_dir or not os.path.exists(self.rose_files_dir):
                self.logger.error("ROSE files directory is not specified or does not exist.")
                print("Error: ROSE files directory is not specified or does not exist.")
                return pd.DataFrame(), pd.DataFrame(), {}
            
            # Call the function with directory and BED filter
            self.logger.info(f"Processing all samples in directory: {self.rose_files_dir}")
            combined_results, samples_df, stats_summary = find_gene_linked_super_enhancers_in_directory(
                self.rose_files_dir, 
                self._gene_bed_filter
            )
            
            # Filter to get only super enhancers (isSuper == 1)
            super_only_results = pd.DataFrame()
            if not combined_results.empty and 'isSuper' in combined_results.columns:
                super_only_results = combined_results[combined_results['isSuper'] == 1].copy()


                # Make sure isSuper is also at the rightmost position in super_only_results
                if 'isSuper' in super_only_results.columns:
                    isSuper_col = super_only_results.pop('isSuper')
                    super_only_results['isSuper'] = isSuper_col

            
            # Store results in internal variables
            self._se_directory_results[gene_symbol] = {
                "gene_symbol": gene_symbol,
                "combined_results": combined_results,
                "super_only_results": super_only_results,
                "samples_df": samples_df,
                "stats_summary": stats_summary,
                "timestamp": datetime.datetime.now().isoformat()
            }
            self._last_searched_gene_directory = gene_symbol
            
            # Save results to files if requested
            if save_results:
                self.save_directory_search_results(
                    gene_symbol=gene_symbol,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    file_format=file_format
                )
            
            # Display summary if verbose is True
            if verbose:
                self._display_directory_search_summary(gene_symbol)
            
            return combined_results, samples_df, stats_summary
            
        except Exception as e:
            self.logger.exception(f"Error in find_gene_linked_super_enhancers_in_directory: {e}")
            print(f"Error searching for super enhancers: {e}")
            return pd.DataFrame(), pd.DataFrame(), {}

    def _display_directory_search_summary(self, gene_symbol: str) -> None:
        """
        Helper method to display a summary of directory search results.
        
        Parameters:
        -----------
        gene_symbol : str
            Gene symbol to display results for
        """
        if gene_symbol not in self._se_directory_results:
            print(f"No directory search results found for gene '{gene_symbol}'.")
            return
        
        result_data = self._se_directory_results[gene_symbol]
        stats = result_data["stats_summary"]
        
        print("\n" + "="*80)
        print(f"SUPER ENHANCER SEARCH SUMMARY FOR GENE: {gene_symbol}")
        print("="*80)
        
        print(f"\nTotal samples processed: {stats['total_sample_count']}")
        print(f"Samples with any enhancers linked to {gene_symbol}: {stats['all_matched']['count']} ({stats['all_matched']['percentage']}%)")
        print(f"Samples with super enhancers linked to {gene_symbol}: {stats['super_matched']['count']} ({stats['super_matched']['percentage']}%)")
        
        # Display sample list if it's short enough
        if len(stats['super_matched']['sample_list']) <= 10:
            print(f"\nSamples with super enhancers: {', '.join(stats['super_matched']['sample_list'])}")
        else:
            print(f"\nSamples with super enhancers: {', '.join(stats['super_matched']['sample_list'][:5])}... (and {len(stats['super_matched']['sample_list'])-5} more)")
        
        # Display sample DataFrame if available
        if not result_data["samples_df"].empty:
            print("\nSAMPLE SUMMARY TABLE:")
            display(result_data["samples_df"])
        
        # Display counts of super enhancers
        super_only = result_data.get("super_only_results", pd.DataFrame())
        if not super_only.empty:
            super_count = len(super_only)
            super_samples = super_only['sample_name'].nunique() if 'sample_name' in super_only.columns else 0
            print(f"\nFound {super_count} super enhancers across {super_samples} samples linked to gene {gene_symbol}")

    def save_directory_search_results(
        self,
        gene_symbol: Optional[str] = None,
        output_dir: str = './output',
        output_prefix: Optional[str] = None,
        file_format: str = "csv"
    ) -> None:
        """
        Save previously searched gene's directory results to files.
        
        Parameters:
        -----------
        gene_symbol : Optional[str]
            Gene symbol to save results for. If None, uses the last searched gene.
        output_dir : str
            Directory to save files.
        output_prefix : Optional[str]
            Prefix for output filenames.
        file_format : str
            File format: "csv" or "tsv" (default: "csv").
        """
        self.logger.info(f"Saving directory search results for gene_symbol={gene_symbol if gene_symbol else 'last searched gene'}.")
        
        try:
            # Validate file format
            if file_format.lower() not in ["csv", "tsv"]:
                self.logger.warning(f"Invalid file format: {file_format}. Using 'csv' instead.")
                file_format = "csv"
            
            # Set delimiter and extension based on file format
            delimiter = "," if file_format.lower() == "csv" else "\t"
            extension = file_format.lower()
            
            if gene_symbol is None:
                if self._last_searched_gene_directory is None:
                    self.logger.error("No gene has been searched. Please specify a gene symbol.")
                    print("Error: No gene has been searched. Please specify a gene symbol.")
                    return
                gene_symbol = self._last_searched_gene_directory
            
            if gene_symbol not in self._se_directory_results:
                self.logger.error(f"No directory results found for gene '{gene_symbol}'. Please search for this gene first.")
                print(f"Error: No directory results found for gene '{gene_symbol}'. Please search for this gene first.")
                return
            
            result_data = self._se_directory_results[gene_symbol]
            
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate timestamp for unique filenames
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Use specified prefix or default to gene symbol
            prefix = output_prefix if output_prefix else gene_symbol
            
            # Save combined results DataFrame
            if not result_data["combined_results"].empty:
                combined_file = os.path.join(output_dir, f"{prefix}_all_enhancers_{timestamp}.{extension}")
                result_data["combined_results"].to_csv(combined_file, sep=delimiter, index=False)
                self.logger.info(f"All enhancers results saved to {combined_file}")
                print(f"All enhancers results saved to {combined_file}")
            
            # Save super enhancers only DataFrame
            if "super_only_results" in result_data and not result_data["super_only_results"].empty:
                super_file = os.path.join(output_dir, f"{prefix}_super_enhancers_only_{timestamp}.{extension}")
                result_data["super_only_results"].to_csv(super_file, sep=delimiter, index=False)
                self.logger.info(f"Super enhancers only results saved to {super_file}")
                print(f"Super enhancers only results saved to {super_file}")
            
            # Save samples summary DataFrame
            if not result_data["samples_df"].empty:
                samples_file = os.path.join(output_dir, f"{prefix}_samples_summary_{timestamp}.{extension}")
                result_data["samples_df"].to_csv(samples_file, sep=delimiter, index=False)
                self.logger.info(f"Samples summary saved to {samples_file}")
                print(f"Samples summary saved to {samples_file}")
            
            # Save statistics summary as JSON
            stats_file = os.path.join(output_dir, f"{prefix}_stats_summary_{timestamp}.json")
            with open(stats_file, 'w') as f:
                # Convert sets to lists for JSON serialization
                json.dump(result_data["stats_summary"], f, indent=2)
            self.logger.info(f"Statistics summary saved to {stats_file}")
            print(f"Statistics summary saved to {stats_file}")
            
        except Exception as e:
            self.logger.exception(f"Error in save_directory_search_results: {e}")
            print(f"Error saving directory search results: {e}")

    def get_directory_super_enhancers(self, gene_symbol: str) -> Optional[pd.DataFrame]:
        """
        Returns the Super Enhancers only (isSuper==1) DataFrame for the specified gene.
        
        Parameters:
        -----------
        gene_symbol : str
            Gene symbol to retrieve Super Enhancers for
            
        Returns:
        --------
        Optional[pd.DataFrame]
            DataFrame containing only Super Enhancers, or None if no search results exist for the gene
        """
        if gene_symbol not in self._se_directory_results:
            self.logger.warning(f"No directory search results found for gene '{gene_symbol}'.")
            return None
        
        result = self._se_directory_results.get(gene_symbol)
        if "super_only_results" in result:
            return result["super_only_results"]
        return None










    ##### find super enhancer only #####


    def find_gene_linked_super_enhancers_only(
        self,
        gene_symbol: str,
        output_dir: Optional[str] = None,
        save_results: bool = False,
        output_prefix: Optional[str] = None,
        file_format: str = "csv",
        verbose: bool = True,
        log_level: str = "INFO"
    ) -> pd.DataFrame:
        """
        Find only super enhancers (isSuper == 1) linked to the specified gene
        across all samples in the ROSE directory.
        The order of display is:
        (1) "P2GL Search Results..." line
        (2) P2GL table (via search_gene_enhancer_links)
        (3) "Super Enhancer Search Results..." line
        (4) super enhancer table
        """
        self.logger.info(f"Finding only super enhancers linked to gene '{gene_symbol}' across all samples in ROSE directory.")
        try:
            from IPython.display import Markdown, display
            
            # (1) Show "P2GL Search Results..." first
            display(Markdown(f"**P2GL Search Results for gene: {gene_symbol}**"))
            
            # (2) Call search_gene_enhancer_links without redirect_stdout
            #     so that the P2GL table is displayed as-is.
            self.search_gene_enhancer_links(
                gene_symbol=gene_symbol,
                display_full_info=False,
                save_tsv=False,
                store_bed_filter=True
            )
            
            # Check if the gene-specific BedTool was created
            if not hasattr(self, '_gene_bed_filter') or self._gene_bed_filter is None:
                self.logger.error(f"Could not create BED filter for gene '{gene_symbol}'")
                if verbose:
                    self.logger.error(f"No enhancers found for gene '{gene_symbol}'")
                return pd.DataFrame()
            
            # Check ROSE directory
            if not self.rose_files_dir or not os.path.exists(self.rose_files_dir):
                self.logger.error("ROSE files directory is not specified or does not exist.")
                if verbose:
                    self.logger.error("ROSE files directory is not specified or does not exist.")
                return pd.DataFrame()
            
            # Temporarily change log level
            original_level = self.logger.level
            try:
                if log_level:
                    level_value = getattr(logging, log_level.upper())
                    self.logger.setLevel(level_value)
                
                # Use the existing function to get combined_results etc.
                combined_results, samples_df, stats_summary = self._find_gene_linked_super_enhancers_in_directory(
                    self.rose_files_dir,
                    self._gene_bed_filter,
                    verbose=verbose
                )
            finally:
                self.logger.setLevel(original_level)
            
            # Filter only super enhancers
            super_only_results = pd.DataFrame()
            if not combined_results.empty and 'isSuper' in combined_results.columns:
                super_only_results = combined_results[combined_results['isSuper'] == 1].copy()
                if 'isSuper' in super_only_results.columns:
                    isSuper_col = super_only_results.pop('isSuper')
                    super_only_results['isSuper'] = isSuper_col
            
            # Store results internally
            self._se_directory_results[gene_symbol] = {
                "gene_symbol": gene_symbol,
                "combined_results": combined_results,
                "super_only_results": super_only_results,
                "samples_df": samples_df,
                "stats_summary": stats_summary,
                "timestamp": datetime.datetime.now().isoformat()
            }
            self._last_searched_gene_directory = gene_symbol
            
            # Save if requested
            if save_results and not super_only_results.empty:
                self.save_super_enhancers_only(
                    gene_symbol=gene_symbol,
                    output_dir=output_dir,
                    output_prefix=output_prefix,
                    file_format=file_format
                )
            
            if verbose:
                self._display_super_enhancers_only_summary(gene_symbol)
            
            # (3) Show "Super Enhancer Search Results..." line
            display(Markdown(f"**Super Enhancer Search Results for gene: {gene_symbol}**"))
            
            # (4) Display super enhancers table with 1-based index
            if not super_only_results.empty:
                display_columns = [
                    'sample_name', 'REGION_ID', 'CHROM', 'START', 'STOP',
                    'NUM_LOCI', 'SIGNAL', 'enhancerRank', 'super_percentile'
                ]
                available_cols = [col for col in display_columns if col in super_only_results.columns]
                display_df = super_only_results[available_cols].reset_index(drop=True)
                display_df.index += 1
                display(display_df)

            return super_only_results
        
        except Exception as e:
            self.logger.exception(f"Error in find_gene_linked_super_enhancers_only: {e}")
            if verbose:
                self.logger.error(f"Error searching for super enhancers: {e}")
            return pd.DataFrame()






    def _find_gene_linked_super_enhancers_in_directory(
        self,
        directory_path: str,
        gene_bed_filter: BedTool,
        verbose: bool = True
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Dict]:
        """
        Internal method that wraps find_gene_linked_super_enhancers_in_directory function 
        but captures its output and redirects to the logger system.
        
        Parameters:
        -----------
        directory_path : str
            Path to directory containing ROSE files (*_AllEnhancers.table.txt)
        gene_bed_filter : BedTool
            BedTool object containing gene filter information
        verbose : bool
            Whether to log detailed information
            
        Returns:
        --------
        tuple of (pd.DataFrame, pd.DataFrame, dict)
            - First DataFrame: Combined results from all samples
            - Second DataFrame: Summary information for all samples
            - Third dict: Simple statistical summary with sample counts and lists
        """
        # Create a StringIO object to capture stdout
        output_buffer = StringIO()
        
        # Use a context manager to redirect stdout temporarily
        with contextlib.redirect_stdout(output_buffer):
            # Call the original function
            combined_results, samples_df, stats_summary = find_gene_linked_super_enhancers_in_directory(
                directory_path, gene_bed_filter
            )
        
        # Get the captured output and split into lines
        captured_output = output_buffer.getvalue()
        output_lines = captured_output.strip().split('\n')
        
        # Log each line that was printed
        for line in output_lines:
            if line.strip():  # Skip empty lines
                # Log different types of lines with appropriate log levels
                if "Error" in line or "Exception" in line:
                    self.logger.error(line)
                elif "Warning" in line:
                    self.logger.warning(line)
                elif verbose:  # Only log detailed info if verbose is True
                    self.logger.info(line)
        
        return combined_results, samples_df, stats_summary

    def _display_super_enhancers_only_summary(self, gene_symbol: str) -> None:
        """
        Helper method to display a summary of super enhancers only for a gene.
        Uses logger instead of print statements.
        
        Parameters:
        -----------
        gene_symbol : str
            Gene symbol to display super enhancers summary for
        """
        if gene_symbol not in self._se_directory_results:
            self.logger.info(f"No directory search results found for gene '{gene_symbol}'.")
            return
        
        result_data = self._se_directory_results[gene_symbol]
        
        super_only = result_data.get("super_only_results", pd.DataFrame())
        if super_only.empty:
            self.logger.info(f"No super enhancers (isSuper==1) found for gene '{gene_symbol}'.")
            return
                
        # Count stats
        super_count = len(super_only)
        super_samples = super_only['sample_name'].nunique() if 'sample_name' in super_only.columns else 0
        
        self.logger.info("\n" + "="*80)
        self.logger.info(f"SUPER ENHANCERS ONLY SUMMARY FOR GENE: {gene_symbol}")
        self.logger.info("="*80)
        
        self.logger.info(f"\nFound {super_count} super enhancers across {super_samples} samples linked to gene {gene_symbol}")
        
        # List samples with super enhancers if not too many
        if super_samples <= 10 and 'sample_name' in super_only.columns:
            sample_list = sorted(super_only['sample_name'].unique())
            self.logger.info(f"Samples with super enhancers: {', '.join(sample_list)}")
        
        # Show summary table by sample if available
        if 'sample_name' in super_only.columns:
            sample_counts = super_only.groupby('sample_name').size().reset_index(name='count')
            sample_counts = sample_counts.sort_values('count', ascending=False)
            self.logger.info("\nSUPER ENHANCERS PER SAMPLE:")
            # For each sample, log the count (as tables can't be directly logged)
            for _, row in sample_counts.iterrows():
                self.logger.info(f"  {row['sample_name']}: {row['count']} super enhancers")
        
        # Show first few super enhancers if available
        if not super_only.empty:
            self.logger.info("\nSAMPLE SUPER ENHANCERS (first 5):")
            display_columns = ['sample_name', 'CHROM', 'START', 'STOP', 'REGION_ID', 'enhancerRank']
            display_cols = [col for col in display_columns if col in super_only.columns]
            
            # Show only first 5 rows
            for i, (_, row) in enumerate(super_only[display_cols].head().iterrows()):
                self.logger.info(f"  {i+1}. {' | '.join([f'{col}: {row[col]}' for col in display_cols])}")






    def save_super_enhancers_only(
        self,
        gene_symbol: Optional[str] = None,
        output_dir: str = './output',
        output_prefix: Optional[str] = None,
        file_format: str = "csv"
    ) -> None:
        """
        Save previously searched gene's super enhancers only (isSuper==1) results to a file.
        
        Parameters:
        -----------
        gene_symbol : Optional[str]
            Gene symbol to save results for. If None, uses the last searched gene.
        output_dir : str
            Directory to save files.
        output_prefix : Optional[str]
            Prefix for output filenames.
        file_format : str
            File format: "csv" or "tsv" (default: "csv").
        """
        self.logger.info(f"Saving super enhancers only for gene_symbol={gene_symbol if gene_symbol else 'last searched gene'}.")
        
        try:
            # Validate file format
            if file_format.lower() not in ["csv", "tsv"]:
                self.logger.warning(f"Invalid file format: {file_format}. Using 'csv' instead.")
                file_format = "csv"
            
            # Set delimiter and extension based on file format
            delimiter = "," if file_format.lower() == "csv" else "\t"
            extension = file_format.lower()
            
            if gene_symbol is None:
                if self._last_searched_gene_directory is None:
                    self.logger.error("No gene has been searched. Please specify a gene symbol.")
                    print("Error: No gene has been searched. Please specify a gene symbol.")
                    return
                gene_symbol = self._last_searched_gene_directory
            
            if gene_symbol not in self._se_directory_results:
                self.logger.error(f"No directory results found for gene '{gene_symbol}'. Please search for this gene first.")
                print(f"Error: No directory results found for gene '{gene_symbol}'. Please search for this gene first.")
                return
            
            result_data = self._se_directory_results[gene_symbol]
            
            # Check if super enhancers only data exists
            if "super_only_results" not in result_data or result_data["super_only_results"].empty:
                self.logger.warning(f"No super enhancers found for gene '{gene_symbol}'.")
                print(f"Warning: No super enhancers found for gene '{gene_symbol}'.")
                return
            
            # Create output directory if it doesn't exist
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate timestamp for unique filenames
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Use specified prefix or default to gene symbol
            prefix = output_prefix if output_prefix else gene_symbol
            
            # Save super enhancers only DataFrame
            super_file = os.path.join(output_dir, f"{prefix}_super_enhancers_only_{timestamp}.{extension}")
            result_data["super_only_results"].to_csv(super_file, sep=delimiter, index=False)
            self.logger.info(f"Super enhancers only results saved to {super_file}")
            print(f"Super enhancers only results saved to {super_file}")
            
        except Exception as e:
            self.logger.exception(f"Error in save_super_enhancers_only: {e}")
            print(f"Error saving super enhancers only results: {e}")





    def plot_super_enhancers_overall_distribution(
        self,
        gene_symbol: Optional[str] = None,
        figsize: Tuple[int, int] = (15, 8),
        colors: List[str] = None,
        title: Optional[str] = None,
        save_fig: bool = False,
        fig_format: str = 'png',
        output_dir: str = './output',
        output_prefix: Optional[str] = None,
        save_data: bool = False,
        data_format: str = 'csv',
        show_plot: bool = True,
        percentile_column: str = 'super_percentile',
        bins: int = 20
    ) -> Optional[Tuple[plt.Figure, plt.Figure, plt.Figure]]:
        """
        Plot distribution of super enhancer percentiles across all samples.
        Creates three separate plots: histogram, boxplot+swarm, and a pie chart showing sample distribution.
        
        Parameters:
        -----------
        gene_symbol : Optional[str]
            Gene symbol to plot. If None, uses the last searched gene.
        figsize : Tuple[int, int]
            Figure size (width, height) in inches
        colors : List[str]
            List of colors for the different plots
        title : Optional[str]
            Custom title for plots. If None, a default title is generated.
        save_fig : bool
            Whether to save the figures
        fig_format : str
            Figure format: 'png', 'svg', or 'pdf' (default: 'png')
        output_dir : str
            Directory to save figure and data
        output_prefix : Optional[str]
            Prefix for output filenames. If None, uses gene symbol.
        save_data : bool
            Whether to save the data used for visualization
        data_format : str
            Data format: 'csv' or 'tsv' (default: 'csv')
        show_plot : bool
            Whether to display the plots in Jupyter
        percentile_column : str
            Column name for percentile data (default: 'super_percentile')
        bins : int
            Number of bins for histogram (default: 20)
                    
        Returns:
        --------
        Optional[Tuple[plt.Figure, plt.Figure, plt.Figure]]
            If show_plot is True, returns the three figures (histogram_fig, boxswarm_fig, pie_fig)
        """
        self.logger.info(f"Plotting super enhancers percentile distribution for gene_symbol={gene_symbol if gene_symbol else 'last searched gene'}.")
        
        try:

            
            # Set default colors if not provided
            if colors is None:
                colors = ['#3498db', '#2ecc71', '#e74c3c', '#9b59b6']
            
            # Get the super enhancers DataFrame
            if gene_symbol is None:
                if self._last_searched_gene_directory is None:
                    self.logger.error("No gene has been searched. Please specify a gene symbol.")
                    print("Error: No gene has been searched. Please specify a gene symbol.")
                    return None
                gene_symbol = self._last_searched_gene_directory
                
            se_df = self.get_directory_super_enhancers(gene_symbol)
            if se_df is None or se_df.empty:
                self.logger.warning(f"No super enhancers found for gene '{gene_symbol}'.")
                print(f"Warning: No super enhancers found for gene '{gene_symbol}'.")
                return None
            
            # Check if percentile column exists - terminate if not found
            if percentile_column not in se_df.columns:
                self.logger.error(f"Column '{percentile_column}' not found in the DataFrame. Available columns: {list(se_df.columns)}")
                print(f"Error: Percentile column '{percentile_column}' not found. This column should be present in the original data.")
                return None
            
            # Extract percentile values
            percentiles = se_df[percentile_column].dropna().values
            
            if len(percentiles) == 0:
                self.logger.warning(f"No valid percentile values found for gene '{gene_symbol}'.")
                print(f"Warning: No valid percentile values found for gene '{gene_symbol}'.")
                return None
            
            # Generate title if not provided
            if title is None:
                title_base = f"Super Enhancer Percentile Distribution for Gene: {gene_symbol}"
            else:
                title_base = title
            
            # Calculate statistics for display
            stats_text = (
                f"Count: {len(percentiles)}\n"
                f"Mean: {np.mean(percentiles):.2f}\n"
                f"Median: {np.median(percentiles):.2f}\n"
                f"Std Dev: {np.std(percentiles):.2f}\n"
                f"Min: {np.min(percentiles):.2f}\n"
                f"Max: {np.max(percentiles):.2f}\n"
                f"25th: {np.percentile(percentiles, 25):.2f}\n"
                f"75th: {np.percentile(percentiles, 75):.2f}"
            )
            
            # ---------- FIGURE 1: HISTOGRAM ----------
            self.logger.info("Creating histogram plot...")
            histogram_fig = plt.figure(figsize=figsize)
            ax = histogram_fig.add_subplot(111)
            
            # Plot histogram
            sns.histplot(percentiles, bins=bins, kde=False, color=colors[0], ax=ax)
            
            # Set labels and title
            ax.set_title(f"{title_base} - Histogram", fontsize=16)
            ax.set_xlabel('Percentile', fontsize=14)
            ax.set_ylabel('Frequency', fontsize=14)
            
            # Add stats text box
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
                    fontsize=10, verticalalignment='top', bbox=props)
            
            # Adjust layout
            histogram_fig.tight_layout()
            
            # Save histogram if requested
            if save_fig:
                os.makedirs(output_dir, exist_ok=True)
                prefix = output_prefix if output_prefix else gene_symbol
                hist_fig_path = os.path.join(output_dir, f"{prefix}_histogram.{fig_format}")
                histogram_fig.savefig(hist_fig_path, format=fig_format, dpi=300, bbox_inches='tight')
                self.logger.info(f"Histogram saved to {hist_fig_path}")
                print(f"Histogram saved to {hist_fig_path}")
            
            # ---------- FIGURE 2: BOXPLOT + SWARM ----------
            self.logger.info("Creating boxplot with swarm plot...")
            boxswarm_fig = plt.figure(figsize=figsize)
            ax = boxswarm_fig.add_subplot(111)
            
            # Create a DataFrame for seaborn
            box_data = pd.DataFrame({percentile_column: percentiles})
            
            # Plot boxplot
            sns.boxplot(y=percentile_column, data=box_data, color=colors[1], ax=ax)
            
            # Add swarm plot on top of the boxplot
            sns.swarmplot(y=percentile_column, data=box_data, color=colors[2], 
                        size=8, alpha=0.7, ax=ax)
            
            # Set labels and title
            ax.set_title(f"{title_base} - Boxplot with Distribution", fontsize=16)
            ax.set_ylabel('Percentile', fontsize=14)
            ax.set_xticks([])  # Remove x-ticks
            
            # Add stats text box
            props = dict(boxstyle='round', facecolor='white', alpha=0.8)
            ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
                    fontsize=10, verticalalignment='top', bbox=props)
            
            # Adjust layout
            boxswarm_fig.tight_layout()
            
            # Save boxplot+swarm if requested
            if save_fig:
                os.makedirs(output_dir, exist_ok=True)
                prefix = output_prefix if output_prefix else gene_symbol
                box_fig_path = os.path.join(output_dir, f"{prefix}_boxswarm.{fig_format}")
                boxswarm_fig.savefig(box_fig_path, format=fig_format, dpi=300, bbox_inches='tight')
                self.logger.info(f"Boxplot+swarm saved to {box_fig_path}")
                print(f"Boxplot+swarm saved to {box_fig_path}")
            
            # ---------- FIGURE 3: PIE CHART FOR SAMPLE DISTRIBUTION ----------
            self.logger.info("Creating pie chart for sample distribution...")
            
            # Access sample statistics
            stats_summary = None
            if gene_symbol in self._se_directory_results:
                if "stats_summary" in self._se_directory_results[gene_symbol]:
                    stats_summary = self._se_directory_results[gene_symbol]["stats_summary"]
            
            if stats_summary is None:
                self.logger.warning("Sample statistics not found. Skipping pie chart.")
                pie_fig = None
            else:
                pie_fig = plt.figure(figsize=(10, 8))
                ax = pie_fig.add_subplot(111)
                
                # Extract data for pie chart
                total_samples = stats_summary.get('total_sample_count', 0)
                super_matched_count = stats_summary.get('super_matched', {}).get('count', 0)
                samples_without_se = total_samples - super_matched_count
                
                # Data for pie chart
                labels = [f'Samples with SE ({super_matched_count})', 
                        f'Samples without SE ({samples_without_se})']
                sizes = [super_matched_count, samples_without_se]
                
                # Create pie chart
                wedges, texts, autotexts = ax.pie(
                    sizes, 
                    labels=None,  # We'll add custom legend instead
                    autopct='%1.1f%%',
                    startangle=90,
                    wedgeprops={'edgecolor': 'w', 'linewidth': 1},
                    textprops={'fontsize': 14}
                )
                
                # Set equal aspect ratio to ensure circular pie
                ax.axis('equal')
                
                # Add title
                ax.set_title(f"SE Distribution Across Samples for Gene: {gene_symbol}", fontsize=16)
                
                # Add legend
                ax.legend(wedges, labels, loc='center left', bbox_to_anchor=(0.9, 0.5), 
                        fontsize=12, frameon=True)
                
                # Add some additional sample statistics as text
                sample_stats_text = (
                    f"Total Samples: {total_samples}\n"
                    f"Samples with SE: {super_matched_count} ({stats_summary.get('super_matched', {}).get('percentage', 0):.1f}%)\n"
                    f"Samples without SE: {samples_without_se} ({100 - stats_summary.get('super_matched', {}).get('percentage', 0):.1f}%)"
                )
                
                # Add text box with sample statistics
                props = dict(boxstyle='round', facecolor='white', alpha=0.8)
                ax.text(-0.2, -0.15, sample_stats_text, transform=ax.transAxes, 
                        fontsize=12, va='top', bbox=props)
                
                # Adjust layout
                pie_fig.tight_layout()
                
                # Save pie chart if requested
                if save_fig and pie_fig is not None:
                    os.makedirs(output_dir, exist_ok=True)
                    prefix = output_prefix if output_prefix else gene_symbol
                    pie_fig_path = os.path.join(output_dir, f"{prefix}_sample_pie.{fig_format}")
                    pie_fig.savefig(pie_fig_path, format=fig_format, dpi=300, bbox_inches='tight')
                    self.logger.info(f"Sample distribution pie chart saved to {pie_fig_path}")
                    print(f"Sample distribution pie chart saved to {pie_fig_path}")
            
            # Display quartile information
            q1 = np.percentile(percentiles, 25)
            q2 = np.percentile(percentiles, 50)
            q3 = np.percentile(percentiles, 75)
            iqr = q3 - q1
            
            quartile_text = (
                f"Quartiles for {gene_symbol} Super Enhancer Percentiles:\n"
                f"Q1 (25%): {q1:.2f}\n"
                f"Q2 (50%, Median): {q2:.2f}\n"
                f"Q3 (75%): {q3:.2f}\n"
                f"IQR: {iqr:.2f}"
            )
            print(quartile_text)
            
            # Save data if requested
            if save_data:
                try:
                    # Get the dataframe with relevant data
                    se_df = self.get_directory_super_enhancers(gene_symbol)
                    if se_df is not None and not se_df.empty:
                        # Validate file format
                        if data_format.lower() not in ["csv", "tsv"]:
                            self.logger.warning(f"Invalid file format: {data_format}. Using 'csv' instead.")
                            data_format = "csv"
                        
                        # Set delimiter and extension based on file format
                        delimiter = "," if data_format.lower() == "csv" else "\t"
                        extension = data_format.lower()
                        
                        # Ensure output directory exists
                        os.makedirs(output_dir, exist_ok=True)
                        
                        # Use specified prefix or default to gene symbol
                        prefix = output_prefix if output_prefix else gene_symbol
                        
                        # Define filename without timestamp
                        filename = f"{prefix}_percentile_data.{extension}"
                        output_path = os.path.join(output_dir, filename)
                        
                        # Save the data
                        se_df.to_csv(output_path, sep=delimiter, index=False)
                        
                        self.logger.info(f"Super enhancers percentile data saved to {output_path}")
                        print(f"Super enhancers percentile data saved to {output_path}")
                    else:
                        self.logger.warning(f"No data to save for gene '{gene_symbol}'")
                        print(f"Warning: No data to save for gene '{gene_symbol}'")
                except Exception as e:
                    self.logger.exception(f"Error saving super enhancers percentile data: {e}")
                    print(f"Error saving super enhancers percentile data: {e}")
            
            # Show plots if requested
            if show_plot:
                plt.show()
            else:
                plt.close(histogram_fig)
                plt.close(boxswarm_fig)
                if pie_fig is not None:
                    plt.close(pie_fig)
            
            # Return figures for further customization
            if show_plot:
                return histogram_fig, boxswarm_fig, pie_fig
            else:
                return None
        
        except Exception as e:
            self.logger.exception(f"Error in plot_super_enhancers_overall_distribution: {e}")
            print(f"Error plotting super enhancers overall distribution: {e}")
            return None





    @property
    def FDR(self) -> float:
        """Return the current FDR threshold."""
        return self._FDR  

    @property
    def r(self) -> float:
        """Return the current r threshold."""
        return self._r  

    @property
    def network_FDR(self) -> float:
        """Return the network-specific FDR threshold."""
        return self._network_FDR  

    @property
    def network_r(self) -> float:
        """Return the network-specific r threshold."""
        return self._network_r  

    @property
    def bed_filter(self) -> BedTool:
        """Return the BedTool object after initial filtering."""
        if not self._bed_filter:
            raise ValueError("Bed filter is not set. Please run filter_p2g_file first.")
        return self._bed_filter

    @property
    def bed_before_threshold(self) -> BedTool:
        """Return the BedTool object before network threshold filtering."""
        if not self._bed_before_threshold_network:
            raise ValueError("Bed before network threshold is not set. Please run create_network first.")
        return self._bed_before_threshold_network  

    @property
    def bed_after_threshold(self) -> BedTool:
        """Return the BedTool object after network threshold filtering."""
        if not self._bed_after_threshold_network:
            raise ValueError("Bed after network threshold is not set. Please run create_network first.")
        return self._bed_after_threshold_network  

    @property
    def bed_df_concat_sort(self) -> BedTool:
        """Return the sorted BedTool object after network threshold filtering."""
        if not self._bed_df_concat_sort_network:
            raise ValueError("Sorted BED data for network is not set. Please run create_network first.")
        return self._bed_df_concat_sort_network  

    @property
    def enhancer_gene_edges(self) -> pd.DataFrame:
        if self._df_bed_mergeSE_to_gene_network is None or self._df_bed_mergeSE_to_gene_network.empty:
            raise ValueError("DataFrame for mergeSE to gene is not set or is empty. Please run create_network first.")
        return self._df_bed_mergeSE_to_gene_network[self._df_bed_mergeSE_to_gene_network['attributes'].apply(lambda x: x.get('edge_type') == 'enhancer-gene')]

    @property
    def mergeSE_enhancer_edges(self) -> pd.DataFrame:
        if self._df_bed_mergeSE_to_gene_network is None or self._df_bed_mergeSE_to_gene_network.empty:
            raise ValueError("DataFrame for mergeSE to gene is not set or is empty. Please run create_network first.")
        return self._df_bed_mergeSE_to_gene_network[self._df_bed_mergeSE_to_gene_network['attributes'].apply(lambda x: x.get('edge_type') == 'mergeSE-enhancer')]
    
    @property
    def df_bed_mergeSE_to_gene(self) -> pd.DataFrame:
        if self._df_bed_mergeSE_to_gene_network is None or self._df_bed_mergeSE_to_gene_network.empty:
            raise ValueError("DataFrame for mergeSE to gene is not set or is empty. Please run create_network first.")
        return self._df_bed_mergeSE_to_gene_network

    @property
    def DG(self) -> nx.DiGraph:
        """Return the created network graph."""
        if not self._DG_network:
            raise ValueError("DiGraph is not set. Please run create_network first.")
        return self._DG_network  

    

    @property
    def selected_se_region(self) -> Optional[str]:
        """Return the currently selected SE region."""
        return self._selected_se_region

    @property
    def selected_se_index(self) -> Optional[int]:
        """Return the currently selected SE region's index."""
        return self._selected_se_index


    @property
    def network_threshold(self) -> int:
        if not hasattr(self, '_threshold_network'):
            raise AttributeError("Network threshold is not set. Please run create_network first.")
        return self._threshold_network



    


    @property
    def subgraph_info_filename(self) -> Optional[str]:

        return getattr(self, '_subgraph_info_filename', None)

    @property
    def subgraph_info_df(self) -> Optional[pd.DataFrame]:

        return getattr(self, '_subgraph_info_df', None)
    
    @property
    def selected_subgraph_id(self) -> Optional[int]:
        """Returns the ID of the selected subgraph."""
        return getattr(self, '_selected_subgraph_id', None)

    @property
    def selected_subgraph(self) -> Optional[nx.DiGraph]:
        """Returns the selected subgraph."""
        return getattr(self, '_selected_subgraph', None)




    @property
    def temp_concat_full(self) -> pd.DataFrame:

        if self._temp_concat_full is None:
            self.logger.error("Concatenated SE data is not available. Please run analyze_merge_SE first.")
            raise ValueError("Concatenated SE data is not available. Please run analyze_merge_SE first.")
        return self._temp_concat_full

    @property
    def temp_full_df_edit(self) -> pd.DataFrame:

        if self._temp_full_df_edit is None:
            self.logger.error("Merged SE data is not available. Please run analyze_merge_SE first.")
            raise ValueError("Merged SE data is not available. Please run analyze_merge_SE first.")
        return self._temp_full_df_edit

    @property
    def sorted_df(self) -> pd.DataFrame:

        if self._sorted_df is None:
            self.logger.error("Sorted SE data is not available. Please run analyze_merge_SE first.")
            raise ValueError("Sorted SE data is not available. Please run analyze_merge_SE first.")
        return self._sorted_df



    @property
    def current_sort_key(self) -> Optional[str]:

        return self._current_sort_key



    @property
    def temp_full_df_info(self) -> pd.DataFrame:
        """
        Returns the DataFrame containing the merged SE information.
        
        Raises:
            ValueError: If analyze_merge_SE_info() has not been executed.
        """
        if not hasattr(self, '_temp_full_df_info') or self._temp_full_df_info is None:
            self.logger.error("Temporary full DataFrame info is not available. Please run analyze_merge_SE_info() first.")
            raise ValueError("Temporary full DataFrame info is not available. Please run analyze_merge_SE_info() first.")
        return self._temp_full_df_info


    @property
    def gene_bed_filter(self) -> Optional[BedTool]:
        """
        Returns the gene-specific BED filter created by search_gene_enhancer_links.
        
        Returns:
            Optional[BedTool]: BedTool object containing enhancers for the last searched gene,
                            or None if no filter has been created.
        """
        if hasattr(self, '_gene_bed_filter'):
            return self._gene_bed_filter
        else:
            self.logger.warning("No gene-specific BED filter available. Run search_gene_enhancer_links first.")
            return None

    @property
    def gene_enhancer_results(self) -> dict:
        """
        Returns a dictionary of all gene enhancer search results.
        
        Returns:
            dict: Dictionary with gene symbols as keys and search results as values
        """
        return self._gene_enhancer_results

    @property
    def last_gene_search_result(self) -> Optional[dict]:
        """
        Returns the search result for the last searched gene.
        
        Returns:
            Optional[dict]: Result dictionary for the last searched gene, or None if no gene has been searched
        """
        if self._last_searched_gene is None:
            self.logger.warning("No gene has been searched yet.")
            return None
        return self._gene_enhancer_results.get(self._last_searched_gene)

    @property
    def last_searched_gene(self) -> Optional[str]:
        """
        Returns the symbol of the last searched gene.
        
        Returns:
            Optional[str]: Last searched gene symbol, or None if no gene has been searched
        """
        return self._last_searched_gene

    @property
    def available_gene_results(self) -> List[str]:
        """
        Returns a list of gene symbols for which search results are available.
        
        Returns:
            List[str]: List of gene symbols with available search results
        """
        return list(self._gene_enhancer_results.keys())


    @current_sort_key.setter
    def current_sort_key(self, new_sort_key: str) -> None:

        self.logger.info(f"Setting current_sort_key to '{new_sort_key}'.")
        self._current_sort_key = new_sort_key









    @property
    def se_directory_results(self) -> dict:
        """
        Returns a dictionary of all directory search results.
        
        Returns:
            dict: Dictionary with gene symbols as keys and search results as values
        """
        return self._se_directory_results

    @property
    def last_directory_search_result(self) -> Optional[dict]:
        """
        Returns the search result for the last searched gene directory.
        
        Returns:
            Optional[dict]: Result dictionary for the last searched gene directory, or None if no search has been performed
        """
        if self._last_searched_gene_directory is None:
            self.logger.warning("No gene directory has been searched yet.")
            return None
        return self._se_directory_results.get(self._last_searched_gene_directory)

    @property
    def last_directory_search_gene(self) -> Optional[str]:
        """
        Returns the symbol of the last searched gene directory.
        
        Returns:
            Optional[str]: Last searched gene directory symbol, or None if no search has been performed
        """
        return self._last_searched_gene_directory

    @property
    def available_directory_results(self) -> List[str]:
        """
        Returns a list of gene symbols for which directory search results are available.
        
        Returns:
            List[str]: List of gene symbols with available directory search results
        """
        return list(self._se_directory_results.keys())

    @property
    def directory_super_enhancers(self) -> Optional[pd.DataFrame]:
        """
        Returns the Super Enhancers only (isSuper==1) DataFrame for the last searched gene.
        
        Returns:
            Optional[pd.DataFrame]: DataFrame containing only Super Enhancers, or None if no search has been performed
        """
        if self._last_searched_gene_directory is None:
            self.logger.warning("No gene directory has been searched yet.")
            return None
        
        result = self._se_directory_results.get(self._last_searched_gene_directory)
        if result and "super_only_results" in result:
            return result["super_only_results"]
        return None
    





    @property
    def all_super_enhancers(self) -> Dict[str, pd.DataFrame]:
        """
        Returns a dictionary of all super enhancers (isSuper==1) DataFrames for all searched genes.
        
        Returns:
            Dict[str, pd.DataFrame]: Dictionary with gene symbols as keys and super enhancers DataFrames as values
        """
        result = {}
        for gene_symbol in self._se_directory_results:
            se_data = self._se_directory_results[gene_symbol].get("super_only_results")
            if se_data is not None and not se_data.empty:
                result[gene_symbol] = se_data
        return result

    @property
    def super_enhancers_count(self) -> Dict[str, int]:
        """
        Returns counts of super enhancers for each gene that has been searched.
        
        Returns:
            Dict[str, int]: Dictionary with gene symbols as keys and counts of super enhancers as values
        """
        result = {}
        for gene_symbol in self._se_directory_results:
            se_data = self._se_directory_results[gene_symbol].get("super_only_results")
            if se_data is not None:
                result[gene_symbol] = len(se_data)
            else:
                result[gene_symbol] = 0
        return result

    @property
    def super_enhancers_sample_counts(self) -> Dict[str, Dict[str, int]]:
        """
        Returns counts of super enhancers per sample for each gene.
        
        Returns:
            Dict[str, Dict[str, int]]: Dictionary with gene symbols as keys and dictionaries of sample counts as values
        """
        result = {}
        for gene_symbol in self._se_directory_results:
            se_data = self._se_directory_results[gene_symbol].get("super_only_results")
            if se_data is not None and not se_data.empty and 'sample_name' in se_data.columns:
                sample_counts = se_data['sample_name'].value_counts().to_dict()
                result[gene_symbol] = sample_counts
            else:
                result[gene_symbol] = {}
        return result