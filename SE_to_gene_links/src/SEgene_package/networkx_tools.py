import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import japanize_matplotlib
import pathlib
import os
import collections
import importlib

from typing import List, Tuple

from natsort import natsorted, order_by_index, index_natsorted

from itertools import islice
from pybedtools import BedTool

import matplotlib.ticker as ticker

import networkx as nx

import copy

from IPython.display import SVG, display, display_svg, Image , clear_output

from pygenomeviz import GenomeViz


from SEgene_package.tools import bed_filter_from_p2g_df_all
from SEgene_package.tools import se_list_to_filtered_region_bed
from SEgene_package.tools import filestamp

from SEgene_package.tools import se_list_to_filtered_region_bed_custom







from SEgene_package.tools import bed_filter_from_p2g_df



from SEgene_package.tools import get_wa_wb_df


def csv_path_to_df_RNA(path_RNA):
    df_RNA_data = pd.read_csv(path_RNA, dtype={"start": int, "end": int})
    return df_RNA_data

def rna_data_to_dic_metadata(df_RNA, gtf=True):

    df_RNA_edit = df_RNA.copy()

    
    if gtf:
        df_RNA_edit["start"] = df_RNA_edit["start"] - 1

    df_RNA_metadata = df_RNA_edit.iloc[: ,0:5]
    df_RNA_metadata_set_index = df_RNA_metadata.set_index("symbol")
    dic_RNA_metadata = df_RNA_metadata_set_index.to_dict(orient='index')

    return dic_RNA_metadata


def make_DG_from_df_p2gl(df_input_p2gl,dic_RNA):

    
    

    
    DG = nx.DiGraph()
    
    for _, row in df_input_p2gl.iterrows():
        enhancer_peak = row['PeakID']
        gene_symbol = row['symbol']

        
        
        assert row["chr"] == dic_RNA[gene_symbol]["chr"], f'row_chr:{row["chr"]}, dic_chr:{dic_RNA[gene_symbol]["chr"]}'
        
        
        DG.add_node(enhancer_peak, 
                   node_type = "enhancer",
                   genome_chr = row["chr"],
                   idxP = row["idxP"],
                   start = row["Start"],
                   end = row["End"]
                  )

        
        DG.add_node(gene_symbol, 
                   node_type = "gene",
                   genome_chr = row["chr"],
                   idxE = row["idxE"],
                   strand = row["strand"],
                   TSS = row["g.start"],
                   start = dic_RNA[gene_symbol]["start"],
                   end = dic_RNA[gene_symbol]["end"],
                   gene_search_region_start = row["tx.s"],
                   gene_search_region_end = row["tx.e"]
                  )

        
        DG.add_edge(enhancer_peak, gene_symbol, 
                   edge_type = "enhancer-gene", 
                   genome_chr = row["chr"],
                   weight = row["r"], 
                   FDR = row["FDR"],
                   distance = row["distance"])
        
    return DG



def make_DG_from_mergeSE_to_enhancer(df_input):
    DG=nx.DiGraph()

    for _, row in df_input.iterrows():
        mergeSE_name = row['name']
        enhancer_name = row['blockCount']

        
        DG.add_node(mergeSE_name, 
                   node_type="mergeSE",
                   chr=row["chrom"],
                   start=row["start"],
                   end=row["end"],
                   SE_count=row["score"]
                  )

        
        DG.add_node(enhancer_name, 
                   node_type="enhancer",
                   chr=row["thickStart"],
                   start=row["thickEnd"],
                   end=row["itemRgb"]
                  )

        
        DG.add_edge(mergeSE_name, enhancer_name, 
                   type="mergeSE-enhancer")
        
    return DG



def DG_from_listpath_and_p2gl(list_path, df_p2gl_filtered, threshold=10):
    
    
    
    
    bed_mergeSE = se_list_to_filtered_region_bed(
                                    list_path, 
                                    df_p2gl_filtered,
                                    threshold=10
                                    )[0]

    
    
    bed_p2gl_filter = bed_filter_from_p2g_df(df_p2gl_filtered)
    
    

    df_mergeSE_wa_wb = get_wa_wb_df(bed_mergeSE, bed_p2gl_filter)
    
    
    
    df_mergeSE_wa_wb_SE_to_E = df_mergeSE_wa_wb.iloc[:, 0:10]
    
    

    df_mergeSE_link = df_mergeSE_wa_wb_SE_to_E.drop_duplicates()

    

    return_DG = make_DG_from_mergeSE_to_enhancer(df_mergeSE_link)
    
    return return_DG



def df_create_for_networkx_mergeSE_to_gene(list_path, df_p2gl_filtered,threshold=10):
    



    bed_p2gl_filter_all = bed_filter_from_p2g_df_all(df_p2gl_filtered)


    bed_mergeSE_region_threshould_filtered = se_list_to_filtered_region_bed(list_path, 
                                   df_p2gl_filtered,
                                   threshold=threshold
                                  )[0]

    bed_mergeSE_region_threshould_filtered_wa_wb_all = bed_mergeSE_region_threshould_filtered.intersect(bed_p2gl_filter_all, 
                                                                                                      wa=True, wb=True)


    
    custom_names = ['chrom', 'start', 'end', 'name', 'all_count', 'sample_count', 
                    'chr_p2gl', 'peak_start', 'peak_end', 'peak_ID', 
                    'gene_symbol',
                    'g.start', 'tx.s', 'tx.e',
                    'gene_strand',  
                    'r', 'distance']

    df_bed_mergeSE_region_threshould_filtered_wa_wb_all = bed_mergeSE_region_threshould_filtered_wa_wb_all.to_dataframe(
                                                                        bed_mergeSE_region_threshould_filtered_wa_wb_all,
                                                                        names=custom_names)
    assert df_bed_mergeSE_region_threshould_filtered_wa_wb_all["chrom"].equals(
        df_bed_mergeSE_region_threshould_filtered_wa_wb_all["chr_p2gl"]), "chr error!"

    return df_bed_mergeSE_region_threshould_filtered_wa_wb_all



def create_MultiDiGraph_mergeSE_to_gene(list_path,
                                       df_p2gl_filtered,
                                       dic_rna,
                                       threshold=10):
    
    
    
    
    df_bed_mergeSE_region_threshould_filtered_wa_wb_all = df_create_for_networkx_mergeSE_to_gene(
                                        list_path=list_path,
                                        df_p2gl_filtered=df_p2gl_filtered,
                                        threshold=threshold
                                      )

    


    M_DG=nx.MultiDiGraph()

    for _, row in df_bed_mergeSE_region_threshould_filtered_wa_wb_all.iterrows():
        mergeSE_name = row['name']
        gene_name = row['gene_symbol']
        edge_name = row['peak_ID']

        
        M_DG.add_node(mergeSE_name,                 
                    node_type="mergeSE",
                    genome_chr=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                    SE_count=row["all_count"],
                    sample_count=row["sample_count"]
                  )

        
        M_DG.add_node(gene_name, 
                   node_type = "gene",
                   genome_chr = row["chr_p2gl"],
                   strand = row["gene_strand"],
                   TSS = row["g.start"],
                   start = dic_rna[gene_name]["start"],
                   end = dic_rna[gene_name]["end"],
                   gene_search_region_start = row["tx.s"],
                   gene_search_region_end = row["tx.e"]
                  )

        
        M_DG.add_edge(mergeSE_name, gene_name,
                      edge_type="mergeSE-gene",
                      edge_ID=edge_name,
                      edge_r=row["r"],
                      distance=row["distance"]
                     )
    return M_DG




def create_3_layer_DiGraph_mergeSE_to_gene(list_path,
                                    df_p2gl_filtered,
                                    dic_rna,
                                    threshold=10):
    
    
    
    df_bed_mergeSE_region_threshould_filtered_wa_wb_all = df_create_for_networkx_mergeSE_to_gene(
                                        list_path=list_path,
                                        df_p2gl_filtered=df_p2gl_filtered,
                                        threshold=threshold
                                      )

    
    
    
    
    temp_DG = make_DG_from_df_p2gl(df_input_p2gl=df_p2gl_filtered,
                                   dic_RNA=dic_rna)
    
    
    
    
    for _, row in df_bed_mergeSE_region_threshould_filtered_wa_wb_all.iterrows():
        
        mergeSE_name = row['name']
        enhancer_name = row['peak_ID']
        
        
        
        assert temp_DG.has_node(enhancer_name),f"enhancer {enhancer_name} not found!"
        
        
        temp_DG.add_node(mergeSE_name,                 
                    node_type="mergeSE",
                    genome_chr=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                    SE_count=row["all_count"],
                    sample_count=row["sample_count"]
                  )
        
        
        temp_DG.add_edge(mergeSE_name, enhancer_name,
                      edge_type="mergeSE-enhancer",
                      genome_chr = row["chrom"])
        
    return temp_DG


def create_SE_to_gene_subgraph_fig(DG, 
                                   components_number, 
                                   output_dir_str="./output",
                                   filename="",
                                   save_fig_file=False
                                  ):
    
    imput_all_SE_DG_weak_components = nx.weakly_connected_components(DG)


    list_imput_all_SE_DG_weak_components_sorted = sorted(imput_all_SE_DG_weak_components, key=len, reverse=True)  

    list_set_of_components = list_imput_all_SE_DG_weak_components_sorted
    
    temp_set_of_nodes = list_set_of_components[components_number - 1]
    
    input_DG = DG.subgraph(temp_set_of_nodes)
    
    
    

    








    
    
    start_position = min([temp_node["start"] for _, temp_node in input_DG.nodes(data=True)]) -1000
    end_position = max([temp_node["end"] for _, temp_node in input_DG.nodes(data=True)]) + 1000

    region_size = end_position - start_position 

    set_chr = list(input_DG.nodes(data="genome_chr"))[0][1]
    
    print(f"pos: {set_chr} {start_position}-{end_position}")


    


    gv = GenomeViz(
        fig_track_height=2,
        feature_track_ratio=1.5,
        link_track_ratio=5.0,
        tick_track_ratio=0.5,
        tick_style="bar",
        align_type="center",
    )

    dic_enhancer= {"name": "enhancer", 
                 "size": region_size,
                 "start_pos": start_position}


    dic_gene= {"name": "gene", 
                "size": region_size,
                "start_pos": start_position}


    track_enhancer = gv.add_feature_track(name=dic_enhancer["name"],
                                 size=dic_enhancer["size"],
                                 start_pos=dic_enhancer["start_pos"],
                                 labelmargin=0.03, 
                                 linecolor="grey", 
                                 linewidth=3)


    track_gene = gv.add_feature_track(name=dic_gene["name"],
                                 size=dic_gene["size"],
                                 start_pos=dic_gene["start_pos"],
                                 labelmargin=0.03, 
                                 linecolor="grey", 
                                 linewidth=3)

    for temp_node, temp_data in input_DG.nodes(data=True):

    
        if temp_data["node_type"] == 'mergeSE':
            track_enhancer.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=1, 
                                      plotstyle="box",
                                      label=f"{temp_node}_sample_{temp_data['sample_count']}", 
                                      labelcolor="black", 
                                      facecolor="yellow", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)

        if temp_data["node_type"] == 'enhancer':
            track_enhancer.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=-1, 
                                      plotstyle="box",
                                      label=temp_node, 
                                      labelcolor="black", 
                                      facecolor="lightgreen", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)

        if temp_data["node_type"] == 'gene':

            if temp_data["strand"] == "+":
                temp_strand = 1
            else:
                temp_strand = -1

            track_gene.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=temp_strand, 
                                      label=temp_node, 
                                      labelcolor="black", 
                                      facecolor="pink", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)



    for node_1, node_2, temp_data in input_DG.edges(data=True):
        if temp_data["edge_type"] == "enhancer-gene":
            gv.add_link(track_link1=(dic_enhancer["name"], 
                                     input_DG.nodes(data=True)[node_1]["start"], 
                                     input_DG.nodes(data=True)[node_1]["end"]),
                        track_link2=(dic_gene["name"], 
                                     input_DG.nodes(data=True)[node_2]["start"], 
                                     input_DG.nodes(data=True)[node_2]["end"]),
                                v=temp_data["weight"], 
                                 vmin=0.4, 
                                 vmax=1, 
                                normal_color="skyblue",
                                curve=True)


    fig = gv.plotfig()
    
    

    gv.set_colorbar(fig,
                    ["skyblue"],
        alpha = 0.8,
        vmin = 0,
        vmax = 1,
        bar_height = 0.5,
        bar_width = 0.05,
        bar_left = 1.02,
        bar_bottom = 0,
        bar_label  = "R",
        bar_labelsize = 15,
        tick_labelsize = 10)

    if save_fig_file==True:
        gv.savefig_html(f"{output_dir_str}/{filestamp()}_{filename}c{str(components_number).zfill(2)}_pos_{set_chr}_{start_position}-{end_position}.html")





def check_and_return_chr(input_Graph):

    dic_get_node_attributes = nx.get_node_attributes(input_Graph, "genome_chr")
    set_chr = set(dic_get_node_attributes.values())
    if len(set_chr)==1:
        TF_value = True
        return_chr = set_chr.pop()
    else:
        TF_value = False
        return_chr = None
    
    return TF_value, return_chr




def create_SE_to_gene_graph_fig(input_DG, 
                                   output_dir_str="/home/shinkai/local_for_download/networkx_pybedtools/output",
                                   filename="",
                                   save_fig_file=False
                                  ):
    


    
    
    
    tf_check, set_chr = check_and_return_chr(input_DG)
    
    assert tf_check == True ,"chr_check error!"
    

    








    
    
    start_position = min([temp_node["start"] for _, temp_node in input_DG.nodes(data=True)]) -1000
    end_position = max([temp_node["end"] for _, temp_node in input_DG.nodes(data=True)]) + 1000

    region_size = end_position - start_position 


    
    print(f"pos: {set_chr} {start_position}-{end_position}")


    


    gv = GenomeViz(
        fig_track_height=2,
        feature_track_ratio=1.5,
        link_track_ratio=5.0,
        tick_track_ratio=0.5,
        tick_style="bar",
        align_type="center",
    )

    dic_enhancer= {"name": "enhancer", 
                 "size": region_size,
                 "start_pos": start_position}


    dic_gene= {"name": "gene", 
                "size": region_size,
                "start_pos": start_position}


    track_enhancer = gv.add_feature_track(name=dic_enhancer["name"],
                                 size=dic_enhancer["size"],
                                 start_pos=dic_enhancer["start_pos"],
                                 labelmargin=0.03, 
                                 linecolor="grey", 
                                 linewidth=3)


    track_gene = gv.add_feature_track(name=dic_gene["name"],
                                 size=dic_gene["size"],
                                 start_pos=dic_gene["start_pos"],
                                 labelmargin=0.03, 
                                 linecolor="grey", 
                                 linewidth=3)

    for temp_node, temp_data in input_DG.nodes(data=True):

    
        if temp_data["node_type"] == 'mergeSE':
            track_enhancer.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=1, 
                                      plotstyle="box",
                                      label=f"{temp_node}_sample_{temp_data['sample_count']}", 
                                      labelcolor="black", 
                                      facecolor="yellow", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)

        if temp_data["node_type"] == 'enhancer':
            track_enhancer.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=-1, 
                                      plotstyle="box",
                                      label=temp_node, 
                                      labelcolor="black", 
                                      facecolor="lightgreen", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)

        if temp_data["node_type"] == 'gene':

            if temp_data["strand"] == "+":
                temp_strand = 1
            else:
                temp_strand = -1

            track_gene.add_feature(start= temp_data["start"],                                    
                                      end=temp_data["end"], 
                                      strand=temp_strand, 
                                      label=temp_node, 
                                      labelcolor="black", 
                                      facecolor="pink", 
                                      linewidth=1, 
    
    
    
    
                                      arrow_shaft_ratio=1.0)



    for node_1, node_2, temp_data in input_DG.edges(data=True):
        if temp_data["edge_type"] == "enhancer-gene":
            gv.add_link(track_link1=(dic_enhancer["name"], 
                                     input_DG.nodes(data=True)[node_1]["start"], 
                                     input_DG.nodes(data=True)[node_1]["end"]),
                        track_link2=(dic_gene["name"], 
                                     input_DG.nodes(data=True)[node_2]["start"], 
                                     input_DG.nodes(data=True)[node_2]["end"]),
                                v=temp_data["weight"], 
                                 vmin=0.4, 
                                 vmax=1, 
                                normal_color="skyblue",
                                curve=True)


    fig = gv.plotfig()
    
    

    gv.set_colorbar(fig,
                    ["skyblue"],
        alpha = 0.8,
        vmin = 0,
        vmax = 1,
        bar_height = 0.5,
        bar_width = 0.05,
        bar_left = 1.02,
        bar_bottom = 0,
        bar_label  = "R",
        bar_labelsize = 15,
        tick_labelsize = 10)

    if save_fig_file==True:
        gv.savefig_html(f"{output_dir_str}/{filestamp()}_{filename}_pos_{set_chr}_{start_position}-{end_position}.html")






def filtering_se_by_sample_threshold(input_Graph, 
                                    sample_count_for_threshold):

    input_Graph_for_node_remove = copy.deepcopy(input_Graph)


    

    list_remove_node_by_se_threshold = []

    for node, data in input_Graph_for_node_remove.nodes(data=True):
        if (data['node_type'] == 'mergeSE'
            and data['sample_count'] < sample_count_for_threshold ):
            list_remove_node_by_se_threshold.append(node)

    for node in list_remove_node_by_se_threshold:
        input_Graph_for_node_remove.remove_node(node)
        
        
    return input_Graph_for_node_remove


def return_list_subgraph_WCC(input_Graph):

    
    generator_weak_components = nx.weakly_connected_components(input_Graph)

    
    list_set_of_components = sorted(generator_weak_components, key=len, reverse=True)  

    return list_set_of_components

def return_subgraph_WCC(input_Graph):

    
    generator_weak_components = nx.weakly_connected_components(input_Graph)

    
    list_set_of_components = sorted(generator_weak_components, key=len, reverse=True)  

    
    subgraphs = []
    for nodes in list_set_of_components:
        subgraph = input_Graph.subgraph(nodes)
        subgraphs.append(subgraph)
    
    return subgraphs




def temp_call_graphviz(G):
    
    A = nx.nx_agraph.to_agraph(G)


    
    for node in A.nodes():
        if node.attr['node_type'] == 'enhancer':
            node.attr['color'] = 'blue'
            node.attr['shape'] = 'box'  
        elif node.attr['node_type'] == 'gene':
            node.attr['color'] = 'red'
            
            node.attr['shape'] = 'ellipse'

    
    A.write('test.dot')

    
    A.draw('test.png', prog='fdp', format='png')

    
    return Image(filename='test.png')


def get_all_merge_network_fig(mergeSE_to_gene_DG, path_output_svg):

    
    
    agraph_mergeSE_to_gene = nx.nx_agraph.to_agraph(mergeSE_to_gene_DG)
    agraph_mergeSE_to_gene.draw(path_output_svg, prog='fdp', format='svg')





    
def get_node_count_in_subgraphs_fig(input_DG,sort_type="all"):

    
    DG=input_DG

    if sort_type=="all":    

        
        components = nx.weakly_connected_components(DG)

        
        node_counts = [len(component) for component in components]

        
        node_counts.sort(reverse=True)

        
        plt.plot(node_counts)
        plt.title('Node counts in subgraphs')
        plt.xlabel('Subgraph index (sorted by total node count)')
        plt.ylabel('Count')
        plt.show()
        






def df_create_for_networkx_mergeSE_to_gene_custom(
    list_path: List[str], 
    df_p2gl_filtered: pd.DataFrame, 
    threshold: int = 10
) -> Tuple[pd.DataFrame, BedTool, BedTool, BedTool]:

    
    bed_after_threshold, bed_before_threshold, bed_df_concat_sort = se_list_to_filtered_region_bed_custom(
        list_se=list_path, 
        df_p2gl=df_p2gl_filtered, 
        threshold=threshold
    )
    
    
    bed_p2gl_filter_all = bed_filter_from_p2g_df_all(df_p2gl_filtered)
    
    
    bed_mergeSE_region_threshould_filtered_wa_wb_all = bed_after_threshold.intersect(
        bed_p2gl_filter_all, 
        wa=True, 
        wb=True
    )
    
    
    custom_names = [
        'chrom', 'start', 'end', 'name', 'all_count', 'sample_count', 
        'chr_p2gl', 'peak_start', 'peak_end', 'peak_ID', 
        'gene_symbol',
        'g_start', 'tx_s', 'tx_e',
        'gene_strand',  
        'r', 'distance'
    ]
    
    
    df_bed_mergeSE_region_threshould_filtered_wa_wb_all = bed_mergeSE_region_threshould_filtered_wa_wb_all.to_dataframe(
        bed_mergeSE_region_threshould_filtered_wa_wb_all,
        names=custom_names
    )
    
    
    assert df_bed_mergeSE_region_threshould_filtered_wa_wb_all["chrom"].equals(
        df_bed_mergeSE_region_threshould_filtered_wa_wb_all["chr_p2gl"]
    ), "chrom chr_p2gl not match"
    
    return df_bed_mergeSE_region_threshould_filtered_wa_wb_all, bed_after_threshold, bed_before_threshold, bed_df_concat_sort




def create_3_layer_DiGraph_mergeSE_to_gene_custom(
    list_path: List[str],
    df_p2gl_filtered: pd.DataFrame,
    dic_rna: dict,
    threshold: int = 10
) -> Tuple[nx.DiGraph, BedTool, BedTool, BedTool]:

    
    df_mergeSE_to_gene, bed_after_threshold, bed_before_threshold, bed_df_concat_sort = df_create_for_networkx_mergeSE_to_gene_custom(
        list_path=list_path,
        df_p2gl_filtered=df_p2gl_filtered,
        threshold=threshold
    )
    
    
    
    
    temp_DG = make_DG_from_df_p2gl(
        df_input_p2gl=df_p2gl_filtered,
        dic_RNA=dic_rna
    )
    
    
    for _, row in df_mergeSE_to_gene.iterrows():
        mergeSE_name = row['name']
        enhancer_name = row['peak_ID']
        
        
        if not temp_DG.has_node(enhancer_name):
            raise AssertionError(f"enhancer {enhancer_name} not found!")
        
        
        temp_DG.add_node(
            mergeSE_name,
            node_type="mergeSE",
            genome_chr=row["chrom"],
            start=row["start"],
            end=row["end"],
            SE_count=row["all_count"],
            sample_count=row["sample_count"]
        )
        
        
        temp_DG.add_edge(
            mergeSE_name,
            enhancer_name,
            edge_type="mergeSE-enhancer",
            genome_chr=row["chrom"]
        )
        
    return temp_DG, bed_after_threshold, bed_before_threshold, bed_df_concat_sort
