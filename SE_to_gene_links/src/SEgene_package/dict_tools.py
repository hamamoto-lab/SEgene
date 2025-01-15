
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


import matplotlib.ticker as ticker


def make_dict_se_region_id_to_data(file_path):
    
    table_new_index=['REGION_ID',
     'CHROM',
     'START',
     'STOP',
     'NUM_LOCI',
     'CONSTITUENT_SIZE',
     'SIGNAL',
     'CONTROL_SIGNAL',
     'enhancerRank',
     'isSuper']
    
    df_temp_table_all=pd.read_table(file_path,comment="#").set_axis(table_new_index, axis='columns')
    
    df_temp_table_all_set_region_id_index=df_temp_table_all.set_index('REGION_ID')
    
    d_index = df_temp_table_all_set_region_id_index.to_dict(orient='index',into=collections.OrderedDict)

    
    
    if len(df_temp_table_all)!=len(d_index):
        raise Exception('dict size Error!')
    
    
    return(d_index)


def make_dict_se_rank_to_data(file_path):
    
    table_new_index=['REGION_ID',
     'CHROM',
     'START',
     'STOP',
     'NUM_LOCI',
     'CONSTITUENT_SIZE',
     'SIGNAL',
     'CONTROL_SIGNAL',
     'enhancerRank',
     'isSuper']
    
    df_temp_table_all=pd.read_table(file_path,comment="#").set_axis(table_new_index, axis='columns')
    
    df_temp_table_all_set_region_id_index=df_temp_table_all.set_index('enhancerRank')
    
    d_index = df_temp_table_all_set_region_id_index.to_dict(orient='index',into=collections.OrderedDict)

    if len(df_temp_table_all)!=len(d_index):
        raise Exception('dict size Error!')

    return(d_index)


def get_simple_sample_id_from_se_table(pathlib_se_file):
    with pathlib_se_file.open() as f:
        sample_name=f.readline()
        output_name=sample_name.lstrip("#").rstrip()
        return output_name



def sort_df_for_bed_by_chr(input_df):


    input_for_sort_df=input_df

    sort_chr_list=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 
                       'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                   'chrX','chrY']

    cols = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    output_df = pd.DataFrame(index=[], columns=cols)

    for temp_chr in sort_chr_list:

        temp_concat_df = input_for_sort_df[input_for_sort_df["chrom"].isin([temp_chr])]

        temp_concat_df_sort_by_chr_start = temp_concat_df.sort_values('chromStart')

        output_df = pd.concat([output_df, temp_concat_df_sort_by_chr_start], axis=0)
            
    return output_df


def return_df_for_bed_from_se_count_list(input_series_se_count,input_path_enhancer_table):


    input_se_list=input_series_se_count.index.tolist()

    input_se_dic=make_dict_se_region_id_to_data(input_path_enhancer_table)


    bed_columns=["chrom",
              "chromStart",
              "chromEnd",
              "name",
              "score",
              "strand"]

    df_for_bed_6 = pd.DataFrame(index=[], columns=bed_columns)


    for se_name in input_se_list:

        add_bed_data=[
        input_se_dic[se_name]["CHROM"],
        input_se_dic[se_name]["START"],
        input_se_dic[se_name]["STOP"],
        se_name,
        input_series_se_count[se_name],
        "."
        ]

        series_add_bed = pd.Series(add_bed_data, index=df_for_bed_6.columns, name=se_name)

        df_for_bed_6=df_for_bed_6.append(series_add_bed)


    return_df=sort_df_for_bed_by_chr(df_for_bed_6)




def return_df_for_bed_from_se_count_list_for_se2gene_file_output_only_rose_se(input_series_se_count,input_path_enhancer_table):



    input_se_list=input_series_se_count.index.tolist()

    input_se_dic=make_dict_se_region_id_to_data(input_path_enhancer_table)


    bed_columns=["chrom",
              "chromStart",
              "chromEnd",
              "name",
              "score",
              "strand"]

    df_for_bed_6 = pd.DataFrame(index=[], columns=bed_columns)


    for se_name in input_se_list:

        add_bed_data=[
        input_se_dic[se_name]["CHROM"],
        input_se_dic[se_name]["START"],
        input_se_dic[se_name]["STOP"],
        se_name,
        input_series_se_count[se_name],
        "."
        ]

        series_add_bed = pd.Series(add_bed_data, index=df_for_bed_6.columns, name=se_name)


        
        df_for_bed_6 = pd.concat([df_for_bed_6, series_add_bed.to_frame().T])



    return df_for_bed_6