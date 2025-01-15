
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import japanize_matplotlib
import pathlib
import os
import collections
import datetime
from natsort import natsorted, order_by_index, index_natsorted

from itertools import islice
from pybedtools import BedTool

from functools import partial
from typing import Optional
from typing import List, Tuple


import matplotlib.ticker as ticker

import random
import string


def sort_bed(bed_object):
    bed_df=bed_object.to_dataframe()
    bed_df.sort_values("end", inplace=True)
    bed_df.sort_values("start", inplace=True)
    index_order = index_natsorted(bed_df.chrom)
    order_by_index_order = order_by_index(bed_df.index, index_order)
    return_df = bed_df.reindex(order_by_index_order)
    return_bed = BedTool.from_dataframe(return_df)
    return return_bed





from SEgene_package.dict_tools import return_df_for_bed_from_se_count_list_for_se2gene_file_output_only_rose_se




def return_bed_obj_tupple_from_rose_all_table(path_all_enhancer_table):

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
    

    df_temp_table_all=pd.read_table(path_all_enhancer_table,comment="#").set_axis(table_new_index, axis='columns')
    

    df_temp_table_all["FIXED_SIGNAL"]=df_temp_table_all["SIGNAL"]-df_temp_table_all["CONTROL_SIGNAL"]

    df_temp_table_all_assign=df_temp_table_all.assign(STRAND=".")
    

    
    df_temp_table_super=df_temp_table_all_assign.query("isSuper == 1")
    df_temp_table_typical=df_temp_table_all_assign.query("isSuper == 0")

    

    df_temp_table_all_to_bed=df_temp_table_all_assign[["CHROM","START","STOP","REGION_ID","FIXED_SIGNAL","STRAND"]]
    
    df_temp_table_super_to_bed=df_temp_table_super[["CHROM","START","STOP","REGION_ID","FIXED_SIGNAL","STRAND"]]
    
    df_temp_table_typical_to_bed=df_temp_table_typical[["CHROM","START","STOP","REGION_ID","FIXED_SIGNAL","STRAND"]]

    

    
    bed_obj_all=BedTool.from_dataframe(df_temp_table_all_to_bed)
    bed_obj_super=BedTool.from_dataframe(df_temp_table_super_to_bed)
    bed_obj_typical=BedTool.from_dataframe(df_temp_table_typical_to_bed)
    
    
    return bed_obj_all,bed_obj_super,bed_obj_typical



def bed_filter_from_p2g_output(path_p2g_output,
                                fdr_filter=0.01,
                                r_filter=0.5):

    df_temp=pd.read_table(path_p2g_output) 
    df_temp_fdr=df_temp.query("FDR < " + str(fdr_filter))
    df_temp_fdr_r=df_temp_fdr.query("r > " + str(r_filter))
    
    filter_query_list=["chr","Start", "End", "PeakID", "symbol", "strand"]
    df_data_for_filter=df_temp_fdr_r[filter_query_list]
    return BedTool.from_dataframe(df_data_for_filter)

def bed_default_from_p2g_output(path_p2g_output):


    df_temp=pd.read_table(path_p2g_output)

    filter_query_list=["chr","Start", "End", "PeakID", "symbol", "strand"]
    df_data_for_filter=df_temp[filter_query_list]
    return BedTool.from_dataframe(df_data_for_filter)

def make_gene_p2g_filter(p2g_filter_bed,key_gene_list):
    df_from_bed=p2g_filter_bed.to_dataframe()
    return_df=df_from_bed[df_from_bed["score"].isin(key_gene_list)]
    return_bed=BedTool.from_dataframe(return_df)
    return return_bed

def bed_filter_from_p2g_df(p2g_df):

    df_temp=p2g_df.copy()

    
    filter_query_list=["chr","Start", "End", "PeakID", "symbol", "strand"]
    df_data_for_filter=df_temp[filter_query_list]
    return_bed = BedTool.from_dataframe(df_data_for_filter)
    return sort_bed(return_bed)


def get_index_and_value(sr_se_count):
    dim2list_sr_se_count=sr_se_count.reset_index().values.tolist()
    temp_se_count_index, temp_se_count_value = zip(*dim2list_sr_se_count)
    return temp_se_count_index,temp_se_count_value



def se_count_from_bedtools_wa_wb(bedtools_wa_wb):

    sr_se_count_bedtools_wa_wb=(bedtools_wa_wb["name"].
                                            value_counts(ascending=False))
    return sr_se_count_bedtools_wa_wb

def peak_merge_se_count_from_bedtools_wa_wb(bedtools_wa_wb):

    peak_merge_wa_wb=bedtools_wa_wb.drop_duplicates(subset="blockCount",keep="last")

    return peak_merge_wa_wb

def gene_merge_se_count_from_bedtools_wa_wb(bedtools_wa_wb):
    gene_merge_wa_wb=bedtools_wa_wb.drop_duplicates(subset=["name","blockSizes"],keep="last")

    return gene_merge_wa_wb



def gene_count_from_wa_wb(df_wa_wb):
    list_gene_count=df_wa_wb["blockSizes"].value_counts(ascending=False)
    return list_gene_count


def test_test(a):
    print(a)




def return_p2g_df(p2g_file):
    df_temp=pd.read_table(p2g_file)
    return df_temp

def filter_p2g_df(p2g_df,
                    fdr_filter=0.01,
                    r_filter=0.5):
    df_temp=p2g_df
    df_temp_fdr=df_temp.query("FDR <= " + str(fdr_filter))
    df_temp_fdr_r=df_temp_fdr.query("r >= " + str(r_filter))

    return df_temp_fdr_r



def print_gene_from_p2g_df(
    p2g_df: pd.DataFrame, 
    output_file: Optional[str] = None, 
    show_item_number: bool = True, 
    show_genes: bool = True, 
    show_gene_count: bool = False, 
    sort_genes: bool = False, 
    save_and_print: bool = False
) -> None:

    gene_counts: pd.Series = p2g_df["symbol"].value_counts()
    gene_array: list[str] = (
        gene_counts.index.tolist() if sort_genes else p2g_df["symbol"].unique().tolist()
    )
    
    output_lines: list[str] = []


    if show_item_number:
        output_lines.append(f"The following {len(gene_array)}  genes were output:")
        output_lines.append("") 


    if show_genes:
        for gene in gene_array:
            line: str = f"{gene} {gene_counts[gene]}" if show_gene_count else gene
            output_lines.append(line)


    if output_file:
        with open(output_file, "w", encoding="utf-8") as f:
            f.write("\n".join(output_lines))
        if save_and_print:
            print("\n".join(output_lines))
        print(f"Output written to {output_file}.")
    else:
        print("\n".join(output_lines))


def count_genes_from_bed(
    bed: BedTool, 
    output_file: Optional[str] = None, 
    show_item_number: bool = True, 
    show_genes: bool = True, 
    show_gene_count: bool = False, 
    sort_genes: bool = False, 
    save_and_print: bool = False
) -> None:


    bed_df: pd.DataFrame = bed.to_dataframe(names=["chrom", "start", "end", "name", "symbol", "strand"])

    print_gene_from_p2g_df(
        p2g_df=bed_df, 
        output_file=output_file, 
        show_item_number=show_item_number, 
        show_genes=show_genes, 
        show_gene_count=show_gene_count, 
        sort_genes=sort_genes, 
        save_and_print=save_and_print
    )


    

def return_super_bed_from_enhancer_table(path_enhancer_table):
    tupple_bed_rose_table=return_bed_obj_tupple_from_rose_all_table(path_enhancer_table) 

    _,bed_table_super,_=tupple_bed_rose_table

    return bed_table_super

def return_bed_form_df(df):
    x = BedTool.from_dataframe(df)
    return x

def return_df_from_bed(pybedtools_bed):
    y = pybedtools_bed.to_dataframe()
    return y

def get_wa_wb_bed(bed_a,bed_b):
    return_bed=bed_a.intersect(bed_b, wa=True, wb=True)
    return return_bed

def get_wa_wb_df(bed_a,bed_b):
    return return_df_from_bed(get_wa_wb_bed(bed_a,bed_b))



def print_gene_list(gene_list):
    print(len(gene_list))
    print()
    for item in gene_list:
        print(item)

def p2g_return_only_in_first(a,b):
    a_list=a["symbol"].unique().tolist()
    b_list=b["symbol"].unique().tolist()
    only_in_a= [x for x in a_list if x not in b_list]
    return only_in_a

def p2g_return_both(a,b):
    a_list=a["symbol"].unique().tolist()
    b_list=b["symbol"].unique().tolist()
    in_both=[x for x in a_list if x in b_list]
    return in_both

def get_same_df(df1, df2, column1, column2):
    common_values = df1.loc[:, [column1, column2]].merge(df2.loc[:, [column1, column2]], on=[column1, column2], how='inner')
    mask1 = df1.apply(lambda x: x[column1] in common_values[column1].values and x[column2] in common_values[column2].values, axis=1)
    df1_filtered = df1[mask1]

    return df1_filtered

def get_dif_df(df1, df2, column1, column2):
    common_values = df1.loc[:, [column1, column2]].merge(df2.loc[:, [column1, column2]], on=[column1, column2], how='inner')
    mask1 = df1.apply(lambda x: x[column1] not in common_values[column1].values or x[column2] not in common_values[column2].values, axis=1)
    df1_filtered = df1[mask1]

    return df1_filtered



def se2gene_file_output_only_rose_se(bed_obj_bed_filter, path_enhancer_table):
    

    
    
    tupple_bed_rose_table=return_bed_obj_tupple_from_rose_all_table(path_enhancer_table)    
    
    bed_table_all,bed_table_super,_=tupple_bed_rose_table
    
    
 
    result_bed_intersect=bed_table_super.intersect(bed_obj_bed_filter, wa=True, wb=True)
    
    df_bed_table_super_coverd_wa_wb=result_bed_intersect.to_dataframe()



    sr_se_count_df_bed_table_super_coverd_wa_wb=se_count_from_bedtools_wa_wb(df_bed_table_super_coverd_wa_wb)

    super_count_df_for_bed=return_df_for_bed_from_se_count_list_for_se2gene_file_output_only_rose_se(
                                                            sr_se_count_df_bed_table_super_coverd_wa_wb,
                                                            path_enhancer_table
    )
    
    
    
    

    gene_merge_df=gene_merge_se_count_from_bedtools_wa_wb(df_bed_table_super_coverd_wa_wb)

    sr_gene_merge_count=se_count_from_bedtools_wa_wb(gene_merge_df)

    super_count_df_for_bed_gene_count=return_df_for_bed_from_se_count_list_for_se2gene_file_output_only_rose_se(
                                                            sr_gene_merge_count,
                                                            path_enhancer_table
    )
    


    peak_merge_df=peak_merge_se_count_from_bedtools_wa_wb(df_bed_table_super_coverd_wa_wb)

    sr_peak_merge_count=se_count_from_bedtools_wa_wb(peak_merge_df)


    super_count_df_for_bed_peak_count=return_df_for_bed_from_se_count_list_for_se2gene_file_output_only_rose_se(
                                                            sr_peak_merge_count,
                                                            path_enhancer_table
    )
    
    
    super_bed_file_set_dic={"link_all":super_count_df_for_bed,
                      "link_gene":super_count_df_for_bed_gene_count,
                      "link_peak":super_count_df_for_bed_peak_count
                     }

    
    

    super_gene_count_list_link_all=gene_count_from_wa_wb(df_bed_table_super_coverd_wa_wb)

    super_gene_count_list_link_gene=gene_count_from_wa_wb(gene_merge_df)

    super_gene_count_list_link_peak=gene_count_from_wa_wb(peak_merge_df)

    super_gene_count_list_dic={"link_all":super_gene_count_list_link_all,
                      "link_gene":super_gene_count_list_link_gene,
                      "link_peak":super_gene_count_list_link_peak
                     }

    return super_bed_file_set_dic,super_gene_count_list_dic





def return_all_se_concat(input_se_file_list, filter_bed):


    cols = ["chrom","chromStart", "chromEnd", "name", "score", "strand"]
    temp_concat_for_bed_df = pd.DataFrame(index=[], columns=cols)
    
    for path_now_processing_sample in input_se_file_list:

        temp_return_bed_dic,temp_return_gene_count_dic=se2gene_file_output_only_rose_se(filter_bed,
                                                                                       path_now_processing_sample)
        temp_concat_for_bed_df = pd.concat([temp_concat_for_bed_df, temp_return_bed_dic["link_gene"]], axis=0)
        
    result_bed_df = sort_df_for_bed_by_chr(temp_concat_for_bed_df)
    
    return result_bed_df


def return_merge_se_count_df(concat_df):
    bed_file_result_bed_df=BedTool.from_dataframe(concat_df)

    test_merge_bed_count=bed_file_result_bed_df.merge(c=1,o="count")
    
    df_test_merge_bed_count=test_merge_bed_count.to_dataframe()
    
    series_se_data = df_test_merge_bed_count["chrom"] \
    + "_" \
    + df_test_merge_bed_count["start"].astype(str) \
    + "_" \
    + df_test_merge_bed_count["end"].astype(str)
    
    df_se_merge_count=pd.DataFrame({"se_data": series_se_data,
             "se_count" : df_test_merge_bed_count["name"] })
    
    return df_se_merge_count

def make_sort_se_merge_count_graph(se_count_df,graph_max_number):
    x_list=se_count_df["se_data"].to_list()
    y_list=se_count_df["se_count"].to_list()
    
    max_number=graph_max_number
    title_str="merge-se number count"

    x_list_graph=x_list[0:max_number]
    y_list_graph=y_list[0:max_number]
    plt.bar(x_list_graph,y_list_graph)
    plt.tight_layout()    

    plt.xticks(rotation=90)

    plt.gca().get_yaxis().set_major_locator(ticker.MaxNLocator(integer=True))

    plt.xlabel('SE')
    plt.ylabel('link number')
    plt.title(title_str)   
    plt.show()    
    plt.clf()
    plt.close()


def get_merge_se_name_series(input_bed_df):
    series_se_data = input_bed_df["chrom"] \
    + "_" \
    + input_bed_df["start"].astype(str) \
    + "_" \
    + input_bed_df["end"].astype(str)
    return series_se_data

def return_merge_se_bed_bed6(merge_se_bed3):
    df_merge_se_bed3=merge_se_bed3.to_dataframe()
    new_df_name_series=get_merge_se_name_series(df_merge_se_bed3)
    new_df_rename = df_merge_se_bed3.rename(columns={"name": "score"})
    new_df_rename.insert(3, "name", new_df_name_series)
    new_df_rename["strand"]="."
    bed_new_df_rename=BedTool.from_dataframe(new_df_rename)
    return bed_new_df_rename

def get_online_bed(chr_txt, start, end):
    temp_str=str(chr_txt)+" "+str(start)+" "+str(end)
    temp_online_bed = BedTool(temp_str, from_string=True)
    return temp_online_bed

def return_bed3(input_bed):
    temp_bed_df=input_bed.to_dataframe()
    temp_bed3_df=temp_bed_df[["chrom", "start", "end"]]
    return BedTool.from_dataframe(temp_bed3_df)

def return_all_count_gene_df(super_enhancer_list, bed_filter):

    series_list=[]


    for item in super_enhancer_list:

        temp_se_file=item


        file_name=item.name

        sample_name=file_name.split("-")[0]

        temp_se_bed=return_super_bed_from_enhancer_table(temp_se_file)

        temp_wa_wb=get_wa_wb_df(temp_se_bed, bed_filter)


        gene_merge_df=gene_merge_se_count_from_bedtools_wa_wb(temp_wa_wb)

        gene_merged_gene_count_series=gene_count_from_wa_wb(gene_merge_df)

        gene_merged_gene_count_series.name=sample_name

        series_list= series_list + [gene_merged_gene_count_series]


    return_df =pd.DataFrame(series_list)
    
    return_df = return_df.fillna(0)

    return_df_T = return_df.T

    return return_df_T



def se_file_to_se_bed(se_file_path):
    tupple_bed_rose_table=return_bed_obj_tupple_from_rose_all_table(se_file_path)    
    
    _,bed_table_super,_=tupple_bed_rose_table
    
    return bed_table_super





def return_filtered_bed_df_from_se_file_and_filter_bed(input_se_file_path, bed_filter):
    temp_se_bed = se_file_to_se_bed(input_se_file_path)
    temp_se_bed_sorted = sort_bed(temp_se_bed)

    sorted_bed_filter=sort_bed(bed_filter)
    temp_se_bed_sorted_wa = temp_se_bed_sorted.intersect(sorted_bed_filter, wa=True)
    temp_se_bed_sorted_wa_df = temp_se_bed_sorted_wa.to_dataframe()


    temp_se_bed_sorted_wa_df_drop = temp_se_bed_sorted_wa_df.drop_duplicates()
    return temp_se_bed_sorted_wa_df_drop


def return_concat_df(input_df_list):
    
    template_columns = input_df_list[0].columns
    if all(temp_df.columns.equals(template_columns) for temp_df in input_df_list):
        return pd.concat(input_df_list, ignore_index=True)
    else :
        raise ValueError("Columns are not uniform.")
    

def p2gl_path_to_filter_df(path_p2g_output, 
                           fdr_filter=0.05, 
                           r_filter=0.5):
    df_p2gl = pd.read_table(path_p2g_output,
                            dtype={"g.start": int, 
                                   'tx.s': int, 
                                   'tx.e': int,
                                    'idxE': int,
                                   'Start': int,
                                    'End': int,
                                    'idxP': int
                                  })
    df_distance = df_p2gl[["Start","End"]].sub(df_p2gl["g.start"], axis=0).abs()
    
    series_distance = df_distance.min(axis=1)
    
    df_p2gl["distance"] = series_distance
    
    df_p2gl_FDR = df_p2gl.query(f'FDR <= {fdr_filter}')

    df_p2gl_FDR_r = df_p2gl_FDR.query(f'r >= {r_filter}')
    return df_p2gl_FDR_r

def get_SE_sample_count(temp_merge_SE_name, mergeSE_wa_wb_sample_add):

    df_temp_count = mergeSE_wa_wb_sample_add.query(f"name == '{temp_merge_SE_name}'")
    
    count_number = df_temp_count['sample_name'].nunique()
    
    return count_number

def make_sample_dict(list_se_path):
    sample_dict={}
    for temp_file_path in list_se_path:
        with open(temp_file_path, 'r') as file:
            header_one_line = file.readline().strip()
            get_sample_name = header_one_line.lstrip('#').split(" ")[0]

        temp_df = pd.read_table(temp_file_path, header=5, usecols=[0])

        list_keys = temp_df.iloc[:,0].tolist()
        temp_dict =dict.fromkeys(list_keys, get_sample_name)
        sample_dict = sample_dict | temp_dict
        
    return sample_dict


def add_sample_count_to_bed6(bed_input, bed_all_sample_concat_SE, temp_se_path):


    temp_allSE = bed_input.intersect(bed_all_sample_concat_SE, wa=True, wb=True).to_dataframe()

    df_bed_input = bed_input.to_dataframe()



    temp_dict_sample_name = make_sample_dict(temp_se_path)

    temp_allSE["sample_name"] = temp_allSE["blockCount"].map(temp_dict_sample_name)

    df_bed_input["sample_count"] = df_bed_input["name"].apply(get_SE_sample_count, args=(temp_allSE,))



    df_bed_input_edit = df_bed_input.drop('strand', axis=1)

    df_bed_input_return = df_bed_input_edit.rename(columns={'sample_count': 'strand'})
    
    bed_return = BedTool.from_dataframe(df_bed_input_return)
    
    return bed_return




def se_list_to_filtered_region_bed(list_se, df_p2gl, threshold=0):
    
    bed_p2gl_filter = bed_filter_from_p2g_df(df_p2gl)
    

    list_processed_df = [
        return_filtered_bed_df_from_se_file_and_filter_bed(path_temp, bed_p2gl_filter)
        for path_temp in list_se]
    
    df_concat = return_concat_df(list_processed_df)

    bed_df_concat_sort = sort_bed(BedTool.from_dataframe(df_concat))


    bed_count = bed_df_concat_sort.merge(c=2,o="count")


    bed_bed6_count= return_merge_se_bed_bed6(bed_count)

    bed_bed6_count_fix = add_sample_count_to_bed6(bed_bed6_count,
                                                 bed_df_concat_sort,
                                                 list_se)

    bed_return = bed_bed6_count_fix


    if threshold >= 1:

        df_threshold_filtered = bed_bed6_count_fix.to_dataframe().query(f"strand >= {threshold}")
        
        bed_return = BedTool.from_dataframe(df_threshold_filtered)
        
    
    return bed_return,bed_df_concat_sort



def bed_filter_from_p2g_df_plus(p2g_df):

    df_temp=p2g_df.copy()

    
    filter_query_list=["chr","Start", "End", "PeakID", "symbol", "strand", "r", "distance"]
    df_data_for_filter=df_temp[filter_query_list]
    return_bed = BedTool.from_dataframe(df_data_for_filter)
    return sort_bed(return_bed)



def bed_filter_from_p2g_df_all(p2g_df):

    df_temp=p2g_df.copy()
    
    filter_query_list=["chr","Start", "End", "PeakID", "symbol", 
                       "g.start","tx.s","tx.e",
                       "strand", "r", "distance"]
    df_data_for_filter=df_temp[filter_query_list]
    
    return_bed = BedTool.from_dataframe(df_data_for_filter)
    return sort_bed(return_bed)



def bed_filter_from_p2g_df_all_check(p2g_df,dic_rna):

    df_temp=p2g_df.copy()

    
    df_temp["gene_start"] =df_temp["symbol"].apply(lambda temp_gene_symbol: dic_rna[temp_gene_symbol]['start'])
    df_temp["gene_end"] =df_temp["symbol"].apply(lambda temp_gene_symbol: dic_rna[temp_gene_symbol]['end'])   
    
    filter_query_list=["chr","Start", "End", "PeakID", "symbol", 
                       "g.start","gene_start","gene_end",
                       "strand", "r", "distance"]
    df_data_for_filter=df_temp[filter_query_list]
    
    return_bed = BedTool.from_dataframe(df_data_for_filter)
    return sort_bed(return_bed)


def filestamp(length=8,sep=None):
    assert 1 <= length <= 20,"timestamp length error" 
      
    stamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')[:length]
    
    if sep is not None:
        assert 4 <=length
        stamp=stamp[:4]+sep+stamp[4:]

    return stamp



def sort_se_merge_count_df(merge_se_count_df):
    sort_by_count_df_se_merge_count=merge_se_count_df.sort_values("se_count",  ascending=False)
    return sort_by_count_df_se_merge_count



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



def get_se_concat_df_and_bed_filter_c4_collapse(input_concat,input_bed_filter,filter_num=0):
    bed_file_result_bed_df=BedTool.from_dataframe(input_concat)
    result=bed_file_result_bed_df.merge(c=4,o="count,collapse").intersect(input_bed_filter, wa=True,wb=True)
    test_result_df=result.to_dataframe()
    test_result_df_filtered=test_result_df.query(f'name >= {filter_num}')
    
    return test_result_df_filtered

def gene_call_from_concat_df(input_concat,input_bed_filter,filter_num=0):
    test_result_df_filtered=get_se_concat_df_and_bed_filter_c4_collapse(input_concat,input_bed_filter,filter_num)
    values = test_result_df_filtered['blockCount'].tolist()
    unique_values = list(set(values))
    for value in unique_values:
        print(value)



def c4_collapse_df_to_count_df(input_c4_df):


    input_c4_df_for_count = input_c4_df.filter(items=["name","score","blockCount"],  axis='columns')


    input_c4_df_for_count["peak_count"]=input_c4_df_for_count.groupby(["score", "blockCount"]).transform("count")


    input_c4_df_for_count_merge=input_c4_df_for_count.drop_duplicates()

    return input_c4_df_for_count_merge




def region_str_to_bed(input_str):
    chr_str,start_str,end_str=input_str.split("_")
    start_int=int(start_str)
    end_int=int(end_str)
    return_bed=get_online_bed(chr_str, start_int, end_int)
    return return_bed


def get_region_p2gl_info(input_chr_str,
                         input_start_pos,
                         input_end_pos,
                         temp_input_p2gl_bed):
    
    bed_for_check=get_online_bed(input_chr_str, input_start_pos, input_end_pos)

    temp_wa_wb=get_wa_wb_df(bed_for_check, temp_input_p2gl_bed)
    
    temp_wa_wb_value_count=temp_wa_wb["thickEnd"].value_counts()
    
    temp_wa_wb_value_count_gene_list=temp_wa_wb_value_count.index.tolist()
    
    temp_wa_wb_gene_count=len(temp_wa_wb_value_count_gene_list)
    
    temp_wa_wb_peak_value_count_series=temp_wa_wb["thickStart"].value_counts()
    
    temp_wa_wb_peak_list=temp_wa_wb_peak_value_count_series.index.tolist()

    temp_wa_wb_peak_count=len(temp_wa_wb_peak_list)
    
    all_line_count=len(temp_wa_wb)
    
    return_info_dic={'gene_count_series': temp_wa_wb_value_count,
                 'gene_list': temp_wa_wb_value_count_gene_list, 
                 'gene_count': temp_wa_wb_gene_count,
                 'peak_count_series': temp_wa_wb_peak_value_count_series,
                 'peak_list': temp_wa_wb_peak_list,
                 'peak_count':temp_wa_wb_peak_count,
                 'all_count': all_line_count}
    
    return return_info_dic


def get_region_p2gl_info_input_SE_str(temp_SE_input, temp_input_p2gl_bed):
    chr_str,start_str,end_str=temp_SE_input.split("_")
    input_chr_str=chr_str
    input_start_pos=int(start_str)
    input_end_pos=int(end_str)
    
    return_dic=get_region_p2gl_info(input_chr_str,
                         input_start_pos,
                         input_end_pos,
                         temp_input_p2gl_bed)
    return return_dic
    




def return_all_se_concat_full_df(
    input_se_file_list: List[str], 
    filter_bed: BedTool
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    cols = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "sample_ID", "sample_name"]
    temp_concat_for_bed_df = pd.DataFrame(columns=cols)
    
    sample_num = 0
    
    for path_now_processing_sample in input_se_file_list:
        temp_return_bed_dic, _ = se2gene_file_output_only_rose_se(
            filter_bed, path_now_processing_sample
        )
        
        temp_file_name = path_now_processing_sample.name
        temp_sample_name = temp_file_name.split("_R1_peaks_")[0]
        
        temp_concat_sample_df = temp_return_bed_dic["link_gene"]
        temp_concat_sample_df["sample_ID"] = sample_num
        temp_concat_sample_df["sample_name"] = temp_sample_name
        
        temp_concat_for_bed_df = pd.concat([temp_concat_for_bed_df, temp_concat_sample_df], axis=0)
        
        sample_num += 1
    
    no_sort_df = temp_concat_for_bed_df.copy()
    
    result_bed_df = sort_df_for_bed_by_chr(temp_concat_for_bed_df)
    
    return result_bed_df, no_sort_df

def return_merge_se_count_df_full(concat_df_full):


    bed_file_result_bed_df=BedTool.from_dataframe(concat_df_full)

    test_merge_bed_count=bed_file_result_bed_df.merge(c=1,o="count")
    
    df_test_merge_bed_count=test_merge_bed_count.to_dataframe()
    
    series_se_data = df_test_merge_bed_count["chrom"] \
    + "_" \
    + df_test_merge_bed_count["start"].astype(str) \
    + "_" \
    + df_test_merge_bed_count["end"].astype(str)
    
    df_se_merge_count=pd.DataFrame({"se_data": series_se_data,
             "se_count" : df_test_merge_bed_count["name"] })
    
    

    
    
    return df_se_merge_count



def sort_se_merge_SAMPLE_count_df(merge_se_count_df):
    sort_by_count_df_se_merge_count=merge_se_count_df.copy().sort_values("sample_count",  ascending=False)
    return sort_by_count_df_se_merge_count


def sort_se_merge_SE_count_df(merge_se_count_df):
    sort_by_count_df_se_merge_count=merge_se_count_df.copy().sort_values("se_count",  ascending=False)
    return sort_by_count_df_se_merge_count



def get_series_se_data(df_SE_merge_bed_count):

    series_se_data = df_SE_merge_bed_count["chrom"] \
                        + "_" \
                        + df_SE_merge_bed_count["start"].astype(str) \
                        + "_" \
                        + df_SE_merge_bed_count["end"].astype(str)

    return series_se_data




def get_temp_SE_concat_intersect_for_sample_count(temp_concat_full):
    
    bed_file_result_bed_df=BedTool.from_dataframe(temp_concat_full)
    
    test_merge_bed_count=bed_file_result_bed_df.merge(c=1,o="count")

    df_test_merge_bed_count=test_merge_bed_count.to_dataframe()
    
    series_se_data = get_series_se_data(df_test_merge_bed_count)
    
    df_test_merge_bed_count_test=df_test_merge_bed_count.copy()
    
    df_test_merge_bed_count_test["name"]=series_se_data
    
    bed_df_test_merge_bed_count_test=return_bed_form_df(df_test_merge_bed_count_test)

    bed_temp_concat_full=return_bed_form_df(temp_concat_full)

    temp_SE_concat_intersect_for_sample_count=get_wa_wb_df(bed_df_test_merge_bed_count_test, bed_temp_concat_full)
    
    
    return temp_SE_concat_intersect_for_sample_count





def return_merge_se_count_df_full_edit(
    concat_df_full: pd.DataFrame, 
    p2gl_bed: BedTool
) -> pd.DataFrame:


    bed_file_result_bed_df = BedTool.from_dataframe(concat_df_full)
    test_merge_bed_count = bed_file_result_bed_df.merge(c=1, o="count")

    df_test_merge_bed_count = test_merge_bed_count.to_dataframe()

    series_se_data = get_series_se_data(df_test_merge_bed_count)

    df_se_merge_count = pd.DataFrame({
        "se_data": series_se_data,
        "se_count": df_test_merge_bed_count["name"]
    })

    temp_se_concat_intersect = get_temp_SE_concat_intersect_for_sample_count(concat_df_full)

    df_se_merge_count["sample_count"] = series_se_data.apply(
        lambda se_name: temp_se_concat_intersect.query("name == @se_name")["blockStarts"].nunique()
    )

    get_region_info = partial(get_region_p2gl_info_input_SE_str, temp_input_p2gl_bed=p2gl_bed)

    dict_series_for_new_add_df = df_se_merge_count["se_data"].apply(get_region_info)

    key_list_for_add_mergeSE_df = [
        "gene_list", 
        "gene_count", 
        "peak_list", 
        "peak_count", 
        "all_count"
    ]


    for key in key_list_for_add_mergeSE_df:
        df_se_merge_count[key] = dict_series_for_new_add_df.apply(lambda x: x[key])

    return df_se_merge_count



def return_merge_se_count_df_full_with_info(concat_df_full, p2gl_bed):


    bed_file_result_bed_df=BedTool.from_dataframe(concat_df_full)


    test_merge_bed_count=bed_file_result_bed_df.merge(c=[1,4,4,4,4,
                                                         8,8,8,8],
                                                      o=["count","collapse","distinct","count","count_distinct",
                                                                   "collapse","distinct","count","count_distinct"])
    

    df_test_merge_bed_count=test_merge_bed_count.to_dataframe()
    
    

    
    series_se_data = get_series_se_data(df_test_merge_bed_count)
    
    df_se_merge_count=pd.DataFrame({"se_data": series_se_data,
             "se_count" : df_test_merge_bed_count["name"] })
    

    sample_count_list=[]
    
    
    temp_SE_concat_intersect_for_sample_count=get_temp_SE_concat_intersect_for_sample_count(concat_df_full)


    for temp_SE_name in series_se_data:

        temp_df_SE_wa_wb_for_count=temp_SE_concat_intersect_for_sample_count.query("name == @temp_SE_name")
        temp_SE_sample_list=temp_df_SE_wa_wb_for_count["blockStarts"].tolist()
        temp_SE_sample_set=set(temp_SE_sample_list)
        temp_SAMPLE_count=len(temp_SE_sample_set)

        sample_count_list.append(temp_SAMPLE_count)

    
    df_se_merge_count["sample_count"]=sample_count_list
    

    
    temp_check_df=df_se_merge_count.copy()
    

    
    get_region_with_fixed_p2gl = partial(get_region_p2gl_info_input_SE_str, temp_input_p2gl_bed=p2gl_bed)
    

    

    dict_series_for_new_add_df = temp_check_df["se_data"].apply(get_region_with_fixed_p2gl)

    key_list_for_add_mergeSE_df=["gene_list", 
                             "gene_count",
                             "peak_list", 
                             "peak_count", 
                             "all_count"]
    
    temp_new_df_for_add =  df_se_merge_count.copy()


    for key in key_list_for_add_mergeSE_df:
        temp_new_df_for_add[key] = dict_series_for_new_add_df.apply(lambda x: x[key])
    

    temp_new_df_for_add["SE_count_distinct"]=df_test_merge_bed_count["thickEnd"]
    temp_new_df_for_add["SAMPLE_count"]=df_test_merge_bed_count["blockSizes"]
    temp_new_df_for_add["SAMPLE_count_distinct"]=df_test_merge_bed_count["blockStarts"]
    temp_new_df_for_add["SE_collapse"]=df_test_merge_bed_count["score"]
    temp_new_df_for_add["SE_distinct"]=df_test_merge_bed_count["strand"]
    temp_new_df_for_add["SAMPLE_collapse"]=df_test_merge_bed_count["itemRgb"]
    temp_new_df_for_add["SAMPLE_distinct"]=df_test_merge_bed_count["blockCount"]

    

    return_df=temp_new_df_for_add.set_index("se_data")

    return return_df




def return_random_filename(name_length=8):
    alphabed_pool=string.ascii_letters + string.digits
    now_time = datetime.datetime.now()
    today_str=now_time.strftime('%Y%m%d')
    temp_str="".join(random.choices(alphabed_pool, k=name_length))
    temp_file_name=today_str+"_"+temp_str
    return temp_file_name




def extract_sample_bed_from_concat_df_full(concat_full, input_sample_name):
    
    sample_df=concat_full.query("sample_name == @input_sample_name")
    
    sample_bed=return_bed_form_df(sample_df)
    
    return sample_bed




def make_dic_of_ChIP_peak_str(tsv_path):
    
    temp_input_ChIP_tsv=tsv_path

    temp_chIP_all_tsv_df_for_dic= pd.read_table(temp_input_ChIP_tsv,
                                        index_col=0,
                                       usecols=[0, 1, 2, 3])
    temp_chIP_all_tsv_df_for_dic.columns = ["chr", "start" , "end"]
    
    temp_ChIP_dict=temp_chIP_all_tsv_df_for_dic.to_dict(orient='index')
    
    return temp_ChIP_dict

def return_bedpe_bed_from_p2gl_tsv(input_p2gl_path, 
                                    input_fdr=0.05,
                                    input_r=0.5):
    temp_p2gl_file= pd.read_table(input_p2gl_path)

    filtered_p2g_file=filter_p2g_df(temp_p2gl_file,
                                    fdr_filter=input_fdr,
                                    r_filter=input_r
                                   )
    
    
    temp_bedpe_df=pd.DataFrame()
    temp_bedpe_df["chrom1"]=filtered_p2g_file["chr"]
    temp_bedpe_df["start1"]=filtered_p2g_file["tx.s"]
    temp_bedpe_df["end1"]=filtered_p2g_file["tx.e"]
    temp_bedpe_df["chrom2"]=filtered_p2g_file["chr"]
    temp_bedpe_df["start2"]=filtered_p2g_file["Start"]
    temp_bedpe_df["end2"]=filtered_p2g_file["End"]
    temp_bedpe_df["name"]=filtered_p2g_file["PeakID"]+ "_" +filtered_p2g_file["symbol"]
    temp_bedpe_df["score"]=filtered_p2g_file["r"]
    temp_bedpe_df["strand1"]=filtered_p2g_file["strand"]
    temp_bedpe_df["strand2"]="."
    temp_bedpe_df["symbol"]=filtered_p2g_file["symbol"]
    temp_bedpe_df["PeakID"]=filtered_p2g_file["PeakID"]
    temp_bedpe_df["FDR"]=filtered_p2g_file["FDR"]

    temp_bedpe_bed=return_bed_form_df(temp_bedpe_df)
    
    return temp_bedpe_bed


def add_region_data_to_info(info_df):
    df=info_df.copy()
    
    df['chr'] = df.index.str.split('_').str[0]
    df['start'] = df.index.str.split('_').str[1].astype(int)
    df['end'] = df.index.str.split('_').str[2].astype(int)

    columns = ['chr', 'start', 'end'] + [col for col in df.columns if col not in ['chr', 'start', 'end']]
    df = df[columns]

    return df

def sort_gtf(gtf_file, output_file):
    gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)


    gtf_sorted = gtf.sort_values(by=[0, 3, 4])  

    gtf_sorted.to_csv(output_file, sep='\t', header=False, index=False)


def modify_gtf_line(line):
    fields = line.strip().split('\t')
    
    fields[2] = 'exon'
    
    chrom = fields[0]
    start = fields[3]
    end = fields[4]
    
    gene_id = f"{chrom}_{start}_{end}"
    
    gene_name = f"{chrom}_{start}"
    
    new_attrs = f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{gene_id}";'
    fields[8] = new_attrs
    
    return '\t'.join(fields)

def create_bed6_from_df(df):

    selected_columns = df.iloc[:, :6]
    
    selected_columns.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand']

    bed = BedTool.from_dataframe(selected_columns)
    return bed






def se_list_to_filtered_region_bed_custom(
    list_se: List[str], 
    df_p2gl: pd.DataFrame, 
    threshold: int = 0
) -> Tuple[BedTool, BedTool, BedTool]:

    bed_p2gl_filter = bed_filter_from_p2g_df(df_p2gl)
    
    list_processed_df = [
        return_filtered_bed_df_from_se_file_and_filter_bed(path_temp, bed_p2gl_filter)
        for path_temp in list_se
    ]
    
    df_concat = return_concat_df(list_processed_df)
    
    bed_df_concat_sort = sort_bed(BedTool.from_dataframe(df_concat))
    
    bed_count = bed_df_concat_sort.merge(c=2, o="count")
    
    bed_bed6_count = return_merge_se_bed_bed6(bed_count)
    
    bed_bed6_count_fix = add_sample_count_to_bed6(
        bed_bed6_count,
        bed_df_concat_sort,
        list_se
    )
    
    bed_before_threshold = bed_bed6_count_fix
    
    if threshold >= 1:
        df_threshold_filtered = bed_bed6_count_fix.to_dataframe().query(f"strand >= {threshold}")
        bed_after_threshold = BedTool.from_dataframe(df_threshold_filtered)
    else:
        bed_after_threshold = bed_before_threshold
    
    return bed_after_threshold, bed_before_threshold, bed_df_concat_sort
