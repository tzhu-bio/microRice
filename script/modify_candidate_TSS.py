#!/usr/bin/env python

import argparse
import pyranges as pr
import pandas as pd
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process miRNA data.")
    parser.add_argument("--distance_threshold", type=int, default=5000, help="Distance threshold for filtering")
    return parser.parse_args()  


if __name__ == "__main__":
    
    args = parse_arguments()
    distance_threshold = args.distance_threshold
    pd.options.mode.chained_assignment = None 
    
    intragenic_same_strand=pd.read_csv("./miRNA_classification.txt",sep="\t",header=None)
    intragenic_same_strand_list=intragenic_same_strand[intragenic_same_strand[1]=="intragenic_same_strand"][0]

    tss_candidate=pd.read_csv("./tss_candidate.csv",sep="\t")
    tss_candidate=pr.PyRanges(tss_candidate)
    intra_tss=tss_candidate[tss_candidate.Name.isin(intragenic_same_strand_list)]
    intra_tss_longer_than_threshold = intra_tss[intra_tss.Distance > distance_threshold]
    tss1=tss_candidate[~tss_candidate.Name.isin(intra_tss_longer_than_threshold.df['Name'])]
    tss1.category='H3K4me3_enriched_zone'

    intragenic_host_tss = pd.read_csv("./intragenic_miRNA_samestrand_hostgene.csv", sep="\t")

    # 根据条件给 tss_midpoint 赋值
    intragenic_host_tss['tss_start'] = np.where(intragenic_host_tss['Strand_b'] == '+', 
                                                intragenic_host_tss['Start_b'], 
                                                intragenic_host_tss['End_b'])
    intragenic_host_tss['tss_end'] = np.where(intragenic_host_tss['Strand_b'] == '+', 
                                              intragenic_host_tss['tss_start'] - 1000, 
                                              intragenic_host_tss['tss_start'] + 1000)

    # 将 tss_start 和 tss_end 进行排序
    intragenic_host_tss[['tss_start', 'tss_end']] = np.sort(intragenic_host_tss[['tss_start', 'tss_end']], axis=1)

    # 将 tss_start 和 tss_end 列转换为整数类型
    intragenic_host_tss['tss_start'] = intragenic_host_tss['tss_start'].astype(int)
    intragenic_host_tss['tss_end'] = intragenic_host_tss['tss_end'].astype(int)

    intra_cotranscribe_with_host=intragenic_host_tss[intragenic_host_tss['Name'].isin(intra_tss_longer_than_threshold.df['Name'])]
    intra_cotranscribe_with_host['sep']="."
    intra_cotranscribe_with_host['Distance'] = np.where(intra_cotranscribe_with_host['Strand_b'] == '+', 
                                              intra_cotranscribe_with_host['Start'] - intra_cotranscribe_with_host['tss_end'], 
                                              intra_cotranscribe_with_host['tss_start'] - intra_cotranscribe_with_host['End'])
    intra_cotranscribe_with_host=intra_cotranscribe_with_host.loc[:,["Chromosome", "tss_start", "tss_end", "Name", "sep", "Strand", "Start", "End", "Distance"]]

    tss1.sep="."
    tss1=tss1.df[["Chromosome","Start_b","End_b","Name","sep","Strand_b","Start","End","Distance","score","category"]]
    intra_cotranscribe_with_host.columns=["Chromosome","Start_b","End_b","Name","sep","Strand_b","Start","End","Distance"]
    intra_cotranscribe_with_host["score"]= np.nan
    intra_cotranscribe_with_host["category"]="hostgene_TSS"
    tss_total=pd.concat([tss1, intra_cotranscribe_with_host], axis=0)
    tss_total['Chromosome']=tss_total['Chromosome'].astype(str)
    tss_total['Start_b']=tss_total['Start_b'].astype(int)
    tss_total['End_b']=tss_total['End_b'].astype(int)
    tss_total_sorted = tss_total.sort_values(by=[tss_total.columns[0], tss_total.columns[1], tss_total.columns[2]])
    tss_total_sorted.to_csv("./miRNA_TSS.bed",sep="\t",header=None,index=False)
    
    print("miRNA_TSS.bed")

