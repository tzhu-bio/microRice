#!/usr/bin/env python

import argparse
import pyranges as pr
import pandas as pd
import numpy as np
import glob
import os

def merge_results(score_cutoff, miRNA_file, *files):
    # 读取 miRNA BED 文件
    miRNA = pr.read_bed(miRNA_file)
    
    # 读取结果文件
    results = []
    for file in files:
        results.append(pd.read_csv(file, sep="\t"))
    
    # 合并结果数据框并选择最大的 score
    merged_df = pd.concat(results, ignore_index=True)
    merged_df = merged_df.loc[merged_df.groupby('sequence')['score'].idxmax()]
    
    # 拆分 sequence 列，并将结果赋值给新的列
    merged_df[['Chromosome', 'Start', 'End']] = merged_df['sequence'].str.split(":|-", expand=True)
    
    # 删除 sequence 列
    merged_df.drop(columns=['sequence'], inplace=True)
    
    # 筛选出 score 大于指定阈值的行
    merged_TSS_candidate = merged_df.loc[merged_df['score'] > score_cutoff]
    
    # 将结果写入新文件
    return(merged_TSS_candidate)

def unstranded_upstream(peaks, miRNA):
    peaks = peaks.copy()
    peaks.ID = np.arange(len(peaks))
    p = pr.concat([peaks, peaks]).sort()
    p.Strand = np.tile(["+", "-"], int(len(p)/2))
    miRNAminus = miRNA["-"].nearest(p["-"], how="upstream")
    miRNAplus = miRNA["+"].nearest(p["+"], how="upstream")
    df = pr.concat([miRNAminus, miRNAplus]).df
    return pr.PyRanges(df)

if __name__ == "__main__":
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="Process files.")    
    
    parser.add_argument("-a", "--miRNAfile", required=True, help="Path to miRNA BED file")
    parser.add_argument("-b", "--resultfiles", nargs='+', required=True, help="Paths to result files")
    parser.add_argument("-cutoff", "--score_cutoff", type=float, default=0.6, required=True, help="Score cutoff value")
    
    args = parser.parse_args()
    
    files = []
    for pattern in args.resultfiles:
        files.extend(glob.glob(pattern))

    # 执行合并结果的操作
    peaks=merge_results(args.score_cutoff, args.miRNAfile, *files)

    
    # 执行 unstranded_upstream 操作
    peaks = pr.PyRanges(peaks)
    miRNA = pr.read_bed(args.miRNAfile)
    tss_candidate = unstranded_upstream(peaks, miRNA)
    tss_candidate=tss_candidate[['Name',"Start_b","End_b","Strand_b","Distance","score"]]
    tss_candidate.df.to_csv("tss_candidate.csv", index=False,sep="\t")
       
    
    print("tss_candidate.csv")
