#!/usr/bin/env python

import sys
import pyranges as pr
import argparse

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Process files.")  
    parser.add_argument("-a", "--miRNAfile", required=True, help="Path to miRNA BED file")
    parser.add_argument("-b", "--genefile", required=True, help="Paths to result files")    
    args = parser.parse_args()

    miRNA = pr.read_bed(args.miRNAfile)
    gene = pr.read_bed(args.genefile)

    intragenic_same_strand = miRNA.overlap(gene, strandedness="same")
    intragenic_same_strand_list = set(intragenic_same_strand.df['Name'])
    
    
    result = {}
    for name in miRNA.df['Name']:
        if name in intragenic_same_strand_list:
            result[name] = 'intragenic_same_strand'
        else:
            result[name] = 'others'

    with open('miRNA_classification.txt', 'w') as f:
        for miRNA_name, category in result.items():
            f.write(f"{miRNA_name}\t{category}\n")   
    
    intra=miRNA.join(gene, how = 'left')
    intra_host_gene=intra[intra.End_b != -1]
    intra_host_gene=intra_host_gene[intra_host_gene.Strand == intra_host_gene.Strand_b]
    intra_host_gene_df=intra_host_gene[['Name',"Start_b","End_b","Strand_b","Name_b"]].df
    intra_host_gene_df.to_csv("intragenic_miRNA_samestrand_hostgene.csv", index=False,sep="\t")
    
    print("miRNA_classification.txt")
    print("intragenic_miRNA_samestrand_hostgene.csv")
