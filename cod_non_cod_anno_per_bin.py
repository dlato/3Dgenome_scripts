#!/usr/local/bin/python3
#
#####################################
#
# go through gff file and create a new column to re-code the classification of the gene to match the GENCODE classification (https://ftp.ebi.ac.uk/pub/databases/gencode/_README_stats.txt)
# used for "annotation" of common interactions
#
# Developed by Daniella F. Lato
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: cod_non_cod_anno_per_bin.py GFF_FILE BIN_SIZE OUTPUT_FILE_PATH_AND_NAME
# input:
#         - gff file
#         - bin size (bp, numeric)
#         - full path and file name for output annotation file
#
#####################################

#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
from gtfparse import read_gtf
#import Bio
import pandas as pd
import numpy as np
from math import floor
#from Bio import SeqIO
#import pprint
#import gffutils
#from BCBio.GFF import GFFExaminer
#from BCBio import GFF
# FOR DEBUGGING ONLY pdb.set_trace()

gff_file = sys.argv[1]
bin_size = sys.argv[2]
outfile = sys.argv[3]
df = read_gtf(gff_file)

#lncRNA_vals <- ["processed_transcript", "lincRNA", "3prime_overlapping_ncrna", "antisense",  "non_coding", "sense_intronic", "sense_overlapping", "TEC", "known_ncrna", "macro_lncRNA",  "bidirectional_promoter_lncrna", "lncRNA"]
# filter DataFrame to gene entries
df_genes = df[df["feature"] == "gene"]
df_genes.loc[:,'broad_class'] = "not_prot"
# protein coding
df_genes.loc[df_genes['gene_type'] == "protein_coding", 'broad_class'] = 'prot'
# lncRNA
df_genes.loc[df_genes['gene_type'] == "processed_transcript", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "lincRNA", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "3prime_overlapping_ncrna", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "antisense", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "non_coding", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "sense_intronic", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "sense_overlapping", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "TEC", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "known_ncrna", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "macro_lncRNA", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "bidirectional_promoter_lncrna", 'broad_class'] = 'lncRNA'
df_genes.loc[df_genes['gene_type'] == "lncRNA", 'broad_class'] = 'lncRNA'
# sncRNA
df_genes.loc[df_genes['gene_type'] == "snRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "snoRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "rRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "Mt_tRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "Mt_rRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "misc_RNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "miRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "ribozyme", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "sRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "scaRNA", 'broad_class'] = 'sncRNA'
df_genes.loc[df_genes['gene_type'] == "vaultRNA", 'broad_class'] = 'sncRNA'
# pseudogene"
df_genes.loc[df_genes['gene_type'].str.contains("pseudogene"), 'broad_class'] = 'pseudogene'
# Immunoglobulin/T-cell receptor gene segments
df_genes.loc[df_genes['gene_type'].str.contains("TR_"), 'broad_class'] = 'Tcell'
df_genes.loc[df_genes['gene_type'].str.contains("IG_"), 'broad_class'] = 'immunoglob'
# dealing with odd scRNA (classify as lncRNA)
df_genes.loc[df_genes['gene_type'] == "scRNA", 'broad_class'] = 'lncRNA'
##print(df_genes['gene_type'].str.contains("pseudogene"))
##print(df[df['gene_type'].str.contains("pseudogene")])
##df_genes[df_genes['gene_type'].str.contains('pseudogene', case=False, na=False)]
#print(df_genes.loc[df_genes['gene_id'] == "ENSG00000186092.6"])
#print(df_genes.loc[df_genes['gene_id'] == "ENSG00000243485.5"])
#print(df_genes.loc[df_genes['gene_id'] == "ENSG00000278267.1"])
#print(df_genes.loc[df_genes['gene_id'] == "ENSG00000233084.2"])
#print(df_genes)
##df_genes['broad_class'] = np.where(df['gene_type'] == "protein_coding", "prot", "not_prot")
##print(df_genes)
##df_genes_chr1 = df_genes[df_genes["seqname"] == "1"]

#add col for what bin the gene is in, listed as START of bin
def rounddown(x,y):
    return int(floor(x / int(y))) * int(y)
df_genes['bin_start'] = df_genes.apply(lambda x: rounddown(x.start, bin_size), 1)
df_genes['bin_end'] = df_genes.apply(lambda x: rounddown(x.end, bin_size), 1)
df_genes = df_genes.replace('', np.nan, regex=True)
print(df_genes)
#print("genes that are between bins")
#print(df_genes[df_genes['bin_start'] != df_genes['bin_end']])
df_genes.to_csv(outfile, sep="\t", index=False, na_rep='NULL')
