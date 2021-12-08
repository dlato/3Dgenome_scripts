#!/usr/local/bin/python3
#
# Program will go through the genbank file and will obtain the locus
# ID for each gene as well as the start and end positions
# run as "get_gene_pos_locus.py *.gbk > *gene_locus_location"
#
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
from gtfparse import read_gtf
#import Bio
import pandas as pd
import numpy as np
#from Bio import SeqIO
#import pprint
#import gffutils
#from BCBio.GFF import GFFExaminer
#from BCBio import GFF
# FOR DEBUGGING ONLY pdb.set_trace()

gff_file = sys.argv[1]
bin_size = sys.argv[3]
print("bin_size")
print(bin_size)
print("gff_file")
print(gff_file)

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



#examiner = GFFExaminer()
#in_handle = open(gff_file)
#print(pprint.pprint(examiner.parent_child_map(in_handle)))
#in_handle.close()
#
#in_handle = open(gff_file)
#for rec in GFF.parse(in_handle):
#    print(rec)
#in_handle.close()
#faa_filename = sys.argv[2]
#output_handle = open(faa_filename, "w")
#output_handle.write("gbk_start\tgbk_end\tgbk_midpoint\tgbk_gene_id\tgbk_locus_tag\tgbk_old_locus_tag\tgbk_strand\n")
#for gb_record in SeqIO.parse(open(sys.argv[1],"r"), "genbank") :
#    genes = gb_record.features
##    print(genes)
#    for seq_feature in gb_record.features :
##        if seq_feature.type=="RNA":
#        m = re.search("[a-z]RNA|CDS",seq_feature.type)
#        if m:
##        output_handle.write(str(type(seq_feature.type)))
##        match_object = re.search(r'*RNA', str(seq_feature.type))
##        output_handle.write(match_object)
##        if (seq_feature.type=="CDS" or match_object):
#            if (seq_feature.qualifiers):
##                if ('pseudo' in seq_feature.qualifiers):
##                    continue
##                #gene is not a pseudogene
##                else:
#                end = seq_feature.location.end
#                start = seq_feature.location.start
#                start = start + 1
#                tmp = end + start
#                midpoint = int(tmp/2)
#                if seq_feature.strand == 1:
#               #leading strand = 0
#                    cod_val = 0
#                else:
#               #lagging strand = 1
#                    cod_val = 1
#                #dealing with the fact that there might not be a
#                #gene id present
#                if ('gene' in seq_feature.qualifiers):
#                    gene_id = seq_feature.qualifiers['gene'][0]
#                #gene is is missing
#                else:
#                    gene_id = 'NA'
#                #dealing with the fact that there might not be a
#                #locus tag present
#                if ('locus_tag' in seq_feature.qualifiers):
#                    locus_tag = seq_feature.qualifiers['locus_tag'][0]
#                #gene is is missing
#                else:
#                    locus_tag = 'NA'
#                #dealing with the fact that there might not be a
#                #old locus tag present
#                if ('old_locus_tag' in seq_feature.qualifiers):
#                    old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
#                #gene is is missing
#                else:
#                    old_locus_tag = 'NA'
#                #print out all the info
#                output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
#                       start,
#                       seq_feature.location.end,
#                       midpoint,
#                       gene_id,
#                       locus_tag,
#                       old_locus_tag,
#                       cod_val))

#faa_filename = sys.argv[2]
#output_handle = open(faa_filename, "w")
#output_handle.write("gbk_start\tgbk_end\tgbk_midpoint\tgbk_gene_id\tgbk_locus_tag\tgbk_old_locus_tag\tgbk_strand\n")
#for gb_record in SeqIO.parse(open(sys.argv[1],"r"), "genbank") :
#    genes = gb_record.features
##    print(genes)
#    for seq_feature in gb_record.features :
##        if seq_feature.type=="RNA":
#        m = re.search("[a-z]RNA|[a-z][a-z]RNA",seq_feature.type)
#        if m:
##        output_handle.write(str(type(seq_feature.type)))
##        match_object = re.search(r'*RNA', str(seq_feature.type))
##        output_handle.write(match_object)
##        if (seq_feature.type=="CDS" or match_object):
#            if (seq_feature.qualifiers):
#                if ('pseudo' in seq_feature.qualifiers):
#                    continue
#                #gene is not a pseudogene
#                else:
#                    end = seq_feature.location.end
#                    start = seq_feature.location.start
#                    start = start + 1
#                    tmp = end + start
#                    midpoint = int(tmp/2)
#                    if seq_feature.strand == 1:
#                        #leading strand = 0
#                        cod_val = 0
#                    else:
#                        #lagging strand = 1
#                        cod_val = 1
#                    #dealing with the fact that there might not be a
#                    #gene id present
#                    if ('gene' in seq_feature.qualifiers):
#                        gene_id = seq_feature.qualifiers['gene'][0]
#                    #gene is is missing
#                    else:
#                        gene_id = 'NA'
#                    #dealing with the fact that there might not be a
#                    #locus tag present
#                    if ('locus_tag' in seq_feature.qualifiers):
#                        locus_tag = seq_feature.qualifiers['locus_tag'][0]
#                    #gene is is missing
#                    else:
#                        locus_tag = 'NA'
#                    #dealing with the fact that there might not be a
#                    #old locus tag present
#                    if ('old_locus_tag' in seq_feature.qualifiers):
#                        old_locus_tag = seq_feature.qualifiers['old_locus_tag'][0]
#                    #gene is is missing
#                    else:
#                        old_locus_tag = 'NA'
#                #print out all the info
#                output_handle.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
#                       start,
#                       seq_feature.location.end,
#                       midpoint,
#                       gene_id,
#                       locus_tag,
#                       old_locus_tag,
#                       cod_val))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#output_handle.write("0\t0\t0\t0\t0\t0\t0\n")
#
#output_handle.close()
#
