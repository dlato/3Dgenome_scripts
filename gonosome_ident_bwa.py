#!/usr/local/bin/python3
#####################################
# 
# Program will calculate mapped read coverage per chromosome
# can be used to identifiy biological sex
# FOR READS ALIGNED WITH BWA
#
# Developed by Daniella F. Lato 
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: gonosome_ident.py READ_LENGTH CHROM_LENGTH_FILE NUMBER_OF_MAPPED_READS_FILE > OUTPUT_FILE
# input:
#         - sequencing read length (number)
#         - chromosome sizes file (full path, tab separated: chromosome length)
#         - *created with script* number of mapped reads per chromosome file (full path, chromosome on one line, corresponding number of mapped reads on next)
#         - name of output file (full path, string)
# standard out:
#         - space separated file with columns: chromosome mapped_coverage
#
#####################################
import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import csv
import Bio
import numpy as np
import pandas as pd
import os #importing operating system commands
import glob #importing glob which helps read files

# read in files
readlen = sys.argv[1] # input read length
chrlen_file = sys.argv[2] # chromosome length file
clen_df = pd.read_csv(chrlen_file, delimiter='\t', names=['chr', 'len'])
mappedfile = sys.argv[3] # mapped reads per chromosome file
chrom_l = []
mapped_l = []
with open(mappedfile) as f:
        for line in f:
            line = line.strip()
            if "chr" in line:
                chrom_l.append(line)
            else:
                mapped_l.append(line)
mapped_df = pd.DataFrame({'chr':chrom_l,
                              'mapped':mapped_l})

#claculate coverage per chromosome
print("chromosome mapped_coverage")
for c in chrom_l:
    N = mapped_df.loc[mapped_df['chr'] == c, 'mapped']
    G = clen_df.loc[clen_df['chr'] == c, 'len']
    cov = int(N)*(int(readlen) / int(G))
    print(c,cov)
