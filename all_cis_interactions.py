#!/usr/local/bin/python3
#####################################
# 
# obtain all possible trans interactions for the genome
#
# Developed by Daniella F. Lato 
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: all_cis_interactions.py BED_FILE MATRIX_FILE > OUTPUT_PATH
# input:
#         - bed file (full path)
# standard output:
#         - tab separated file with columns: chrA startA endA chrB startB endB
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
bed = sys.argv[1]
#bed file
beddf = pd.read_csv(bed, delimiter='\t', names=['chrA', 'st', 'end'])
uchrs = beddf.chrA.unique()
ldf = len(beddf)
#vector with starts of each chrom
chrn = []
pos = []
for A in range(0,ldf):
    if A == 0:
        pos.append(A)
        chrn.append(beddf.chrA[A])
    if A+1 != ldf:
        if beddf.chrA[A] != beddf.chrA[A+1]:
           pos.append(A+1)
           chrn.append(beddf.chrA[A+1])
#append last line to vec
pos.append(ldf)

#loop through each row and print all cis interactions with that row
lastchr = beddf.chrA[0]
for A in range(0,len(pos)):
    if A == len(pos)-1:
        break
    chrloop = beddf.chrA[pos[A]:pos[A+1]]
    chrloop = chrloop.to_numpy()
    stloop = beddf.st[pos[A]:pos[A+1]]
    stloop = stloop.to_numpy()
    endloop = beddf.end[pos[A]:pos[A+1]]
    endloop = endloop.to_numpy()
    ll = len(chrloop)
    c = 0
    for B in range(0,ll):
        for C in range(c,ll):
            print(chrloop[B], stloop[B], endloop[B], chrloop[C], stloop[C], endloop[C])
        c += 1
