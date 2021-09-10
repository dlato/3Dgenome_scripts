#!/usr/local/bin/python3
#####################################
# 
# Program will convert HiC-Pro output to 3Dflow input
#
# Developed by Daniella F. Lato 
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: hicproTo3Dflow.py BED_FILE MATRIX_FILE > OUTPUT_PATH
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
chrn = []
pos = []
beddf = pd.read_csv(bed, delimiter='\t', names=['chrA', 'st', 'end'])
uchrs = beddf.chrA.unique()
ldf = len(beddf)
#vector with starts of each chrom
for A in range(0,ldf):
    if A+1 != ldf:
        if beddf.chrA[A] != beddf.chrA[A+1]:
           pos.append(A+1)
           chrn.append(beddf.chrA[A+1])

c = 0
#loop through each row and print all interactions with that row
lastchr = beddf.chrA[0]
for A in range(0,ldf):
    if beddf.chrA[A] == beddf.chrA[ldf-1]:
        break
    chrloop = beddf.chrA[pos[c]:]
    chrloop = chrloop.to_numpy()
    stloop = beddf.st[pos[c]:]
    stloop = stloop.to_numpy()
    endloop = beddf.end[pos[c]:]
    endloop = endloop.to_numpy()
    ll = len(chrloop)
    for B in range(0,ll):
        print(beddf.chrA[A], beddf.st[A], beddf.end[A], chrloop[B], stloop[B], endloop[B])
    if lastchr != beddf.chrA[A]:
        c += 1
    lastchr = beddf.chrA[A]

