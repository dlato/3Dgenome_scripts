#!/usr/local/bin/python3
#####################################
# 
# Program will convert HiC-Pro output to 3Dflow input
#
# Developed by Daniella F. Lato 
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: hicproTo3Dflow.py BED_FILE MATRIX_FILE OUTPUT_PATH RESOLUTION
# input:
#         - bed file (full path)
#         - matrix file (full path)
#         - output path
#         - resolution
# output: two files:
#         - trans.RESOLUTION_iced.sorted.txt
#         - cis.RESOLUTION_iced.sorted.txt
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
#mat = sys.argv[2]
#out = sys.argv[3]
#res = sys.argv[4]
##opening files to write to
#out_path = out
#out_cis = out + "/cis." + res + "_iced.sorted.txt"
#out_trans = out + "/trans." + res + "_iced.sorted.txt"
#out_fcis = open(out_cis,"w+")
#out_ftrans = open(out_trans,"w+")

#bed file
chrn = []
pos = []
beddf = pd.read_csv(bed, delimiter='\t', names=['chrA', 'st', 'end'])
print(beddf.head())
print(beddf.tail())

uchrs = beddf.chrA.unique()
print(uchrs)
ldf = len(beddf)
print(ldf)
print(beddf.chrA[0])
print(beddf.chrA[1])

#vector with starts of each chrom
for A in range(0,ldf):
    if A+1 != ldf:
        if beddf.chrA[A] != beddf.chrA[A+1]:
           pos.append(A+1)
           chrn.append(beddf.chrA[A+1])

c = 0
#print(beddf.chrA[pos[0]:])
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

