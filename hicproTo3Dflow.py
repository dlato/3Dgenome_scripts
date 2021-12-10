#!/usr/local/bin/python3
#####################################
# 
# Program will convert HiC-Pro output to 3Dflow input
#
# Developed by Daniella F. Lato 
# email: daniellalato@gmail.com
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
import os #importing operating system commands
import glob #importing glob which helps read files

# read in files
bed = sys.argv[1]
mat = sys.argv[2]
out = sys.argv[3]
res = sys.argv[4]
#opening files to write to
out_path = out
out_cis = out + "/cis." + res + "_iced.sorted.txt"
out_trans = out + "/trans." + res + "_iced.sorted.txt"
out_fcis = open(out_cis,"w+")
out_ftrans = open(out_trans,"w+")

#bed file
pos = []
ID = 0
with open(bed, 'r') as f:
    for row in f:
        sub_row = re.split("\t",row)
        tmp_pos = sub_row[0] + "/" + sub_row[1] + "/" + sub_row[2]
        pos.append(tmp_pos)
        ID += 1

#matrix file
with open(mat, 'r') as f:
    for row in f:
        sub_row = re.split("\t",row)
        idA = int(sub_row[0]) -1 
        idB = int(sub_row[1]) -1 
#link matrix with bed
        A = pos[idA]
        B = pos[idB]
        out_id = "anchor_" + A + "__target_" + B + "\t" + sub_row[2]
#print cis to cis and trans to trans
        A_spl = re.split("/",A)
        B_spl = re.split("/",B)
        if (A_spl[0] == B_spl[0]):
            out_fcis.write(str(out_id))
        else:
            out_ftrans.write(str(out_id))
