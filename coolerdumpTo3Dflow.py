#!/usr/local/bin/python3
#####################################
# 
# Program will convert a cooler dumped file (plain text version of cooler) to 3Dflow input format
#
# Developed by Daniella F. Lato 
# email: daniellalato@gmail.com
# GitHub: https://github.com/dlato
#
# to run: coolerdumpTo3Dflow.py COOLER_DUMPED_FILE OUTPUT_PATH RESOLUTION
# input:
#         - cooler dumped tab separated file (full path and file name)
#         - output path
#         - resolution (bp, numeric)
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
mat = sys.argv[1]
out = sys.argv[2]
res = sys.argv[3]
#opening files to write to
out_path = out
out_cis = out + "/cis." + res + "_iced.sorted.txt"
out_trans = out + "/trans." + res + "_iced.sorted.txt"
out_fcis = open(out_cis,"w+")
out_ftrans = open(out_trans,"w+")

#read in cooler dump file
with open(mat, 'r') as f:
    for row in f:
        sub_row = re.split("\t",row)
# remove rows with no balanced interaction value
        if (sub_row[7] != "\n"):
            A = sub_row[0] 
            B = sub_row[3] 
            out_id = "anchor_" + A + "/" + sub_row[1] + "/" + sub_row[2] + "__target_" + B + "/" + sub_row[4] + "/" + sub_row[5] + "\t" + f'{float(sub_row[7]):.8f}' + "\n"
#        print(out_id)
#print cis to cis and trans to trans
            A_spl = re.split("/",A)
            B_spl = re.split("/",B)
            if (A_spl[0] == B_spl[0]):
                out_fcis.write(str(out_id))
            else:
                out_ftrans.write(str(out_id))
