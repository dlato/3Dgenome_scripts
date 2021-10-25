#!/usr/local/bin/python3
#####################################
#
# Program will read in clustalW alignment file and annotate CpG methylation changes in samples compared to reference
#
# Developed by Daniella F. Lato
# email: DaniellaLato@GMail.com
# GitHub: https://github.com/dlato
#
# to run: cpg_anno.py CLUSTALW_ALN_FILE
# input:
#         - clustalW alignment file, reference sequence MUST be first!
#         - chromosome sizes file (full path, tab separated: chromosome length)
#         - *created with script* number of mapped reads per chromosome file (full path, chromosome on one line, corresponding number of mapped reads on next)
#         - name of output file (full path, string)
#         - strand (options: "top" or "bottom"
# out to file:
#         - comma separated file with columns:
#                 NUCLEOTIDE_INDEX	REFERENCE_SEQ	REFERENCE_METHYLATION	SAMPLE1_METHYLATION	SAMPLE2_METHYLATION	...
#         one column per sequence in alignment (reference seq will always be the first column)
#
#####################################

import sys, re
import pdb # FOR DEBUGGING ONLY import pdb
import pandas as pd #importing pandas
import numpy as np #importing Numpy
import csv
from numpy import nan
from Bio import SeqIO
from Bio import AlignIO
from functools import partial
# FOR DEBUGGING ONLY pdb.set_trace()

#below will open the mafft file and get the aln info
align = AlignIO.read(sys.argv[1], "clustal")
# ALIGNMENT STARTS AT POS 0 IN THE ARRAY
seqnum = len(align)
#aln_len = len(align[4].seq)
#print(len(align[4].seq))
#print(align[:,661])

# annotate reference seq
# top/forward strand
ref = align[0].seq
refAnno = []
rseq = []
for i in range(len(ref)):
    nuc = ref[i]
    if i in [len(ref)-1, 0]:
        # dealing with gap in ref
        if nuc in ["-"]:
            refAnno.append("NA")
            rseq.append(nuc)
        #case when at beginning of aln
        if i == 0:
            nexnuc = ref[i+1]
            if nuc in ["C","c"]:
                if nexnuc in ["G","g"]:
                    refAnno.append("1")
                    rseq.append(nuc)
                else:
                    refAnno.append("NA")
                    rseq.append(nuc)
            else:
                refAnno.append("NA")
                rseq.append(nuc)
        #case when at end of aln
        if i == len(ref)-1:
            prenuc = ref[i-1]
            if nuc in ["G","g"]:
                if prenuc in ["c","C"]:
                    # changing G's to NA because only care about the C's and if they have changed
                    refAnno.append("NA")
#                    refAnno.append("1")
                    rseq.append(nuc)
                else:
                    refAnno.append("NA")
                    rseq.append(nuc)
            else:
                refAnno.append("NA")
                rseq.append(nuc)
    else:
        # dealing with gap in ref
        if nuc in ["-"]:
            refAnno.append("NA")
            rseq.append(nuc)
        nexnuc = ref[i+1]
        prenuc = ref[i-1]
        if nuc in ["C","c"]:
            if nexnuc in ["G","g"]:
                refAnno.append("1")
                rseq.append(nuc)
            else:
                refAnno.append("NA")
                rseq.append(nuc)
        if nuc in ["G","g"]:
            if prenuc in ["C","c"]:
                # changing G's to NA because only care about the C's and if they have changed
#                refAnno.append("1")
                refAnno.append("NA")
                rseq.append(nuc)
            else:
                refAnno.append("NA")
                rseq.append(nuc)
        if nuc in ["T","t", "A","a"]:
            refAnno.append("NA")
            rseq.append(nuc)
              
allseqAnno = []
#annotate other samples
for s in range(1,seqnum):
    seq = align[s].seq
    seqAnno = []
    for n in range(len(seq)):
        nuc = seq[n]
        # not a CpG site
        if refAnno[n] in ["NA"]:
            seqAnno.append("NA")
        else:
            #matches ref
            if nuc in ref[n]:
                seqAnno.append("1")
            #does not match ref
            else:
                # gap in aln (no seq info)
                if nuc in ["-"]:
                    seqAnno.append("gap")
                # mismatch
                else:
                    # ambig nuc (K,N..etc)
                    if nuc not in ["A", "a","C","c","G","g","T","t"]:
                        seqAnno.append("oddNuc")
                    else:
                        seqAnno.append("0")
    allseqAnno.append(seqAnno)



# print out CpG ident for all sequences (including the ref)
# get sequence IDs from alignment
IDs = []
for i in range(0,seqnum):
    tmp_id = align[i].id
    IDs.append(tmp_id)
#add ref annotation to new list
anno = allseqAnno
anno.insert(0, refAnno)

#save info as pandas df
df = pd.DataFrame(anno)
df = df.transpose()
df.columns = IDs

#write to file
out = sys.argv[2]
#strand info
strand = sys.argv[3]
if strand in ["bottom"]:
    #add ref seq as column
    colname = IDs[0] + "_seq"
    df = df.assign(tmpcol = rseq)
    first_column = df.pop("tmpcol")
    df.insert(0,"tmpcol",rseq)
    df.columns = [colname if x=='tmpcol' else x for x in df.columns]
    # reverse order of df
    df = df.iloc[::-1]
    #add in nuc index starting at 1
    ind = [i for i in range(1,len(df.index)+1)]
    df = df.assign(NucPos = ind)
    first_column = df.pop('NucPos')
    df.insert(0,'NucPos',ind)
    #print to csv
    df.to_csv(out, index=False)
else:
    #add ref seq as column
    colname = IDs[0] + "_seq"
    df = df.assign(tmpcol = rseq)
    first_column = df.pop("tmpcol")
    df.insert(0,"tmpcol",rseq)
    df.columns = [colname if x=='tmpcol' else x for x in df.columns]
    #add in nuc index starting at 1
    ind = [i for i in range(1,len(df.index)+1)]
    df = df.assign(NucPos = ind)
    first_column = df.pop('NucPos')
    df.insert(0,'NucPos',ind)
    #print to csv
    df.to_csv(out, index=False)
