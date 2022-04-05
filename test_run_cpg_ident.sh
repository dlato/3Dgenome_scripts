#!/bin/bash
#PBS -l nodes=1:ppn=40,mem=20gb,walltime=1:00:00 #specify memory requirements
#PBS -j eo
#PBS -N Trans_all_inters_df #name your job
#PBS -V
#PBS -m e #sends you an email when your job is done

module load python/3.8.1
module load perl/5.28.0
module load clustalw/2.1
module load emboss/6.6.0

#!/bin/sh

#######################################
# DO NOT EDIT
cod=/hpf/largeprojects/pmaass/Daniella/Katerina/CpG_ident
read_len=100
chrom_len_file=/hpf/largeprojects/pmaass/resources/human/GRCh38.p10_v26/hg38.chrom.sizes
chrList=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts/human_chr_list.txt
#######################################


#######################################
# EDIT PATHS AND FILE NAMES HERE
inpath=/hpf/largeprojects/pmaass/Daniella/Katerina/CpG_ident
reffile=CISTR_bisulfite.fasta
#samplefiles=top_file_list.txt # list of full path to files and full file name for all samples from same strand. on file name per line
samplefiles=bottom_file_list.txt # list of full path to files and full file name for all samples from same strand. on file name per line
#strand=top # which strand are your samples sequenced on? options: "top" or "bottom"
strand=bottom # which strand are your samples sequenced on? options: "top" or "bottom"
mapped_reads_file=/hpf/largeprojects/pmaass/Sergio/19_10_PM001/analysis/HiC-Pro/human/Astrocytes/03_11_18_Astro/hic_results/data/Astro_human/Astro_human.allValidPairs
outpath=/hpf/largeprojects/pmaass/Daniella/Katerina/CpG_ident
outfile=CpG_test.csv
#######################################

## get cat commands for samples
#perl -p -i -e "s/\n/ /g" $inpath/$samplefiles

#for samples aligning to top (forward) strand
if [ "$strand" = "top" ]
then
    echo $strand
    #get one input multi-fasta file for alignment
    touch $outpath/aln_in.txt
    cat $inpath/$reffile >> $outpath/aln_in.txt
    while IFS= read -r file
            do
            echo $file
            cat $file >> $outpath/aln_in.txt
    done < "$inpath/$samplefiles"
fi
#for samples aligning to bottom (reverse) strand
if [ "$strand" = "bottom" ]
then
    echo $strand
    #reverse complement of the sequence
   revseq $inpath/$reffile -outseq $outpath/rev.fasta 
    #get one input multi-fasta file for alignment
    touch $outpath/aln_in.txt
    cat $outpath/rev.fasta >> $outpath/aln_in.txt
    while IFS= read -r file
            do
            echo $file
            cat $file >> $outpath/aln_in.txt
    done < "$inpath/$samplefiles"
fi
# align sequences with clustal w (default parameters)
clstcmd=$(echo "clustalw2 -INFILE=$outpath/aln_in.txt -OUTORDER=INPUT")
echo $clstcmd
eval $clstcmd

# identify CpG meth
echo "identifying CpG methylation ..."
python3 $cod/cpg_anno.py $outpath/aln_in.aln $outpath/$outfile $strand


echo "DONE :)"
#
#
#
#
#
##convert spaces to tabs in output file
#perl -p -i -e "s/ /\t/g" $outpath/$outfile
