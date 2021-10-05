#!/bin/bash
#PBS -l nodes=1:ppn=40,mem=20gb,walltime=1:00:00 #specify memory requirements
#PBS -j eo
#PBS -N Trans_all_inters_df #name your job
#PBS -V
#PBS -m e #sends you an email when your job is done

module load python/3.8.1
module load perl/5.28.0

#!/bin/sh

#######################################
# DO NOT EDIT
cod=/hpf/largeprojects/pmaass/programs/
read_len=100
chrom_len_file=/hpf/largeprojects/pmaass/resources/human/GRCh38.p10_v26/hg38.chrom.sizes
chrList=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts/human_chr_list.txt
#######################################


#######################################
# EDIT PATHS AND FILE NAMES HERE
mapped_reads_file=/hpf/largeprojects/pmaass/Sergio/19_10_PM001/analysis/HiC-Pro/human/Astrocytes/03_11_18_Astro/hic_results/data/Astro_human/Astro_human.allValidPairs
outpath=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
outfile=chrom_mapped_coverage.txt
#######################################

#read in each line of chromosome file
#touch $outpath/num_mapped_per_chrom.txt
#while IFS= read -r file
#        do
#        echo $file >> $outpath/num_mapped_per_chrom.txt
## grab total number of mapped reads per chrom, save to file
#        grep -P "$file\t" $mapped_reads_file | wc -l >> $outpath/num_mapped_per_chrom.txt
#done < "$chrList"

python3 $cod/gonosome_ident.py $read_len $chrom_len_file $outpath/num_mapped_per_chrom.txt > $outpath/$outfile

#convert spaces to tabs in output file
perl -p -i -e "s/ /\t/g" $outpath/$outfile
