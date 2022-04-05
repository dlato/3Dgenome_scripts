#!/bin/bash
#PBS -l nodes=1:ppn=40,mem=20gb,walltime=1:00:00 #specify memory requirements
#PBS -j eo
#PBS -N Trans_all_inters_df #name your job
#PBS -V
#PBS -m e #sends you an email when your job is done

module load python/3.8.1

#!/bin/sh

#######################################
# DO NOT EDIT
cod=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
#######################################


#######################################
# EDIT PATHS AND FILE NAMES HERE
input_bed_path=/hpf/largeprojects/pmaass/3D-flow/3D-flow2/scripts/src
input_bed=human.1000000.bed
outpath=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
outfile=all_trans_interactions_1Mb.txt
#outfile=TEST_all_trans_interactions_1Mb.txt
#######################################

python3 $cod/all_trans_interactions.py $input_bed_path/$input_bed > $outpath/$outfile
