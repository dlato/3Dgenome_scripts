#!/bin/bash
#PBS -l nodes=1:ppn=40,mem=20gb,walltime=1:00:00 #specify memory requirements
#PBS -j eo
#PBS -N Trans_all_inters_df #name your job
#PBS -V
#PBS -m e #sends you an email when your job is done

module load python/3.8.1
#module load perl/5.28.0
#module load clustalw/2.1
#module load emboss/6.6.0

#!/bin/sh

#######################################
# DO NOT EDIT
cod=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
#######################################


#######################################
# EDIT PATHS AND FILE NAMES HERE
#inpath=/hpf/largeprojects/pmaass/resources/human/GRCh38.p13_v32
#gfffile=gencode.v32.annotation.gtf
inpath=/hpf/largeprojects/pmaass/Daniella/juicer_test
cooldumpfile=4dn.bsorted.chr21_22_only.pairs.cool.txt
res=1000000 # resolution size in bp
outpath=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
#######################################

python3 $cod/coolerdumpTo3Dflow.py $inpath/$cooldumpfile $outpath $res

echo "DONE :)"
