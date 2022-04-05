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
inpath=/hpf/largeprojects/pmaass/resources/human/GRCh38.p13_v32
gfffile=gencode.v32.annotation.gtf
#inpath=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
#gfffile=gencode.v32.annotation_subset.gtf
bin=1000000 # bin size in bp
outpath=/hpf/largeprojects/pmaass/Daniella/SickKids_scripts
outfile=hg38_p13_v32_annotation.txt
#######################################
#echo "python path"
#echo $PYTHONPATH
export PYTHONPATH="/hpf/largeprojects/pmaass/programs/gffutils_python_3.8.1"
export PYTHONPATH="/hpf/largeprojects/pmaass/programs/python_v3.8.1/gtfparse"
#export PYTHONPATH="/hpf/largeprojects/pmaass/programs/python_v3.8.1/bcbio_gff"
#echo "python path"
#echo $PYTHONPATH
echo "sum annotation per bin ..."
python3 $cod/cod_non_cod_anno_per_bin.py $inpath/$gfffile $bin $outpath/$outfile

echo "DONE :)"
