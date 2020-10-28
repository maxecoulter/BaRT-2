
#####Get CCS and FLNC reads using SMRT-tools 5.1
#!/bin/bash

#$ -cwd
#$ -j yes
#$-l hostname="n17*"


#Script for running TAMA merge on all TAMA collapse samples for BaRTv2.0. This is to create the initial .bed file for filtering. Script contains hard paths. 

#############run TAMA merge#################

###generate filelist for TAMA merge
cd /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples #This is the directory where sample folders are contained

rm filelist.txt filelist1.txt



conda activate oldpython #Need python 2 to run TAMA

ls */*new_flnc_collapsed.bed | awk '{split($1,a,".");print $1"\tcapped\t1,1,1\t"a[1];}'>> filelist.txt

cat filelist.txt | awk '{print $1"\t"$2"\t"$3"\t"$4}'> filelist1.txt

###run TAMA merge
###5' and 3' threshold 0 nt



python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_merge.py -f filelist1.txt -p merged -m 0 -a 0 -z 0 -d merge_dup 
conda deactivate


#Get gene and transcript numbers from bed file
cut -f4 merged.bed | awk '{split($1, a, ";");print a[1]}'| sort | uniq | wc -l
cut -f4 merged.bed | awk '{split($1, a, ";");print a[2]}'| sort | uniq | wc -l








