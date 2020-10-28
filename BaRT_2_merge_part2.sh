
#####
#!/bin/bash

#$ -cwd
#$ -j yes
#$-l hostname="n17*"


#Script for merging Iso-Seq based transcriptome after filtering using TAMA merge. This removes redundant transcripts by merging transcripts with same splice junctions and similar ends. 


prefix=$1 #The file prefix
path=/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples
path2=/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/Splicejunction


#Run final tama merge
echo "Running final tama merge"
conda activate oldpython

cd $path
python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_merge.py -f ${path}/file_list_filtered.txt -p ${path2}/${prefix}_final_merge100 -m 0 -a 100 -z 100 -d merge_dup 
conda deactivate


echo "Number of genes:"
cut -f4 ${path2}/${prefix}_final_merge.bed  | awk '{split($1, a, ";");print a[1]}'| sort | uniq | wc -l
echo "Number of transcripts:"
cut -f4 ${path2}/${prefix}_final_merge.bed  | awk '{split($1, a, ";");print a[2]}'| sort | uniq | wc -l






###statisitics


cut -f4 ${path2}/${prefix}_final_merge.bed | awk '{split($1, a, ";");print a[1]}'| sort | uniq | wc -l
cut -f4 ${path2}/${prefix}_final_merge.bed | awk '{split($1, a, ";");print a[2]}'| sort | uniq | wc -l

###statisitics for barley cv. Morex 1 ONT sample
## genes transcripts





