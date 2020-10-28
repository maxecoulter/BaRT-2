

#!/bin/bash

#$ -cwd
#$ -j yes
#$-l hostname="n17*"
#$ -pe mpi 28

#Script for merging short read and long read based transriptomes. For BaRTv2.0 the short read transcriptome was in split chromosomes, so long read assembly also needed to be split. After this, transcriptomes are merged and filtered based on method described in BaRTv2.0 paper. FASTA file is generated for resulting transcriptome, and duplicate transcripts are removed

iso_file=$1
short_read=$2
version=$3
path=/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD #Path to merged transcriptome folder
genome=/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/barke_split.fasta
RTDfasta=/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT_2_${version}.fasta #Output fasta


#First convert low confidence genes to split chromosomes
echo "splitting chromosomes for low confidence"
conda activate ggplot_python
python convert_bed_to_split_chromosomes.py ${iso_file}_low_confidence.bed ${iso_file}_low_confidence_split.bed
conda deactivate

#Get list of files for merging: BaRT_Iso, BaRT_Illumina, low confidence Iso-Seq
rm /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelist.txt
ls ${iso_file}_final_merge100_split.bed | awk '{split($1,a,".");print $1"\tcapped\t1,1,1\tIso4";}'>> /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelist.txt
ls $short_read | awk '{split($1,a,".");print $1"\tno_cap\t2,1,2\tIllumina2";}'>> /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelist.txt
ls ${iso_file}_low_confidence_split.bed | awk '{split($1,a,".");print $1"\tcapped\t2,1,2\tIsolowconfidence";}'>> /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelist.txt

echo "Merging BaRT Iso, BaRT Illumina, and low confidence Iso"
conda activate oldpython #TAMA requires python 2
#NOw merge with short reads

python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_merge.py -f ${path}/filelist.txt -p ${path}/merged_BaRT2_w100_merge${version} -m 0 -a 100 -z 100 -d merge_dup
conda deactivate


#Now run filtering
echo "Filtering merged file"
conda activate ggplot_python
python short_long_RTD_filter.py -m /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT2_w100_merge${version}_merge.txt -b /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT2_w100_merge${version}.bed -gtf /mnt/shared/scratch/je42879/BaRTD_03July20_annotation_only/BaRT_July20_TS_transfix.gtf -monoexon T -o /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT_filtered${version}.bed
conda deactivate

echo "Running tama merge one last time!"
#Now run tama merge again
rm /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelistfiltered.txt
ls /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT_filtered${version}.bed | awk '{split($1,a,".");print $1"\tcapped\t1,1,1\tIsomergedfiltered";}'>> /mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/filelistfiltered.txt

conda activate oldpython
path=/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD
python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_merge.py -f ${path}/filelistfiltered.txt -p ${path}/BaRT_2_${version}_with_duplicates -m 0 -a 0 -z 0 -d merge_dup
conda deactivate



#Get stats
echo "Getting stats"
echo "Stats for BaRT Iso:"
conda activate ggplot_python
#First BaRT Iso:
python merge_file_analysis.py -i ${iso_file}_final_merge100_split.bed -o BaRT_Iso${version}

#Now BaRT Illumina:
#python merge_file_analysis.py -i $short_read -o BaRT_Illumina_2
conda deactivate



#Get fasta file
echo "Get fasta file for BaRT 2"

RTDfasta2=${RTDfasta%%.fasta}_withduplicates.fasta

conda activate bedtools

bedtools getfasta -fi $genome -bed ${path}/BaRT_2_${version}_with_duplicates.bed -split -name -s -fo $RTDfasta2 | fold -w 60

conda deactivate

echo "remove duplicates"

conda activate seqkit

cat $RTDfasta2 | seqkit rmdup -s -o $RTDfasta -d duplicated.fasta -D duplicated.detail.txt
conda deactivate



#Remove duplicates from .bed file as well

python filter_bed_duplicates.py -di /mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/transcriptomes/duplicated.detail.txt -b ${path}/BaRT_2_${version}_with_duplicates.bed -o ${path}/BaRT_2_${version}.bed -f $RTDfasta -fd $RTDfasta2

#Now BaRT 2
conda activate ggplot_python
python merge_file_analysis.py -i ${path}/BaRT_2_${version}.bed -o BaRT_2_${version}
conda deactivate

#Problem identified - fix fasta file so it matches .bed file
#python duplicate_fix.py -f $RTDfasta2 -b ${path}/BaRT_2_${version}.bed -o $RTDfasta



echo "Script complete!"