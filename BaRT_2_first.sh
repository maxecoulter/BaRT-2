
#!/bin/bash

#$ -cwd
#$ -j yes
#$-l hostname="n17*"
#$ -pe mpi 28 
#$ -t 1-20

#Script for processing raw iso-seq data fo BaRTv2.0. Script will find raw subreads, extract name and create a seperate folder for each sample. Each sample will then be worked on as an array job (for BaRTv2.0 this was 20 samples). In brief, CCS reads are created, mapped to the Barke genome, and used to create a basic unfiltered transcriptome annotation using TAMA collapse.

#Script contains hard paths. 


#Find all isoseq subread files, create folders for results if needed
#find . -name "*m*.subreads.bam" > readfile

#find . -name "m54203_191120_152842.subreads.bam" > readfile
#find . -name "m54203_191121_114918.subreads.bam" >> readfile
find . -wholename "./*/*/*/*/*/*/*/m*.subreads.bam" > extrafolders
for i in `cat extrafolders`
do
FName=${i##*/};#<bam file name>.bam
FName2=${FName%%.*};#bam file name
mkdir $FName2
done

#for i in `cat readfile`
#do
#FName=${i##*/};#<bam file name>.bam
#FName2=${FName%%.*};#bam file name
#mkdir $FName2
#done


#for i in `cat readfile` 
#do

INFILE=`awk "NR==$SGE_TASK_ID" readfile`
#INFILE=./m54203_190507_063222/m54203_190507_063222.subreads.bam
FName=${INFILE##*/};#<bam file name>.bam
FName1=${INFILE##./};#bam input (*/<bam file name>.bam
FName2=${FName%%.*};#bam file name
bamdirectory=/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/$FName1
#mv $INFILE ${FName2}/
#FName=${i%/*};#./<bam file name>
#FName1=${i##*/};#bam input (<bam file name>.bam
#FName2=${FName1%%.*};#bam file name 

#mv $i ${FName2}/ #Move .bam file to 

cd ${FName2}/



#mkdir ${FName2}

echo "working on ccs for" $FName2
#conda activate isoseqold
conda activate isoseq3
###step1: Circular Consensus Sequence calling ccs 4.0.0 (commit v4.0.0)
ccs --min-rq 0.9 -j 28 $bamdirectory ${FName2}.ccs.bam
#Use older version for now
#ccs --noPolish --minLength=300 --minPasses=1 --minZScore=-999 --maxDropFraction=0.8 --minPredictedAccuracy=0.8 --minSnr=4 $bamdirectory ${FName2}new.ccs.bam

echo "working on lima for" $FName2
###step 2 Primer removal and demultiplexing lima 1.10.0 (commit v1.10.0)
lima ${FName2}.ccs.bam /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/teloprimer_lima.fasta ${FName2}new.fl.bam --isoseq --peek-guess
#singularity run --app pbtranscript /mnt/apps/singularity/images/smrt-tools/5.1.0/smrt-tools.img classify  --flnc ${FName2}_isoseq_flnc.fasta --nfl ${FName2}_isoseq_nfl.fasta -p /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/teloprimer.fasta --cpus 15 ${FName2}.ccs.bam ${FName2}_draft.fasta

echo "working on refine for" $FName2
###step 3: refine: 1)Trimming of poly(A) tails; 2) Rapid concatmer identification and removal isoseq3 3.2.2 (commit v3.2.2)
isoseq3 refine ${FName2}new.fl.F0_5p--R0_3p.bam /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/teloprimer_lima.fasta ${FName2}new.flnc.bam
conda deactivate

####convert bam files to fasta files
conda activate samtoolsoldSSL
echo "Converting bam to fasta..." $FName2
samtools fasta ${FName2}new.flnc.bam > ${FName2}new.flnc.fasta
conda deactivate

conda activate biopython
####remove poly-A tail remnants for Telo prime
python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_go/sequence_cleanup/tama_flnc_polya_cleanup.py -f ${FName2}new.flnc.fasta -p ${FName2}new.flnc1
conda deactivate


conda activate minimap2
echo "Mapping" $FName2 "reads to Barke genome"
####minimap2 2.17-r941
###output sam format
###--splice-flank=no no assumptions of the nucleotide after GT and AG
#minimap2 -ax splice -uf -G 6000 --secondary=no  ${FName2}.flnc.fasta > ${FName2}_isoseq_flnc.fasta.sam 2> ${FName2}_isoseq_flnc.fasta.sam.log
minimap2 -ax splice:hq -uf -G 6000 --secondary=no /mnt/shared/projects/barley/201903_RTD2/BarkeAssembly/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta ${FName2}new.flnc1.fa > ${FName2}new_flnc.fasta.sam 2> ${FName2}new_flnc.fasta.sam.log
conda deactivate




##sort the sam file: samtools 1.9 Using htslib 1.9
#convert sam to bam
conda activate samtoolsoldSSL
echo "convert SAM to BAM"
samtools view \
-b \
-T /mnt/shared/projects/barley/201903_RTD2/BarkeAssembly/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta \
${FName2}new_flnc.fasta.sam \
> ${FName2}new_flnc.fasta.unsorted.bam


echo "sorting sam file"

samtools sort -o ${FName2}new_flnc.fasta.sorted.bam ${FName2}new_flnc.fasta.unsorted.bam
echo "filtering reads with no mapping"
##filtered out reads with no mapping -F 4
samtools view -h -F 4 ${FName2}new_flnc.fasta.sorted.bam > ${FName2}new_flnc.fasta.sorted.filtered.sam

conda deactivate

###Run TAMA collapse version 2019.09.13
###5' and 3' threshold 10nt
###-m Exon/Splice junction threshold
echo "running tama for sample" $FName2
#conda activate oldpython
conda activate oldpython
#conda activate
python /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/tama/tama_collapse.py \
-s ${FName2}new_flnc.fasta.sorted.filtered.sam \
-f /mnt/shared/projects/barley/201903_RTD2/BarkeAssembly/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta \
-p ${FName2}new_flnc_collapsed -d merge_dup -x capped -m 0 -a 0 -z 0 -sj sj_priority -lde 30 -sjt 30
conda deactivate

cd ..

#done

###########use files generated by CCS, classify and TAMA
##/mnt/shared/scratch/rz41242/2017Isoseq/sam

awk 'BEGIN{print "fileName\t#ZMW\t#CCS"}'> ReadStat

for i in `cat readfile` 
do
#i=./m54203_190605_093609/m54203_190605_093609.subreads.bam
#i=./m54203_190508_154014/m54203_190508_154014.subreads.bam
FName=${i%/*};
FName1=${i##*/};
FName2=${FName1%%.*};

####statistics up to CCS
cat /mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/${FName2}/ccs_report.txt | awk -v fileN="$FName2" '{if(NR==1){zwm=$5}; if(NR==2){zwm1=$6};}END{print fileN"\t"zwm"\t"zwm1>> "ReadStat"}'

done

####FLNC
awk 'BEGIN{print "#FL"}'> FLStat

for i in `cat readfile` 
do 

FName=${i%/*};
FName1=${i##*/};
FName2=${FName1%%.*};

cat ./${FName2}/${FName2}new.fl.lima.summary | awk '{ if(NR==2){ print $7 >> "FLStat"}}'

done

#####mapped reads, collapsed transcripts/genes

awk 'BEGIN{print "#mappedLocations"}'> mappedLocationStat
awk 'BEGIN{print "#UniquemappedFLNCs"}'> UniquemappedFLNCs
awk 'BEGIN{print "#transcripts"}'> TranscriptStat
awk 'BEGIN{print "#genes"}'> GeneStat

for i in `cat readfile` 
do 

FName=${i%/*};
FName1=${i##*/};
FName2=${FName1%%.*};
conda activate samtoolsoldSSL
samtools view -c -F 4 ${FName2}/${FName2}new_flnc.fasta.sorted.filtered.sam >> mappedLocationStat
conda deactivate
cat ${FName2}/${FName2}new_flnc_collapsed_read.txt | cut -f1 | sed '1d' | sort | uniq | wc -l >> UniquemappedFLNCs

cat ${FName2}/${FName2}new_flnc_collapsed_trans_report.txt | cut -f1 | sed '1d'| wc -l >> TranscriptStat

cat ${FName2}/${FName2}new_flnc_collapsed_trans_report.txt | cut -f1 | sed '1d'| tr "." "\t" | cut -f1 | sort | uniq | wc -l >> GeneStat

done
######
paste ReadStat FLStat UniquemappedFLNCs mappedLocationStat  TranscriptStat GeneStat > ReadStat1.txt







