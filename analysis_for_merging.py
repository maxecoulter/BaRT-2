"""Obtain important information from merge report of merging short and long read RTDs. Create venn diagrams based on this. Getting an idea of how both datasets compare so that filtering can be decided."""

#from python_modules import TSS_functions as f
from upsetplot import from_memberships
from upsetplot import plot
from matplotlib import pyplot
pyplot.rc('font', size=8)
import pandas as pd
from plotnine import *
import matplotlib
from matplotlib import pyplot as plt
from matplotlib_venn import venn2,venn2_circles
from matplotlib_venn import venn3,venn3_circles

import statistics

mergefile = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/merged_BaRT2_w100_merge10_merge.txt"

class Bedline:
	def __init__(self,line):
		line1 = line.strip("\n").split("\t")
		self.chromosome = line1[0]
		transcriptinfo = line1[3].split(";")#transcript id is transcriptinfo[1]
		self.gene = transcriptinfo[0]
		self.transcript = transcriptinfo[1]
		self.strand = line1[5]
		self.start = int(line1[6])
		self.end = int(line1[7])
		self.exon_no = int(line1[9])
		self.list_exon_lengths = line1[10].split(",")
		self.list_exon_starts = line1[11].split(",")
	def sjcoordinates(self):#Create sj coordinate for each sj in transcript
		list_sj_coordinates = []
		for i in range(0,int(self.exon_no)-1):#Go through each intron in .bed file (0 - n)
			intron_start=int(self.start)+int(self.list_exon_starts[i])+int(self.list_exon_lengths[i])#Extract start of intron from .bed file
			intron_end=int(self.start)+int(self.list_exon_starts[i+1])#Extract end of intron from .bed file
			list_sj_coordinates.append(self.chromosome+"_"+str(intron_start)+"_"+str(intron_end)+"_"+self.strand)
		return list_sj_coordinates
	def sequence_length(self):#Generate sequence length for each transcript
		length = 0
		for exon_length in self.list_exon_lengths:
			length = length + int(exon_length)
		return length

#rc('font',**{'family':'serif','serif':['Times']})
font = {'family' : 'sans-serif', 'weight' : 'normal','size':16}

matplotlib.rc('font', **font)

def venn_diagram(filename,listofsubsets,labels):
	plt.figure(figsize=(7,7))
	v = venn2(listofsubsets, set_labels = labels, set_colors = ('purple', 'green'), alpha = 0.7)
	c = venn2_circles(listofsubsets, linestyle='-', linewidth=0.4, color="black")
	for text in v.set_labels:
		text.set_fontsize(16)
	for text in v.subset_labels:
		text.set_fontsize(15)
	plt.title("")
	plt.savefig(filename, dpi=1000)

def venn_diagram3(filename,listofsubsets,labels):
	plt.figure(figsize=(8,8))
	v = venn3(listofsubsets, set_labels = labels)
	c = venn3_circles(listofsubsets, linestyle='-', linewidth=1, color="black")
	plt.title("")
	plt.savefig(filename, dpi=1000)





transcripts = [set(),set(),set()] #long, short, long low confidence
genes = [set(),set(),set()]
sjs = [set(),set(),set()]




gene_exon_numbers_short = {}
gene_exon_numbers_long = {}
gene_exon_numbers_long_low = {}

short_lengths = []
long_lengths = []

short_exon_lengths = []
long_exon_lengths = []
long_exonlow_lengths = []

short_exon_numbers = []
long_exon_numbers = []

gene_transcripts = {}
transcript_objects = {}

for line in open(mergefile):
	bed = Bedline(line)#This class is for bedfiles, but the input is slightly different. In this case transcript = bed.gene and source_transcript = bed.transcript
	transcript = bed.gene
	gene_transcripts.setdefault(transcript.split(".")[0],set()).add(bed)
	transcript_objects[bed.gene] = bed
	if bed.transcript.split("_")[0] == "Illumina2":
		short_exon_numbers.append(bed.exon_no)
		short_lengths.append(bed.end - bed.start)
		short_exon_lengths.append(bed.sequence_length())
		transcripts[1].add(bed.gene)
		genes[1].add(bed.gene.split(".")[0])
		gene_exon_numbers_short.setdefault(bed.gene.split(".")[0],set()).add(bed.exon_no)
		for sj in bed.sjcoordinates():
			sjs[1].add(sj)
	elif bed.transcript.split("_")[0] == "Iso4":
		long_exon_numbers.append(bed.exon_no)
		long_lengths.append(bed.end - bed.start)
		long_exon_lengths.append(bed.sequence_length())
		transcripts[0].add(bed.gene)
		genes[0].add(bed.gene.split(".")[0])
		gene_exon_numbers_long.setdefault(bed.gene.split(".")[0],set()).add(bed.exon_no)
		for sj in bed.sjcoordinates():
			sjs[0].add(sj)
	elif bed.transcript.split("_")[0] == "Isolowconfidence":
		transcripts[2].add(bed.gene)
		genes[2].add(bed.gene.split(".")[0])
		long_exonlow_lengths.append(bed.sequence_length())
		gene_exon_numbers_long_low.setdefault(bed.gene.split(".")[0],set()).add(bed.exon_no)
		for sj in bed.sjcoordinates():
			sjs[2].add(sj)
	else:
		print("Error!")
		print(bed.transcript.split("_")[1])
		break




single_exon = [set()] * 3
multi_exon = [set()] * 3


for gene in gene_exon_numbers_short.keys():
	exon_numbers = gene_exon_numbers_short[gene]
	if exon_numbers == {1}:
		single_exon[1].add(gene)
	else:
		multi_exon[1].add(gene)

for gene in gene_exon_numbers_long.keys():
	exon_numbers = gene_exon_numbers_long[gene]
	if exon_numbers == {1}:
		single_exon[0].add(gene)
	else:
		multi_exon[0].add(gene)

for gene in gene_exon_numbers_long_low.keys():
	exon_numbers = gene_exon_numbers_long_low[gene]
	if exon_numbers == {1}:
		single_exon[2].add(gene)
	else:
		multi_exon[2].add(gene)

#Use f.venn_diagram for two way comparisons
#First compare Iso-Seq and Illumina datasets
venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_long_transcripts",transcripts[:-1],("Iso-Seq","Illumina"))

venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_long_genes",genes[:-1],("Iso-Seq","Illumina"))

venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_long_sjs",sjs[:-1],("Iso-Seq","Illumina"))
#Next compare Illumina and lc Iso_Seq datasets

venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_lowc_transcripts",transcripts[1:],("Illumina","Iso-Seq no TSS"))

venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_lowc_genes",genes[1:],("Illumina","Iso-Seq no TSS"))

venn_diagram("/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/merged_RTD/short_lowc_sjs",sjs[1:],("Illumina","Iso-Seq no TSS"))


#Venn diagrams show low overlap between transcripts. Why is this? Added 18/08/21

illumina_categories = {"illumina_only":0,"ends":0,"no_phase":0}#categories of illumina transcripts wihtout iso-seq support, Illumina_only: trancsripts that are from genes without Iso-seq support, ends: Transcripts that have the same set of sjs as an iso-seq transcript - the difference is therefore due to ends, no_phase: There trancsripts are from genes with Iso-seq transcripts but it is out of phase from all of them
isoseq_only_transcripts = set([transcript for transcript in transcripts[0] if transcript not in transcripts[1]])
illumina_only_transcripts = set([transcript for transcript in transcripts[1] if transcript not in transcripts[0]])


for transcript in illumina_only_transcripts:#This transcript is not shared between Illumina and Iso-seq data bases. WHy not?
	gene = transcript.split(".")[0]
	if gene not in genes[0]:#Gene is not an Iso-seq gene
		illumina_categories["illumina_only"] += 1
	else:
		gene_t = gene_transcripts[gene]
		#Find a transcript with the same SJ coordinates
		transcript_object = transcript_objects[transcript]
		ends = False
		for gene_transcript in gene_t:
			if gene_transcript.gene in isoseq_only_transcripts:#If transcript being compared is only Iso-seq based
				if set(transcript_object.sjcoordinates()) == set(gene_transcript.sjcoordinates()):#If transcript being compared has the same set of sjs
					illumina_categories["ends"] += 1
					ends = True
					break #This transript must difer from an Iso-seq just by ends
		if not ends:
			illumina_categories["no_phase"] += 1

assert(sum(illumina_categories.values()) == len(illumina_only_transcripts))#Each illumina only transcript should be in one category

print(f"Different categories of Illumina only transcripts: {illumina_categories}")
			


			










"""
venn_diagram3("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/short_long_low_genes_single",single_exon,("Iso-Seq single exon genes","Illumina single exon genes","Iso-Seq low confidence single exon genes"))


venn_diagram3("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/short_long_genes_multi",multi_exon,("Iso-Seq multi exon genes","Illumina multi exon genes","Iso-Seq low confidence multi exon genes"))



venn_diagram3("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/short_long_low_transcripts",transcripts,("Iso-Seq transcripts","Illumina transcripts","Iso-Seq low confidence transcripts"))

venn_diagram3("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/short_long_low_genes",genes,("Iso-Seq genes","Illumina genes","Iso-Seq low confidence genes"))

venn_diagram3("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/short_long_low_sjs",sjs,("Iso-Seq splice junctions","Illumina splice junctions","Iso-Seq low confidence splice junctions"))

#Make upset and length plots:


#transcripts = [set(),set(),set()] #long, short, long low confidence

for sets,name in zip([transcripts,genes,sjs],["Transcripts","Genes","Splice junctions"]):
	upset1 = from_memberships([['Iso-Seq LC'],['Iso-Seq HC'],['Iso-Seq LC', 'Iso-Seq HC'],['Illumina'],['Illumina', 'Iso-Seq LC'],['Illumina', 'Iso-Seq HC'],['Illumina', 'Iso-Seq HC', 'Iso-Seq LC'],],data=[len(sets[2] - sets[1] - sets[0]), len(sets[0] - sets[1] - sets[2]), len(sets[0] & sets[2] - sets[1]), len(sets[1] - sets[2] - sets[0]), len(sets[1] & sets[2] - sets[0]), len(sets[1] & sets[0] - sets[2]), len(sets[1] & sets[0] & sets[2])])
	plot(upset1,show_counts='%d')
	pyplot.savefig(f"/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/example_upset_{name}")

length_type = ["Iso-Seq HC"] * len(long_exon_lengths) + ["Iso-Seq LC"] * len(long_exonlow_lengths) + ["Illumina"] * len(short_exon_lengths)
length = long_exon_lengths + long_exonlow_lengths + short_exon_lengths
pandadict = {"Source":length_type,"Length":length}
length_plot_data = pd.DataFrame(pandadict)

t = ggplot(aes(x='Source', y='Length'), data = length_plot_data) + geom_boxplot() + xlab('') + ylab('Sequence length') + theme_bw() + theme(text=element_text(size=16))
t.save("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/lengths_boxplot")


statistics.mean(short_lengths)
statistics.mean(long_lengths)

statistics.mean(short_exon_lengths)
statistics.mean(long_exon_lengths)

plt.hist(short_lengths,bins=50)
plt.xlabel('length (nt)')
plt.ylabel('frequency')
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Illumina transcript lengths including intron.png", dpi=150)
plt.clf()

plt.hist(short_exon_lengths,bins=50)
plt.xlabel('length (nt)')
plt.ylabel('frequency')
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Illumina transcript lengths.png", dpi=150)
plt.clf()






plt.hist(long_lengths,bins=50)
plt.xlabel('length (nt)')
plt.ylabel('frequency')
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Iso3 transcript lengths including introns.png", dpi=150)
plt.clf()

plt.hist(long_exon_lengths,bins=50)
plt.xlabel('length (nt)')
plt.ylabel('frequency')
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Iso3 transcript lengths.png", dpi=150)
plt.clf()

plt.hist(short_exon_numbers,bins=50)
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Illumina transcript exon numbers.png", dpi=150)
plt.clf()

plt.hist(long_exon_numbers,bins=50)
plt.savefig("/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/histogram BaRT Iso3 transcript exon numbers.png", dpi=150)
plt.clf()

"""