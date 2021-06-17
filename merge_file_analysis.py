
import sys
import pandas as pd
from plotnine import *
import argparse

infile3 = "/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/Splicejunction/20_samples_NoNs_sj_st_end_20_60_binomial_simple_lowexpressedgene_rescuemostsjs_method6_multi_map_fix_final_merge.bed"

def bar_plot(input_dictionary,filename,x_axis_title,y_axis_title,bar_colour,axis_angle_90):#Creates barplot in ggplot style, from dictionary with keys as x axis and values as y axis
	histokeys = list(input_dictionary.keys())
	histovalues = []
	for genes in input_dictionary.values():
		if type(genes) == set or type(genes) == list:
			histovalues.append(len(genes))
		elif type(genes) == int or type(genes) == float:
			histovalues.append(genes)
		else:
			print("Type error: dictionary values must be integer, set or list")
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	pandadict = {x_axis_title:histokeys,y_axis_title:histovalues}
	histodata = pd.DataFrame(pandadict)
	t = ggplot(aes(x = x_axis_title, weight = y_axis_title), data = histodata) + geom_bar(fill = bar_colour,color='black')
	p = t + theme_bw() + labs(title="", x = x_axis_title,y = y_axis_title)
	if type(histokeys[-1]) == int:
		if max(histokeys) >= 100:
			if max(histokeys) >= 5000:
				maxlabel = int(max(histokeys)/1000)*1000
				breaksize = maxlabel/1000
				breaks1 = []
				for i in range(0,int(breaksize)):
					breaks1.append(i*1000)
				breaks1.append(maxlabel)
				#breaks1.append(max(histokeys))
				p = p + scale_x_continuous(breaks = breaks1)
			else:
				p = p + scale_x_continuous(breaks = [0,20,40,60,80,100])
	if axis_angle_90==True:
		p = p + theme(axis_text_x  = element_text(size = 14, angle = 90))
	p.save(filename)


class Bedline:
	def __init__(self,line):
		line1 = line.strip("\n").split("\t")
		self.chromosome = line1[0]
		transcriptinfo=line1[3].split(";")#transcript id is transcriptinfo[1]
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
			intron_start = int(self.start)+int(self.list_exon_starts[i])+int(self.list_exon_lengths[i])#Extract start of intron from .bed file
			intron_end = int(self.start)+int(self.list_exon_starts[i+1])#Extract end of intron from .bed file
			sjnumber = i#sj no is same as intron number
			list_sj_coordinates.append(self.chromosome+"_"+str(intron_start)+"_"+str(intron_end)+"_"+self.strand)
		return list_sj_coordinates
	def sequence_length(self):#Generate sequence length for each transcript
		length = 0
		for exon_length in self.list_exon_lengths:
			length = length + int(exon_length)
		return length

parser = argparse.ArgumentParser(description='Analyse .bed file')
parser.add_argument('-i', dest = 'bed', type = str, help = 'The bed input file')
parser.add_argument('-o', dest = 'output',type = str, help = 'Name of output prefix')
args = parser.parse_args()

startsdict = {}
endsdict = {}
genes = set()
transcripts = set()
gene_transcript_dict_merged = {}

long_sjs = set()
gene_exon_no_dict = {}
single_exon_transcripts = set()

with open(args.bed) as bedfile:
	for line in bedfile:
		bed = Bedline(line)
		genes.add(bed.gene)
		transcripts.add(bed.transcript)
		long_sjs.update(set(bed.sjcoordinates()))
		gene_transcript_dict_merged.setdefault(bed.gene,set()).add(bed.transcript)
		gene_exon_no_dict.setdefault(bed.gene,set()).add(bed.exon_no)
		if bed.exon_no == 1:
			single_exon_transcripts.add(bed.transcript)
		if bed.strand == "+":
			startsdict.setdefault(bed.chromosome,set()).add(bed.start)
			endsdict.setdefault(bed.chromosome,set()).add(bed.end)
		elif bed.strand == "-":
			startsdict.setdefault(bed.chromosome,set()).add(bed.end)
			endsdict.setdefault(bed.chromosome,set()).add(bed.start)
		else:
			print(f'Warning! {bed.gene} does not have strand and so cannot be used for TSS/ TES analysis')
			print(f'{bed.strand}')
			#print("Failed")
			#sys.exit()

#Get number of single exon genes
single_exon_genes = set()
for gene in gene_exon_no_dict.keys():
	if gene_exon_no_dict[gene] == {1}:
		single_exon_genes.add(gene)




#Count starts
startcounts = 0
endscounts = 0

for chrom in startsdict.keys():
	startcounts += len(startsdict[chrom])
	endscounts += len(endsdict[chrom])


#Histogram of transcript numbers
transcript_numbers = {}
for i in range(1,10):
	transcript_numbers[i] = 0
transcript_numbers["10+"] = 0
for gene in gene_transcript_dict_merged.keys():
	number = len(gene_transcript_dict_merged[gene])
	if number in transcript_numbers.keys():
		transcript_numbers[number] += 1
	else:
		if number >= 10:
			transcript_numbers["10+"] += 1



bar_plot(transcript_numbers,"/mnt/shared/scratch/mc42302/201903_RTD2/Figures/Num of transcripts per gene histogram" + args.output,"Number of transcripts per gene","Frequency","blue",False)



print(f'Stats for {args.output}')
print(f'# genes: {len(genes)}')
print(f'# transcripts: {len(transcripts)}')
print(f'# Unique splice junctions: {len(long_sjs)}')
print(f'# Single exon genes: {len(single_exon_genes)} ({(len(single_exon_genes)/len(genes)) * 100}%)')
print(f'# Multi exon genes: {len(genes) - len(single_exon_genes)}')
print(f'# Single transcript genes: {transcript_numbers[1]} ({(transcript_numbers[1]/len(genes)) * 100}%)')
print(f'# Single exon transcripts: {len(single_exon_transcripts)} ({(len(single_exon_transcripts)/len(transcripts)) * 100}%)')
print(f'# Multi exon transcripts: {len(transcripts) - len(single_exon_transcripts)}')

