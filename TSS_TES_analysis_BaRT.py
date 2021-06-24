"""Identify motifs and check if enriched"""
import re
import random
import matplotlib
import pandas as pd
import numpy
import Bio
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from importlib import reload


#import scipy
#import statsmodels
#import ggplot

from plotnine import *

outpath = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/TSS_TES_validation/"
bart1_gtf = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT.1.0_Unsplit_unpadded.gtf"
#bart2_gtf = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT_2_10.gtf"
bart2_gtf = "/mnt/shared/projects/jhi/barley/201903_RTD2/results/transcriptomes/BaRT2v18/BaRT2v18.gtf" #whole chromosomes
bartiso_gtf = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT_Iso4_split.gtf"
bart_illumina_gtf = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/transcriptomes/BaRT_Illumina3.gtf"


barke_genome = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/barke_split.fasta"
morex_genome = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/150831_barley_pseudomolecules.fasta"
barke_genome_whole = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/Benchmarking/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta"

bart1TSS = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bart1_TSS.txt"
bart2TSS = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bart2_TSS.txt"

bart1TES = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bart1_TES.txt"
bart2TES = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bart2_TES.txt"

bartisoTSS = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bartiso_TSS.txt"
bartisoTES = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bartiso_TES.txt"

bartilluminaTSS = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bartillumina_TSS.txt"
bartilluminaTES = "/mnt/shared/scratch/mcoulter/mc42302/201903_RTD2/scripts/bartillumina_TES.txt"

def get_motif_data(motif,TSS,chrom_index,strand_index,TSS_index,genome):
	print(genome.keys())
	print(TSS)
	histo = [0] * 1101
	histokeys = list(range(-550,551))
	for n,line in enumerate(open(TSS)):
		if not n:
			continue
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index]
		strand = line1[strand_index]
		TSS_position = int(line1[TSS_index])
		if strand == "+":
			sequence = genome[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome[chromosome][TSS_position - 550:TSS_position + 551].upper())
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo[i] += 1
	histo_dict = dict(zip(histokeys,histo))
	histo_dict2 = {a:z for a,z in histo_dict.items() if a >= -500 and a <= 500}
	return histo_dict2

def generate_motif_plot_allcompare(motif,genome1_comparisons,genome2_comparisons,data_type,chrom_index,strand_index,TSS_index,genome1,genome2,filename,colours,other_genome_comparisons,other_genome):
	"""Create plot with enriched sites, not enriched sites and random sites"""
	outputs = []
	totals = []
	if other_genome_comparisons:
		for TSS in other_genome_comparisons:
			histo_dict = get_motif_data(motif,TSS,chrom_index,strand_index,TSS_index,other_genome)
			outputs.append(histo_dict)
			totals.append(len(open(TSS).readlines()))
	for TSS in genome1_comparisons:
		histo_dict = get_motif_data(motif,TSS,chrom_index,strand_index,TSS_index,genome1)
		outputs.append(histo_dict)
		totals.append(len(open(TSS).readlines()))
	for TSS in genome2_comparisons:
		histo_dict = get_motif_data(motif,TSS,chrom_index,strand_index,TSS_index,genome2)
		outputs.append(histo_dict)
		totals.append(len(open(TSS).readlines()))
	#####
	#NOw do control run
	histo_control = [0] * 1101
	strands = ["+","-"]
	for i in range(0,len(open(genome1_comparisons[0]).readlines())):
		#First gent random chromosome
		random.shuffle(strands)
		strand = strands[0]
		chromosomes = list(genome1.keys())
		random.shuffle(chromosomes)
		chromosome = chromosomes[0]
		TSS_position = random.randint(550,len(genome1[chromosome]) - 551)
		if strand == "+":
			sequence = genome1[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome1[chromosome][TSS_position - 550:TSS_position + 551].upper())
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_control[i] += 1
	histokeys = list(range(-550,551))
	histo_dict_control = dict(zip(histokeys,histo_control))
	histo_dict_control2 = {a:z for a,z in histo_dict_control.items() if a >= -500 and a <= 500}
	outputs.append(histo_dict_control2)
	totals.append(len(open(genome1_comparisons[0]).readlines()))
	line_plot_3(outputs,filename,f"position relative to predicted {data_type} (nt)","Instances",False,totals,colours)
	line_plot_3(outputs,filename + "_frequency",f"position relative to predicted {data_type} (nt)",f"Motif instances per {data_type}",True,totals,colours)


def generate_motif_plot_bothcompare(motif,enriched_TSS,enriched_TSS2,chrom_index,strand_index,TSS_index,chrom_index2,strand_index2,TSS_index2,genome1,genome2,filename,colour,colour2):
	"""Create plot with enriched sites, not enriched sites and random sites"""
	histo = [0] * 1101
	histokeys = list(range(-550,551))
	for n,line in enumerate(open(enriched_TSS)):
		if not n:
			continue
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index]
		strand = line1[strand_index]
		TSS_position = int(line1[TSS_index])
		if strand == "+":
			sequence = genome1[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome1[chromosome][TSS_position - 550:TSS_position + 551].upper())
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo[i] += 1
	histo_dict = dict(zip(histokeys,histo))
	histo_dict2 = {a:z for a,z in histo_dict.items() if a >= -500 and a <= 500}
	########
	histo_no_enriched = [0] * 1101
	for n,line in enumerate(open(enriched_TSS2)):#gene, chromosome,TSS,strand
		if not n:
			continue
		line1 = line.rstrip("\n").split("\t")
		chromosome = line1[chrom_index2]
		strand = line1[strand_index2]
		TSS_position = int(line1[TSS_index2])
		if strand == "+":
			sequence = genome2[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome2[chromosome][TSS_position - 550:TSS_position + 551].upper())
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		#now find if motif is in region +/- 500
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_no_enriched[i] += 1
	histo_dict_notenriched_full = dict(zip(histokeys,histo_no_enriched))
	histo_dict_notenriched = {a:z for a,z in histo_dict_notenriched_full.items() if a >= -500 and a <= 500}
	#####
	
	#NOw do control run
	histo_control = [0] * 1101
	strands = ["+","-"]
	for i in range(0,max([len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines())])):
		#First gent random chromosome
		random.shuffle(strands)
		strand = strands[0]
		chromosomes = list(genome1.keys())
		random.shuffle(chromosomes)
		chromosome = chromosomes[0]
		TSS_position = random.randint(550,len(genome1[chromosome]) - 551)
		if strand == "+":
			sequence = genome1[chromosome][TSS_position - 550:TSS_position + 551].upper()
		elif strand == "-":
			sequence_reverse = Seq(genome1[chromosome][TSS_position - 550:TSS_position + 551].upper())
			sequence = str(sequence_reverse.reverse_complement())
		else:
			print("Input error")
		motif_positions = motif_pos(motif,sequence)
		if motif_positions:
			#assert len(motif_position) == 2
			for motif_position in motif_positions:
				for i in range(motif_position[0],motif_position[1]):
					histo_control[i] += 1
	histo_dict_control = dict(zip(histokeys,histo_control))
	histo_dict_control2 = {a:z for a,z in histo_dict_control.items() if a >= -500 and a <= 500}
	line_plot_2(histo_dict2,histo_dict_notenriched,histo_dict_control2,filename,"position relative to predicted TSS (nt)","Instances",colour,colour2,False,len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines()))
	line_plot_2(histo_dict2,histo_dict_notenriched,histo_dict_control2,filename + "_frequency","position relative to predicted TSS (nt)","Motif instances per TSS",colour,colour2,True,len(open(enriched_TSS).readlines()),len(open(enriched_TSS2).readlines()))


def get_bart_ends(input_name,outname):
	"""Produces outfiles with TSS and TES information that can be used for downstream analysis"""
	TSS_output = open(f"{outname}_TSS.txt","w")
	TES_output = open(f"{outname}_TES.txt","w")
	TSS_set = set()
	TES_set = set()
	for n,line in enumerate(open(input_name)):
		line1 = line.rstrip("\n").split("\t")
		if line1[2] == "transcript":
			if line1[6] == "+":
				TSS = str(int(line1[3]) - 1)
				TES = line1[4]
			elif line1[6] == "-":
				TES = str(int(line1[3]) - 1)
				TSS = line1[4]
			else:
				print(f"strand error in line {n}")#many transcripts in bart1 do not have strands. Ignore these
				continue
			chromosome = line1[0]
			strand = line1[6]
			if (chromosome,strand,TSS) not in TSS_set:
				TSS_output.write(chromosome +"\t" + strand +"\t" + TSS + "\n")
				TSS_set.add((chromosome,strand,TSS))
			if (chromosome,strand,TES) not in TES_set:
				TES_output.write(chromosome +"\t" + strand +"\t" + TES + "\n")
				TES_set.add((chromosome,strand,TES))
	TSS_output.close()
	TES_output.close()

def line_plot_2(input_dictionary1,input_dictionary2,input_dictionary3,filename,x_axis_title,y_axis_title,bar_colour,bar_colour2,is_frequency,total1,total2):
	histokeys = list(input_dictionary1.keys())
	histovalues1 = list(input_dictionary1.values())
	histovalues2 = list(input_dictionary2.values())
	histovalues3 = list(input_dictionary3.values())
	if is_frequency:
		histovalues1 = [i/total1 for i in histovalues1]
		histovalues2 = [i/total2 for i in histovalues2]
		histovalues3 = [i/total2 for i in histovalues3]
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	pandadict = {x_axis_title:histokeys,y_axis_title:histovalues1,"second":histovalues2,"random":histovalues3}
	histodata = pd.DataFrame(pandadict)
	t = ggplot(aes(x = x_axis_title), data = histodata) + geom_line(aes(y = y_axis_title), color = bar_colour) + geom_line(aes(y = "second"), color = bar_colour2) + geom_line(aes(y = "random"), color = "grey")
	p = t + theme_bw() + labs(title="", x = x_axis_title,y = y_axis_title) + theme(axis_text_x=element_text(size=15), axis_text_y=element_text(size=15),axis_title=element_text(size=15))
	p.save(filename, dpi=1000)

def line_plot_3(inputs,filename,x_axis_title,y_axis_title,is_frequency,totals,colours):
	name = 0
	pandadict = {}
	for i in range(0,len(inputs)):
		input_dict = inputs[i]
		total = totals[i]
		pandadict[x_axis_title] = list(input_dict.keys())
		name += 1
		if is_frequency:
			pandadict[str(name)] = [i/total for i in list(input_dict.values())]
		else:
			pandadict[str(name)] = list(input_dict.values())
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	histodata = pd.DataFrame(pandadict)
	t = ggplot(aes(x = x_axis_title), data = histodata)
	i = 0
	for name in pandadict.keys():
		if name != x_axis_title:
			t += geom_line(aes(y = name), color = colours[i])
			i += 1
	p = t + theme_bw() + labs(title="", x = x_axis_title,y = y_axis_title)
	p.save(filename, dpi=1000)

def motif_pos(motif,sequence):
	p = re.compile(motif)
	if not p.search(sequence):
		return []
	else:
		coordinates = []
		iterator = p.finditer(sequence)
		for match in iterator:
			coordinates.append(match.span())
		return coordinates

def parse_genome(genome_input):
	genome = open(genome_input).read()
	#Split genome into chromosome dictionary with key = chr id and item chromosome sequence
	chromosomes = genome.split(">")
	del chromosomes[0]
	chromosomedict = {}
	for chromosome in chromosomes:
		splitchrom = chromosome.split("\n")
		chromid = splitchrom[0]
		chromsequence = "".join(splitchrom[1:])
		chromosomedict[chromid] = chromsequence #Genome sequence now in dictionary with chromosomes as key
	return chromosomedict

#First get end information for both bart 1 and bart2 transcriptomes
#this may not work - many Bart 2 ends are illumina based

#First parse genomes:
barke_chromosomes = parse_genome(barke_genome)
morex_chromosomes = parse_genome(morex_genome)
barke_chromosomes_whole = parse_genome(barke_genome_whole)

get_bart_ends(bart1_gtf,"bart1")
get_bart_ends(bart2_gtf,"bart2")
get_bart_ends(bartiso_gtf,"bartiso")
get_bart_ends(bart_illumina_gtf,"bartillumina")


print(f"Number of BaRT 1 transcript 5' ends: {len(open(bart1TSS).readlines())}")
print(f"Number of BaRT 2 transcript 5' ends: {len(open(bart2TSS).readlines())}")
print(f"Number of BaRT 1 transcript 3' ends: {len(open(bart1TES).readlines())}")
print(f"Number of BaRT 2 transcript 3' ends: {len(open(bart2TES).readlines())}")


#Compare all transcriptome TSS and TES

genome1_comparisons = [bartisoTSS,bartilluminaTSS]
genome2_comparisons = [bart1TSS]
other_genome_comparisons = [bart2TSS]
colours = ["blue","green","purple","red","grey"]

#TSS first:
tatamotif = "TATA[AT]A[AT]"
generate_motif_plot_allcompare(tatamotif,genome1_comparisons,genome2_comparisons,"TSS",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "TATA TSS BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)

Inr = "[CT]TCA[ATCG]T[CT][CT]"
generate_motif_plot_allcompare(Inr,genome1_comparisons,genome2_comparisons,"TSS",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "Inr TSS BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)

ypatch3 = "C[CT]TC[CT][CT]CC[CT]C"

generate_motif_plot_allcompare(ypatch3,genome1_comparisons,genome2_comparisons,"TSS",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "ypatch3 TSS BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)

kozakmotif5 = "[ATCG][AC][AG][AC][ATCG]ATGGCG" 
generate_motif_plot_allcompare(kozakmotif5,genome1_comparisons,genome2_comparisons,"TSS",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "kozak TSS BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)

#TES
genome1_comparisons = [bartisoTES,bartilluminaTES]
genome2_comparisons = [bart1TES]
other_genome_comparisons = [bart2TES]
colours = ["blue","green","purple","red","grey"]


CFlm = "TGTA"
generate_motif_plot_allcompare(CFlm,genome1_comparisons,genome2_comparisons,"TES",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "CFlm TES BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)

PAS = "AATAAA"
generate_motif_plot_allcompare(PAS,genome1_comparisons,genome2_comparisons,"TES",0,1,2,barke_chromosomes,morex_chromosomes,outpath + "PAS TES BaRT1 BaRT2 all compared",colours,other_genome_comparisons,barke_chromosomes_whole)




#TSS first:
tatamotif = "TATA[AT]A[AT]"
generate_motif_plot_bothcompare(tatamotif,bart2TSS,bart1TSS,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "TATA TSS BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1

Inr = "[CT]TCA[ATCG]T[CT][CT]"
generate_motif_plot_bothcompare(Inr,bart2TSS,bart1TSS,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "Inr TSS BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1


ypatch3 = "C[CT]TC[CT][CT]CC[CT]C"

generate_motif_plot_bothcompare(ypatch3,bart2TSS,bart1TSS,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "ypatch3 TSS BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1

kozakmotif5 = "[ATCG][AC][AG][AC][ATCG]ATGGCG" 
generate_motif_plot_bothcompare(kozakmotif5,bart2TSS,bart1TSS,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "kozakmotif TSS BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1

#Now TES:

CFlm = "TGTA"
generate_motif_plot_bothcompare(CFlm,bart2TES,bart1TES,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "CFlm TES BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1

PAS = "AATAAA"

generate_motif_plot_bothcompare(PAS,bart2TES,bart1TES,0,1,2,0,1,2,barke_chromosomes_whole,morex_chromosomes,outpath + "PAS TES BaRT1 BaRT2 compared","blue","red")#blue is bart2, red is bart1