
"""Splice-centric based Iso-Seq Filtering for Transcripts in Transcriptome (SIFT)
Author: Max Coulter, based on algorithms created by Runxuan Zhang and John Brown and Max Coulter
For creating RTD based on Iso-Seq and Illumina data, using output from TAMA
"""


import time
import sys
import os
import logging
import re
import math
from operator import add
from itertools import repeat
import matplotlib
import pandas as pd
import numpy
from plotnine import *

#For venn diagram
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import statsmodels
from statsmodels import stats
from statsmodels.stats import multitest
import argparse


import scipy.stats as stats
from scipy.stats import fisher_exact as fishers_exact_test


start1 = time.time()

conversiondict = {"chr1H":205415724,"chr2H":295953730,"chr3H":267034671,"chr4H":273224149,"chr5H":209770093,"chr6H":246027230,"chr7H":316571039}#For converting split Barke file to whole
readsupportcutoff = 5 #Min read support for splice junction
min_overhang = 10 #min overhang for splice junction


splicejunctionfilter = True
transcript_start_end_filter = True
lowexpressedgene_rescue_most_sjs = False #As above, but just rescues the longest transcripts with the most sjs
keep_proximal = False #For high confidence genes, keep the transcript with most proximal start(despite low read support)
low_expression_cutoff = "NA" #Expression of gene at or below which gene is considered lowly expressed. These genes are treated slighly differently in TSS filtering, with higher p value required for first TSS and possibility for TSS to be supported by 2 read starts within the startwindow
cluster_predicted_ends = False

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
			intron_start=int(self.start)+int(self.list_exon_starts[i])+int(self.list_exon_lengths[i])#Extract start of intron from .bed file
			intron_end=int(self.start)+int(self.list_exon_starts[i+1])#Extract end of intron from .bed file
			list_sj_coordinates.append(self.chromosome+"_"+str(intron_start)+"_"+str(intron_end)+"_"+self.strand)
		return list_sj_coordinates
	def sequence_length(self):#Generate sequence length for each transcript
		length = 0
		for exon_length in self.list_exon_lengths:
			length = length + int(exon_length)
		return length

class Genes:
	Instances = {}
	reads_dict = {}
	reads_Ns_dict = {}
	transcripts_dict = {}
	def __init__(self,bed,reads,is_Ns):
		if is_Ns:
			self.gene = bed.gene
			self.reads_Ns_dict.setdefault(self.gene,set()).update(reads)
		else:
			self.gene = bed.gene
			self.chromosome = bed.chromosome
			self.strand = bed.strand
			self.reads_dict.setdefault(self.gene,set()).update(reads)
			self.transcripts_dict.setdefault(self.gene,set()).add(bed.transcript)
			self.Instances[self.gene] = self
	def get(self,gene,transcripts):
		self.gene = gene
		self.transcripts = self.transcripts_dict[gene]
		self.chromosome = self.Instances[gene].chromosome
		self.strand = self.Instances[gene].strand
		self.reads = self.reads_dict[gene]
		if self.strand =="+":
			self.startindex = 1
			self.endindex = 2
		elif self.strand =="-":
			self.startindex = 2
			self.endindex = 1
		else:
			print("Strand error!")
		if gene in self.reads_Ns_dict.keys():
			self.reads_Ns = self.reads_Ns_dict[gene]
			self.expression = len(self.reads) + len(self.reads_Ns)
		else:
			self.expression = len(self.reads)
		self.five_start_transcript = furthest_start2(self.transcripts,transcripts) #Find most proximal start trasncript
		self.three_end_transcript = furthest_end2(self.transcripts,transcripts)
		self.five_start = transcripts.get(self.five_start_transcript).start
		self.three_end = transcripts.get(self.three_end_transcript).end
		if self.strand == "+":
			self.total_length = self.three_end - self.five_start
		else:
			self.total_length = self.five_start - self.three_end
		return self

class Reads:
	Instances = {}
	sjs_dict = {}
	transcripts_dict = {}
	def __init__(self,line,type_file):
		line1 = line.strip("\n").split("\t")
		if type_file == "multi":
			self.read = line1[0]
			self.sjs_dict.setdefault(self.read,set()).add(line1[4])
			self.chromosome = line1[4].split("_")[0]
			self.start = int(line1[6])
			self.end = int(line1[7])
			self.strand = line1[4].split("_")[-1]
		elif type_file == "single":
			self.read = line1[0]
			self.chromosome = line1[2]
			self.start = int(line1[3])
			self.end = int(line1[4])
			self.strand = line1[5]
			self.sjs_dict[self.read] = set()
		else:
			print("Specify file type correctly for Read class please")
			sys.exit()
		if self.strand == "+":
				self.start_s = self.start
				self.end_s = self.end
		elif self.strand == "-":
			self.start_s = self.end
			self.end_s = self.start
		else:
			print("strand error in Reads class!")
			sys.exit()
		self.Instances[self.read] = self
	def get(self,read):
		self.read = read
		self.sjs = self.sjs_dict[read]
		self.chromosome = self.Instances[read].chromosome
		self.start = self.Instances[read].start
		self.end = self.Instances[read].end
		self.strand = self.Instances[read].strand
		return self

class Splice_Junctions:
	Instances = {}
	reads_dict = {}
	readkeys_dict = {}
	errors_dict = {}
	def __init__(self,line):
		line1 = line.strip("\n").split("\t")
		self.sj = line1[4]#
		self.read = line1[0]
		self.sjno = line1[2]
		self.error = line1[3]
		self.readkey = self.read + "_" + self.sjno
		self.sequence = line1[5].upper()
		self.reads_dict.setdefault(self.sj,set()).add(line1[0])
		self.readkeys_dict.setdefault(self.sj,{})[self.readkey] = self.error
		self.errors_dict.setdefault(self.sj,set()).add(self.error)
		self.Instances[self.sj] = self
	def get(self,sj):
		self.sj = sj
		self.reads = self.reads_dict[sj]
		self.readkeys = self.readkeys_dict[sj].keys()
		self.errors = self.errors_dict[sj]
		self.sequence = self.Instances[sj].sequence
		split = sj.split("_")
		self.chromosome = split[0]
		self.intronstart = int(split[1])
		self.intronend = int(split[2])
		self.strand = split[-1]
		return self
	def get_error(self,sj,readkey):
		return self.readkeys_dict[sj][readkey]

class Transcripts:
	Instances = {}
	def __init__(self,bed,reads):
		self.transcript = bed.transcript
		self.gene = bed.gene
		self.chromosome = bed.chromosome
		self.strand = bed.strand
		start = bed.start
		end = bed.end
		if self.strand == "+":
			self.start = start
			self.end = end
		elif self.strand =="-":
			self.start = end
			self.end = start
		self.sjs = set(bed.sjcoordinates())
		self.reads = reads
		self.Instances[self.transcript] = self
	def get(self,transcript):
		self.transcript = transcript
		self.gene = self.Instances[transcript].gene
		self.chromosome = self.Instances[transcript].chromosome
		self.strand = self.Instances[transcript].strand
		self.start = self.Instances[transcript].start
		self.end = self.Instances[transcript].end
		if self.strand == "+":
			self.total_length = self.end - self.start
		elif self.strand == "-":
			self.total_length = self.start - self.end
		self.sjs = self.Instances[transcript].sjs
		self.reads = self.Instances[transcript].reads
		self.read = list(self.reads)[0]
		return self



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

def bionomial_method(n,h,t): #n is expression, h is number of reads supporting TSS being tested and #t is total number of TSS
	try:
		return math.factorial(n)/(math.factorial(n-h)*math.factorial(h))*(1/t)**h*((t-1)/t)**(n-h)
	except OverflowError:
		try:
			return math.factorial(n)//int(math.factorial(n-h)*math.factorial(h))*(1/t)**h*((t-1)/t)**(n-h)
		except OverflowError:
			return 0

def combine_starts(maxstartdict,startwindow,p_dict):#Combine starts with start support in starts dictionary
	newdict = {}
	lost_starts = set()
	for start in maxstartdict.keys():
		if start in lost_starts:
			continue
		inwindowstarts = set()
		readsupport = maxstartdict[start]
		for i in range(start - startwindow,start + startwindow):
			if i != start:
				if i in maxstartdict.keys():
					inwindowstarts.add(i)
					readsupport += maxstartdict[i]
					lost_starts.add(i)
		newstart = int((start + sum(inwindowstarts))/(1 + len(inwindowstarts)))
		p = p_dict[start][1]
		p_dict[newstart] = readsupport, p
		newdict[newstart] = readsupport
	return newdict

def create_end_dict(gene_object,whichend,transcripts):
	startsdict = {} #All possible starts in gene
	if whichend == 5:
		if gene_object.strand == "+":
			for i in range(0,gene_object.total_length):
				startsdict[gene_object.five_start + i] = 0
		else:
			for i in range(0,gene_object.total_length):
				startsdict[gene_object.five_start - i] = 0
	elif whichend == 3:
		if gene_object.strand == "+":
			for i in range(0,gene_object.total_length):
				startsdict[gene_object.three_end - i] = 0
		else:
			for i in range(0,gene_object.total_length):
				startsdict[gene_object.three_end + i] = 0
	totalreadcounts = 0
	for transcript in gene_object.transcripts:#Go through every transcript in gene
		transcript_obj = transcripts.get(transcript)#Transcript(transcript)
		if whichend == 5:
			transcript_start = transcript_obj.start
			#transcript_start = int(trans_start_end_dict[transcript][gene_object.startindex])
		elif whichend == 3:
			transcript_start = transcript_obj.end
			#transcript_start = int(trans_start_end_dict[transcript][gene_object.endindex])
		else:
			print("Fatal error")
			sys.exit()
		if transcript_start in startsdict.keys():
			try:
				readcount = len(transcript_obj.reads)
				#readcount = len(trans_read_dict[transcript])
				startsdict[transcript_start] += readcount
				totalreadcounts += readcount
			except KeyError:
				print("Error: No expression for transcript "+ transcript)
				#try:
					#readcount = len(trans_read_single_exon[transcript])
					#startsdict[transcript_start] += readcount
					#totalreadcounts += readcount
				#except KeyError:
					#print("Error: No expression for transcript "+ transcript)
		else:
			print("Error! transcript start not in histo")
			f"{transcript_start}"
	return startsdict

def binomial_simple(self,mid_confidence_start_transcripts,trueTSS_dict,trueTSS_dict2,trueTES_dict2,degredation_dict,degraded_genes,startwindow,endwindow,transcripts,lowconfidence_genes):
	startsdict = create_end_dict(self,5,transcripts) #All possible starts in gene
	enddict = create_end_dict(self,3,transcripts) #All possible starts in gene
	if self.expression == "NA":
		print("Error!")
		return "Error!"
	elif self.expression == 1:
		lowconfidence_genes.add(self.gene)
		return
	elif low_expression_cutoff != "NA" and self.expression < low_expression_cutoff:
		mid_gene = False
		#Here a low expressed gene does not have multiple starts supporting one location, yet may still have 2 reads with starts within a window. in this case, find the locations where this is the case, and take only the location closest to the proximal TSS
		for transcript in self.transcripts:
			transcript_object = transcripts.get(transcript)#Transcript(transcript)
			if five_prime_support_hunter2(transcript_object.start,startwindow,startsdict,False) and three_prime_support_hunter3(transcript_object,enddict,endwindow,False):
			#if five_prime_support_hunter(transcript,startwindow,self.chromosome,self.strand,False,transcript_start,2):#Find if this transcript is supported by normal means
				mid_confidence_start_transcripts.add(transcript)
				mid_gene = True
		if not mid_gene:
			lowconfidence_genes.add(self.gene)
		return
	else:
		p_threshold = 0.01
		n = int(self.expression)
		p_dict = {}
		p_dict_end = {}
		total_maxstartposition_dict = {} #Position is key, read support is value
		#bar_plot(startsdict, gene + " TSS","Start position nt","Read support","blue",True)
		total_maxendposition_dict = {}
		totalstarts = sum(startsdict.values())
		possible_starts = [i for i in startsdict.keys() if startsdict[i] > 0]
		possible_ends = [i for i in enddict.keys() if enddict[i] > 0]
		support_start_dict, start_ordered_support = order_end_dict(possible_starts,startsdict)
		support_end_dict, end_ordered_support = order_end_dict(possible_ends,enddict)
		for start_support in start_ordered_support:
			test_end = True
			for start in support_start_dict[start_support]:
				if not significance_binomial(start,startsdict,possible_starts,n,self.gene,lowconfidence_genes,total_maxstartposition_dict,p_dict,p_threshold):
					test_end = False
					break
			if not test_end:
				break
		for end_support in end_ordered_support:
			test_end = True
			for end in support_end_dict[end_support]:
				if not significance_binomial(end,enddict,possible_ends,n,self.gene,lowconfidence_genes,total_maxendposition_dict,p_dict_end,p_threshold):
					test_end = False
					break
			if not test_end:
				break
		if len(total_maxstartposition_dict.keys()) == 0 or len(total_maxendposition_dict.keys()) == 0:
			if low_expression_cutoff == "NA":
				#Here a low expressed gene does not have multiple starts supporting one location, yet may still have 2 reads with starts within a window. in this case, find the locations where this is the case, and take only the location closest to the proximal TSS
				mid_gene = False
				for transcript in self.transcripts:
					transcript_object = transcripts.get(transcript)#Transcript(transcript)
					if five_prime_support_hunter2(transcript_object.start,startwindow,startsdict,False) and three_prime_support_hunter3(transcript_object,enddict,endwindow,False):
						mid_confidence_start_transcripts.add(transcript)
						mid_gene = True
				if not mid_gene:
					lowconfidence_genes.add(self.gene)
			else:
				lowconfidence_genes.add(self.gene)
			return
		if keep_proximal:
			five_start = self.five_start
			if five_start not in total_maxstartposition_dict.keys():
				total_maxstartposition_dict[five_start] = startsdict[five_start]
				p_dict[five_start] = startsdict[five_start], "NA"
				
		#Now combine all the starts/ends iteratively within window
		if cluster_predicted_ends:
			total_maxstartposition_dict_combined = combine_starts(total_maxstartposition_dict,startwindow,p_dict)
			total_maxendposition_dict_combined = combine_starts(total_maxendposition_dict,endwindow,p_dict_end)
		else:
			total_maxstartposition_dict_combined = total_maxstartposition_dict
			total_maxendposition_dict_combined = total_maxendposition_dict
		for start in total_maxstartposition_dict_combined.keys():
			trueTSS_dict2.setdefault(self.gene,dict())[start] = p_dict[start]
			trueTSS_dict.setdefault(self.chromosome,set()).add(start)
		for end in total_maxendposition_dict_combined.keys():
			trueTES_dict2.setdefault(self.gene,dict())[end] = p_dict_end[end]
		if self.expression >= 100:
			degredation_dict[self.gene] = (sum(startsdict.values())/totalstarts)*100
			if degredation_dict[self.gene] > 90:
				degraded_genes.add(self.gene)

def filter_bed(bedinput,genes,transcripts,fulltranscriptlist,fullgenelist,filteredtranscriptlist,filteredgenelist,gene_rescue_set,transcript_rescue_set,polyareaddict,mid_genes,mid_confidence_start_transcripts,lowconfidence_genes,trueTSS_dict2,trueTES_dict2,all_good_sjs,path,outputfile_prefix,startwindow,endwindow):
	mainoutput = open(path + outputfile_prefix + "_merged_filtered.bed","w")
	transcriptsremoved = open(path + outputfile_prefix + "_removed.bed","w")
	lowconfidence_good = open(path + outputfile_prefix + "_low_confidence.bed","w")
	transcript_info_file = open(path + outputfile_prefix + "_transcript_info","w")
	transcript_info_file.write("Transcript\tPolyA Error\tlow expressed\trescued\thigh confidence\t5' support\t3' support\tsplice junctions supported\n")
	gene_transcript_dict_binomial = {}
	for n,line in enumerate(open(bedinput)):
		#print(str(n))
		bed = Bedline(line)
		fulltranscriptlist.add(bed.transcript)
		fullgenelist.add(bed.gene)
		transcript_info_file.write(bed.transcript + "\t")
		if bed.gene not in genes.Instances.keys():
			print("Gene " + bed.gene + " not in database!")
			transcriptsremoved.write(line)
			continue
		try:
			transcript_object = transcripts.get(bed.transcript)
		except KeyError:
			print("Transcript " + bed.transcript + " not in database!")
			transcriptsremoved.write(line)
			continue
		#Transcript must pass parameters before it is considered high confidence
		transcript_support = {"no 3 prime polyA":False,"high confidence":False,"5 prime support":False,"3 prime support":False,"sj support":False}
		transcript_info = {"polyA error": False,"low expressed":False,"rescued":False}
		#Check to see if transcript has polyA error:
		if transcript_object.read not in polyareaddict.keys():
			transcript_support["no 3 prime polyA"] = True
			transcript_info_file.write("False\t")
		else:
			transcript_info_file.write("True\t")
		#Check to see if transcript has high confident starts and ends
		if bed.gene in mid_genes:
			if bed.transcript in mid_confidence_start_transcripts:
				transcript_info["low expressed"] = True
				#transcript_info_file.write("True\t")
				transcript_support["high confidence"] = True
				transcript_support["5 prime support"] = True
				transcript_support["3 prime support"] = True
			else:
				pass
		elif bed.gene in lowconfidence_genes:
			transcript_info["low expressed"] = True
			if bed.transcript in transcript_rescue_set: #IS it rescued?
				transcript_support["high confidence"] = True
				transcript_info["rescued"] = True
			if set(bed.sjcoordinates()) < all_good_sjs.keys() and transcript_support["no 3 prime polyA"]:#Check if all sjs in hc sj list
				lowconfidence_good.write(line) #These are transcripts which may get support from Illumina transcripts later
		elif bed.gene in trueTSS_dict2.keys() and bed.gene in trueTES_dict2.keys():
			#transcript_info_file.write("False\t")
			if five_prime_support_hunter2(transcript_object.start,int(startwindow/2),trueTSS_dict2[bed.gene],True):
				transcript_support["5 prime support"] = True
			if three_prime_support_hunter3(transcript_object,trueTES_dict2[bed.gene],int(endwindow/2),True):
				transcript_support["3 prime support"] = True
			if transcript_support["3 prime support"] and transcript_support["5 prime support"]:
				transcript_support["high confidence"] = True
		else:
			pass
		#Now check splicejunctions
		if not splicejunctionfilter:
			transcript_support["sj support"] = True
		else:
			if set(bed.sjcoordinates()) < all_good_sjs.keys():#Check if all sjs in hc sj list
				transcript_support["sj support"] = True
		#Now see if transcript passes or fails
		if all(transcript_support.values()):
			mainoutput.write(line)
			filteredtranscriptlist.add(bed.transcript)
			filteredgenelist.add(bed.gene)
			gene_transcript_dict_binomial.setdefault(bed.gene,set()).add(bed.transcript)
		else:
			if transcript_info["rescued"]:
				mainoutput.write(line)
				filteredtranscriptlist.add(bed.transcript)
				filteredgenelist.add(bed.gene)
				gene_transcript_dict_binomial.setdefault(bed.gene,set()).add(bed.transcript)
			else:
				if bed.gene not in gene_transcript_dict_binomial.keys():
					gene_transcript_dict_binomial[bed.gene] = set()
				transcriptsremoved.write(line)
		#print(str(transcript_support))
		#print(str(transcript_info))
		transcript_info_file.write(str(transcript_info["low expressed"]) + "\t" + str(transcript_info["rescued"]) + "\t" + str(transcript_support["high confidence"]) + "\t" + str(transcript_support["5 prime support"]) + "\t" + str(transcript_support["3 prime support"]) + "\t" + str(transcript_support["sj support"]) + "\n")
	mainoutput.close()
	transcriptsremoved.close()
	lowconfidence_good.close()
	transcript_info_file.close()
	return gene_transcript_dict_binomial

def filter_sjs(splice_junctions,chromosomedict,hammingthreshold):
	all_good_sjs = {}
	list_of_sj_fails = {}
	for sj_coordinate in splice_junctions.Instances.keys():
		sj = splice_junctions.get(sj_coordinate)
		if sj.strand == "+":
			if sj.sequence == "GTAG" or sj.sequence == "GCAG" or sj.sequence == "GTAT":
				pass
			else:
				list_of_sj_fails.setdefault(sj_coordinate,{"non_canonical":False,"error":False,"RT switch":False})["non_canonical"] = True
		elif sj.strand == "-":
			if sj.sequence == "CTAC" or sj.sequence == "CTGC" or sj.sequence == "ATAC":
				pass
			else:
				list_of_sj_fails.setdefault(sj_coordinate,{"non_canonical":False,"error":False,"RT switch":False})["non_canonical"] = True
		else:
			sys.exit("Strand error")
		#Find SJs with zero mismatches with 10nt around the SJs
		#First generate dictionary of key to error information
		#Next use sj_coordinate_based_dict to see if sj has support from other reads
		#If yes, keep read
		#Then in tama merge table find out if sj has alternative support. If yes keep read if no throw it out
		for key in sj.readkeys:
			errorprofile = splice_junctions.get_error(sj_coordinate,key)
			error = sj_error(errorprofile)
			if not error:
				break
		if error:
			list_of_sj_fails.setdefault(sj_coordinate,{"non_canonical":False,"error":False,"RT switch":False})["error"] = True
		###get the sequences near the SJs for template switching analysis (edit distance/hamming distance)
		###calculate the similarities between end of exon and end of next intron(field 3 and 5), or beginning of intron and the beginning of next exon, there has to >8 matches between two sequences to count as RT switching/repeat
		if template_switch(sj,chromosomedict,hammingthreshold):
			list_of_sj_fails.setdefault(sj_coordinate,{"non_canonical":False,"error":False,"RT switch":False})["RT switch"] = True
		#if sj not in error list, add to hc list
		if not sj_coordinate in list_of_sj_fails.keys():
			all_good_sjs[sj_coordinate] = sj.sequence
	return all_good_sjs, list_of_sj_fails

def five_prime_support_hunter2(start,startwindow,startsdict,highexpressed):
	try:
		if startsdict[start] > 1:
			return True
	except KeyError:
		pass
	except TypeError:
		pass
	for i in range(start - startwindow,start + startwindow + 1):
		if i in startsdict.keys():
			if not highexpressed:
				if startsdict[i] >= 1 and i != start:
					return True
			else:
				return True
	return False


#New 02/06/20

def furthest_end2(gene_transcripts,transcripts):#Finds transcript with most proximal 3' start
	firstend = -1
	for t in gene_transcripts:
		transcript = transcripts.get(t)
		if firstend == -1:
			firstend = transcript.end
			furthest_end_transcript = t
		elif firstend < transcript.end and transcript.strand == "+":
			firstend = transcript.end
			furthest_end_transcript = t
		elif firstend > transcript.end and transcript.strand == "-":
			firstend = transcript.end
			furthest_end_transcript = t
		else:
			pass
	return furthest_end_transcript

#New 3/6/20
def furthest_start2(gene_transcripts,transcripts):#gene_class_dict,transcript_class_dict):
	firststart = -1
	for t in gene_transcripts:
		transcript = transcripts.get(t)
		if firststart == -1:
			firststart = transcript.start
			furthest_start_transcript = t
		elif firststart > transcript.start and transcript.strand == "+":
			firststart = transcript.start
			furthest_start_transcript = t
		elif firststart < transcript.start and transcript.strand == "-":
			firststart = transcript.start
			furthest_start_transcript = t
		else:
			pass
	return furthest_start_transcript

def generate_t_read_dict(reads):#read_dict,readchromdict,start_end_dict):
	t_read_dict = {}
	for r in reads.Instances.keys():
		read = reads.get(r)
		t_read_dict.setdefault(read.chromosome + "_" + str(read.start) + "_" + str(read.end) + "_" + "_".join(sorted(list(read.sjs))),set()).add(r)
	return t_read_dict

def hamming_distance(s1, s2):
	assert len(s1) == len(s2)
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def histogram_from_list(list1,no_breaks,filename,x_axis_title,colour): #Takes list of floats, creates histogram. Way better then numpy.histogram!
	maximum = max(list1)
	range1 = maximum
	breaksize = range1/no_breaks
	histo = {}
	histo[0] = 0
	complete = 0
	while complete <= maximum:
		complete += breaksize
		histo[complete] = 0
	for key in histo.keys():
		next1 = key + breaksize
		for item in list1:
			if item >= key and item < next1:
				histo[key] +=1
			else:
				pass
	bar_plot(histo,filename,x_axis_title,"frequency",colour,False)

def inv(t):
	if t < 0:
		n = abs(t)
	elif t > 0:
		n =- t
	elif t == 0:
		n = t
	else:
		print("inverse error!")
	return n

def in_window(first_start,second_start,windowsize):
	difference = first_start - second_start
	if difference < 0:
		difference = inv(difference)
	return difference <= windowsize

def is_proximal(starttest,total_maxstartposition_dict,strand):
	if strand == "+":
		start = min(total_maxstartposition_dict.keys())
		return starttest < start
	else:
		start = max(total_maxstartposition_dict.keys())
		return starttest > start

def order_end_dict(possible_ends,enddict):
	support_end = {}
	for end in possible_ends:
		support_end.setdefault(enddict[end],set()).add(end)
	return support_end, sorted(list(support_end.keys()),reverse = True)

def parse_bed_input(bedinput,t_read_dict,reads,all_matched_reads,transcript_no_read_support,matched_coordinates,notmatched_coordinates,Ns):
	if not Ns:
		bedinput1 = open(bedinput)
	else:
		bedinput1 = bedinput
	for line in bedinput1:
		#print(line)
		bed = Bedline(line) #create object bed of class Bedline
		t_readkey = bed.chromosome + "_" + str(bed.start + 1) + "_" + str(bed.end + 1) + "_" + "_".join(sorted(bed.sjcoordinates()))
		if bed.transcript == "G68.264":
			print(line)
			print(t_readkey)
			print(str(t_readkey in t_read_dict.keys()))
			if "m54203_180927_185924/51511488/ccs" in reads.Instances.keys():
				read = reads.get("m54203_180927_185924/51511488/ccs")
				print(str(read.start))
				print(str(read.end))
				print(str(read.sjs))
			else:
				print("Error! Read m54203_180927_185924/51511488/ccs not in read Instances")
			#sys.exit()
		try:
			read_list = t_read_dict[t_readkey]#Extracted readlist is set
			transcripts = Transcripts(bed,read_list)
			genes = Genes(bed,read_list,Ns)
			all_matched_reads.update(read_list)
			matched_coordinates.add(t_readkey)
		except KeyError:
			transcript_no_read_support.add(bed.transcript)
			print("Error - transcript coordinates not matching any read coordinates")
			print(bed.transcript)
			print(t_readkey)
			print(str(bed.sjcoordinates()))
			notmatched_coordinates.add(t_readkey)
			continue
	return genes, transcripts

def parse_genome(genome_input):
	genome = open(genome_input).read().replace("\n","")
	#Split genome into chromosome dictionary with key = chr id and item chromosome sequence
	chromosomes = genome.split(">")
	del chromosomes[0]
	chromosomedict = {}
	for chromosome in chromosomes:
		chromid=chromosome[0:5]
		chromsequence=chromosome[5:]
		chromosomedict[chromid]=chromsequence #Genome sequence now in dictionary with chromosomes as key
	return chromosomedict

def parse_polyA(input_file,polyathreshold):
	#Dictionary of reads with possible false ends, read id is key with list as items, list[0] is strand, list[1] is percentage A
	polyareaddict = {}
	for line in open(input_file):
		if line.startswith("cluster"):
			continue
		else:
			line1 = line.split("\t")
			if float(line1[3]) >= polyathreshold:
				polyareaddict[line1[0]] = [line1[2],line1[3]]
			else:
				pass
	return polyareaddict

def parse_short_read_STAR(input1,conversiondict,readsupportcutoff,min_overhang):
	short_read_sj = set()
	for line in open(input1):
		line1 = line.strip("\n").split("\t")
		#Convert split chromosomes to whole chromosomes
		chromosome = line1[0].split("_")
		chrom = chromosome[0]
		if chromosome[-1] == "2":
			start = conversiondict[chrom] + int(line1[1])
			end = conversiondict[chrom] + int(line1[2])
		elif chromosome[-1] == "1":
			start = int(line1[1])
			end = int(line1[2])
		else:#Chromosome Un
			start = int(line1[1])
			end = int(line1[2])
		strand = line1[3]
		intronmotif = line1[4]
		read_support_single = int(line1[6])
		read_support_multi = int(line1[7])
		max_overhang = int(line1[8])
		if max_overhang < min_overhang:#Filter sjs with max overhang smaller than threshold
			continue
		if read_support_single <= readsupportcutoff and read_support_multi <= readsupportcutoff:#Filter sjs with low read support
			continue
		if intronmotif == 0:#Filter non canonical sj
			continue
		if strand == "1":
			short_read_sj.add(chrom + "_" + str(start-1) + "_" + str(end) + "_" + "+")
		elif strand == "2":
			short_read_sj.add(chrom + "_" + str(start - 1) + "_" + str(end) + "_" + "-")
		else:
			print("Read error!")
	return short_read_sj

def parse_sj_table(sj_input):
	starttime = time.time()
	for line in open(sj_input):
		reads = Reads(line,"multi")
		splice_junctions = Splice_Junctions(line)
	print(str(time.time() - starttime))
	return reads, splice_junctions

def rescue_low_expressed(lowconfidence_genes,genes,transcripts,all_good_sjs,short_read_sj):
	gene_rescue_set = set()
	transcript_rescue_set = set()
	if lowexpressedgene_rescue_most_sjs:
		for gene_name in lowconfidence_genes:
			gene = genes.get(gene_name,transcripts)
			sj_count = 0
			maxlen = 0
			maxtranscript = ""
			for transcript_name in gene.transcripts:
				transcript = transcripts.get(transcript_name)
				if len(transcript.sjs) != 0:#Some genes have single exon fragments which cannot be rescued by sjs
					if transcript.sjs < all_good_sjs.keys() and transcript.sjs < short_read_sj and len(transcript.sjs) > sj_count and transcript.total_length > maxlen:#If a splice junction is not supported, transcript cannot be rescued
						sj_count = len(transcript.sjs)
						maxlen = transcript.total_length
						maxtranscript = transcript_name
			if maxtranscript == "":
				pass
			else:
				gene_rescue_set.add(gene_name)
				transcript_rescue_set.add(maxtranscript)
	return gene_rescue_set, transcript_rescue_set

def scatter_plot(inputx,inputy,filename,x_axis_title,y_axis_title,bar_colour,axis_angle_90):#inputs must be lists of integers
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	pandadict = {x_axis_title:inputx,y_axis_title:inputy}
	scatterdata=pd.DataFrame(pandadict)
	t = ggplot(aes(x=x_axis_title, y=y_axis_title), scatterdata) + geom_point()
	p = t + theme_bw() + labs(title="", x=x_axis_title,y=y_axis_title)
	p.save(filename)

def short_sj_rescue(splice_junctions,all_good_sjs,short_read_sj,list_of_sj_fails):
	sj_rescue_set = set()
	for splice_junction in splice_junctions.Instances.keys():
		sj = splice_junctions.get(splice_junction)
		if splice_junction not in all_good_sjs.keys():
			if splice_junction in short_read_sj and not list_of_sj_fails[splice_junction]["RT switch"]:
				all_good_sjs[splice_junction] = sj.sequence
				sj_rescue_set.add(splice_junction)
	return all_good_sjs, sj_rescue_set

def significance_binomial(end,enddict,possible_ends,n,gene,lowconfidence_genes,total_maxposition_dict,p_dict,p_threshold):
	endsupport = enddict[end]
	totalends = len(possible_ends)
	expected = n / totalends
	if endsupport > expected and totalends > 1:
		try:
			p = bionomial_method(n, endsupport, totalends)
		except ZeroDivisionError:
			print(gene)
			print("maxstart:" + str(endsupport))
			print("expected: " + str(expected))
			print("n: " + str(n))
			return
		#p_test = True
	#elif len(possible_ends) == 1 and endsupport > 1:
		#p  = 0
		#p_test = False
	else:
		return False
	if p < p_threshold:
		total_maxposition_dict[end] = endsupport
		#if p_test:
		p_dict[end] = endsupport, p
		return True
	else:
		return False

def single_exon_parser(single_exon_input):#read_class_dict,single_exon_pos):#,single_exon_read,,chrom_start_read_dict,chrom_end_read_dict):
	for line in open(single_exon_input):
		reads = Reads(line,"single")
	return reads

def sj_error(errorprofile):
	"""Find characters that are not _"""
	errorset = set(errorprofile.split(">")[0][21:31])
	errorset.update(set(errorprofile.split(">")[1][0:10]))
	return not {"_"} == errorset

def template_switch(sj,chromosomedict,hammingthreshold):
	#Ends
	exonendsequence = chromosomedict[sj.chromosome][sj.intronstart - 8:sj.intronstart].upper()
	intronendsequence = chromosomedict[sj.chromosome][sj.intronend - 7:sj.intronend + 1].upper()
	endhamming = hamming_distance(exonendsequence,intronendsequence)
	#Starts
	exonstartsequence = chromosomedict[sj.chromosome][sj.intronend + 1:sj.intronend + 9].upper()
	intronstartsequence = chromosomedict[sj.chromosome][sj.intronstart:sj.intronstart + 8].upper()
	starthamming = hamming_distance(exonstartsequence,intronstartsequence)
	return endhamming <= hammingthreshold or starthamming <= hammingthreshold

def three_prime_support_hunter3(transcript,enddict,endwindow,highexpressed):#Finds support for TES in a given window
	if not highexpressed:
		if enddict[transcript.end] > 1:
			return True
		#if len(transcript.readcoverage_end) > 1:
			#return True
	for e in range(transcript.end - endwindow,transcript.end + endwindow + 1):
		if e in enddict.keys():
			if not highexpressed:
				if enddict[e] >= 1 and e != transcript.end:
					return True
			else:
				return True
	return False

def unique_values_in_dict_sets(inputdict):#For each value in set per key in dictionary, find out if that value is unique to that key. if so add value to outputdictionary
	outputdict = {}
	for sample in inputdict.keys():
		outputdict[sample] = set()
		for gene in inputdict[sample]:
			notunique = False
			for sample1 in inputdict.keys():
				if sample1 == sample:
					pass
				else:
					if gene in inputdict[sample1]:
						notunique = True
					else:
						pass
			if notunique == False:
				outputdict[sample].add(gene)
			else:
				pass
	return outputdict
		
def venn_diagram(filename,listofsubsets,labels):
	from matplotlib_venn import venn2,venn2_circles
	plt.figure(figsize=(7,7))
	v = venn2(listofsubsets, set_labels = labels)
	c = venn2_circles(listofsubsets, linestyle='-', linewidth=1, color="black")
	plt.title("")
	plt.savefig(filename)

def write_outputs(path,path2,outputfile_prefix,lowconfidence_genes,tama_merge_input,fulltranscriptlist,filteredtranscriptlist,fullgenelist,filteredgenelist,mid_genes,gene_rescue_set,trueTSS_dict2,reads,genes,trueTES_dict2,transcripts,all_good_sjs,list_of_sj_fails,allreadlengths):
	#Now write low confidence genes to file
	low_confidence_output = open(path + outputfile_prefix + "_low_confidence_genes.txt","w")
	for gene in lowconfidence_genes:
		low_confidence_output.write(gene + "\n")
	low_confidence_output.close()
	reportout=open(path2 + "filtering_report.txt","w")
	filelist=open(path2 + "file_list_filtered.txt","w")#Input for tama merge
	filelist.write(tama_merge_input + "\t" + "capped" + "\t" + "1,1,1" + "\t" + tama_merge_input + "\n")
	filelist.close()
	##Print number of transcripts in, number out and number removed, same with genes
	#print ('{:,}'.format(value))
	reportout.write("Total number of reads in input: " + '{:,}'.format(len(reads.Instances.keys())) + "\n")
	reportout.write(f"Longest mapped read length: {max(allreadlengths)}")
	reportout.write("Total number of trancripts in input: " + '{:,}'.format(len(fulltranscriptlist)) + "\n")
	reportout.write("Total number of trancripts that passed filtering: "+'{:,}'.format(len(filteredtranscriptlist)) + "\n")
	reportout.write("Total number of genes in input: " + '{:,}'.format(len(fullgenelist)) + "\n")
	reportout.write("Total number of genes that passed filtering: " + '{:,}'.format(len(filteredgenelist)) + "\n")
	reportout.write("Total number of genes with enriched TSS and TES: " + '{:,}'.format(len(trueTSS_dict2.keys())) + "\n")
	reportout.write("Total number of genes with starts and ends predicted from fixed window: " + '{:,}'.format(len(mid_genes)) + "\n")
	reportout.write("Total number of low confidence genes: " + '{:,}'.format(len(lowconfidence_genes)) + "\n")
	reportout.write("Total number of low confidence genes rescued (these are in the filtered RTD): " + '{:,}'.format(len(gene_rescue_set)) + "\n")
	reportout.write("Total number of high confidence splice junctions: " + '{:,}'.format(len(all_good_sjs)) + "\n")
	reportout.write("Total number of low confidence splice junctions (removed): " + '{:,}'.format(len(list_of_sj_fails.keys())) + "\n")
	reportout.write("Total number of low confidence splice junctions (removed) non canonical: " + '{:,}'.format(len([sj for sj in list_of_sj_fails.keys() if list_of_sj_fails[sj]["non_canonical"]])) + "\n")
	reportout.write("Total number of low confidence splice junctions (removed) sequence error: " + '{:,}'.format(len([sj for sj in list_of_sj_fails.keys() if list_of_sj_fails[sj]["error"]])) + "\n")
	reportout.write("Total number of low confidence splice junctions (removed) potential template switch: " + '{:,}'.format(len([sj for sj in list_of_sj_fails.keys() if list_of_sj_fails[sj]["RT switch"]])) + "\n")
	
	#reportout.write("Total number of low confidence genes recued by short read splice junctions: " + str(len(gene_rescue_set_again)))
	#print("Number of transcripts with single exon removed: "+'{:,}'.format(len(single_exon_transcripts_low_read)))
	reportout.close()
	starts_out = open(path2 + "TSS_information.txt","w")#Columns: gene, TSS start pos, read support, p value
	ends_out = open(path2 + "TES_information.txt","w")#Columns: gene, TES start pos, read support, p value
	for gene in trueTSS_dict2.keys():
		chromosome = genes.get(gene,transcripts).chromosome
		strand = genes.get(gene,transcripts).strand
		for start in trueTSS_dict2[gene].keys():
			starts_out.write(gene + "\t" + chromosome + "\t" + str(start) + "\t" + str(trueTSS_dict2[gene][start][0]) + "\t" + str(trueTSS_dict2[gene][start][1]) + "\t" + strand + "\n")
	for gene in trueTES_dict2.keys():
		chromosome = genes.get(gene,transcripts).chromosome
		strand = genes.get(gene,transcripts).strand
		for end in trueTES_dict2[gene].keys():
			ends_out.write(gene + "\t" + chromosome + "\t" + str(end) + "\t" + str(trueTES_dict2[gene][end][0]) + "\t" + str(trueTES_dict2[gene][end][1]) + "\t" + strand + "\n")
	starts_out.close()
	all_starts = open(path2 + "all_TSS_with_support.txt","w")
	for transcript in transcripts.Instances.keys():
		t = transcripts.Instances[transcript]
		all_starts.write(t.chromosome + "\t" + str(t.start) + "\t" + str(len(t.reads)) + "\t" + t.strand + "\t" + t.gene + "\n")
	ends_out.close()
	all_starts.close()
	
	



#Input - Combined file of all samples with following columns:#Order: readid,transid,sjno,error_profile,sj_coordinates,sj_sequence

def main():
	#Parse arguments
	parser = argparse.ArgumentParser(description='Filtering of Iso-Seq based bed file using high confidence sjs and TSS/TES sites. Outputs a filtered bed file as well as other figures')
	parser.add_argument('-i', dest = 'infile', type = str, help = 'Splice junction table produced by generate_filter_information.py. Include full path')
	parser.add_argument('-bed', dest = 'bedinput', type = str, help = 'Tama merge output bedfile (with transcripts overlapping Ns removed). Include full path', default = "")
	parser.add_argument('-bedn', dest = 'bedinput2', type = str, help = 'Tama merge output bedfile (before N removal).If N removal not required put bedinput here. Include full path')
	parser.add_argument('-sr', dest = 'short_read_input', type = str, help = 'STAR input of short read sjs')
	
	parser.add_argument('-g', dest = 'genome_input', type = str, help = 'Genome used for mapping')
	parser.add_argument('-pA', dest = 'polyainput', type = str, help = 'PolyA input from tama collapse')
	parser.add_argument('-s', dest = 'single_exon_input', type = str, help = 'Single exon input file from generate_filter_information.py')
	parser.add_argument('-o', dest = 'outputfile_prefix', type = str, help = 'Prefix to be used in outputs')
	parser.add_argument('--hamming', dest = 'hammingthreshold', type = int, help = 'For template switching,threshold hamming distance below which sj considered RT switching. FOr example 2 means a difference of 2 bases in 8. Default = 1', default = 1)
	parser.add_argument('--polyA', dest = 'polyathreshold', type = int, help = "threshold for percentage of As at 3' end of gene, above which read is removed", default = 80)
	parser.add_argument('--st_window', dest = 'startwindow', type = int, help = "Size of the window for removing unsupported 5' (+/- n)", default = 20)
	parser.add_argument('--end_window', dest = 'endwindow', type = int, help = "Size of the window for removing unsupported 3' (+/- n)", default = 60)
	args = parser.parse_args()
	path = "/".join(args.infile.split("/")[:-1]) + "/"
	#path = "/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/Splicejunction/"#Path to where input files are kept (except tama merge outputs and inputs)
	path2 = "/".join(args.bedinput2.split("/")[:-1]) + "/"
	#path2 = "/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/"
	tama_merge_input = path + args.outputfile_prefix + "_merged_filtered.bed"
	if not args.bedinput:
		bedinput = args.bedinput2
	else:
		bedinput = args.bedinput
	start1 = time.time()
	print("parsing input file")
	reads, splice_junctions = parse_sj_table(args.infile)
	stop2 = time.time() - start1

	print("Parsing single exon input file")
	reads = single_exon_parser(args.single_exon_input)

	print("Finding reads for merged transcripts...")
	#t_read_dict = {}#chromosome[0], start[1], end[2], set of sjs[3]seperated by "_" is key, set of reads is item
	t_read_dict = generate_t_read_dict(reads)
	##
	print("Matching transcripts to reads...")
	all_matched_reads = set()
	transcript_no_read_support = set()
	matched_coordinates = set()
	notmatched_coordinates = set()
	genes, transcripts = parse_bed_input(bedinput,t_read_dict,reads,all_matched_reads,transcript_no_read_support,matched_coordinates,notmatched_coordinates,False)


	#Create histogram of mapped read lengths
	allreadlengths = []
	for line in open(bedinput):
		bed = Bedline(line)
		try:
			transcript = transcripts.get(bed.transcript)
		except:
			continue
		allreadlengths += [bed.sequence_length()] * len(transcript.reads)
	histogram_from_list(allreadlengths,50,"histogram_of_mapped_exonreadlengths_breaks50","Read length","blue")
	print(f"max read length (exon only) = {max(allreadlengths)}")

	allreadlengths = []
	for t in transcripts.Instances.keys():
		transcript = transcripts.get(t)
		allreadlengths += [transcript.total_length] * len(transcript.reads)

	histogram_from_list(allreadlengths,50,"histogram_of_mapped_readlengths_breaks50","Read length","blue")
	#Problem idenfified: Genes with transcripts removed can have scewed fragment profiles. These transcripts have been rmeoved due to overlap with Ns. Therefore find these transcripts. Then remove these genes
	bedlines = set(open(bedinput).readlines())
	bedlines2 = set(open(args.bedinput2).readlines()) #This bedfile has transcripts overlapping polyNs

	bedlines_overlapNs = []
	#Find which transcripts overlap Ns
	for line in bedlines2:
		if line not in bedlines:
			bedlines_overlapNs.append(line)
		else:
			pass
	genes, transcripts = parse_bed_input(bedlines_overlapNs,t_read_dict,reads,all_matched_reads,transcript_no_read_support,matched_coordinates,notmatched_coordinates,True)
	not_matched_read_co = t_read_dict.keys() - matched_coordinates
	print("transcripts no read support " + str(len(transcript_no_read_support)))
	reads_not_matched = reads.Instances.keys() - all_matched_reads
	print("reads with no transcript support " + str(len(reads_not_matched)))
	print(list(reads_not_matched)[0:10])
	print(list(transcript_no_read_support)[0:10])


	########
	print("Creating database of high confidence TSS and TES...")
	trueTSS_dict = {}
	trueTSS_dict2 = {}
	trueTES_dict2 = {}
	lowconfidence_genes = set()
	degredation_dict = {}#Gene is key, percent degredation is value
	degraded_genes = set()
	mid_confidence_start_transcripts = set()

	for gene in genes.Instances.keys():
		#print(gene)
		gene_process = genes.get(gene,transcripts)
		binomial_simple(gene_process,mid_confidence_start_transcripts,trueTSS_dict,trueTSS_dict2,trueTES_dict2,degredation_dict,degraded_genes,args.startwindow,args.endwindow,transcripts,lowconfidence_genes)
	mid_genes = set()
	for t in mid_confidence_start_transcripts:
		mid_genes.add(t.split(".")[0])

	assert len(mid_genes) + len(trueTSS_dict2.keys()) + len(lowconfidence_genes) == len(genes.Instances.keys())

	"""for gene in degredation_dict.keys():
		self = Gene(gene)
		degradation = degredation_dict[gene]
		total_degradation += degradation * self.expression
		total_reads += self.expression

	print(str(total_degradation/total_reads))
	#Degradation on 23/03/20 was 4.9 percent"""

	print("Parsing polyA input...")
	#Aim: Remove transcripts with polyAs 20bp at 3' end of gene
	#Need combined polyA.txt file from tama collapse outputs
	polyareaddict = parse_polyA(args.polyainput,args.polyathreshold)

	#Parse genome for template switching analysis
	print("Opening genome...")
	chromosomedict = parse_genome(args.genome_input)
	print("Filter splice junctions...")
	all_good_sjs, list_of_sj_fails = filter_sjs(splice_junctions,chromosomedict,args.hammingthreshold)
	####
	#Now rescue more sjs with short read information
	#Parse short read sj information (in tabulated STAR format)
	short_read_sj = set()
	if args.short_read_input:
		###
		short_read_sj = parse_short_read_STAR(args.short_read_input,conversiondict,readsupportcutoff,min_overhang)
		#Check whether sjs match up
	all_sj = set(all_good_sjs.keys())#Just use high confident sjs, 162756 30/12/19
		####
	in_both = set()
	for sj in short_read_sj:
		if sj in all_sj:
			in_both.add(sj)
	vennlist=[short_read_sj,all_sj]
	venn_names=("Illumina","Pacbio")
	venn_diagram(path2 + "Splice junctions.png",vennlist,venn_names)
	####
	all_good_sjs, sj_rescue_set = short_sj_rescue(splice_junctions,all_good_sjs,short_read_sj,list_of_sj_fails)
	#First identify genes with all transcripts thrown out
	#Iterate through these, rescue based on double sj support alone
	# rescue transcript with greatest number of supported sjs
	gene_rescue_set, transcript_rescue_set = rescue_low_expressed(lowconfidence_genes,genes,transcripts,all_good_sjs,short_read_sj)
	#######
	###Open merge.bed 
	print("writing output to file")
	fulltranscriptlist = set()
	filteredtranscriptlist = set()
	fullgenelist = set()
	filteredgenelist = set()
	#####
	gene_transcript_dict_binomial = filter_bed(bedinput,genes,transcripts,fulltranscriptlist,fullgenelist,filteredtranscriptlist,filteredgenelist,gene_rescue_set, transcript_rescue_set,polyareaddict,mid_genes,mid_confidence_start_transcripts,lowconfidence_genes,trueTSS_dict2,trueTES_dict2,all_good_sjs,path,args.outputfile_prefix,args.startwindow,args.endwindow)
	#####
	write_outputs(path,path2,args.outputfile_prefix,lowconfidence_genes,tama_merge_input,fulltranscriptlist,filteredtranscriptlist,fullgenelist,filteredgenelist,mid_genes,gene_rescue_set,trueTSS_dict2,reads,genes,trueTES_dict2,transcripts,all_good_sjs,list_of_sj_fails,allreadlengths)
	#Create more barplots!
	sjfails_non_canonical = []
	sj_error_only = []
	sjfails_RTswitch = []
	for sj in list_of_sj_fails.keys():
		if list_of_sj_fails[sj]["non_canonical"]:
			sjfails_non_canonical.append(sj)
		if list_of_sj_fails[sj]["error"]:
			sj_error_only.append(sj)
		if list_of_sj_fails[sj]["RT switch"]:
			sjfails_RTswitch.append(sj)
	##
	error_dict={"Non canonical":sjfails_non_canonical,"Sequence errors":sj_error_only,"RT switching":sjfails_RTswitch}
	bar_plot(error_dict,path2+'SJ errors.png',"","Number of splice junctions","blue",False)
	readfails = open(path + args.outputfile_prefix + "_reads_no_match","w")
	for read in reads_not_matched:
		readfails.write(read + "\n")
	readfails.close()
	################################
	print("Done. Time of first loop: "+str(stop2))
	################################
	stop_time = time.time() - start1
	###################
	print("Done. Total time: "+str(stop_time))

if __name__ == "__main__":
	main()

