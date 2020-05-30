
#Script to predict likely PCR product sizes and match to RT-PCR data

#Rules:
#Products >600bp and < 90bp are a bit dodge, these will not be counted
#Primers will amplify where they match best. The hit does not necessarily need to be 100% however - if a primer has a different base at its proximal end, it will bind fine (even if two bases). Primers will bind if a match is missing in the middle, but not as well, so if there is a better binding site, this will dominate.
#ACGTC
#|||||
#ACGTC
#TGCAG
#First work out whether primer is a good hit according to the citeria above. If it is, store the information for that primer in 
#Furthermore, once likely products from different transcripts are generated, these need to be clustered. As with Paulo's script, products will be clustered if they are within a window of +/- 6, and will be counted as the corresponding RT-PCR product using the same criteria.
#Also possible 2 products may be amplified from same transcript, these will ahve to be flagged

#See below for BLAST output format
#"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen salltitles qcovs qcovhsp"
#53F	G3395;G3395.23(-)	100.000	18	0	0	1	18	421	438	0.010	36.2	18	3092	G3395;G3395.23(-)	100	100
#Parse BLAST output
import sys
import os
from plotnine import *
import pandas as pd
import numpy
import scipy.stats
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import argparse


#window = 4

#BLAST_input = "/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/Primers_vs_AllTranscripts.txt"

#RT_PCR_input = "/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/rtPCR_productsAndProportions.txt"

#salmon_quants_folder = "/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/salmon_quants/quants/"

#scatterplot_output = "/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/RT_PCR/BaRT_2_Iso1_correlation_complete_matches_spuriousRT_PCR_removed"

class Parse_BLAST:
	def __init__(self,line):
		line1 = line.rstrip("\n").split("\t")
		self.primer = line1[0]
		self.primer_direction = self.primer[-1]
		self.primer_name = self.primer[:-1]
		self.transcript_id = line1[1]
		self.gene = self.transcript_id.split(";")[0]
		self.pident = float(line1[2])
		self.hit_length = int(line1[3])
		self.mismatches = int(line1[4])
		self.qstart = int(line1[6])
		self.qend = int(line1[7])
		self.sstart = int(line1[8])
		self.send = int(line1[9])
		if self.sstart > self.send and self.primer_direction == "F" or self.sstart < self.send and self.primer_direction == "R":
			self.wrong_way = True
		else:
			self.wrong_way = False
		self.bitscore = float(line1[11])
		self.qlen = int(line1[12])
		self.slen = int(line1[13])
		self.qcovs = int(line1[15])
		self.qcovhsp = int(line1[16])
	def process(self,primer_hits_F,primer_hits_R):
		if self.primer_direction == "F":
			if self.pident == 100 and self.qcovs == 100 and self.qcovhsp == 100: #Perfect match
				try:
					primerdict = primer_hits_F[self.primer_name]
					primerdict.setdefault("perfect",set()).add(self)
					primer_hits_F[self.primer_name] = primerdict
				except KeyError:
					primer_hits_F[self.primer_name] = {"perfect":{self}}
			elif self.mismatches == 0 and self.qstart in range(2,4) and self.qend == self.qlen: #Not quite as good but pretty good
				try:
					primerdict = primer_hits_F[self.primer_name]
					primerdict.setdefault("second",set()).add(self)
					primer_hits_F[self.primer_name] = primerdict
				except KeyError:
					primer_hits_F[self.primer_name] = {"second":{self}}
			elif self.mismatches == 1 and self.qcovs == 100 and self.qcovhsp == 100: #Least good
				try:
					primerdict = primer_hits_F[self.primer_name]
					primerdict.setdefault("third",set()).add(self)
					primer_hits_F[self.primer_name] = primerdict
				except KeyError:
					primer_hits_F[self.primer_name] = {"third":{self}}
			else:
				pass
		elif self.primer_direction == "R":
			if self.pident == 100 and self.qcovs == 100 and self.qcovhsp == 100: #Perfect match
				try:
					primerdict = primer_hits_R[self.primer_name]
					primerdict.setdefault("perfect",set()).add(self)
					primer_hits_R[self.primer_name] = primerdict
				except KeyError:
					primer_hits_R[self.primer_name] = {"perfect":{self}}
			elif self.mismatches == 0 and self.qstart in range(2,4) and self.qend == self.qlen: #Not quite as good but pretty good
				try:
					primerdict = primer_hits_R[self.primer_name]
					primerdict.setdefault("second",set()).add(self)
					primer_hits_R[self.primer_name] = primerdict
				except KeyError:
					primer_hits_R[self.primer_name] = {"second":{self}}
			elif self.mismatches == 1 and self.qcovs == 100 and self.qcovhsp == 100: #Least good
				try:
					primerdict = primer_hits_R[self.primer_name]
					primerdict.setdefault("third",set()).add(self)
					primer_hits_R[self.primer_name] = primerdict
				except KeyError:
					primer_hits_R[self.primer_name] = {"third":{self}}
			else:
				pass

#Cluster similar products from different transcripts (+/- window)
# #Need to have this format: {1:{1,2,3},23:{23,24,25}
def cluster_products(product_size,product_sizes,window):
	for product in product_sizes:
		for i in range(product[0] - window, product[0] + window +1):
			if i == product_size:
				return True
	return False

def cluster_transcriptome_products(product_dict,products_clustered):
	for primer_name in product_dict.keys():
		products_clustered[primer_name] = {}
		clusters = {}
		clustersinfo = {}
		for product_info in product_dict[primer_name]:
			if clusters == {}:
				clustersinfo[product_info[0]] = product_info
				clusters[product_info[0]] = [product_info]
			else:
				clustername = False
				for cluster1 in clusters.keys():
					if cluster_products(product_info[0],clusters[cluster1],6):
						clustername = cluster1
						break
				if not clustername:
					clusters[product_info[0]] = [product_info]
				else:
					clusters[clustername].append(product_info)
		#Now get averages of clusters
		average_clusters = {}
		for cluster1 in clusters.keys():
			size_list = []
			for product_info in clusters[cluster1]:
				size_list.append(product_info[0])
			average_size = int(sum(size_list)/len(size_list))
			average_clusters[average_size] = clusters[cluster1]
		products_clustered[primer_name] = average_clusters

def combine_salmon_reps(salmon_quants_dict,repinfo):
	salmon_quants_dict2 = {}
	for transcript in salmon_quants_dict:
		samples = salmon_quants_dict[transcript].keys()
		sample_rep_dict = {} #sample name is key, rep names are set
		for sample in samples:
			sample_rep_dict.setdefault(sample[:inv(repinfo)],set()).add(sample)
		sample_average_quants = {}#sample is key, average quant is value
		quants = []
		for sample in sample_rep_dict.keys():
			for rep in sample_rep_dict[sample]:
				quants.append(salmon_quants_dict[transcript][rep])
			sample_average_quants[sample] = sum(quants)/len(quants)
		salmon_quants_dict2[transcript] = sample_average_quants
	return salmon_quants_dict2

def convert_names(salmon_quants_dict2):
	#Now convert names (just for bart1)
	#Note from Micha No, these are Paulo's original names (from Stringtie). We converted them to BART IDs because the original names were awful and useless. The ID numbers have stayed the same though -- so converting is easy. So for example MSTRG.22 becomes BART1_0-p00022 (or BART1_0-u00022 in the unpadded version of BART). Syntax is BART1_0-[p|u] followed by the five digit ID number, left-padded with leading zeroes.
	#My gene name: MSTRG.1.1gene=MSTRG.1
	#BART1_0-p26461.009
	#BART1_0-p42835.011
	salmon_quants_dict3 = {}
	for transcript in salmon_quants_dict2:
		if transcript.startswith("MSTRG"):
			t_split = transcript.split("gene")[0].split(".")
			bartnumber = ""
			decimal = ""
			for i in range(0,(5 - len(t_split[1]))):
				bartnumber += "0"
			for i in range(0,(3 - len(t_split[-1]))):
				decimal += "0"
			bartnumber += str(t_split[1])
			decimal += str(t_split[-1])
			salmon_quants_dict3["BART1_0-u" + bartnumber + "." + decimal] = salmon_quants_dict2[transcript]
		else:
			return salmon_quants_dict2
	print(list(salmon_quants_dict3.keys())[0:10])
	return salmon_quants_dict3

#Identify all possible products <600bp and >90bp, rank in order of likelyhood
#Go through perfect hits first
#find possible corresponding hits in R
#Perfect hitting primers will hit multiple transcripts. Only take these products
def get_primer_products(primer_name,set_F,set_R,rank,product_dict):
	for primer_info_F in set_F:
		for primer_info_R in set_R:
			if primer_info_F.transcript_id == primer_info_R.transcript_id:
				if not primer_info_F.wrong_way or not primer_info_R.wrong_way:
					product_size = primer_info_R.sstart - primer_info_F.sstart
					if product_size < 600 and product_size > 90:
						product_dict.setdefault(primer_name,[]).append([product_size,primer_info_F.transcript_id,rank])

def get_proportions(complete_matches,matched_products,products_clustered,salmon_quants_dict,RT_PCR_proportions_dict,all_RT_PCR_proportions,all_transcriptome_proportions,outtable):
	out = open(outtable,"w")
	out.write("primer name\tSample name\tRT PCR product\tMatched clustered product\tBLAST product\tPrimer match information\tTranscript matched\tTPM of transcript\tRT PCR product proportion\tTotal transcript quant proportion\n")
	for primer_name in complete_matches:
		RT_PCR_products = matched_products[primer_name]
		combined_tpms_dict = {}
		transcriptome_proportions = {}
		for RT_product in RT_PCR_products:
			transcriptome_products = matched_products[primer_name][RT_product]
			sample_quants_dict = {}
			for t_clustered_product in transcriptome_products:
				#sample is key, list of TPM values is value
				for t_product in products_clustered[primer_name][t_clustered_product]:
						transcript = t_product[1]
						try:
							transcript_quants = salmon_quants_dict[transcript]
							for sample in transcript_quants.keys():
								sample_quants_dict.setdefault(sample,[]).append(transcript_quants[sample])
						except KeyError:
							print("Warning: Transcript " + transcript + " was not found in salmon quants")
			combined_tpms_dict[RT_product] = {}
			for sample in sample_quants_dict.keys():
				combined_tpms_dict[RT_product][sample] = sum(sample_quants_dict[sample])
		for RT_product in combined_tpms_dict.keys():
			for sample in combined_tpms_dict[RT_product]:
				transcriptome_proportions[sample] = {}
		for sample in transcriptome_proportions.keys():
			total_tpm = 0
			for RT_product in combined_tpms_dict.keys():
				total_tpm += combined_tpms_dict[RT_product][sample]
			for RT_product in combined_tpms_dict.keys():
				if total_tpm != 0:
					proportion = combined_tpms_dict[RT_product][sample] / total_tpm
				else:
					proportion = 0
				transcriptome_proportions[sample][RT_product] = proportion
		for sample in transcriptome_proportions.keys():
			for RT_product in transcriptome_proportions[sample].keys():
				for matched_product in matched_products[primer_name][RT_product]:
					for match in products_clustered[primer_name][matched_product]:
						transcript = match[1]
						#product_dict[product_size,primer_info_F.transcript_id,rank]
						out.write(primer_name + "\t")
						out.write(sample + "\t")
						out.write(str(RT_product) +"\t")
						out.write(str(matched_product) + "\t")
						out.write(str(match[0]) + "\t")
						out.write(",".join(match[2]) + "\t")
						out.write(transcript + "\t")
						out.write(str(salmon_quants_dict[transcript][sample]) + "\t")
						try:
							out.write(str(RT_PCR_proportions_dict[primer_name][RT_product][sample]) + "\t")
						except KeyError:
							out.write("\t")
						out.write(str(transcriptome_proportions[sample][RT_product]) + "\n")
				try:#If missing value skip
					all_RT_PCR_proportions.append(RT_PCR_proportions_dict[primer_name][RT_product][sample])
					all_transcriptome_proportions.append(transcriptome_proportions[sample][RT_product])
				except KeyError:
					pass
	out.close()

def get_primer_proportions(complete_matches,matched_products,products_clustered,salmon_quants_dict,RT_PCR_proportions_dict,all_RT_PCR_proportions,all_transcriptome_proportions,scatterplot_output):
	for primer_name in complete_matches:
		all_RT_PCR_proportions = []
		all_transcriptome_proportions = []
		RT_PCR_products = matched_products[primer_name]
		combined_tpms_dict = {}
		transcriptome_proportions = {}
		for RT_product in RT_PCR_products:
			transcriptome_products = matched_products[primer_name][RT_product]
			sample_quants_dict = {}
			for t_clustered_product in transcriptome_products:
				#sample is key, list of TPM values is value
				for t_product in products_clustered[primer_name][t_clustered_product]:
						transcript = t_product[1]
						try:
							transcript_quants = salmon_quants_dict[transcript]
							for sample in transcript_quants.keys():
								sample_quants_dict.setdefault(sample,[]).append(transcript_quants[sample])
						except KeyError:
							print("Warning: Transcript " + transcript + " was not found in salmon quants")
			combined_tpms_dict[RT_product] = {}
			for sample in sample_quants_dict.keys():
				combined_tpms_dict[RT_product][sample] = sum(sample_quants_dict[sample])
		for RT_product in combined_tpms_dict.keys():
			for sample in combined_tpms_dict[RT_product]:
				transcriptome_proportions[sample] = {}
		for sample in transcriptome_proportions.keys():
			total_tpm = 0
			for RT_product in combined_tpms_dict.keys():
				total_tpm += combined_tpms_dict[RT_product][sample]
			for RT_product in combined_tpms_dict.keys():
				proportion = combined_tpms_dict[RT_product][sample] / total_tpm
				transcriptome_proportions[sample][RT_product] = proportion
		for sample in transcriptome_proportions.keys():
			for RT_product in transcriptome_proportions[sample].keys():	
				try:#If missing value skip
					all_RT_PCR_proportions.append(RT_PCR_proportions_dict[primer_name][RT_product][sample])
					all_transcriptome_proportions.append(transcriptome_proportions[sample][RT_product])
				except KeyError:
					pass
		scatter_plot(all_transcriptome_proportions,all_RT_PCR_proportions,scatterplot_output + primer_name,"transcriptome proportions","RT PCR proportions","black",False)
		print("Pearson correlation for " + primer_name + ": " + str(scipy.stats.pearsonr(all_transcriptome_proportions,all_RT_PCR_proportions)))
		print("Spearman ranked correlation for " + primer_name + ": " + str(scipy.stats.spearmanr(all_transcriptome_proportions,all_RT_PCR_proportions)))

def inv(t):
	if t < 0:
		n=abs(t)
	elif t > 0:
		n=-t
	elif t==0:
		n=t
	else:
		print("inverse error!")
	return n

def in_window(number1,number2,window):
	for i in range(number1 - window, number1 + window + 1):
		if i == number2:
			return True
	return False

def match_to_RT_PCR_products(products_clustered,matched_products,complete_matches,perfect_matches,one_match_only,RT_PCR_size_dict,window):
	for primer_name in products_clustered.keys():
		try:
			RT_PCR_products1 = RT_PCR_size_dict[primer_name]
		except KeyError:
			print("Warning! Primer set " + primer_name + " not in RT PCR dataset")
			continue
		product_dict = {}
		if len(products_clustered[primer_name].keys()) < len(RT_PCR_products1) and len(RT_PCR_products1) == 2:
			one_match_only.add(primer_name)
			continue
		RT_products_matched = set()
		for product in products_clustered[primer_name].keys():
			for product_main in RT_PCR_size_dict[primer_name]:
				if product in RT_products_matched:
					continue
				if in_window(product,product_main,window - 1):
					product_dict.setdefault(product_main,set()).add(product)
					RT_products_matched.add(product)
		matched_products[primer_name] = product_dict
		if product_dict.keys() == RT_PCR_size_dict[primer_name]:
			if len(product_dict.keys()) == len(products_clustered[primer_name].keys()):
				perfect_matches.add(primer_name) #Same number of products in both RT-PCR and predicted from transcriptome, and products are the same
				complete_matches.add(primer_name)
			else:
				complete_matches.add(primer_name) #Products from RT-PCR are all present in transcriptome data, but some predicted products are not in RT-PCR data (this is entirely possible as they could be found on a sample not in RT-PCR data)

def parse_RT_PCR_table(RT_PCR_input,RT_PCR_size_dict,RT_PCR_proportions_dict):
	for line in open(RT_PCR_input):
		if line.startswith("Primer"):
			samples = line.rstrip("\n").split("\t")
			continue
		line1 = line.rstrip("\n").split("\t")
		primer_name = line1[0]
		size = int(float(line1[1]))
		try:
			sizedict = RT_PCR_proportions_dict[primer_name]
		except KeyError:
			sizedict = {}
		proportions = {}
		for i in range(2,len(line1)):
			try:
				proportions[samples[i]] = float(line1[i])
			except ValueError:
				print("Warning! Primer set has missing values")
		if sum(proportions.values()) < 0.5 or sum(proportions.values()) == 12:
			continue#Some weird stuff in this RT_PCR table!!Remove	
		sizedict[size] = proportions
		RT_PCR_proportions_dict[primer_name] = sizedict
		RT_PCR_size_dict.setdefault(primer_name,set()).add(size)

def parse_salmon_input(salmon_quants_folder,salmon_quants_dict):
	bigdict = {}
	listofsamples = []
	for sample in os.listdir(salmon_quants_folder):
		if "." in sample:#Only iterate through folders, not files
			continue
		listofsamples.append(sample)
		#Extract name  and TPM values for each sample, have sample as key with (transcript,TPM) as values
		try:
			for line in open(salmon_quants_folder + sample + "/quant.sf"):
				if line.startswith("Name"):
					continue
				line1 = line.split("\t")
				bigdict.setdefault(line1[0],[]).append((sample,line1[3]))
		except FileNotFoundError:#If folder does not have a quant.sf file, just ignore
			continue
	for transcript in bigdict.keys():
		sample_dict = {}
		for sample_tpm in bigdict[transcript]:
			sample_dict[sample_tpm[0]] = float(sample_tpm[1])
		salmon_quants_dict[transcript] = sample_dict

def products_in_order(primer_hits_F,primer_hits_R,all_primer_pairs_with_hits,product_dict,failed_primers,category_dict):
	categories = ["perfect","second"] #category third is parsed, but at the moment is not included in the analysis. This is because primers of the third category are likely to be outcompeted by primer matches of other categories. The below algorithm would need to be adjusted to take this into account. 
	for primer_name in all_primer_pairs_with_hits:
		f_category_list = []
		r_category_list = []
		possible_combinations = []
		for category in categories:
			if category in primer_hits_F[primer_name].keys():
				f_category_list.append(category)
		for category in categories:
			if category in primer_hits_R[primer_name].keys():
				r_category_list.append(category)
		#Now create list of possible combinations in order of rank
		for i in range(0,len(f_category_list)):
			for t in range(0,len(r_category_list)):
				possible_combinations.append((f_category_list[i],r_category_list[t]))
		category_dict[primer_name] = possible_combinations
		#Now generate products
		combinations_completed = set()
		for combination in possible_combinations:
			if combination in combinations_completed:
				continue
			if combination[0] == combination[1]:
				get_primer_products(primer_name,primer_hits_F[primer_name][combination[0]],primer_hits_R[primer_name][combination[1]],combination,product_dict)
			else:
				try:
					get_primer_products(primer_name,primer_hits_F[primer_name][combination[0]],primer_hits_R[primer_name][combination[1]],combination,product_dict)
					combinations_completed.add(combination)
				except KeyError:
					get_primer_products(primer_name,primer_hits_F[primer_name][combination[1]],primer_hits_R[primer_name][combination[0]],combination,product_dict)
					new_combination = (combination[1],combination[0])
					combinations_completed.add(new_combination)
				try:#Try the other way
					get_primer_products(primer_name,primer_hits_F[primer_name][combination[1]],primer_hits_R[primer_name][combination[0]],combination,product_dict)
					new_combination = (combination[1],combination[0])
					combinations_completed.add(new_combination)
				except KeyError:
					pass
		if primer_name not in product_dict.keys():
			print("Primer pair " + primer_name + " did not align properly to transcriptome")
			failed_primers.add(primer_name)

def scatter_plot(inputx,inputy,filename,x_axis_title,y_axis_title,bar_colour,axis_angle_90):#inputs must be lists of integers
	#histokeys,histovalues=histokeys[1:],histovalues[1:]
	pandadict = {x_axis_title:inputx,y_axis_title:inputy}
	scatterdata=pd.DataFrame(pandadict)
	t = ggplot(aes(x=x_axis_title, y=y_axis_title), scatterdata) + geom_point()
	p = t + theme_bw() + labs(title="", x=x_axis_title,y=y_axis_title)
	p.save(filename)

def write_primers_to_fasta(primers,outputfile):
	out = open(outputfile,"w")
	for primer in primers:
		out.write(primer +"F\n")
		out.write(primer + "R\n")
	out.close()


def main():
	parser = argparse.ArgumentParser(description='Compare RT-PCR and Salmon quantifiactions and produce correlations')
	parser.add_argument('-b', dest = 'BLAST_input', type = str, help = 'BLASTn input')
	parser.add_argument('-p', dest = 'RT_PCR_input', type = str, help = 'RT_PCR input file')
	parser.add_argument('-s', dest = 'salmon_quants_folder',type = str, help = 'salmon quants folder path')
	parser.add_argument('-w', dest = 'window',type = int, help = 'Size of window for combining predicted primer products')
	parser.add_argument('-o', dest = 'scatterplot_output',type = str, help = 'Name of scatterplot output file')
	parser.add_argument('-reps', dest = 'reps',type = bool, help = 'Are there salmon quant reps, True or False?',default = False)
	parser.add_argument('-rep_info', dest = 'rep_info',type = int, help = 'Length of rep suffix. E.g <sample_name>_a,<sample_name>_b')
	parser.add_argument('-bart_names', dest = 'bart_names',type = bool, help = 'Particular to Bart1. Need to change names',default = False)
	#parser.parse_args(['--foo', '2'])
	args = parser.parse_args()
	destination_path = "/".join(args.scatterplot_output.split("/")[:-1]) + "/"

	salmon_quants_dict = {} #transcript id is key, value is dictionary with samples as keys and tpms as value
	print("Parsing salmon quant files...")
	parse_salmon_input(args.salmon_quants_folder,salmon_quants_dict)
	if args.reps: #Combines salmon reps together
		salmon_quants_dict = combine_salmon_reps(salmon_quants_dict,args.rep_info)
	if args.bart_names:
		salmon_quants_dict = convert_names(salmon_quants_dict)
	
	primer_hits_F = {}#key is self.primer_name, value is dictionary with key as perfect, second, third, values as primer_info object
	primer_hits_R = {}#key is self.primer_name, value is dictionary with key as perfect, second, third, values as primer_info object
	for line in open(args.BLAST_input):
		primer_info = Parse_BLAST(line)
		primer_info.process(primer_hits_F,primer_hits_R)
	#Generate set of paired primers:
	all_primer_pairs_with_hits = set()
	for primer_name in primer_hits_F.keys():
		if primer_name in primer_hits_R.keys():
			all_primer_pairs_with_hits.add(primer_name)
	print(str(len(all_primer_pairs_with_hits)))
	#Parse RT-PCR table
	#Primer	Size	B1	B4	B5	B6	B7	B46	B48	B53	B71	B74	B79	B87
	#Proportions for each sample
	RT_PCR_size_dict = {} #primer_name is key, set of sizes is value
	RT_PCR_proportions_dict = {}#primer_name is key, dict with paired sample name and proportion is value
	parse_RT_PCR_table(args.RT_PCR_input,RT_PCR_size_dict,RT_PCR_proportions_dict)
	#Cleanup - remove primers that have products within 6bp of each other
	product_dict = {}#primer_name is key, value is list of possible primer products, with information on transcript, and rank (p = perfect, s = second, t = third, order = FR. e.g an F primer with a perfect hit and an R primer with a second rate hit would give the rank ps)
	failed_primers = set() # list of all the failed primers that did not produce any good alignments
	category_dict = {} #primer_name is key, list of categories are value
	products_in_order(primer_hits_F,primer_hits_R,all_primer_pairs_with_hits,product_dict,failed_primers,category_dict)
	products_clustered = {} #key is primer_name, value is dict, with keys as products and values as list of transcripts which 
	cluster_transcriptome_products(product_dict,products_clustered)
	print(str(len(products_clustered.keys())))
	#Now compare this list to the list of predicted products from RT-PCR data
	matched_products = {}#primer_name is key, dictionary (with RT_PCR product as key and list of matched transcriptome products as values) as value
	complete_matches = set()
	perfect_matches = set()
	one_match_only = set()
	match_to_RT_PCR_products(products_clustered,matched_products,complete_matches,perfect_matches,one_match_only,RT_PCR_size_dict,args.window)
	#Now create correlations with RT_PCR data and transcriptome data
	#First parse quants file
	#Transcript	B46	B5	B74	B7	B1	B6	B87	B71	B48	B79	B4	B53
	all_RT_PCR_proportions = []
	all_transcriptome_proportions = []
	get_proportions(complete_matches,matched_products,products_clustered,salmon_quants_dict,RT_PCR_proportions_dict,all_RT_PCR_proportions,all_transcriptome_proportions,destination_path + "results.txt")
	scatter_plot(all_transcriptome_proportions,all_RT_PCR_proportions,args.scatterplot_output,"transcriptome proportions","RT PCR proportions","black",False)
	print("Data points: " + str(len(all_RT_PCR_proportions)))
	print("Pearson correlation: " + str(scipy.stats.pearsonr(all_transcriptome_proportions,all_RT_PCR_proportions)))
	print("Spearman ranked correlation: " + str(scipy.stats.spearmanr(all_transcriptome_proportions,all_RT_PCR_proportions)))
	print("Perfect matches: " + str(len(perfect_matches)))
	print("Complete matches: " + str(len(complete_matches)))
	print("All matches: " + str(len(matched_products.keys())))
	write_primers_to_fasta(complete_matches, destination_path + "primers_complete_match.txt")
	#get_primer_proportions(complete_matches,matched_products,products_clustered,salmon_quants_dict,RT_PCR_proportions_dict,all_RT_PCR_proportions,all_transcriptome_proportions,scatterplot_output)

if __name__ == "__main__":
	main()
