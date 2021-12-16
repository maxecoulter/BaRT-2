"""Script to filter monoexon genes based on TPM values"""
from matplotlib import pyplot as plt
import argparse




class Gtf:
	gtfentries = []
	def __init__(self,line):
		fields = line.strip().split("\t")
		self.chromosome = fields[0]
		self.assembler = fields[1]
		self.type = fields[2]
		self.start = int(fields[3])
		self.end = int(fields[4])
		self.strand = fields[6]
		self.info = fields[8]
		self.info_split = self.info.split(";")
		if self.type == "exon":
			self.exon_number = int(self.info_split[2].split('"')[1])
		for i in range(0,len(self.info_split)):
			if "transcript" in self.info_split[i]:
				self.transcript = self.info_split[i].split('"')[-2]
				break
		self.gene = self.transcript.split(".")[0]
		self.line = line
		self.gtfentries.append(self)
	def summarise(self):
		genes = set()
		transcripts = set()
		for gtf in self.gtfentries:
			genes.add(gtf.gene)
			transcripts.add(gtf.transcript)
		print(f"GTF file has {len(genes)} genes and {len(transcripts)} transcripts")
	def getmonoexons(self):#returns a list of genes that have monoexon only transcripts
		gene_exons = {}
		for gtf in self.gtfentries:
			if gtf.type == "exon":
				gene_exons.setdefault(gtf.gene,set()).add(gtf.exon_number)
		monoexons = [gene for gene, exon_number in gene_exons.items() if exon_number == {1}]
		print(f"Number of monoexon genes in dataset: {len(monoexons)}")
		return monoexons
	def filter(self, genes_for_removal,filtered_out,removed_out):
		"""Filter gtf file by removing genes that are in set provided"""
		with open(filtered_out,"w") as out1, open(removed_out, "w") as out2:
			[out2.write(gtf.line) if gtf.gene in genes_for_removal else out1.write(gtf.line) for gtf in self.gtfentries]


def filter_genes_by_tpm(gene_list, gene_transcripts, tpm_dict, threshold, monoexon_antisense,monoexon_intronic,monoexon_non_coding):
	genes_filtered = set()
	for gene in gene_list:
		gene_tpms = [0] * 20
		for transcript in gene_transcripts[gene]:
			tpms = tpm_dict[transcript]
			gene_tpms = [gene_tpms[i] + tpms[i] for i in range(0,len(tpms))] #Convert transcript tpms to gene tpms
		if sum(num >= threshold for num in gene_tpms) > 1:#If the gene is expressed in at least 2 samples
			genes_filtered.add(gene)
	if monoexon_non_coding:
		print(f"Number of non_coding monoexons: {len(monoexon_non_coding)}")
		print(f"Number of non_coding monoexons removed: {len(monoexon_non_coding - genes_filtered)}")
	if monoexon_antisense:
		print(f"Number of antisense monoexons: {len(monoexon_antisense)}")
		print(f"Number of antisense monoexons removed: {len(monoexon_antisense - genes_filtered)}")
	if monoexon_intronic:
		print(f"Number of intronic monoexons: {len(monoexon_intronic)}")
		print(f"Number of intronic monoexons removed: {len(monoexon_intronic - genes_filtered)}")

	print(f"{len(gene_list)} monoexons of various categories were identified for filtering based on their expression, {len(set(gene_list) - genes_filtered)} were removed due to low expression leaving {len(genes_filtered)}")
	return set(gene_list) - genes_filtered #Output the genes that need to be removed


def get_intronic(new_ids,gene_list,infile):
	#/mnt/shared/scratch/jentizne/BaRTv2_transcriptomes/BaRT_2_13_01March21/BaRT_2_13_01March21_RTDmaker_output/report/RTD_tables/BaRT_2_13_01March21_other_categories_IDs_intronic.tsv
	#List of intronic transcripts. Has header
	gene_set = set(gene_list)
	monoexon_intronic = set()
	for n,line in enumerate(open(infile)):
		if not n:
			continue
		g_old = new_ids[line.strip()].split(".")[0]
		if g_old in gene_set:
			monoexon_intronic.add(g_old)
	return monoexon_intronic #Set of genes

def get_antisense(new_ids,gene_list,infile):
	#/mnt/shared/scratch/jentizne/BaRTv2_transcriptomes/BaRT_2_13_01March21/BaRT_2_13_01March21_RTDmaker_output/report/RTD_tables/BaRT_2_13_01March21_other_categories_IDs_monoexonic_antisense.tsv
	#List of monoexon antisense transcripts. Has header
	gene_set = set(gene_list)
	monoexon_antisense = set()
	for n, line in enumerate(open(infile)):
		if not n:
			continue
		g_old = new_ids[line.strip()].split(".")[0]
		if g_old in gene_set:
			monoexon_antisense.add(g_old)
	return monoexon_antisense

def get_genes(tpm_dict):
	gene_transcripts = {} #First need mapping of gene and transcripts
	for transcript in tpm_dict.keys():
		gene_transcripts.setdefault(transcript.split(".")[0],set()).add(transcript)
	return gene_transcripts


def get_new_transcript_names(infile):
	#/mnt/shared/scratch/jentizne/BaRTv2_transcriptomes/BaRT_2_13_01March21/BaRT_2_13_01March21_RTDmaker_output/report/ID_lookup_tables/BaRT_2_13_01March21_table_ID_lookup_Transcripts.csv
	#First column is original,second is novel. Has header
	new_ids = {}
	for n, line in enumerate(open(infile)):
		if not n:
			continue
		old, new = line.strip().split(",")
		new_ids[new] = old #To get back to the old transcript names which we need
	return new_ids

def get_non_coding(new_ids,gene_list,infile):
	#/mnt/shared/scratch/jentizne/BaRTv2_transcriptomes/BaRT_2_13_01March21/BaRT_2_13_01March21_TS_output/BaRT_2_13_01March21_transfeat/BaRT_2_13_01March21_transfeat.csv
	#This is the main transuite table. From JC: The first column is the Gene_ID and third column is the Coding_potentiality, so you just have to extract transcripts for which the Gene_ID has one and only one Non_Coding flag in the third column.
	# look for "Non_Coding" in column 3
	gene_set = set(gene_list)
	non_coding_monoexon_all = {} #gene is key, value is set of all trnascript coding potentialities
	non_coding_monoexon = set()
	for n, line in enumerate(open(infile)):
		if not n:
			continue
		transcript_id, coding_potentiality = line.split(",")[1:3]
		old_g = new_ids[transcript_id].split(".")[0]
		if old_g in gene_set:
			non_coding_monoexon_all.setdefault(old_g,set()).add(coding_potentiality)
	for old_g, potentialities in non_coding_monoexon_all.items():
		if potentialities == {"Non_Coding"}:#Possible that some monoexon genes will ahve coding and non-coding transcripts. This will just take genes with only non-coding transcripts
			non_coding_monoexon.add(old_g)
	return non_coding_monoexon

def get_transcript_tpms(infile):
	tpm_dict = {}
	for n,line in enumerate(open(infile)):
		if not n:
			continue
		fields = line.strip().split("\t")
		tpm_dict[fields[0]] = [float(tpm) for tpm in fields[1:]] #Transcript id is key, list of tpms is value
	return tpm_dict

def mono_exon_tpm_histogram(tpm_dict, monoexons, gene_transcripts, threshold, output, intronic=None, antisense=None, non_coding=None):#threshold is value at which transcript is considered expressed in sample
	#Each histogram shows per sample support for each gene
	per_gene_histogram(monoexons,gene_transcripts,tpm_dict,threshold,output)
	per_sample_histogram(monoexons,gene_transcripts,tpm_dict,threshold,output + "_per_sample")
	if intronic:
		per_gene_histogram(intronic,gene_transcripts,tpm_dict,threshold,output + "_intronic")
		per_sample_histogram(intronic,gene_transcripts,tpm_dict,threshold,output + "_intronic_per_sample")
	if antisense:
		per_gene_histogram(antisense,gene_transcripts,tpm_dict,threshold,output + "_antisense")
		per_sample_histogram(antisense,gene_transcripts,tpm_dict,threshold,output + "_antisense_per_sample")
	if non_coding:
		per_gene_histogram(non_coding,gene_transcripts,tpm_dict,threshold,output + "_non_coding")
		per_sample_histogram(non_coding,gene_transcripts,tpm_dict,threshold,output + "_noncoding_per_sample")

def per_gene_histogram(gene_list,gene_transcripts,tpm_dict,threshold,output):
	number_of_samples = []
	one_tpms = []
	for gene in gene_list:
		gene_tpms = [0] * 20
		for transcript in gene_transcripts[gene]:
			tpms = tpm_dict[transcript]
			gene_tpms = [gene_tpms[i] + tpms[i] for i in range(0,len(tpms))] #Convert transcript tpms to gene tpms
		number_of_samples.append(sum(float(num) >= threshold for num in gene_tpms))#This is a bit cryptic - it is the number of samples with expression > threshold
		if sum(float(num) >= threshold for num in gene_tpms) == 1:
			one_tpms.append(gene)
	print(f"All numbers:{len(set(number_of_samples))}")
	#Print out histogram numbers
	for i in range(0,21):
		print(f"Number of genes with expression in {i} samples: {number_of_samples.count(i)}")
	
	print(f"Genes with expression in just one sample: {one_tpms[0:10]}")
	plt.figure()
	plt.hist(number_of_samples, bins=21)
	plt.savefig(output + ".png")
	plt.clf()

def per_sample_histogram(gene_list,gene_transcripts,tpm_dict,threshold,output):
	"""Number of genes in each sample that are above the theshold of expression"""
	number_of_samples = []
	for gene in gene_list:
		gene_tpms = [0] * 20
		for transcript in gene_transcripts[gene]:
			tpms = tpm_dict[transcript]
			gene_tpms = [gene_tpms[i] + tpms[i] for i in range(0,len(tpms))] #Convert transcript tpms to gene tpms
		for i in range(0,len(gene_tpms)):
			if gene_tpms[i] >= threshold:
				number_of_samples.append(i + 1)
	plt.figure()
	plt.hist(number_of_samples,bins=20)
	plt.savefig(output + ".png")
	plt.clf()
				





def main():
	parser = argparse.ArgumentParser(description='Filter GFF file based on coordinates')
	parser.add_argument('-gtf', dest = 'gtf_input', type = str, help = 'input file')
	parser.add_argument('-t', dest = 'threshold', type = int, help = 'This is the threshold under which a gene is not considered expressed',default=1)
	parser.add_argument('-tpm', dest = 'tpm_input', type = str, help = 'input file')
	parser.add_argument('-o', dest = 'output',type = str, help = 'output file name')
	parser.add_argument('-new', dest = 'new_transcript_names',type = str, help = '', default="")
	parser.add_argument('-a', dest = 'antisense',type = str, help = '', default="")
	parser.add_argument('-i', dest = 'intronic',type = str, help = '', default="")
	parser.add_argument('-nc', dest = 'non_coding',type = str, help = '', default="")
	parser.add_argument('--f', dest = 'filter', help = 'Filter transcriptome based on criteria', action='store_true')
	args = parser.parse_args()
	for line in open(args.gtf_input):
		gtf = Gtf(line)
	gtf.summarise()

	tpm_dict = get_transcript_tpms(args.tpm_input)
	gene_transcripts = get_genes(tpm_dict)
	monoexons = gtf.getmonoexons()

	if args.new_transcript_names:
		new_ids = get_new_transcript_names(args.new_transcript_names)
	if args.antisense:
		monoexon_antisense = get_antisense(new_ids,monoexons,args.antisense)
	else:
		monoexon_antisense = None
	if args.intronic:
		monoexon_intronic = get_intronic(new_ids,monoexons,args.intronic)
	else:
		monoexon_intronic = None
	if args.non_coding:
		monoexon_non_coding = get_non_coding(new_ids,monoexons,args.non_coding)
	else:
		monoexon_non_coding = None

	mono_exon_tpm_histogram(tpm_dict, monoexons, gene_transcripts, args.threshold, args.output, monoexon_intronic, monoexon_antisense, monoexon_non_coding)

	#Now filter
	if args.filter:
		genes_for_filter = monoexons #Keep it simple
		
		genes_for_removal = filter_genes_by_tpm(genes_for_filter, gene_transcripts, tpm_dict, args.threshold, monoexon_antisense,monoexon_intronic, monoexon_non_coding)
		print(f"Number of genes that will be removed: {len(genes_for_removal)}")

		gtf.filter(genes_for_removal, args.output + ".gtf", args.output + "_removed.gtf")

if __name__ == "__main__":
	main()
	









