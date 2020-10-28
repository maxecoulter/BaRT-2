
#Notes on merge_merge format:
#Column 3 (0 based) merged transcript name ; input transcript name (which can be got by line[3].split(";")[1].split("_")[-1]
#chr1H   72168   73280   G1.1;Splicejunction/20_samples_NoNs_sj_st_end_18_52_newmethod4_correctedp_lowexpressedgene_rescuemostsjs_merged_filtered.bed_G1.20	40      +       72168   73280   255,0,0 1       1112    0

#gene as key, with transcript objects as set
#each transcript object should have atributes: how many transcripts support transcript, and from long, short, or both
#Need to find out for each gene whether it has transcripts from, long, short or both
#If it has mutliple long read based transcripts, discard short read transcripts
#If it only has short read transcripts, keep.
# If it only has 1 long read transcript, check if gene is in low confidence list

#So first parse merge_merged file (bed112 format)
import argparse
import statistics
low_confidence_coverage_threshold = 0.8
#merge_input = "/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT_merge.txt" #The merge_merge info fle
#low_confidence_gene_list = "/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/Splicejunction/20_samples_NoNs_sj_st_end_18_52_newmethod4_correctedp_lowexpressedgene_rescuemostsjs_test_low_confidence_genes.txt"
#bed_input = "/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT.bed" #The tanscriptome to be filtered
#output = "/mnt/shared/scratch/mc42302/201903_RTD2/merged_RTD/merged_BaRT_filtered.bed"


class Bedline:
	def __init__(self,line,merge):
		line1 = line.split("\t")
		#line1 = line.rstrip("\n").split("\t")
		self.chromosome = line1[0]
		transcriptinfo = line1[3].split(";")#original source and transcript id is transcriptinfo[1]
		self.transcript = transcriptinfo[0]
		self.bed_transcript = transcriptinfo[1]#Using same class for different file types. Use this attribute for transcript id if file is a proper RTD bed file
		self.gene = self.transcript.split(".")[0]
		if merge:
			self.source = transcriptinfo[1].split("_")[0]
			self.source_transcript = "_".join(transcriptinfo[1].split("_")[1:])
			self.source_gene = self.source_transcript.split(".")[0]
		self.start = int(line1[1])
		self.end = int(line1[2])
		self.strand = line1[5]
		self.exon_no = int(line1[9])
		self.list_exon_lengths = line1[10].split(",")
		self.list_exon_starts = line1[11].split(",")
	def __repr__(self):
		return f'{self.transcript} from source {self.source_transcript}'
	def process(self,gene_dict):#gene_dict: gene is key, [0] is list of long read transcript objects and [1] is list of short read transcript objects
		if self.source == "Iso4": #if transcript is from IsoSeq library only
			gene_dict.setdefault(self.gene,[[],[]])[0].append(self)
		elif self.source == "Isolowconfidence":
			gene_dict.setdefault(self.gene,[[],[]])[0].append(self)
		elif self.source == "Illumina2":
			gene_dict.setdefault(self.gene,[[],[]])[1].append(self)
		else:
			print("Error lib fail")
	def sjcoordinates(self):#Create sj coordinate for each sj in transcript
		list_sj_coordinates = []
		for i in range(0,int(self.exon_no) - 1):#Go through each intron in .bed file (0 - n)
			intron_start = int(self.start) + int(self.list_exon_starts[i]) + int(self.list_exon_lengths[i])#Extract start of intron from .bed file
			intron_end = int(self.start) + int(self.list_exon_starts[i + 1])#Extract end of intron from .bed file
			list_sj_coordinates.append(self.chromosome + "_" + str(intron_start) + "_" + str(intron_end) + "_" + self.strand)
		return list_sj_coordinates
	def sequence_length(self):#Generate sequence length for each transcript
		length = 0
		for exon_length in self.list_exon_lengths:
			length = length + int(exon_length)
		return length

class Gtfline:
	coding_transcripts = set()
	length_dict = {} #CDS length
	def __init__(self,line):
		line1 = line.rstrip("\n").split("\t")
		self.chromosome = line1[0]
		self.assembler = line1[1]
		self.type = line1[2]
		self.start = int(line1[3])
		self.end = int(line1[4])
		self.score = line1[5]
		self.strand = line1[6]
		self.information = line1[8]
	def process(self):
		if self.type == "CDS":
			self.transcript_id = self.information.split(";")[1].replace('"','').split(" ")[-1]
			self.coding_transcripts.add(self.transcript_id)
			self.length_dict.setdefault(self.transcript_id,[]).append(self.end - self.start)


#gtfinput = "/mnt/shared/scratch/je42879/BaRTD_03July20_annotation_only/BaRT_July20_TS_transfix.gtf"




def filter_transcripts(gene_dict,low_confidence_coverage_threshold,monoexon,illumina_coding_transcripts,length_dict,remove_chimeras,iso_sjs):
	hc_transcript_set = set()
	removed_illumina = set()
	removed_pacbio = set()
	kept_low_confidence = set()
	CDS_lengths = []
	for gene in gene_dict.keys():
		if gene == "G9252":
			print(f'gene {gene} has a problem')
			print(str(gene_dict[gene]))
			for transcript_object in gene_dict[gene][0]:
				print(f'{transcript_object} is from {transcript_object.source}')
		chimera = False
		if len(gene_dict[gene][0]) > 0:#First scenario:If Pacbio transcript, 
			source_genes = set([transcript.source_gene for transcript in gene_dict[gene][0]])
			sources = {transcript.source_gene:transcript.source for transcript in gene_dict[gene][0]}
			sources_hc = set([source_gene for source_gene in sources.keys() if sources[source_gene] == "Iso4"])
			if remove_chimeras:
				if len(sources_hc) >= 3: #If 3 or more hc ISo-Seq transcripts, solve chimera by removing Illumina
					chimera = True
			if len(gene_dict[gene][1]) > 0:
				longest_short = 0
				mean_short = statistics.mean([transcript_object.sequence_length() for transcript_object in gene_dict[gene][1]])
				for transcript_object in gene_dict[gene][1]:#First get longest illumina transcript
					#hc_transcript_set.add(transcript_object.transcript) #If low confidence, add illumina transcript
					if transcript_object.sequence_length() > longest_short:
						longest_short = transcript_object.sequence_length()
				for transcript_object in gene_dict[gene][0]:
					source_gene = transcript_object.source_transcript.split(".")[0]
					source = sources[source_gene]
					if source == "Isolowconfidence":
						if transcript_object.sequence_length()/mean_short >= low_confidence_coverage_threshold:
							hc_transcript_set.add(transcript_object.transcript)# Take lc iso-seq transcripts that are > threshold length of longest illumina transcript
							kept_low_confidence.add(transcript_object.transcript)
						else:
							removed_pacbio.add(transcript_object.transcript)
					else:
						hc_transcript_set.add(transcript_object.transcript)
				for transcript_object in gene_dict[gene][1]:
					if "Iso4" not in set([sources[source_gene] for source_gene in source_genes]):
						hc_transcript_set.add(transcript_object.transcript)
					else:
						if not set(transcript_object.sjcoordinates()) <= iso_sjs and not chimera: #If short read sjs are not subset of long read sjs:
							hc_transcript_set.add(transcript_object.transcript)
							#sj_positions_short = sorted([int(sj.split("_")[-3]) for sj in transcript_object.sjcoordinates()] + [int(sj.split("_")[-2]) for sj in transcript_object.sjcoordinates()])
							#if sj_positions_short[0] > sj_positions_long[0] and sj_positions_short[-1] < sj_positions_long[-1]: #only keep illumina transcript if within bounds of hc isoseq transcript but has novel sjs
								#hc_transcript_set.add(transcript_object.transcript)
							#else:
								#removed_illumina.add(transcript_object.transcript)
						else:
							removed_illumina.add(transcript_object.transcript)
			else:
				if not set(sources.values()) == {"Isolowconfidence"}:#If gene is high confidence. If no Illumina it will just be one source gene
					for transcript_object in gene_dict[gene][0]:
						hc_transcript_set.add(transcript_object.transcript)
				else:
					for transcript_object in gene_dict[gene][0]:
						removed_pacbio.add(transcript_object.transcript)
		else:# No pacbio transcripts, just add short read transcripts
			for transcript_object in gene_dict[gene][1]:
				hc_transcript_set.add(transcript_object.transcript)
				#if set([transcript_object.exon_no for transcript_object in gene_dict[gene][1]]) == {1}:
					#if transcript_object.source_transcript in illumina_coding_transcripts:
						#CDS_lengths.append(max(list(length_dict[transcript_object.source_transcript])))
						#if max(list(length_dict[transcript_object.source_transcript])) == 89:
							#print(transcript_object.source_transcript)
						#print(f'Mono-exon coding transcript identified: {transcript_object.transcript}')
						#hc_transcript_set.add(transcript_object.transcript)
					#else:
						#removed_illumina.add(transcript_object.transcript)
				#else:
					#hc_transcript_set.add(transcript_object.transcript)
	#print(f'Minimum monoexon CDS Illumina transcript length: {min(CDS_lengths)}')
	return hc_transcript_set, removed_illumina, removed_pacbio, kept_low_confidence



"""
if len(gene_dict[gene][0]) > 1:#First scenario:If more than 1 Pacbio transcript, remove short reads, except if low confidence
	for transcript_object in gene_dict[gene][0]:
		hc_transcript_set.add(transcript_object.transcript)#Just add pacbio transcripts, don't add short read ones
	
elif len(gene_dict[gene][0]) == 0: #Second scenario: No Pacbio transcripts in gene. In this case, just use short read transcripts
	for transcript_object in gene_dict[gene][1]:
		hc_transcript_set.add(transcript_object.transcript)
elif len(gene_dict[gene][0]) == 1 and len(gene_dict[gene][1]) > 0: #Third scenario: 1 pacbio transcript, multiple short read transcripts
	source_gene = gene_dict[gene][0][0].source_transcript.split(".")[0] #The gene that the pacbio transcript originally came from
	for transcript_object in gene_dict[gene][1]:
		if set(transcript_object.sjcoordinates()) <= set(gene_dict[gene][0][0].sjcoordinates()): #See if every sj in short read transcript is in long read transcript
			hc_transcript_set.add(transcript_object.transcript)
	if source_gene not in low_confidence_genes: #If True, this is a gene that did not have good start/end support in pacbio data
		hc_transcript_set.add(gene_dict[gene][0][0].transcript)
elif len(gene_dict[gene][0]) == 1 and len(gene_dict[gene][1]) == 0: #Just 1 pacbio transcript
	hc_transcript_set.add(gene_dict[gene][0][0].transcript)
else:#There may be other criteria that need to be considered.
	print("Gene " + gene + " does not fit any criteria and so will be removed")"""

#def parse_low_confidence_genes(low_confidence_gene_list):
	#return set(open(low_confidence_gene_list).read().split("\n")[:-1])

def get_Illumina_coding_transcripts(gtfinput):
	for line in open(gtfinput):
		coding_info = Gtfline(line)
		coding_info.process()
	print(f'Number of Illumina transcripts with CDS: {len(coding_info.coding_transcripts)}')
	return coding_info.coding_transcripts, coding_info.length_dict

def write_bed(bedinput,hc_transcript_set,output):
	destination_path = "/".join(output.split("/")[:-1]) + "/"
	outfile = open(output, "w")
	transcripts_removed = open(destination_path + "transcripts_removed.bed", "w")
	for line in open(bedinput):
		bed = Bedline(line, False)
		if bed.bed_transcript in hc_transcript_set:
			outfile.write(line)
		else:
			transcripts_removed.write(line)
	outfile.close()

def main():
	#First parse arguments
	parser = argparse.ArgumentParser(description='Filter long and short read transcriptome')
	parser.add_argument('-m', dest = 'merge_input', type = str, help = 'The merge_merge input file')
	parser.add_argument('-b', dest = 'bed_input',type = str, help = 'The RTD bed file to be filtered')
	parser.add_argument('-monoexon', dest = 'monoexon',type = str, help = 'Keep mono-exon genes with only illumina support',default="T")
	parser.add_argument('-gtf', dest = 'gtfinput',type = str, help = 'The Illumina gtf file from Transfix')
	parser.add_argument('-rm_chimeras', dest = 'rm_chimeras',type = str, help = 'Whether to remove Illumina data from potentially chimeric genes on the bases of Iso-Seq data, y or n', default="n")
	parser.add_argument('-o', dest = 'output',type = str, help = 'Name of output file')
	args = parser.parse_args()
	remove_chimeras = True if args.rm_chimeras == "y" else False
	print("Parsing merge file...")
	gene_dict = {}
	iso_sjs = set()
	#Extract information from merge_merge file
	#mergeinput_corrected = open(args.merge_input).read().replace("\n_","M").split("\n")[:-1]
	#print(args.merge_input)
	#print(open(args.merge_input).readlines()[0])
	for n,line in enumerate(open(args.merge_input)):
		#print(str(n))
		#print(line)
		bed = Bedline(line,True)
		bed.process(gene_dict)
		if bed.source == "Iso4":
			iso_sjs.update(set(bed.sjcoordinates()))

	illumina_coding_transcripts, length_dict = get_Illumina_coding_transcripts(args.gtfinput)
	#Now filter transcripts based on certain criteria. Generate a high confidence set of transcripts
	print("Filtering transcripts...")
	hc_transcript_set, removed_illumina, removed_pacbio, kept_low_confidence = filter_transcripts(gene_dict,low_confidence_coverage_threshold,args.monoexon,illumina_coding_transcripts,length_dict,remove_chimeras,iso_sjs)
	#Use this set to filter RTD .bed file
	print("writing to file...")
	write_bed(args.bed_input,hc_transcript_set,args.output)
	print(f'Number of removed illumina transcripts: {len(removed_illumina)}')
	print(f'Number of removed pacbio transcripts: {len(removed_pacbio)}')
	print(f'Number of kept pacbio low confidence transcripts: {len(kept_low_confidence)}')
	print(f'Number of kept pacbio low confidence genes: {len(set([t.split(".")[0] for t in kept_low_confidence]))}')
	




if __name__ == "__main__":
	main()
