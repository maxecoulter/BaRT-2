

import argparse

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
	def sequence_length(self):#Generate sequence length for each transcript
		length = 0
		for exon_length in self.list_exon_lengths:
			length = length + int(exon_length)
		return length
	def sjcoordinates(self):#Create sj coordinate for each sj in transcript
		list_sj_coordinates = []
		for i in range(0,int(self.exon_no)-1):#Go through each intron in .bed file (0 - n)
			intron_start=int(self.start)+int(self.list_exon_starts[i])+int(self.list_exon_lengths[i])#Extract start of intron from .bed file
			intron_end=int(self.start)+int(self.list_exon_starts[i+1])#Extract end of intron from .bed file
			list_sj_coordinates.append(self.chromosome+"_"+str(intron_start)+"_"+str(intron_end)+"_"+self.strand)
		return list_sj_coordinates



parser = argparse.ArgumentParser(description='')
parser.add_argument('-di', dest = 'duplicates_input', type = str,)
parser.add_argument('-b', dest = 'bedinput', type = str,)
parser.add_argument('-o', dest = 'bedoutput',type = str,)
parser.add_argument('-f', dest = 'fasta_input',type = str,)
parser.add_argument('-fd', dest = 'fasta_input_with_duplicates',type = str,)

args = parser.parse_args()

fasta1 = set([line for line in open(args.fasta_input) if line.startswith(">")])
fasta_withduplicates = set([line for line in open(args.fasta_input_with_duplicates) if line.startswith(">")])

fasta_transcripts = set([line.split(";")[1].split("(")[0] for line in fasta1])

duplicate_transcripts1 = set([line.split(";")[1].split("(")[0] for line in list(fasta_withduplicates - fasta1)])
duplicate_transcripts = set()
for line in open(args.duplicates_input):
	duplicates = line.rstrip("\n").split("\t")[-1].split(",")[1:]
	for duplicate in duplicates:
		transcript = duplicate.split(";")[1].split("(")[0]
		duplicate_transcripts.add(transcript)

print(f"It is {duplicate_transcripts1 == duplicate_transcripts} that both sets are the same")

print(duplicate_transcripts)
print(f"Actual duplicate transcripts: {duplicate_transcripts1}")

bedoutput = open(args.bedoutput,"w")
bed_transcripts = set()
for line in open(args.bedinput):
	bed = Bedline(line)
	if bed.transcript not in duplicate_transcripts1:
		bed_transcripts.add(bed.transcript)
		bedoutput.write(line)

bedoutput.close()
try:
	assert fasta_transcripts == bed_transcripts
except AssertionError:
	print(f"Assertion error! The following transcripts are in bed that should not be: {bed_transcripts - fasta_transcripts}")