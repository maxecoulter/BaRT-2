
"""Script for generating table of splice junctions and monoexon reads 

This script is used for generating the above from TAMA collapse output files, outputting a table of splice junctions with information on read, position (chromosome, start end), strand and donor/acceptor sequences. 

A table of monoexon reads are also generated, with information on chromosome, start end, strand.

If you have multiple sequence libraries, both tables from all libraries will need to be joined using the cat command to make one of each table for next filtering step.

Requires just one argument to run, which is the file prefix used for TAMA collapse (including the path)."""

import optparse
from optparse import OptionParser



genome_file = "/mnt/shared/projects/barley/201903_RTD2/BarkeAssembly/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta"

def density_error_parser(fileprefix,readdict,posdict,sj_no_dict,sj_error_profile_dict,linedict,dictofmultimaps,single_exon_dict):
	"""Parses the density error file output from TAMA collapse.
	
	Parameters
	----------
	fileprefix : str
		The file prefix from commandline argument
	readdict : dict
		After parsing, readkey is key, read id is value
	posdict : dict
		After parsing, readkey is key, value is list [<read start>,<read end>]
	sj_no_dict : dict
		After parsing, readkey is key, sj no is value
	sj_error_profile_dict : dict
		After parsing, readkey is key, value is profile of errors 30 nt either side of splice junction
	linedict : dict
		After parsing, readkey is key, value is line of density error file
	dictofmultimaps : dict
		After parsing, readkey is key, value is line of density error file where multimap has occured
	single_exon_dict : dict
		After parsing, read id is key, value is list [chromosome, start, end,strand]
	
	For the above, readkey = <read id>_<sj#>
	"""

	for line in open(fileprefix+"new_flnc_collapsed_local_density_error.txt"):
		if line.startswith("cluster_id"):
			continue
		line1 = line.split("\t")
		if line1[1] == "lde_fail":
			continue
		splicesplit = line1[11].split(";") #Column 11 is sj error profile
		count = 0
		if line1[6] == "1":#Column 6 is number of exons
			if line1[0] in single_exon_dict.keys():
				print("single exon read " + line1[0] + " already has map position")#This is in case of multimapping, just take first position, in BaRT second is always chrUn
				continue
			else:
				single_exon_dict[line1[0]] = [line1[2], line1[3], line1[4], line1[5]]#Chromosome, start, end,strand
				continue
		for sp in splicesplit:
			if sp == "na":
				pass
			else:
				if line1[0] == "m54203_190621_052144/65012252/ccs" or line1[0] == "m54203_190621_052144/26280086/ccs":
					print(str(line1))
				sjkey = line1[0] + "_" + str(count)#Key for all future dictionaries, format <read id>_<sj#>
				if sjkey in linedict.keys():
					dictofmultimaps[sjkey] = line1
					if line1[2] == "chrUn":
						pass
					else:
						print(str(line1))
					#listofmultimaps.append(sjkey)
				else:
					readdict[sjkey] = line1[0]#link key to read
					linedict[sjkey] = str(line1)#link key to whole line
					sj_no_dict[sjkey] = count#integer sj count dict
					sj_error_profile_dict[sjkey] = sp #sj error profile dict
					posdict[sjkey] = [line1[3],line1[4],]#Dictionary with readkey as keys and start and ends recorded in list [0] is start and [1] is end
					count += 1

def generate_sj_coordinate_dict(readdict,transiddict,fail_list,sj_errorlist,transcriptnalist,sj_coordinate_dict):
	"""Link readkey to splice junction coordinate.

	Parameters
	----------
	readdict : dict
		readkey is key, read id is value
	transiddict : dict
		read id is key, transcript id is value
	fail_list : list
		List of reads with no transcript id
	sj_errorlist : list
		List of splice junction keys (<transcript id>_<splice junction #>) without splice junction coordinate
	transcriptnalist : list
		List of transcript ids that are nas
	sj_coordinate_dict : dict
		<transcript id>_<splice junction #> is key, splice junction coordinate is value 
	
	Returns
	-------
	dict
		key is readkey, value is splice junction coordinate
	"""
	
	final_sj_coordinate_dict = {}
	for key in readdict.keys():
		read = readdict[key]
		#print(read)
		try:
			transcript = transiddict[read]
			if transcript == "na":
				transcriptnalist.append(key)#Should be 0
			else:
				pass
		except KeyError:
			#print("no transcript id")
			#transcript="na"
			#transiddict[read]="na"
			fail_list.append(key)#9715 fails - reads not in final transcript list 
			continue
		#print(transcript)
		keysplit = key.split("_")
		number = keysplit[-1]
		#if transcript=="G16147.1":
			#print(key)
		if int(number) >= 100:
			print("number error!")
			break
		else:
			#print(str(number))
			sj_coordinate_key = transcript + "_" + str(number)
			try:
				sj_coordinate = sj_coordinate_dict[sj_coordinate_key]
				final_sj_coordinate_dict[key] = sj_coordinate #Final length 1397819, 1388104 once na's are removed
			except KeyError:
				#print("no sj coordinate")
				sj_coordinate = "na"
				print("Error: No splice junction coordinates identified for " + key)
				print("read not used")
				sj_errorlist.append(sj_coordinate_key)#Mostly na_<number>
	return final_sj_coordinate_dict

def make_sj_sequence_dict(final_sj_coordinate_dict,chromosomedict,sj_fails):
	"""Create dictionary of donor acceptor pairs.
	
	Parameters
	----------
	final_sj_coordinate_dict : dict
		key is readkey, value is splice junction coordinate
	chromosomedict : dict
		key is chromosome, value is chromosome sequence
	sj_fails : list
		Becomes a list of readkeys with no splice junction coordinates
	
	Returns
	-------
	dict
		key is readkey, value is donor acceptor pairs
	"""
	
	sj_sequence_dict = {}
	for key in final_sj_coordinate_dict.keys():
		sj_coordinatelist = final_sj_coordinate_dict[key].split("_")
		if sj_coordinatelist[0] == "na":
			print("na error!")
			#print(key)
			print(str(sj_coordinatelist))
			sj_fails.append(key)
			#break
			#pairs="na"
			#sj_sequence_dict[key]=pairs
		else:
			chromosome = sj_coordinatelist[0]
			start = int(sj_coordinatelist[1])
			end = int(sj_coordinatelist[2])
			startpair = chromosomedict[chromosome][start:start + 2]
			endpair = chromosomedict[chromosome][end - 2:end]
			pairs = startpair + endpair
			sj_sequence_dict[key] = pairs
	return sj_sequence_dict

def parse_genome(genome_file):
	"""Parses genome into dictionary
	
	Parameters
	----------
	genome_file : str
		path to genome file
		
	Returns
	-------
	dict
		key is chromosome, value is chromosome sequence"""
	
	genome = open(genome_file).read().replace("\n","")
	#Split genome into chromosome dictionary with key = chr id and item chromosome sequence
	chromosomes = genome.split(">")
	del chromosomes[0]
	chromosomedict = {}
	for chromosome in chromosomes:
		chromid = chromosome[0:5]
		chromsequence = chromosome[5:]
		chromosomedict[chromid] = chromsequence
	return chromosomedict

def parse_trans_read(fileprefix,transiddict,transposdict,transmulti):
	"""Parse new_flnc_collapsed_trans_read.bed output from TAMA, to provide transcript information for reads.
	
	Parameters
	----------
	fileprefix : str
		The file prefix from commandline argument
	transiddict : dict
		After parsing, read id is key, transcript id is value
	transposdict : dict
		After parsing, transcript id is key, list [start position,end position] is value
	transmulti : dict
		After parsing, read id with multimapping is key, list of transcripts is value
	"""

	for line in open(fileprefix+"new_flnc_collapsed_trans_read.bed"):
		line1 = line.split("\t")
		startpos = line1[1]
		endpos = line1[2]
		trans_read = line1[3].split(";")	#trans id and read in list
		#if trans_read[1]=="m54203_190621_052144/65012252/ccs" or trans_read[1]=="m54203_190621_052144/26280086/ccs":
			#print(trans_read[0])
		#else:
			#pass
		transposdict[trans_read[0]] = [startpos,endpos]
		if trans_read[1] in transiddict.keys():
			transmulti.setdefault(trans_read[1],list()).append(trans_read[0])
		else:
			transiddict[trans_read[1]] = trans_read[0]#Key = read, item = trans id

def single_exon_transcript(single_exon_dict,transiddict):
	"""Give single exon reads a transcript.
	
	Parameters
	----------
	single_exon_dict : dict
		read id is key, value is list [chromosome,start,end,strand]
	transiddict : dict
		read id is key, transcript id is value
	
	Returns
	-------
	dict
		a dictionary of single exon reads, where read id is key, transcript id is value
	"""

	single_exon_reads_transcripts = {}
	for read in single_exon_dict.keys():
		if read in single_exon_reads_transcripts.keys():
			print("Read " + read + " has multimapping!")
			continue
		else:
			single_exon_reads_transcripts[read] = transiddict[read]
	return single_exon_reads_transcripts

def sjs_from_collapse(fileprefix,list_sj_coordinates,list_transcripts,sj_coordinate_dict):
	"""Generate splice junction coordinates for transcripts from collapse.bed file.
	
	Parameters
	----------
	fileprefix : str
		The file prefix from commandline argument
	list_sj_coordinates : list
		After parsing, list of all sj coordinates in format <start>_<end>_<strand>
	list_transcripts : list
		After parsing, List of transcript ids
	sj_coordinate_dict : dict
		After parsing, <transcript id>_<splice junction #> is key, splice junction coordinate is value 
	
	splice junction coordinate is in format <chromosome>_<intron start>_<intron end>_<strand>
	"""

	for line in open(fileprefix + "new_flnc_collapsed.bed"):
		line1 = line.replace("\n","").split("\t")
		chromosome = line1[0]
		transcriptinfo = line1[3].split(";")#transcript id is transcriptinfo[1]
		#if transcriptinfo[1]=="G16147.1" or transcriptinfo[1]=="G233.3_0":
			#print(str(line1))
		strand = line1[5]
		start = line1[6]
		#end = line1[7]
		exon_no = line1[9]
		list_exon_lengths = line1[10].split(",")
		list_exon_starts = line1[11].split(",")
		for i in range(0,int(exon_no)-1):
			intron_start = int(start) + int(list_exon_starts[i]) + int(list_exon_lengths[i])#Start of intron. May not need +1
			intron_end = int(start) + int(list_exon_starts[i + 1])
			sjnumber = i
			list_sj_coordinates.append(str(intron_start) + "_" + str(intron_end) + "_" + strand) 
			list_transcripts.append(transcriptinfo[1])
			sj_coordinate_dict[transcriptinfo[1] + "_" + str(sjnumber)] = chromosome + "_" + str(intron_start) + "_" + str(intron_end) + "_" + strand #link sj coordinates with read

def write_to_file(fileprefix,single_exon_reads_transcripts,single_exon_dict,final_sj_coordinate_dict,transiddict,sj_no_dict,sj_error_profile_dict,sj_sequence_dict,posdict,readdict):
	"""Write data to two output files, a table of splice junction data and a table of single exon read information.

	Parameters
	----------
	fileprefix : str
		The file prefix from commandline argument
	single_exon_reads_transcripts : dict
		read id is key, transcript id is value
	single_exon_dict : dict
		read id is key, value is list [chromosome, start, end,strand]
	final_sj_coordinate_dict : dict
		key is readkey, value is splice junction coordinate
	transiddict : dict
		read id is key, transcript id is value
	sj_no_dict : dict
		readkey is key, sj no is value
	sj_error_profile_dict : dict
		readkey is key, value is profile of errors 30 nt either side of splice junction
	sj_sequence_dict : dict
		key is readkey, value is donor acceptor pairs
	posdict : dict
		readkey is key, value is list [<read start>,<read end>]
	readdict : dict
		readkey is key, read id is value
	
	Format of output files:
	Splicejuntion table columns: read, transcript id, splice junction number, splice junction errors,splice junction coordinate, splice junction donor acceptor pairs, read start, read end

	Single exon table columns: read, transcript id, chromosome, read start, read end, read strand
	"""

	outfile = open(fileprefix + "_splice_junction_table.txt","w")#Order: readid,transid,sjno,error_profile,sj_coordinates,sj_sequence, read start position, read end position
	outfile2 = open(fileprefix + "_single_exon_reads.txt","w")#Order: Read, transcript,start,end,strand
	for read in single_exon_reads_transcripts.keys():
		outfile2.write(read + "\t" + single_exon_reads_transcripts[read] + "\t" + single_exon_dict[read][0] + "\t" + single_exon_dict[read][1] + "\t" + single_exon_dict[read][2] + "\t" + single_exon_dict[read][3] + "\n")
	for key in final_sj_coordinate_dict.keys():
		outfile.write(str(readdict[key]) + "\t" + str(transiddict[readdict[key]]) + "\t" + str(sj_no_dict[key]) + "\t" + str(sj_error_profile_dict[key]) + "\t" + str(final_sj_coordinate_dict[key]) + "\t" + str(sj_sequence_dict[key]) + "\t" + str(posdict[key][0]) + "\t" + str(posdict[key][1]) + "\n")
	outfile.close()
	outfile2.close()

def main(genome_file):
	"""Run the program.
	
	Parameters
	----------
	fileprefix : str
		The file prefix from commandline argument
	genome_file : str
		path to genome file
	"""
	parser = OptionParser()
	args = parser.parse_args()
	fileprefix = args[0]

	readdict = {}
	posdict = {}#Dictionary with readkey as keys and start and ends recorded in list [0] is start and [1] is end
	sj_no_dict = {}
	sj_error_profile_dict = {}
	linedict = {}
	dictofmultimaps = {}
	single_exon_dict = {}#Read is key with start and end as list
	print("Parsing density error file to find splice junctions")
	density_error_parser(fileprefix,readdict,posdict,sj_no_dict,sj_error_profile_dict,linedict,dictofmultimaps,single_exon_dict)

	#245 reads multimapping, all to chrUn
	#Now create transcript id dictionary
	transiddict = {}
	transposdict = {}
	transmulti = {}
	print("Parsing trans_read.bed to link reads to transcripts")
	parse_trans_read(fileprefix,transiddict,transposdict,transmulti)

	single_exon_reads_transcripts = single_exon_transcript(single_exon_dict,transiddict)

	print("Generating splice junction coordinates from collapse.bed file")
	list_sj_coordinates = [] 
	list_transcripts = []
	sj_coordinate_dict = {}
	#List of transcript ids (G1.1 etc.) that match sj coordinates
	#Now extract sj coordinates from .bed file
	sjs_from_collapse(fileprefix,list_sj_coordinates,list_transcripts,sj_coordinate_dict)

	print("linking main key to sj coordinates")
	fail_list = []
	sj_errorlist = []
	transcriptnalist = []
	final_sj_coordinate_dict = generate_sj_coordinate_dict(readdict,transiddict,fail_list,sj_errorlist,transcriptnalist,sj_coordinate_dict)

	print("Opening genome")
	chromosomedict = parse_genome(genome_file)

	print("obtaining splice junction sequences")
	sj_fails = []
	sj_sequence_dict = make_sj_sequence_dict(final_sj_coordinate_dict,chromosomedict,sj_fails)

	print("Writing to file")
	write_to_file(fileprefix,single_exon_reads_transcripts,single_exon_dict,final_sj_coordinate_dict,transiddict,sj_no_dict,sj_error_profile_dict,sj_sequence_dict,posdict,readdict)

if __name__ == "__main__":
	main(genome_file)







