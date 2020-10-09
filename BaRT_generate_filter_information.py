
"""Script for generating table of splice junctions and monoexon reads 

This script is used for generating the above from TAMA collapse output files, outputting a table of splice junctions with information on read, position (chromosome, start end), strand and donor/acceptor sequences. 

A table of monoexon reads are also generated, with information on chromosome, start end, strand.

If you have multiple sequence libraries, both tables from all libraries will need to be joined using the cat command to make one of each table for next filtering step.

Requires just one argument to run, which is the file prefix used for TAMA collapse (including the path)."""

import optparse
from optparse import OptionParser
import sys



def check_both_lists(final_sj_coordinate_dict,single_exon_reads_transcripts,transcript_coordinate_dict,transiddict,posdict,single_exon_dict):
	"""Check if read is present in both lists. if it is, it is a hidden case of multimapping"""
	not_single = set() #these reads are not to be outputted as they already have a map position
	not_multi = set() #these keys are to be ignored as well
	extra_multimaps = set() #Combined list of all extra multimapped reads
	for key in final_sj_coordinate_dict.keys():
		read = "_".join(key.split("_")[:-1])
		if read in single_exon_reads_transcripts.keys():
			extra_multimaps.add(read)
			first_match = posdict[key]
			second_match = single_exon_dict[read]
			actual_match = transcript_coordinate_dict[transiddict[read]]
			if first_match == actual_match:
				not_single.add(read)
			elif second_match == actual_match:
				not_multi.add(key)
			else:
				print("This read/key has more than two matches or there is a fatal error!")
				print(read)
	return not_single, not_multi, extra_multimaps
			

def check_multimap(key,readdict_m,posdict,posdict_m,transcript_coordinate_dict,transiddict):
	"""Check if read is multimapped, and if it is whether it is the first or second match. IF True, the second match (in multimap dataset is correct"""
	if key in readdict_m.keys():
		read = "_".join(key.split("_")[:-1])
		first_match = posdict[key]
		second_match = posdict_m[key]
		actual_match = transcript_coordinate_dict[transiddict[read]]
		if first_match == actual_match:
			return False
		elif second_match == actual_match:
			return True
		else:
			print("This read has more than two matches or there is a fatal error!")
			print(read)
			#sys.exit()
	else:
		return False

def check_multimap_single(read,multimapped_reads_single,single_exon_dict,transcript_coordinate_dict,single_exon_reads_transcripts):
	if read in multimapped_reads_single.keys():
		first_info = single_exon_dict[read]
		second_info = multimapped_reads_single[read]
		actual_info = transcript_coordinate_dict[single_exon_reads_transcripts[read]]
		if first_info == actual_info:
			return False
		elif second_info == actual_info:
			return True
		else:
			print("This read has more than two matches or there is a fatal error!")
			print(read)
			#sys.exit()
	else:
		return False

def density_error_parser(fileprefix,readdict,readdict_m,posdict,posdict_m,sj_no_dict,sj_no_dict_m,sj_error_profile_dict,sj_error_profile_dict_m,linedict,linedict_m,dictofmultimaps,single_exon_dict):
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
	multimapped_reads_single = {}
	multimapped_reads_double = {}
	Ide_fail_reads = set()
	reads_na_sjs = set()
	for line in open(fileprefix + "new_flnc_collapsed_local_density_error.txt"):
		if line.startswith("cluster_id"):
			continue
		line1 = line.split("\t")
		if line1[1] == "lde_fail":
			Ide_fail_reads.add(line1[0] + "\n")
			continue
		splicesplit = line1[11].split(";") #Column 11 is sj error profile
		count = 0
		if line1[6] == "1":#Column 6 is number of exons
			if line1[0] in single_exon_dict.keys():
				multimapped_reads_single[line1[0]] = [line1[2], line1[3], line1[4], line1[5]]
				continue
			else:
				single_exon_dict[line1[0]] = [line1[2], line1[3], line1[4], line1[5]]#Chromosome, start, end,strand
				continue
		for sp in splicesplit:
			if sp == "na":
				reads_na_sjs.add(line1[0] + "\n")
			else:
				sjkey = line1[0] + "_" + str(count)#Key for all future dictionaries, format <read id>_<sj#>
				if sjkey in linedict.keys():
					dictofmultimaps[sjkey] = line1
					multimapped_reads_double[line1[0]] = [line1[2], line1[3], line1[4], line1[5]]
					readdict_m[sjkey] = line1[0]#link key to read
					linedict_m[sjkey] = str(line1)#link key to whole line
					sj_no_dict_m[sjkey] = count#integer sj count dict
					sj_error_profile_dict_m[sjkey] = sp #sj error profile dict
					posdict_m[sjkey] = [line1[2],line1[3],line1[4],line1[5]]#Dictionary with readkey as keys and start and ends recorded in list [0] is start and [1] is end
					count += 1
				else:
					readdict[sjkey] = line1[0]#link key to read
					linedict[sjkey] = str(line1)#link key to whole line
					sj_no_dict[sjkey] = count#integer sj count dict
					sj_error_profile_dict[sjkey] = sp #sj error profile dict
					posdict[sjkey] = [line1[2],line1[3],line1[4],line1[5]]#Dictionary with readkey as keys and start and ends recorded in list [0] is start and [1] is end
					count += 1
	return multimapped_reads_single, multimapped_reads_double, Ide_fail_reads, reads_na_sjs

def generate_sj_coordinate_dict(readdict,transiddict,fail_list,sj_errorlist,transcriptnalist,sj_coordinate_dict,readdict_m,single_exon_dict,posdict,transcript_coordinate_dict):
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
			print("no transcript id")
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
				if read not in single_exon_dict.keys() and key not in readdict_m.keys():#Some reads are multimapped, to both single exon and sj lists. This produces an error if the transcript it matches does not have sjs. If this isn't the case though,
					#if posdict[read + "_0"] == posdict[key]:#This arises if the second mapping has an extra sj
					if posdict[key] == transcript_coordinate_dict[transiddict[read]]: #If 
						print("Error: No splice junction coordinates identified for " + key)
						sj_errorlist.add(read)
	#Now for second multimap reads
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
		transposdict[trans_read[0]] = [startpos,endpos]
		if trans_read[1] in transiddict.keys():#This should not occur
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
	single_exon_reads_transcripts_multi = {}
	for read in single_exon_dict.keys():
		if read in single_exon_reads_transcripts.keys():
			print("Read " + read + " has multimapping!")
			single_exon_reads_transcripts_multi.setdefault(read,set()).add(transiddict[read])
			continue
		else:
			single_exon_reads_transcripts[read] = transiddict[read]
	return single_exon_reads_transcripts, single_exon_reads_transcripts_multi

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
	transcript_coordinate_dict = {}#transcript is key, [chr,start,end,strand] is value
	for line in open(fileprefix + "new_flnc_collapsed.bed"):
		line1 = line.replace("\n","").split("\t")
		chromosome = line1[0]
		transcriptinfo = line1[3].split(";")#transcript id is transcriptinfo[1]
		#if transcriptinfo[1]=="G16147.1" or transcriptinfo[1]=="G233.3_0":
			#print(str(line1))
		strand = line1[5]
		start = line1[6]
		end = line1[7]
		exon_no = line1[9]
		list_exon_lengths = line1[10].split(",")
		list_exon_starts = line1[11].split(",")
		transcript_coordinate_dict[transcriptinfo[1]] = [chromosome, str(int(start) + 1), str(int(end) + 1), strand]
		for i in range(0,int(exon_no)-1):
			intron_start = int(start) + int(list_exon_starts[i]) + int(list_exon_lengths[i])#Start of intron. May not need +1
			intron_end = int(start) + int(list_exon_starts[i + 1])
			sjnumber = i
			list_sj_coordinates.append(str(intron_start) + "_" + str(intron_end) + "_" + strand) 
			list_transcripts.append(transcriptinfo[1])
			sj_coordinate_dict[transcriptinfo[1] + "_" + str(sjnumber)] = chromosome + "_" + str(intron_start) + "_" + str(intron_end) + "_" + strand #link sj coordinates with read
	return transcript_coordinate_dict

def write_to_file(fileprefix,single_exon_reads_transcripts,single_exon_dict,final_sj_coordinate_dict,transiddict,sj_no_dict,sj_error_profile_dict,sj_sequence_dict,posdict,readdict,multimapped_reads_single,transcript_coordinate_dict,readdict_m,posdict_m,sj_no_dict_m,sj_error_profile_dict_m):
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
	#print("NUmber of single exon reads: " + str(len(single_exon_reads_transcripts.keys())))
	#print("Number of splice junctions: " + str(len(final_sj_coordinate_dict.keys())))
	not_single, not_multi, extra_multimaps = check_both_lists(final_sj_coordinate_dict,single_exon_reads_transcripts,transcript_coordinate_dict,transiddict,posdict,single_exon_dict)
	print("Not single: " + str(len(not_single)))
	print("Not multi keys: " + str(len(not_multi)))
	for read in single_exon_reads_transcripts.keys():
		if read in not_single:
			continue
		if check_multimap_single(read,multimapped_reads_single,single_exon_dict,transcript_coordinate_dict,single_exon_reads_transcripts):
			outfile2.write(read + "\t" + single_exon_reads_transcripts[read] + "\t" + multimapped_reads_single[read][0] + "\t" + multimapped_reads_single[read][1] + "\t" + multimapped_reads_single[read][2] + "\t" + multimapped_reads_single[read][3] + "\n")
		else:
			outfile2.write(read + "\t" + single_exon_reads_transcripts[read] + "\t" + single_exon_dict[read][0] + "\t" + single_exon_dict[read][1] + "\t" + single_exon_dict[read][2] + "\t" + single_exon_dict[read][3] + "\n")
	for key in final_sj_coordinate_dict.keys():
		if key in not_multi:
			continue
		if check_multimap(key,readdict_m,posdict,posdict_m,transcript_coordinate_dict,transiddict):
			outfile.write(str(readdict[key]) + "\t" + str(transiddict[readdict[key]]) + "\t" + str(sj_no_dict_m[key]) + "\t" + str(sj_error_profile_dict_m[key]) + "\t" + str(final_sj_coordinate_dict[key]) + "\t" + str(sj_sequence_dict[key]) + "\t" + str(posdict_m[key][1]) + "\t" + str(posdict_m[key][2]) + "\n")
		else:
			outfile.write(str(readdict[key]) + "\t" + str(transiddict[readdict[key]]) + "\t" + str(sj_no_dict[key]) + "\t" + str(sj_error_profile_dict[key]) + "\t" + str(final_sj_coordinate_dict[key]) + "\t" + str(sj_sequence_dict[key]) + "\t" + str(posdict[key][1]) + "\t" + str(posdict[key][2]) + "\n")
	outfile.close()
	outfile2.close()
	return extra_multimaps

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
	options, args = parser.parse_args()
	print(str(args))
	fileprefix = args[0]
	genome_file = args[1] #"/mnt/shared/scratch/mc42302/201903_RTD2/Benchmarking/180903_Barke_Unfiltered_chloro_clean_pseudomolecules_v1.fasta"
	print(fileprefix)

	readdict = {}
	posdict = {}#Dictionary with readkey as keys and start and ends recorded in list [0] is start and [1] is end
	sj_no_dict = {}
	sj_error_profile_dict = {}
	linedict = {}
	dictofmultimaps = {}
	single_exon_dict = {}#Read is key with start and end as list
	readdict_m = {}
	posdict_m = {}
	sj_no_dict_m = {}
	sj_error_profile_dict_m = {}
	linedict_m = {}
	print("Parsing density error file to find splice junctions")
	#
	multimapped_reads_single, multimapped_reads_double, Ide_fail_reads, reads_na_sjs = density_error_parser(fileprefix,readdict,readdict_m,posdict,posdict_m,sj_no_dict,sj_no_dict_m,sj_error_profile_dict,sj_error_profile_dict_m,linedict,linedict_m,dictofmultimaps,single_exon_dict)

	#245 reads multimapping, all to chrUn
	#Now create transcript id dictionary
	transiddict = {}
	transposdict = {}
	transmulti = {}
	print("Parsing trans_read.bed to link reads to transcripts")
	parse_trans_read(fileprefix,transiddict,transposdict,transmulti)
	assert not len(transmulti.keys())#should be 0 reads with multiple transcripts

	single_exon_reads_transcripts, single_exon_reads_transcripts_multi = single_exon_transcript(single_exon_dict,transiddict)
	assert not len(single_exon_reads_transcripts_multi.keys())#should be 0 reads with multiple transcripts

	print("Generating splice junction coordinates from collapse.bed file")
	list_sj_coordinates = [] 
	list_transcripts = []
	sj_coordinate_dict = {}
	#List of transcript ids (G1.1 etc.) that match sj coordinates
	#Now extract sj coordinates from .bed file
	transcript_coordinate_dict = sjs_from_collapse(fileprefix,list_sj_coordinates,list_transcripts,sj_coordinate_dict)

	print("linking main key to sj coordinates")
	fail_list = []
	sj_errorlist = set()
	transcriptnalist = []
	final_sj_coordinate_dict = generate_sj_coordinate_dict(readdict,transiddict,fail_list,sj_errorlist,transcriptnalist,sj_coordinate_dict,readdict_m,single_exon_dict,posdict,transcript_coordinate_dict)

	print("Opening genome")
	chromosomedict = parse_genome(genome_file)

	print("obtaining splice junction sequences")
	sj_fails = []
	sj_sequence_dict = make_sj_sequence_dict(final_sj_coordinate_dict,chromosomedict,sj_fails)

	print("Writing to file")
	extra_multimaps = write_to_file(fileprefix,single_exon_reads_transcripts,single_exon_dict,final_sj_coordinate_dict,transiddict,sj_no_dict,sj_error_profile_dict,sj_sequence_dict,posdict,readdict,multimapped_reads_single,transcript_coordinate_dict,readdict_m,posdict_m,sj_no_dict_m,sj_error_profile_dict_m)
	error_out = open(fileprefix + "_Ide_read_errors","w")
	multimap_out = open(fileprefix + "_multimapped_reads","w")
	sj_fails = open(fileprefix + "_reads_sjna","w")
	for read in Ide_fail_reads:
		error_out.write(read)
	for reads in multimapped_reads_single, multimapped_reads_double, extra_multimaps:
		for read in reads:
			multimap_out.write(read + "\n")
	for readkey in sj_errorlist:
		sj_fails.write(readkey + "\n")
	#for read in multimapped_reads:
		#multimap_out.write(read)
	error_out.close()
	multimap_out.close()
	sj_fails.close()

if __name__ == "__main__":
	main(genome_file)







