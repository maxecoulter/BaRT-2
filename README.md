# BaRTv2.0


Author: Max Coulter

Contact: mecoulter@dundee.ac.uk

These scripts were used for filtering the deep Iso-Seq data for construction of BaRTv2.0,  a barley cv. Barke reference transcriptome. Below is a brief description of how these scripts can be run. A detailed explanation of the methodology used in the scripts is described in the BaRTv2.0 paper (Coulter et al., unpublished).

The filtering is based on the output from TAMA collapse and TAMA merge (https://github.com/GenomeRIK/tama/wiki), a program used to create gene and transcript coordinates from Iso-Seq mappings. The outputs from these programs are required for running of filtering.

The filtering of the Iso-Seq has two modules: **BaRT_generate_filter_information.py** and **BaRT_2_filter_binomial.py**.

### Dependencies:
Both scripts for filtering Iso-seq require python >= 3.6. They also require the following:
- Pandas
- Plotnine
- matplotlib


### BaRT_generate_filter_information.py

This is the first module that needs to be run. It needs to be run on each Iso-seq library seperately, and produces the inputs required for **BaRT_2_filter_binomial.py**. Before running, you will need to run Tama collapse and Tama merge (see https://github.com/GenomeRIK/tama).

The program takes various tama collapse output files (<prefix>_local_density_error.txt, <prefix>_trans_read.bed, <prefix>_collapsed.bed) and parses information to generate information on transcript coordinates, read coordinates, which reads support which transcripts and splice junction information, including types of errors sorrounding splice junctions. These are then written into output tables that can be used by **BaRT_2_filter_binomial.py**.

#### Usage
Usage: python BaRT_generate_filter_information.py input genome file
  
**input** The path to the folder where tama collapse output is, and the file prefix used by tama, e.g /path/to/sample_folder/prefix

**genome file** The path and fasta file of reference genome. e.g /path/to/genome.fasta

#### Output

The following outputs are poduced:

1. **\<prefix>\_splice_junction_table.txt**
 
 This a tab delimited text file in this format:
 
    m54203_191129_023846/20775159/ccs	G1.2	0	______________________________>______________________________	chr1H_105062_105167_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	1	______________________________>_______________XX_____________	chr1H_105209_105713_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	2	______________________________>______________________________	chr1H_105781_105861_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	3	______________________________>______________________________	chr1H_105904_106005_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	4	______________________________>______________________________	chr1H_106037_106589_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	5	______________________________>______________________________	chr1H_106760_106938_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	6	______________________________>_________I____________________	chr1H_107116_107201_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	7	________________________I_____>______________________________	chr1H_107402_107509_-	CTAC	104372	107925
    m54203_191129_023846/20775159/ccs	G1.2	8	____________________________D_>______________________________	chr1H_107700_107776_-	CTAC	104372	107925

Where:

  a) Column 1, read id;

  b) Column 2, transcript id;

  c) Column 3, splice junction number;

  d) Column 4, splice junction error profile +/- 30 bp either side of the splice junction. The > symbol represents the location of the splice junction, with various characters indicating errors (for further details see **prefix_local_density_error.txt** output information in  https://github.com/GenomeRIK/tama/wiki/Tama-Collapse;
  
  e) Column 5, unique splice junction location id, in format \<chromosome>_<start>\_<end>_\<strand>\;
  
  f) Column 6, The donor and acceptor dinucleotides for the splice junction;
  
  g) Column 7, Read map position (start);
  
  h) Column 8, Read map position (end)
  
  
  2. **\<prefix>\_single_exon_reads.txt**
  
  This is a tab delimited text file in this format:
  
    m54203_191129_023846/5243497/ccs	G8.1	chr1H	205480	206287	-
    m54203_191129_023846/53609273/ccs	G26.40 chr1H	2898678	2898891	+
    m54203_191129_023846/51118505/ccs	G27.1	chr1H	2913273	2913832	+
    m54203_191129_023846/22610633/ccs	G36.1	chr1H	4336426	4337540	-
    m54203_191129_023846/39321962/ccs	G37.1	chr1H	4597766	4599137	-
    m54203_191129_023846/24576989/ccs	G37.2	chr1H	4597768	4599133	-
    m54203_191129_023846/29163745/ccs	G42.4	chr1H	5045380	5046038	+
    m54203_191129_023846/30867754/ccs	G46.4	chr1H	5095629	5098167	-
  
  
  a) Column 1, read id;
    
  b) Column 2, transcript id;
  
  c) Column 3, chromosome
  
  d) Read map position (start);
  
  e) Read map position (end)
    
  This is a dataset of all the reads without splice junctions.
  
  3. **\<prefix>\_Ide_read_errors**
  
  This is a list of all the reads that have an Ide failed flag (see https://github.com/GenomeRIK/tama/wiki/Tama-Collapse for further details).
  
  4. **\<prefix>\_multimapped_reads**
  
  This is a list of all reads that map to multiple locations. TAMA will only use one map position to support a transcript, and so the program will discard the alternatve mapping which is not used in the annotation.
  
  5. **\<prefix>\_reads_sjna**
  
  Any splice junctions with errors are reported here. This file should be empty, if not it means there is a problem with the mapping, or the input. The file is a list of splice junction ids (\<read id>_<splice junction number>\). 
  

### BaRT_2_filter_binomial.py

Before running **BaRT_2_filter_binomial.py** the **\<prefix>\_splice_junction_table.txt** outputs from all samples will need to be joined together. This can be done using the linux cat command, e.g:

    cat *_splice_junction_table.txt > all_splice_junction_table.txt

The same needs to be done for the **\<prefix>\_single_exon_reads.txt** file:

    cat *_single_exon_reads.txt > all_single_exon_reads.txt
    
These new tables will form part of the input for **BaRT_2_filter_binomial.py**. It is also assumed that different sample collapse.bed files have been merged using TAMA merge (see https://github.com/GenomeRIK/tama/wiki/Tama-Merge for further details).

Finally poly(A) information from TAMA collapse from each sample needs to be joined together in one file. This can be done as above:

    cat *_polya.txt > all_collapsed_polya.txt

#### Usage

python BaRT_2_filter_binomial.py [-i] [-bed] [-bedn] [-sr] [-g] [-pA] [-s] [-o]

**Optional arguments:**

    -i          splice junction table
    -bed        merged transcriptome .bed file without N containing transcripts (Optional)
    -bedn       merged transcriptome .bed file
    -sr         Short read input (Optional). Splice junction output from STAR mapping of short reads. Currently only works with Barke genome
    -g          Genome reference
    -pA         poly(A) information collated from tama collapse
    -s          single exon input
    -o          output prefix
    --hamming   For template switching,threshold hamming distance below which sj considered RT switching. FOr example 2 means a difference of 2 bases in 8 (default = 1)
    --polyA     Threshold for percentage of As at 3' end of gene, above which read is removed (default = 80)
    --st_window Size of the window for removing unsupported 5' ends (+/- n), (default = 20)
    --end_window Size of the window for removing unsupported 3' ends (+/- n), (default = 60)
    
  Default usage would look like this:
  
  python BaRT_2_filter_binomial.py -i all_splice_junction_table.txt -bedn merge.bed  -g genome.fasta -pA all_collapsed_polya.txt -s all_single_exon_reads.txt -o filtered_transcriptome1
  
  Detailed explanation of arguments:
  
  **-i all_splice_junction_table.txt** This input is the collated per sample output from **BaRT_generate_filter_information.py**. It is described in more detail above. This provides the program with splice junction, read and trnascript information.
  
  **-bedn merge.bed** This is the merged transcriptome from tama merge. 
  
  **-bed merge_noNs.bed** This is an optional argument. With the barley genome there are stretches of Ns caused by gaps in the genome. When reads overlap these regions this can cause downstream annotation problems. It is best to remove all reads overlapping Ns (using bedtools intersect) and use the resulting filtered bed file as an input here. PLEASE NOTE: You still need the original unfiltered merged transcriptome from tama merge as input if you have a file as input here.
  
  **-sr all_short_sjs.tab** This is an optional argument. This allows splice junctions from short reads to be used to support long read based splice junctions. Currently only works for the barley cv. Barke genome.
  
  **-g genome.fasta** The genome reference file used for mapping of Iso-Seq.
  
  **-pA all_collapsed_polya.txt** This is the collated poly(A) output from tama collapse, from \<prefix>\_polya.txt files.
  
  **-s all_single_exon_reads.txt** The collated single exon reads file from **BaRT_generate_filter_information.py**. This provides information on reads without splice junctions.
  
  **-o output_prefix** This is the prefix used in outfiles
  
#### Outputs

**BaRT_2_filter_binomial.py** produces the following outputs:

**\<prefix>\_merged_filtered.bed** This is the filtered .bed file output that is the main result of filtering process. In order to produce an updated Iso-Seq transcriptome, TAMA merge will need to be run again using this as input.

**\<prefix>\_removed.bed** This is a .bed file of transcripts that have been removed from the input .bed file in the filtering process.

**\<prefix>\_low_confidence.bed** This is a .bed file of transcripts from genes with low read support, and so are low confidence with regard to 5' and 3' ends. They will have high confidence splice junctions. In BaRTv2.0 these transcripts were included if they had support from Illumina transcripts.

**\<prefix>\_transcript_info** This is a tab delimited text file showing information on each transcript. The format is as follows:

    Transcript	PolyA Error	low expressed	rescued	high confidence	5' support	3' support	splice junctions supported
    G1.1	False	False	False	False	False	True	True
    G1.2	False	False	False	False	False	True	True
    G1.3	False	False	False	False	False	True	True
    G1.4	False	False	False	False	False	True	True
    G1.5	False	False	False	False	False	True	True
    G1.6	False	False	False	False	False	True	True
    G1.7	False	False	False	False	False	True	True
    G1.8	False	False	False	False	False	True	True
    G1.9	False	False	False	True	True	True	True
    G1.10	False	False	False	True	True	True	True
    
    
    
  a) Column 1, Transcript id;
  
  b) Column 2, True or False, whether trancript has PolyAs at 3' end in genome (possible mispriming). If true, transcript will have been removed;
  
  c) Column 3, True or False, whether transcript is from a low expressed gene;
  
  d) Column 4, True or False, whether the longest transcript with most splice junctions from low expressed genes with no support for 5' and 3' ends has been "rescued" (added to the main .bed file). This is currently not used and is not recommended, so all transcripts will be False here.
  
  e) Column 5, True or False, whether transcript is high confidence. If True, the transcript has good support for 5' and 3' ends;
  
  f) Column 6, True or False, whether transcript has support for 5' end, either by multiple reads within a window or by a high confidence TSS site;
  
  g) Column 7, True or False, whether transcript has support for 3' end, either by multiple reads within a window or by a high confidence TES site;
  
  h) Column 8, True or False, whether transcript has support from high confidence splice junctions. If all transcript splice junctions are present in the high confidence set of splice junctions, this will be True
  
  
  
  


  
  
  
  
  
    



  
  
  



