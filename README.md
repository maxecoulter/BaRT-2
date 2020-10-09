# BaRTv2.0


Author: Max Coulter

Contact: mecoulter@dundee.ac.uk

These scripts were used for filtering the deep Iso-Seq data for construction of BaRTv2.0,  a barley cv. Barke reference transcriptome. Below is a brief description of how these scripts can be run.

The filtering is based on the output from TAMA collapse and TAMA merge (https://github.com/GenomeRIK/tama/wiki), a program used to create gene and transcript coordinates from Iso-Seq mappings. The outputs from these programs are required for running of filtering.

The filtering of the Iso-Seq has two modules: **BaRT_generate_filter_information.py** and **BaRT_2_filter_binomial.py**.

### BaRT_generate_filter_information.py

This is the first module that needs to be run. It needs to be run on each sample seperately, and produces the inputs required for **BaRT_2_filter_binomial.py**. 

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
  
  e) Column 5, unique splice junction location id, in format <chromosome>_<start>_<end>_<strand>;
  
  f) Column 6, The donor and acceptor dinucleotides for the splice junction;
  
  g) Column 7, Read map position (start);
  
  h) Column 8, Read map position (end)
  
  
  2. **\<prefix>\_single_exon_reads.txt**
  
  This is a tab delimited text file in this format:
  
      m54203_191129_023846/5243497/ccs	G8.1	chr1H	205480	206287	-
    m54203_191129_023846/53609273/ccs	G26.40	chr1H	2898678	2898891	+
    m54203_191129_023846/51118505/ccs	G27.1	chr1H	2913273	2913832	+
    m54203_191129_023846/22610633/ccs	G36.1	chr1H	4336426	4337540	-
    m54203_191129_023846/39321962/ccs	G37.1	chr1H	4597766	4599137	-
    m54203_191129_023846/24576989/ccs	G37.2	chr1H	4597768	4599133	-
    m54203_191129_023846/29163745/ccs	G42.4	chr1H	5045380	5046038	+
    m54203_191129_023846/30867754/ccs	G46.4	chr1H	5095629	5098167	-
  
  
  
  
  
  



