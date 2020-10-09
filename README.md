# BaRTv2.0


Author: Max Coulter

Contact: mecoulter@dundee.ac.uk

These scripts were used for filtering the deep Iso-Seq data for construction of BaRTv2.0,  a barley cv. Barke reference transcriptome. Below is a brief description of how these scripts can be run.

The filtering is based on the output from TAMA collapse and TAMA merge (https://github.com/GenomeRIK/tama/wiki), a program used to create gene and transcript coordinates from Iso-Seq mappings. The outputs from these programs are required for running of filtering.

The filtering of the Iso-Seq has two modules: **BaRT_generate_filter_information.py** and **BaRT_2_filter_binomial.py**.

### BaRT_generate_filter_information.py

This is the first module that needs to be run. It needs to be run on each sample seperately, and produces the inputs required for **BaRT_2_filter_binomial.py**. 


