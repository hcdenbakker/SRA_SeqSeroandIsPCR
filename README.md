# SRA_SeqSeroandIsPCR
pipeline for download of Illumina SRA data, assembly and in silico PCR for primer testing in python

Script to (i) download data from SRA, (2) translate .sra files to paired end .fastq.gz,
run SeqSero (https://github.com/denglab/SeqSero), (3) trim reads using Trimmomatic, 
(4) de novo assemble reads, and (5) perfrom in silico PCR 
input: 1) two column tab-delimited file, first column SRA accession, second column sample
        name
       2) file with primer-pairs formatted for IsPCR (space delimited)
       e.g., primer_pair_name FWD_primer Rev_primer  
output: a tab delimeted matrix with the SeqSero results and results of IsPCR (i.e., primer-
pair(length predicted product))
Assumes gnu parallel, sra toolkit, SeqSero, megahit and isPcr are installed in path, while
path to trimmomatic and aspera ascp are given by user (a little bi more complex for aspera)
