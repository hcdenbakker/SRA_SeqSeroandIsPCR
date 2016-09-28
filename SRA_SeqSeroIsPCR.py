'''
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
     
'''
import sys
import os
from sys import argv
import subprocess
import argparse
from collections import defaultdict

def aspera_download(infile):
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sra= accession[0]
        sample=accession[1]
        subprocess.call(["/home/henk/.aspera/connect/bin/ascp -i "
                         "/home/henk/.aspera/connect/etc/asperaweb_id_dsa.openssh -k1 -Tr -l50m "
                         "anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/"
                         +str(sra[:3])+'/'+str(sra[:6])+'/'+str(sra)+'/'+str(sra)+'.sra '
                         +str(sample)+'.sra'],stdout=subprocess.PIPE, shell=True)
    fp.close()

    subprocess.call(["parallel  -j 30 'fastq-dump --split-3 --gzip {}' ::: *.sra"],stdout=subprocess.PIPE, shell=True )


def RunSeqSero(infile):
    #seqsero_results=[]
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        subprocess.call(['SeqSero.py -m 2 -i  '+str(sample)+'_*.fastq.gz > '+str(sample)+'_seqsero.out'],stdout=subprocess.PIPE, shell=True )


def make_SeqSero_dict(infile):
    seqsero_results = defaultdict(list)
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        seqsero_file = str(sample) + '_seqsero.out'
        SeqSero = open(seqsero_file, 'r')
        for i, line in enumerate(SeqSero):
            line = line.strip('\n')
            if line.split('\t')[0] == str('Predicted serotype(s):'):
                predicted_serotype = line.split('\t')[1]
                seqsero_results[sample].append(predicted_serotype)
            else:
                pass
    return seqsero_results



def clean_and_assemble(infile):
    cwd = os.getcwd()
    subprocess.call(['mkdir contigs'],stdout=subprocess.PIPE, shell=True )
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        subprocess.call(['java -jar /home/henk/Trimmomatic-0.36/trimmomatic-0.36.jar '+
            'PE -threads 32 -phred33 '+str(sample)+'_1.fastq.gz '+str(sample)+'_2.fastq.gz '+
            str(sample)+'_1.trimmedP.fastq.gz '+str(sample)+'_1.trimmedS.fastq.gz '+str(sample)+
            '_2.trimmedP.fastq.gz '+str(sample)+'_2.trimmedS.fastq.gz '+
            'ILLUMINACLIP:/home/henk/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 '+
            'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'],stdout=subprocess.PIPE, shell=True )
        subprocess.call(['rm *.trimmedS.fastq.gz'], stdout=subprocess.PIPE, shell=True)
        subprocess.call(['megahit --presets bulk -1 '+str(sample)+'_1.trimmedP.fastq.gz -2 '+
            str(sample)+'_2.trimmedP.fastq.gz -o '+str(sample)+'_megahit'], stdout=subprocess.PIPE, shell=True)
        subprocess.call(['cd '+str(sample)+'_megahit; mv final.contigs.fa '+str(sample)+'.contigs.fa; mv '
                         + str(sample) + '.contigs.fa '+str(cwd)+'/contigs'],stdout=subprocess.PIPE, shell=True)
        subprocess.call(['rm -r '+str(sample)+'_megahit'],stdout=subprocess.PIPE, shell=True )


def isPCR(infile, primer_file):
    cwd = os.getcwd()
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        subprocess.call(['isPcr '+str(cwd)+'/contigs/'+str(sample)+'.contigs.fa '+str(primer_file)+' '+str(sample)+'_is.out']
                        ,stdout=subprocess.PIPE, shell=True )

def make_isPCR_list(infile):
    isPCR_results=[]
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        isPcr_file = str(sample)+'_is.out'
        SeqSero = open(isPcr_file, 'r')
        for i, line in enumerate(SeqSero):
            if line[0] == '>':
                primer_pair = line.split(' ')[1]
                product_length = line.split(' ')[2]
                isPCR_results.append([sample,[primer_pair, product_length]])
            else:
                pass
    return isPCR_results


def create_table(infile, predictions, isPCR_list, primerfile):
    header=['\t'+'SeqSero']
    primer_list = []
    fp = open(primerfile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        primers= line.split(' ')
        header.append(primers[0])
        primer_list.append(primers[0])
    fp.close()
    result_table = open('results.out', mode='a')
    result_table.write('\t'.join(map(str, header))+'\n')
    fp = open(infile, 'r')
    for i, line in enumerate(fp):
        line = line.strip('\n')
        accession = line.split('\t')
        sample = accession[1]
        observations=[]
        observations.append(sample)
        if predictions[sample][0][:3] == 'N/A':
            observations.append('N/A')
        else:
            observations.append(predictions[sample][0])
        amplifications=[]
        amplifications_short=[]
        ordered_amplifications=[]        
        for s in isPCR_list:
            if str(s[0]) == str(sample):
                amplifications.append(s[1])
                amplifications_short.append(s[1][0])
        for p in primer_list:
            if p in amplifications_short:            
                for a in amplifications:
                    if str(p) == str(a[0]):
                        ordered_amplifications.append(str(a[0])+'('+str(a[1])+')')
                        
            else:
                ordered_amplifications.append('none')
                    
        for p in ordered_amplifications:
            observations.append(p)   
        result_table.write('\t'.join(map(str, observations))+'\n')
    fp.close()
    result_table.close()



def main():
    infile = sys.argv[1]
    primerfile = sys.argv[2]
    #aspera_download(infile)
    #RunSeqSero(infile)
    predictions = make_SeqSero_dict(infile)
    print(predictions)
    #clean_and_assemble(infile)
    isPCR(infile, primerfile)
    isPCR_list = make_isPCR_list(infile)
    print(isPCR_list)
    create_table(infile, predictions, isPCR_list, primerfile)





if __name__ == '__main__':
    main()
