#!/bin/bash

#Written by Z. Chikwambi 22/03/2022 for the codeathon project
#Retrieving sequences from SRA database
#Requires SRS-toolkit to be installed first and made available to environment

#Data retrieval from SRA:
# download file: prefetch will download and save SRA file related to SRR accession in 
# the current directory under newly created SRA accession directory

#GWS
prefetch SRR2589044 SRR2584863 SRR2584866  # for a single file

#RNA-seq

#prefetch  SRR5790106 SRR5790104  # multiple files

# convert to FASTQ: fasterq-dump will convert SRR.sra to SRR.fastq

# now you can also replace fastq-dump with fasterq-dump which is much faster 
# and efficient for large datasets
# by default it will use 6 threads (-e option)
# for paired-end data use --split-files (fastq-dump) and -S or --split-files (fasterq-dump) option
# use --split-3 option, if the paired-end reads are not properly arranged (e.g. some reads has lack of mate pair)


#GWS
fasterq-dump  -S  SRR2584866 SRR2584863 SRR2589044  #SRR2589044 #SRR5790106  # single file

#RNA-seq

#fasterq-dump  -S SRR5790106  SRR5790104 # multiple files

