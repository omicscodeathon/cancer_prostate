#!/bin/bash
# This script performs 

# we are going to use trimmomatic v0.39 for quality control
# Trimmomatic filter poor quality reads and trim poor quality bases from NGS sample reads
# Trimmomatic can operate one sample at a time, to iterate through many samples we use a for loop

# The option PE specifies we will be taking paired end file as input, 
# In paired end mode, Trimmomatic expects the two input files, and then the names of the output files
# The first two file names after the ${infile} flag represents our two input files and the next four are our output files
# The option SLIDINGWINDOW:4:20 allows us to use a sliding window of size 4 that will remove bases if their phred score is below 20
# The option ILLUMINACLIP:SRR_adapters.fa allow us to o clip the Illumina adapters from the input file using the adapter sequences
# listed in SRR_adapters.fa


for infile in *_1.fastq
do
base=$(basename ${infile} _1.fastq)
trimmomatic PE ${infile} ${base}_2.fastq ${base}_1.trim.fastq ${base}_1un.trim.fastq ${base}_2.trim.fastq ${base}_2un.trim.fastq SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

done

#making a new directory trimmed_fastq
mkdir -p ./trimmed_fastq

#moving our trimmed fastq files to trimmed_fastq directory
mv *.trim* ./trimmed_fastq

#The command below will remind the user to check the quality of the trimmed reads
echo 'You may consider checking the quality of trimmed reads with fastqc....'
