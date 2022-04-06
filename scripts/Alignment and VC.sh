#!/bin/bash/ 

# This pipeline will perform read alignment to the reference genome and variant calling

# This line below will ensure that our script will exit if an error occurs
set -e

#The line below makes directory results/
mkdir -p ~/codatheon/genome_workflow/Escherichia/data/results

#The line below navigates us to the results directory.
cd ~/codatheon/genome_workflow/Escherichia/data/results/

#The line below create a variable for path to our reference genome
genome=~/codatheon/genome_workflow/Escherichia/data/refGenome/ecoli_rel606.fasta

#The line below indexes the reference genome for use by bwa and samtools
bwa index $genome

#The command below creates a directory structure to store our results for alignment and variant calling processes
mkdir -p ./sam ./bam ./bcf ./vcf

# The for loop allows us to run the variant calling workflow on each of our fastq file.
#The full list of commands below will be executed once for each of the files 
for fq1 in ~/codatheon/genome_workflow/Escherichia/data/trimmed_fastq/*_1.trim.fastq #for all files with forward reads in trimmed directory
    do
    echo "working with file $fq1"  
    
    # Creating base names for each file
    base=$(basename $fq1 _1.trim.fastq) 
    echo "base name is $base"

    fq1=~/codatheon/genome_workflow/Escherichia/data/trimmed_fastq/${base}_1.trim.fastq
    fq2=~/codatheon/genome_workflow/Escherichia/data/trimmed_fastq/${base}_2.trim.fastq
    sam=~/codatheon/genome_workflow/Escherichia/data/results/sam/${base}.aligned.sam
    bam=~/codatheon/genome_workflow/Escherichia/data/results/bam/${base}.aligned.bam
    sorted_bam=~/codatheon/genome_workflow/Escherichia/data/results/bam/${base}.aligned.sorted.bam
    raw_bcf=~/codatheon/genome_workflow/Escherichia/data/results/bcf/${base}_raw.bcf
    variants=~/codatheon/genome_workflow/Escherichia/data/results/vcf/${base}_variants.vcf
    final_variants=~/codatheon/genome_workflow/Escherichia/data/results/vcf/${base}_final_variants.vcf 


    #The command below aligns the reads to the reference genome and output a .sam file
    #Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment
    #Indexing the reference only has to be done once, unless you chose another reference genome   
    bwa mem $genome $fq1 $fq2 > $sam
    #The command below converts the sam file to bam format using the view command and tell this command 
    #that the input is in SAM format (-S) and to output BAM format (-b)
    samtools view -S -b $sam > $bam
    #The command below sorts the BAM file using the sort command from samtools. -o tells the command where to write the output.
    samtools sort -o $sorted_bam $bam
    #The command below indexes the BAM file for display purposes
    samtools index $sorted_bam
    #The command below calculates the read coverage of positions in the genome using bcftools.
    #The flag -O b tells bcftools to generate a bcf format output file, -o specifies where to write the output file, 
    #and -f flags the path to the reference genome.
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    #The command below calls SNVs with bcftools call. We have to specify ploidy with the flag --ploidy, which is one for the haploid E. coli.
    #-m allows for multiallelic and rare-variant calling, -v tells the program to output variant sites only (not every site in the genome), 
    #and -o specifies where to write the output file
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
    #The command below filters the SNVs using vcfutils.pl and reports the SNVs in variant calling format 
    vcfutils.pl varFilter $variants > $final_variants
    
    done
