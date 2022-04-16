!/bin/bash/ 

# This script will perform read alignment to the reference genome and variant calling

# This line below will ensure that our script will exit if an error occurs
set -e

#-------This part performs quality check of untrimmed reads-----------
#changing directory to wgs data directory
cd ~/my_shared_data_folder/prostatecancer/wgs

#creating variables to locate  the REF genome file and gene annotation file
#uncomment to download the reference genome and gtf file
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz 

genome=~//my_shared_data_folder/prostatecancer/wgs/reference/GRCh38.p13_genome.fa
gtf=~/my_shared_data_folder/prostatecancer/referenceGenome/GRCh38.p13.gtf

#creating the output directories
#the -p opton allows to create the whole path of the directories if they dont exist without complaining if they exist
mkdir -p ./fastqc_output ./multiqc_output ./trimmed_output  ./outputs ./trimmed_output/fastqc_reports ./trimmed_output/multiqc_reports

# Run FastQC and move output to appropriate location
# These next two lines will give us a status message to tell us that we are currently running FastQC,
# then will run FastQC on all of the files in our current directory with a .fastq extension

#echo 'Running FastQC....'
fastqc -t 10 ./*.fastq

#Moving the .zip and .html file into fastqc_o
mv *.zip *.html ./fastqc_output

#Running Multiqc on the fastqc files
multiqc -d . -o ./multiqc_output

#------------------- we are going to use trimmomatic v0.39 for quality control-------------------------------------------
 
# Trimmomatic filter poor quality reads and trim poor quality bases from NGS sample reads
# Trimmomatic can operate one sample at a time, to iterate through many samples we use a for loop
# The option PE specifies we will be taking paired end file as input,
# In paired end mode, Trimmomatic expects the two input files, and then the names of the output files
# The first two file names after the ${infile} flag represents our two input files and the next four are our output files
# The option SLIDINGWINDOW:4:20 allows us to use a sliding window of size 4 that will remove bases if their phred score is below 20
# The option ILLUMINACLIP:SRR_adapters.fa allow us to o clip the Illumina adapters from the input file using the adapter sequences
# listed in SRR_adapters.fa

for infile in ./*_1.fastq
do
base=$(basename ${infile} _1.fastq)
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar  PE -threads 6 -summary trimmomatic_summary.txt ${infile} ${base}_2.fastq ${base}_1.trim.fastq ${base}_1un.trim.fastq ${base}_2.trim.fastq ${base}_2un.trim.fastq SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:~/my_shared_data_folder/prostatecancer/rna_seq/rnaseq_data/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:40:15

done

#moving our trimmed fastq files to trimmed_output directory
mv *.trim* ./trimmed_output

#changing directory to trimmed_output
cd trimmed_output

#Running fastqc
#These next two lines will give us a status message that we are currently running FastQC,
echo 'Running FastQC....'
fastqc -t 10  *.fastq

#Moving the .zip and .html files into the fastqc_reports directory
mv ./*.zip fastqc_reports
mv ./*.html fastqc_reports

#Running multiqc
multiqc . -o ./multiqc_report

#Indexing the reference genome for reads alignment / mapping
bwa index $genome

#The command below creates a directory structure to store our results for alignment and variant calling processes
mkdir -p ./sam ./bam ./bcf ./vcf

# The for loop allows us to run the variant calling workflow on each of our fastq file.
#The full list of commands below will be executed once for each of the files 
for fq1 in ./*_1.fastq #for all files with forward reads in trimmed directory
    do
    echo "working with file: " $fq1  
    
    # Creating base names for each file
    base=$(basename $fq1 _1.fastq) 
    echo "base name is: " $base

    fq1=./${base}_1.fastq
    fq2=./${base}_2.fastq
    sam=./sam/${base}.aligned.sam
    bam=./bam/${base}.aligned.bam
    sorted_bam=./bam/${base}.aligned.sorted.bam
    raw_bcf=./bcf/${base}_raw.bcf
    variants=./vcf/${base}_variants.vcf
    final_variants=./vcf/${base}_final_variants.vcf 


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
    
    #The command below calls SNVs and INDELs with bcftools call. We have to specify ploidy with the flag --ploidy, which is one for the haploid E. coli.
    #-m allows for multiallelic and rare-variant calling, -v tells the program to output variant sites only (not every site in the genome), 
    #and -o specifies where to write the output file
    bcftools call  -m -v -o $variants $raw_bcf
    
    #The command below filters the SNVs using vcfutils.pl and reports the SNVs in variant calling format 
    vcfutils.pl varFilter $variants > $final_variants
    
    done
