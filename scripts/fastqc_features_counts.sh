#!/bin/bash

## This script will perform quality control (QC) of reads  upto generating a count matrix
## The input is an rna-seq fastq file, and output a counts matrix with rows as genes and columns as samples

## first give executable rights (chmod +x <name of script_file>)
## run with bash <name of script_file>

#changing directory to rna_seq data directory, change directory to suite your directory
cd ~/my_shared_data_folder/prostatecancer/rna_seq/rnaseq_data

#creating variables to directory containing the reference genome file, index files and gene annotation file
#uncomment to download the reference genome, the gtf file and copy them to the directory
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gzt 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf.gz

ref_genome=~/my_shared_data_folder/prostatecancer/referenceGenome/GRCh38.p13_ref_genome.fna
gtf=~/my_shared_data_folder/prostatecancer/referenceGenome/GRCh38.p13.gtf
indexed_ref=~/my_shared_data_folder/prostatecancer/rna_seq/rnaseq_data/GRCh38.p13_ref_genome

#creating the output directories 
#the -p opton allows to create the whole path of the directories if they dont exist without complaining if they exist
#mkdir -p ./fastqc_output ./multiqc_output ./trimmed_output  ./outputs ./trimmed_output/fastqc_reports ./trimmed_output/multiqc_report

### -------------------------------Run FastQC and move output to appropriate location-------------------------------------------------------###
# These next two lines will give us a status message to tell us that we are currently running FastQC,
# then it will run FastQC on all of the files in our current directory with a .fastq extension
echo 'Running FastQC....'
fastqc -t 10 ./*.fastq

#Moving the .zip and .html file into fastqc_ouput directory
mv *.zip *.html ./fastqc_output 

###----------------------------------------Running Multiqc on the fastqc files--------------------------------------------------------------###
multiqc -d . -o ./multiqc_output

# The next two lines tells us we are unzipping and then we use for loop to unzips the files
# which are in  fastqc_output directory

cd fastqc_output
echo 'Unzipping......'
for filename in ./*.zip 
do 
unzip $filename
done

# Next we concatenate all our summary files into a single output
echo 'Saving summary....'
cat */*summary.txt > fastqc_summaries.txt
cd ..

###------------------------------- we are going to use trimmomatic v0.39 for quality control----------------------------------------------###

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

#moving our trimmed fastq files to trimmed_ directory
mv *.trim* ./trimmed_output

#changing directory to trimmed_output
cd trimmed_output

#-------------------------------------------------Running fastqc on trimmed reads----------------------------------------------------------###
#These next two lines will give us a status message that we are currently running FastQC,
echo 'Running FastQC....'
fastqc -t 10  *.fastq

#Moving the .zip and .html files into the fastqc_reports directory
mv ./*.zip ./*.html fastqc_reports


#Unzipping files in fastqc_reports
#The next two commands notifies us that we are unzipping the files using a for loop
cd fastqc_reports
echo 'Unzipping.....'
for filename in ./*.zip
do
unzip $filename
done

#Next we concatenate all our summary files into a single output, with a status message notifying us
echo 'Saving summary....'
cat */*summary.txt > trimmed_fastqc_summaries.txt

cd ..
#--------------------------------------------------Running multiqc on trimmed reads---------------------------------------------------------###
multiqc . -o ./multiqc_report

###-----------------------indexing the reference genome using Hisat2 build. This process can only be done once optionally-------------------###
hisat2-build $ref_genome $indexed_ref

##aligning all samples to indexed reference genome to generate a sam file as output using hisat2
##-----------------------------------------------------Alignment to ref genome--------------------------------------------------------------###
for r2 in ./*_1.trim.fastq
do
echo $r2
sample=$(basename $r2)
sample=${sample%_1.trim.fastq}
echo "Processing sample: " $sample

hisat2 -p 10 -x $indexed_ref \
-1 ./${sample}_1.trim.fastq \
-2 ./${sample}_2.trim.fastq | \
samtools view -S -b - | \
samtools sort -o ~/my_shared_data_folder/prostatecancer/rna_seq/rnaseq_data/outputs/${sample}_trimmed_sorted.bam - &&
samtools index ~/my_shared_data_folder/prostatecancer/rna_seq/rnaseq_data/outputs/${sample}_trimmed_sorted.bam

done

###----------------------------------------Count mapped reads as a measure of gene expression--------------------------------------------###
#The input will be bam files and a gtf file
#The output will be a cont matrix, with genes as rows and samples as column
#We are using the featureCounts tool because it is accurate, fast and  relatively easy to use

featureCounts -T 12  -a $gtf -o ./outputs/feature_counts.txt ./outputs/*.bam



