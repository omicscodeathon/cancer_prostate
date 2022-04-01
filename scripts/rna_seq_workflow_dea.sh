#!/bin/bash

#This script performs RNA_seq  analysis from transcriptome alignemnt
#visualisation, simple quality control checks and profile transcriptomic
#differences by identifying differentially expressed genes.

#The command below notifies us of any errors, when running the script
set -e


#Changing directory to the data used for analysis otherwise create them by uncommenting the command below
#mv -p course_data/rna_seq_human
cd ~/course_data/rna_seq_human


#alternatively annotation data can be downloaded using the commands below
#wget https://www.dropbox.com/s/nfuea7ik6wlidum/hsapiens_chr21_transcript_to_gene.csv?dl=1 -O data/hsapiens_chr21_transcript_to_gene.csv

#wget https://www.dropbox.com/s/8xt8q1o0aej1ry1/hsapiens_chr21_transcripts.fa?dl=1 -O data/hsapiens_chr21_transcripts.fa


#creating a directory to store the outputs

mkdir -p outputs

#----------------we will start with QC of reads.................
#The next step is to map the RNA-seq reads to the genome using HISAT2
#HISAT2 is a fast and sensitive splice aware aligner. HISAT2 compresses
#the genome using an indexing scheme BWT and FM index to reduce the amount
#of space needed to store the genome.

#Note that we have restricted the number of reads in each sample to those mapping largely
#to a single chromosome to reduce the mapping time

#The command below builds a HISAT2 index for chromosome 21 of the human reference genome using hisat2_build.
hisat2-build data/hsapien_grch38_chr21.fa outputs/hsapien_grch38_chr21_hisat2.idx

#Performing reads alignment for the samples using HISAT2. -x tells us that we chose our index files with hisat2-build
#-1 and -2 are used to indicate the left and right read files for the samples
#-S tells HISAT2 we want the output in sam format
#hisat2 -x outputs/hsapien_grch38_chr21_hisat2.idx -1 data/PT2_1.fastq.gz -2 data/PT2_2.fastq.gz -S outputs/PT2.sam

#Converting the SAM file to BAM format using samtools, BAM file is smaller 
#samtools view -S -o outputs/PT2.bam -b outputs/PT2.sam

#Sorting the BAM file.sorting aligment by the coordinates for each chromosoe to be able to index it
#samtools sort outputs/PT2_sorted.bam outputs/PT2.bam

#Indexing our BAM file so that it can be read efficiently by IGV
#samtools index outputs/PT2_sorted.bam

for r1 in data/*_1.fastq.gz
do
echo $r1
sample=$(basename $r1)
sample=${sample%_1.fastq.gz}
echo "Processing sample: "$sample

hisat2 -x outputs/hsapien_grch38_chr21_hisat2.idx \
-1 data/${sample}_1.fastq.gz \
-2 data/${sample}_2.fastq.gz | \
samtools view -S -b - | \
samtools sort -o outputs/${sample}_sorted.bam - &&
samtools index outputs/${sample}_sorted.bam

done
#The alignments generated here can be used for transcriptome quantification using tools
#such as feautureCounts or htseq-count and a reference transcriptome to quantify read
#counts at the gene or transcript level

#Visualisation can be done using IGV
#uncoomment the lines below to start IGV
#igv &

#Quantifying Transcript Expression with Kallisto
#Kallisto uses a fast alignment-free method for transcript quantification, that is Rather than looking at where the
#reads map, Kallisto uses the compatibility between the reads and transcripts to estimate transcript
#abundance.

#Step 1: building a kallisto index called GRCh38_kallisto from transcript sequences
#To generate the index, Kallisto first
#builds a transcriptome de Bruijn Graph (TBDG) from all of the k-mers (short sequences of k nucleotides) that it finds in the transcriptome.
kallisto index -i outputs/GRCh38_ch21_kallisto data/hsapiens_chr21_transcripts.fa

#Quantifying the transcript expression levels for our samples with 100 bootstrap samples and
#store the results in the output directory of each sample
#i is our index prefix, -o our outout directory and -b option for bootsrapping where by kallisto
#quantifies the uncertainty in its abundance estimates using random resampling and replacement

for r1 in data/*_1.fastq.qz
do
echo $r1
sample=$(basename $r1)
sample=${sample%_1.fastq.gz}
echo "Processing sample: "$sample

kallisto quant -i outputs/GRCh38_ch21_kallisto -o outputs/${sample} -b 100 data/${sample}_1.fastq.gz data/${sample}_2.fastq.gz

done
#Running kallisto quant generated three output files in the PT6 directory
#---*abundance.h5 --> file conataining run info, abundances estimates, bootsrap estimates and transcript length info length
#---*abundance.tsv--> Plain text file containing abundance estimates (doesnt contain bootstrap estimates) 
#---*run_info.json--> JSON file containing information about the run.

#View abundance estimates for each gene for the PT6 sample
head outputs/PT6/abundance.tsv

#Differential Expression Analysis with Sleuth

#The goal of differential expression analysis is to identify genes whose expression levels differ between experimental conditions.
#sleuth is a companion tool for Kallisto. Unlike most other tools, sleuth can utilize the technical variation information
# generated by Kallisto so that you can look at both the technical and biological variation in your dataset.


#Configure the option for a web browser for R at the command line
export R_BROWSER='firefox'

#Navigating to the correct working directory within R:
setwd("~/course_data/rna_seq_human")
#loading the sleuth package
library(sleuth)
library(dplyr)

# load sample metadata
sample_info <- read.table(file="data/sample_info.txt", header = T, sep = "\t")
sample_info

#The command below will create a vector of file paths pointing to kallisto quantification results
kallisto_result_directory <-sapply(X = sample_info$sample,
				function(id) file.path('outputs', id))
kallisto_result_directory

#The command below will configure sample_to_covaraites data frame
s2c <- dplyr::select(sample_info,
	sample = sample, sample_type, ethnicity, individualID=patientID)

kallisto_result_directory <- sapply(sample_info$sample, function(id)
file.path('outputs', id)) # path to kallisto results
s2c <- dplyr::mutate(s2c, path = kallisto_result_directory)
s2c

# transcript to gene annotations
filename <-"data/hsapiens_chr21_transcript_to_gene.csv"
t2g <- read.table(filename, header = T, sep = ',')
names(t2g) <- c('target_id', 'ensembl_gene_id', 'gene_symbol','gene_biotype', 'gene_descrition')
t2g <- t2g[,1:4]
head(t2g)
dim(t2g)

#Creating a sleuth object
#Here it is important to specify the model design where you enumerate the covariates you want to model
#These can include both technical and biological variables under study
#Here we consider the following covariates: * individualID: unique identifier for each individual
#* sample_type: normal vs. cancer (tumor) * ethnicity: African American (AA) vs/ European American (EA)

# create design matrix specification
design <- ~ individualID + sample_type + ethnicity

# Creating a  sleuth object (a group of kallistos) for analysis
so <- sleuth_prep(sample_to_covariates=s2c, # sample_to_covariates data frame
	full_model=design, # model design matrix
	target_mapping = t2g, # transcript to gene annotations
	aggregation_column = "ensembl_gene_id",
	extra_bootstrap_summary = TRUE,
	read_bootstrap_tpm = TRUE,
	transformation_function = function(x) log2(x + 0.5),
	num_cores=1)

#Step 2: Fit the sleuth model , sleuth's measurement error model
# (full - consider all covariates simultaneoulsy)
so <- sleuth_fit(so, formula=design, fit_name="full")

# fit a reduced model - fit a model without the final factor
so <- sleuth_fit(so, formula=~ individualID + sample_type , fit_name="reduced")

#Step 3: Statistical testing between conditions
# likelihood ratio test between the two models
# tests for ancestry-related differences
# --- difference between the full model with 3 coveriates----
# versus the reduced model with 2 covariates
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table_tx <- sleuth_results(obj = so,
				test = 'reduced:full',
				test_type ='lrt',
				show_all = FALSE,
				pval_aggregate = FALSE)
head(sleuth_table_tx, 5)

#From our results we get a table with the top 5 transcripts in the differential expression analysis results table.
#The key columns are: * pval - p-value of the chosen model * qval - false discovery rate adjusted p-value

#Visualize transcript abundance for top hit from differential expression analysis:
topDE_hit <- sleuth_table_tx[1,"target_id"]
# group view
plot_bootstrap(so, topDE_hit, units = "est_counts", color_by = "ethnicity")

# paired sample view
library(ggplot2)
df <- get_bootstrap_summary(so , topDE_hit)

ggplot(data = df, aes(x = sample_type, ymin = min, lower = lower,
	middle = mid, upper = upper, ymax = max)) +
	geom_boxplot(stat="identity", aes(fill=ethnicity)) +
	facet_wrap(~individualID)

#You can look at your sleuth object to see what models you have fit and their design matrix specification like this:
models(so)

#Wald's Testing: testing for significant differences between conditions
#Wald test for individual coefficients (betas) tumor vs. normal( while controlling for
#inter-individual variation and ethnicity specific differences
so <- sleuth_wt(obj = so, which_beta = "sample_typetumor", which_model = 'full')

#summary table
de_sampletype <- sleuth_results(so, test='sample_typetumor',
			test_type = "wt", which_model = "full", pval_aggregate = F)
head(de_sampletype, 5)

#Visaulize the results

df <- get_bootstrap_summary(so , de_sampletype[1,'target_id'])

ggplot(data = df, aes(x = sample_type, ymin = min, lower = lower,
	middle = mid, upper = upper, ymax = max)) +
	geom_boxplot(stat="identity", aes(fill=ethnicity)) +
	facet_wrap(~individualID) 

#Saving the sleuth  result object
save(so, file = "outputs/sleuthObj.RData")

#Our model above generated few ancestry-specific differentially expressed genes after false discovery rate correction 
# launch interactive exploration of sleuth differential expression object
sleuth_live(so)
# if you have issues connecting to the browswer try this alternative:
# sleuth_live(so, options=list(launch.broswer=FALSE))


