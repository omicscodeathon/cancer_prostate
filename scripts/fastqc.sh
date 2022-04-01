#!/bin/bash

# This script will check the quality of our sequencing reads using fastqc

# The command below ensures that our script will exit if an error occurs
set -e

#-------This part performs quality check of untrimmed reads-----------
#Uncomment to run this part

# changing directory to untrimmed_fastqc
#cd ~/codatheon/genome_workflow/Escherichia/data/untrimmed_fastqc
# These next two lines will give us a status message to tell us that we are currently running FastQC,
# then will run FastQC on all of the files in our current directory with a .fastq extension

#echo 'Running FastQC....'
#fastqc ../*.fastq*

#Moving the .zip and .html file into the untrimmed_fastqc directory
#mv ../*.zip . 
#mv ../*.html .

# The next two lines tells us we are unzipping and then we use for loop to unzips the files

#echo 'Unzipping......'
#for filename in *.zip 
#do 
#unzip $filename
#done
# Next we concatenate all our summary files into a single output, with a status message reminding ourselves
#echo 'Saving summary....'
#cat */summary.txt > fastqc_summaries.txt

#.........This part performs quality check of trimmed reads.........................

#Uncomment this part to perform QC on trimmed reads

#Changing directory to trimmed_fastqc
cd ~/codatheon/genome_workflow/Escherichia/data/trimmed_fastq
#These next two lines will give us a status message that we are currently running FastQC,
echo 'Running FastQC....'
fastqc *.fastq

#Moving the .zip and .html files into the trimmed_fastqc directory
mv ./*.zip ~/codatheon/genome_workflow/Escherichia/data/trimmed_fastqc
mv ./*.html ~/codatheon/genome_workflow/Escherichia/data/trimmed_fastqc

#changing directory to trimmed_fastqc 
cd ~/codatheon/genome_workflow/Escherichia/data/trimmed_fastqc

#The next two commands notifies us that we are unzipping the files using a for loop
echo 'Unzipping.....'
for filename in *.zip
do
unzip $filename
done

#Next we concatenate all our summary files into a single output, with a status message notifying us
echo 'Saving summary....'
cat */summary.txt > trimmed_fastqc_summaries.txt






