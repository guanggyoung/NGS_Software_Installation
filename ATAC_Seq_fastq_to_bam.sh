#!/bin/bash

#######################################################################################
####################### Real analysis starts here ####################################


## *********** Definition of several parameters used by all scripts below ********##
## NOTE: 
##      Refine these parameters each time before running the code.
##      The parameter here are set as the current situation. 
##      They may be modified if system changed(Most of the time: different folders used).
##
# Folder contains all folders of samples
# Use full path.
## Get into the folder with the fastq files first

######### Parameters related to sample for analysis #####
sample_name=$1
species=$2

######### Parameters related to Bowtie2 ###############
# Human genome bowtie2 index files
human_bowtie2_index="/data/guang/genome_index/bowtie2/human_h38_release97/GRCh38.97"
# Mouse genome bowtie2 index files
mouse_bowtie2_index="/data/guang/genome_index/bowtie2/mouse_m38_release97/GRCm38.97"
####################################################

######### Parameters related to Trimmomatic ##########
# How to run Trimmomatic (inlcude path of Trimmomatic)
execute_Trimmomatic="java -jar /home/guang/bio_softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"
# Paired-end Trimmomatic adapter file
PE_Trimmomatic_adapter="/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"
# Single-end Trimmomatic adapter file
SE_Trimmomatic_adapter="/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-SE.fa"
######################################################

######## Parameter related to picard ################
execute_picard="java -XX:ParallelGCThreads=4 -jar /home/guang/bio_softwares/picard.jar"
## To use multi-thread
## java -XX:ParallelGCThreads=<num of  threads> -jar <program_name>.jar
#####################################################

# Determine which species is provided so that to determine the right bowtie2 genome inex
if [ $species = "human" ]
	then
		bowtie2_index=$human_bowtie2_index
elif [ $species = "mouse" ]
	then	bowtie2_index=$mouse_bowtie2_index
else
	echo "The script only accept <species> of mouse or human, NO others"
	exit 1
fi
## *************************************************************** ##



## ********** Part 1 Sequencing data quality control ***************##
## This part is shared by all analyzing pipelines related to Next-
## generation sequencing.
## FastQC and MultiQC are used to do this job
##********************** Description of Part 1 ****************************##
## In this part, fastqc program will simply be called to produce quality report of 
## the fastq.gz files named as <sample_name>_1.fastq.gz and <sample_name>_2.fastq.gz.
## -t 4 : using 4 threads
fastqc -t 4 $sample_name*

## *****************************************************************##


## ************ Part 2 Trim fastq.gz files using Trimmomatic ************* ##
##********************** Description of Part 2 ****************************##
## In this part, Trimmomatic will be called to trim out adapters and low quality
## end bases of fastq reads.
##
## This script expect pair-end reads for ATAC-seq, 
## so Trimmomatic will be called using PE mode.
## After Trimmomatic trimming, output files will have names:
## trimmed_<sample_name>_*.fastq.gz
###############################################################################

# Do pair-end trimming using Trimmomatic
# $sample_name*.fastq.gz are the 2 fastq.gz files related to this sample

$execute_Trimmomatic PE -phred33 \
$sample_name*.fastq.gz \
./trimmed\_$sample_name\_paired.R1.fastq.gz \
./trimmed\_$sample_name\_unpaired.R1.fastq.gz \
./trimmed\_$sample_name\_paired.R2.fastq.gz \
./trimmed\_$sample_name\_unpaired.R2.fastq.gz \
ILLUMINACLIP:$PE_Trimmomatic_adapter:2:30:10 LEADING:3 \
TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#Note: 'ILLUMINACLIP:/home/guang/bio_softwares/Trimmomatic-0.39/adapters/TruSeq3-PE.fa' is used
#      This is default set since newer experiments use TruSeq3 kits.

## ********************************************************************* ##


## ******* Part 3 Align trimmed reads to genome using Bowtie2 ******** ##
################# Description of Part 3 ###########################
## Bowtie2 is used to align trimmed pair-end reads to reference genome.
## The species argument while calling this script is used to define the
## the reference genome index.
## The reference genome should be generated before all the analysis.
## bowtie2-build --thread 4 /data/guang/human_genome/Homo_sapiens.GRCh38.97.dna.primary_assembly.fa ./human_h38_release97/GRCh38.97
## 
## Only pair-end reads passed the trimming process will be used in the analysis
###################################################################

# better alignment results are frequently achieved with --very-sensitive
# use -X 2000 to allow larger fragment size (default is 500)
# >2 is used to redirect output.
bowtie2 --very-sensitive -X 2000 -p 6 \
-x $bowtie2_index \
-1 trimmed*$sample_name*_paired.R1.fastq.gz \
-2 trimmed*$sample_name*_paired.R2.fastq.gz \
-S ./$sample_name\_trimmed_pair_end_aln_unsorted.sam \
2> $sample_name.bowtie2.log 
samtools sort -@ 6 -O bam -o $sample_name.sorted.bam ./$sample_name\_trimmed_pair_end_aln_unsorted.sam

# bowtie2 is first called to use the 2 paired-end read files as 
# input to produce a  sam file.
# Then samtools is used to sort the output sam file from bowtie2 
# and sort it into bam file(to save disk space).

samtools index -@ 6 $sample_name.sorted.bam
# index the sorted bam file for IGV usage in the future

## ******* Remove mitochondrial reads from aligned bam file ******** ##
## NOTE: in 1st version of ATAC-Seq, a large percentage of reads came from
## mitochondria.
## $ samtools view SRR891269.sorted.bam | grep MT | wc -l
## 66229582
## samtools view SRR891269.sorted.bam | wc -l
## 109450758


### samtools -h : Include the header in the output.

## Optional
## Check chromosome names in bam file
## samtools idxstats $sample_name.sorted.bam | head -30
## For reference genome downloaded from ENSEMBL, mitochondia is named as "MT"

### For grep: 
### -v 
### --invert-match
### Invert the sense of matching, to select non-matching lines. (-v is specified by POSIX.)
samtools view -@ 6 -h $sample_name.sorted.bam | grep -v MT | samtools sort -@ 6 -O bam -o $sample_name.rmMT.bam

## **************************************************************** ##

##*** Label PCR duplicates in bam file using "picard MarkDuplicates" ***##
$execute_picard MarkDuplicates \
QUIET=true INPUT=$sample_name.rmMT.bam OUTPUT=$sample_name.duplicate.marked.bam METRICS_FILE=$sample_name.rmMT.sorted.metrics \
REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp
# REMOVE_DUPLICATES=false: mark duplicate reads, not remove.
# Change it to true to remove duplicate reads.

## ********* Remove duplicates in bam file using samtools ************* ##
# MarkDuplicates will add a FALG 1024 to duplicate reads, we can remove them using samtools
samtools view -@ 4 -h -b -F 1024 $sample_name.duplicate.marked.bam > $sample_name.rmMT.rmDup.bam
## ******************************************************************** ##

## ************ Process multiple-mapped(non-unique alignments) reads in bam file ************** ##
## Non-uniquely aligned reads can be checked by the -q parameter of samtools view
## Different genome aligners have varied implementation of mapping quality (MAPQ)
## For Bowtie2, people usually use MAPQ > 30 for non-unique alignments
### Remove multi-mapped reads (i.e. those with MAPQ < 30, using -q in SAMtools)
samtools view -@ 4 -h -q 30 $sample_name.rmMT.rmDup.bam > $sample_name.rmMT.rmPCRDup.rmMutimap.bam
## ******************************************************************************************* ##

## ************* Remove other low quality reads ************************* ##
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain only high-quality, properly paired reads -f 2
samtools view -@ 4 -h -b -F 1804 -f 2 $sample_name.rmMT.rmPCRDup.rmMutimap.bam > $sample_name.cleaned.high.quality.bam
## ********************************************************************** ##



