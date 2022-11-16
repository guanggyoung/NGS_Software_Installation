#!/bin/bash
# bash scNMT_DNAseq_QC_trimming_hisat3n_alignment_PE_get_foldername_as_parameter.sh folder_name
set -e
# load modules (if run on HPC)
module load fastqc
module load trimgalore
module load samtools
module load picard
# module load fastqc
# [+] Loading fastqc  0.11.9
# module load trimgalore
# [+] Loading singularity  3.8.5-1  on cn0858
# [+] Loading cutadapt  3.4
# [+] Loading eigen 3.3.9  ...
# [+] Loading trimgalore  0.6.7  ...
#  module load samtools
# [+] Loading samtools 1.14  ...
# module load picard
# [+] Loading picard  2.26.9


##******************* Parameter Definition Finished Here *********************##

##*************** The actual alignment script starts here. ********************##
################################################################################
# The top_level_folder contains folders of all samples.
# Use full directory for the top_level_folder
top_level_folder="/data/guang/DNAseq_data"
## Initialize the sample folder. Must use full path!!

## This sample_folder should contain the paired-end fatq.gz files for one sample.
sample_folder=$1 # This is the folder for one sample (e.g. E1G_1) which contains paired-end fastq files
cd $top_level_folder/$sample_folder # Get into the sample_folder with fastq.gz file(s)
mkdir -p $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
# Make a new folder in sample_folder to store FastQC result; Full path used here.
fastqc -t 8 *.fq.gz -O ./$sample_folder\_fastqc_results/
# -t 8: use 8 threads
# relative path is used. The full path is 
# $top_level_folder/$sample_folder/$sample_folder\_fastqc_results
# Two files generate after calling fastqc: 
#                (fastq.gz filename)_fastqc.html and (fastq.gz filename)_fastqc.zip



########################################################################################
################# Trimming with Trim_Galore
cd $top_level_folder/$sample_folder
# Get into the sample_folder with fastq.gz file
mkdir -p $top_level_folder/$sample_folder/$sample_folder\_trim_galore_trimmed
# Make a new folder in sample_folder to store trimmed result; Full path used here.
fastq_gz_files="$(ls *.fq.gz)"
# Nextera Transposase Adapters were used for scNMT-Seq DNA part
trim_galore --paired --retain_unpaired --phred33 --nextera \
--length 20 \
--output_dir ./$sample_folder\_trim_galore_trimmed \
$fastq_gz_files


######### hisat3n alignment ###############
cd $top_level_folder/$sample_folder
#get into sample_folder
mkdir -p $top_level_folder/$sample_folder/$sample_folder\_hisat3n_alignment
	
# trimmed_fastq_gz_files="$(ls -d $PWD/$sample_folder\_*\_trimmed/*.fq.gz)"
##$PWD is current path
##This command get the full path of trimmed fast.gz files
	
# num_trimmed_fastq="$(echo "$trimmed_fastq_gz_files" | wc -l)"
# if [ $num_trimmed_fastq -ge 2 ] 
## Check numer of trimmed fastq.gz files. If great than or equals to 2, they must come from 
## paired-end read files.
## Treat as trimmed fastq files from paired-end fastq files
# then
paired_trimmed_read1="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_val_1.fq.gz)"
paired_trimmed_read1_basename="$(basename $paired_trimmed_read1)"
paired_trimmed_read2="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_val_2.fq.gz)"
paired_trimmed_read2_basename="$(basename $paired_trimmed_read2)"
unpaired_trimmed_reads="$(ls -d $PWD/$sample_folder\_*\_trimmed/*_unpaired*.fq.gz)"
# use unpaired_reads to store the filenames of the 2 unpaired read files trimmed out by trimmomatic
hisat-3n -p 6 -x /data/guang/mouse_genome_index/GRCm38.plus.lambda.hisat3n_index/GRCm38.plus.lambda \
		 -q -1 $paired_trimmed_read1 -2 $paired_trimmed_read2 \
		 -S $top_level_folder/$sample_folder/$sample_folder\_hisat3n_alignment/$sample_folder\_PE_hisat3n_unique_alignment_results.sam \
		 --base-change C,T --no-spliced-alignment --no-repeat-index --unique-only
		 
cd $top_level_folder/$sample_folder/$sample_folder\_hisat3n_alignment

# Sort SAM file
samtools sort -@ 8 -m 4G $sample_folder\_PE_hisat3n_unique_alignment_results.sam -o $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam -O sam

# SAM file can not be indexed or not needed
# samtools index $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam

# Transfer the SAM file into bam file
samtools view -bS $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam > $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.bam
samtools index $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.bam
# fi

#### hisat3n-table
cd $top_level_folder/$sample_folder/$sample_folder\_hisat3n_alignment
# Deduplicate
# java -Xmx4g -XX:ParallelGCThreads=6 -jar $PICARDJARPATH/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
# 	METRICS_FILE=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam.dup.txt \
# 	VALIDATION_STRINGENCY=LENIENT \
# 	INPUT=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam \
# 	OUTPUT=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.deduplicated.sam

# Make it into single line to make sure picard works
java -Xmx4g -XX:ParallelGCThreads=6 -jar $PICARDJARPATH/picard.jar MarkDuplicates REMOVE_DUPLICATES=true METRICS_FILE=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam.dup.txt VALIDATION_STRINGENCY=LENIENT INPUT=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam OUTPUT=$sample_folder\_PE_hisat3n_unique_alignment_results_sorted.deduplicated.sam

# hisat-3n-table
hisat-3n-table --threads 6 --alignments $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.deduplicated.sam --ref /data/guang/mouse_genome_index/GRCm38.plus.lambda.hisat3n_index/Mus_musculus.GRCm38.dna.GRCm38.dna.primary_assembly_plus_lambda_DNA.fa --output-name $sample_folder\_CtoT.tsv --base-change C,T

# Remove SAM files to save space on disk
rm $sample_folder\_PE_hisat3n_unique_alignment_results.sam # save space of disk
rm $sample_folder\_PE_hisat3n_unique_alignment_results_sorted.sam # save space of disk
