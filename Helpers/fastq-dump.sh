#!/bin/bash

# module load sratoolkit # (if using HPC)
cd /data/guang/Seq_RAW_data
# --gzip to Compress output fastq
# --skip-technical Dump only biological reads
# --dumpbase Formats sequence using base space(default for other machines than SOLiD. SOLiD use color space)
# --split-files Dump each read into separate file. Create two files for pair read
# --readids. Do NOT use it if file will be analyzed using BWA.
fastq-dump --gzip --skip-technical --dumpbase --split-files  SRR689241

