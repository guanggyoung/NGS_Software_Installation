#!/bin/bash
set -e
cat SRA_acc_List.txt | xargs fastq-dump --gzip --skip-technical --dumpbase --split-files $1
