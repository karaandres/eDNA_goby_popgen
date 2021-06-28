#!/bin/bash
for f1 in *_R1.fastq.gz
do
    f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
    java -jar /programs/trimmomatic/trimmomatic-0.39.jar PE -threads 24 -phred33 $f1 $f2 ${f1}_paired.fastq.gz ${f1}_unpaired.fastq.gz ${f2}_paired.fastq.gz ${f2}_unpaired.fastq.gz ILLUMINACLIP:/programs/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 MINLEN:35
done
