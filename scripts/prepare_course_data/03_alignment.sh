#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
READDIR="${HOME}/project/course_data/reads"


mkdir -p "$ALIGNDIR"

for sample in tumor normal
do
    bwa mem \
    "$REFDIR"/ref_genome.fa \
    "$READDIR"/"$sample"_R1.fastq.gz \
    "$READDIR"/"$sample"_R2.fastq.gz \
    2> "$ALIGNDIR"/$sample.bwa.log \
    | samtools sort -o "$ALIGNDIR"/"$sample".bam -
done