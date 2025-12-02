#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"

# decompress the reference
if [ ! -f "$REFDIR"/ref_genome.fa ]; then
    gunzip -c "$REFDIR"/ref_genome.fa.gz > "$REFDIR"/ref_genome.fa
fi
samtools faidx "$REFDIR"/ref_genome.fa

mkdir -p "$ALIGNDIR"

bwa index "$REFDIR"/ref_genome.fa
