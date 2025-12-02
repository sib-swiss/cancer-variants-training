#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"

for sample in tumor normal
do
    gatk MarkDuplicates \
    --INPUT "$ALIGNDIR"/"$sample".rg.bam \
    --OUTPUT "$ALIGNDIR"/"$sample".rg.md.bam \
    --METRICS_FILE "$ALIGNDIR"/marked_dup_metrics_"$sample".txt \
    --CREATE_INDEX true

done