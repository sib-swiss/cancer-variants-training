#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments

for sample in tumor normal
do
    gatk MarkDuplicates \
    --INPUT "$ALIGNDIR"/"$sample".rg.bam \
    --OUTPUT "$ALIGNDIR"/"$sample".rg.md.bam \
    --METRICS_FILE "$ALIGNDIR"/marked_dup_metrics_"$sample".txt 

    samtools index "$ALIGNDIR"/"$sample".rg.md.bam
done