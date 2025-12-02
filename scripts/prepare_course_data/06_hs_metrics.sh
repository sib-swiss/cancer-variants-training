#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
VARIANTDIR="$HOME/project/course_data/variants"
RESOURCEDIR="${HOME}/project/course_data/resources"

for sample in tumor normal
do
    gatk CollectHsMetrics \
    -I "$ALIGNDIR"/"$sample".recal.bam \
    -O "$ALIGNDIR"/"$sample".recal.bam_hs_metrics.txt \
    -R "$REFDIR"/ref_genome.fa \
    --BAIT_INTERVALS "$REFDIR"/exome_regions.bed.interval_list \
    --TARGET_INTERVALS "$REFDIR"/exome_regions.bed.interval_list
done
