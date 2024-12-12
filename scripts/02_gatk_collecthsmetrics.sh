#!/usr/bin/env bash

ALIGNDIR=~/project/course_data/alignments
REFDIR=~/project/course_data/reference
RESOURCEDIR=~/project/course_data/resources
VARIANTDIR=~/project/course_data/variants

for sample in tumor normal
do
    gatk CollectHsMetrics \
    -I "$ALIGNDIR"/"$sample".recal.bam \
    -O "$ALIGNDIR"/"$sample".recal.bam_hs_metrics.txt \
    -R "$REFDIR"/ref_genome.fa \
    --BAIT_INTERVALS "$REFDIR"/exome_regions.bed.interval_list \
    --TARGET_INTERVALS "$REFDIR"/exome_regions.bed.interval_list
done
