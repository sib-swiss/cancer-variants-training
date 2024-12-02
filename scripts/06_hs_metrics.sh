#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
REFDIR=/config/data/reference
RESOURCEDIR=/config/data/resources
VARIANTDIR=/config/data/variants

for sample in tumor normal
do
    gatk CollectHsMetrics \
    -I "$ALIGNDIR"/"$sample".rg.md.bam \
    -O "$ALIGNDIR"/"$sample".rg.md.bam_hs_metrics.txt \
    -R "$REFDIR"/ref_genome.fa \
    --BAIT_INTERVALS "$REFDIR"/exome_regions.bed.interval_list \
    --TARGET_INTERVALS "$REFDIR"/exome_regions.bed.interval_list
done
