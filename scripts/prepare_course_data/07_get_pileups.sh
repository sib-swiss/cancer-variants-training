#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
VARIANTDIR="$HOME/project/course_data/variants"
RESOURCEDIR="${HOME}/project/course_data/resources"

for sample in tumor normal
do
    gatk GetPileupSummaries \
    -I "$ALIGNDIR"/"$sample".rg.md.bam \
    -V "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
    -L "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
    -O "$VARIANTDIR"/"$sample".pileups.table
done
