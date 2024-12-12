#!/usr/bin/env bash

ALIGNDIR=~/project/data/alignments
RESOURCEDIR=~/project/data/resources
VARIANTDIR=~/project/data/variants

for sample in tumor normal
do
    gatk GetPileupSummaries \
    -I "$ALIGNDIR"/"$sample".rg.md.bam \
    -V "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
    -L "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
    -O "$VARIANTDIR"/"$sample".pileups.table
done
