#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
RESOURCEDIR=/config/data/resources
VARIANTDIR=/config/data/variants

for sample in tumor normal
do
    gatk GetPileupSummaries \
    -I "$ALIGNDIR"/"$sample".rg.md.bam \
    -V "$RESOURCEDIR"/af-only-gnomad.hg38.vcf.gz \
    -L "$RESOURCEDIR"/af-only-gnomad.hg38.vcf.gz \
    -O "$VARIANTDIR"/"$sample".pileups.table
done
