#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
REFDIR=/config/data/reference
RESOURCEDIR=/config/data/resources
VARIANTDIR=/config/data/variants

mkdir -p $VARIANTDIR

gatk Mutect2 \
-R "$REFDIR"/ref_genome.fa \
--intervals "$REFDIR"/exome_regions.bed.interval_list \
-I "$ALIGNDIR"/tumor.rg.md.bam \
-I "$ALIGNDIR"/normal.rg.md.bam \
-normal normal \
--germline-resource "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
--panel-of-normals "$RESOURCEDIR"/1000g_pon.hg38.subset.vcf.gz \
-O "$VARIANTDIR"/somatic.vcf.gz