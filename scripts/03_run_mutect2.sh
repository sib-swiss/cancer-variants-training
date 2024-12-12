#!/usr/bin/env bash

ALIGNDIR=~/project/course_data/alignments
REFDIR=~/project/course_data/reference
RESOURCEDIR=~/project/course_data/resources
VARIANTDIR=~/project/course_data/variants

mkdir -p $VARIANTDIR

gatk Mutect2 \
-R "$REFDIR"/ref_genome.fa \
--intervals "$REFDIR"/exome_regions.bed.interval_list \
-I "$ALIGNDIR"/tumor.recal.bam \
-I "$ALIGNDIR"/normal.recal.bam \
-normal normal \
--germline-resource "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
--panel-of-normals "$RESOURCEDIR"/1000g_pon.hg38.subset.vcf.gz \
-O "$VARIANTDIR"/somatic.vcf.gz
