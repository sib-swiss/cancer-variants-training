#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
VARIANTDIR="$HOME/project/course_data/variants"
RESOURCEDIR="${HOME}/project/course_data/resources"

mkdir -p $VARIANTDIR

# takes about 12 minutes to run
gatk Mutect2 \
-R "$REFDIR"/ref_genome.fa \
--intervals "$REFDIR"/exome_regions.bed.interval_list \
-I "$ALIGNDIR"/tumor.recal.bam \
-I "$ALIGNDIR"/normal.recal.bam \
-normal normal \
--germline-resource "$RESOURCEDIR"/af-only-gnomad.hg38.subset.vcf.gz \
--panel-of-normals "$RESOURCEDIR"/1000g_pon.hg38.subset.vcf.gz \
-O "$VARIANTDIR"/somatic.vcf.gz
