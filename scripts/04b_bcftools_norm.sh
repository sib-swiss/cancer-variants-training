#!/usr/bin/env bash

# Define directories
REFDIR="${HOME}/project/course_data/reference"
VARIANTDIR="${HOME}/project/course_data/variants"

bcftools norm -m -any \
-f "$REFDIR"/ref_genome.fa \
-Oz -o "$VARIANTDIR"/somatic.filtered.norm.vcf.gz \
"$VARIANTDIR"/somatic.filtered.vcf.gz

bcftools index --tbi "$VARIANTDIR"/somatic.filtered.norm.vcf.gz
