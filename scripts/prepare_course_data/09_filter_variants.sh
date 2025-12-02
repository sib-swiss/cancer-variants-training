#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
VARIANTDIR="$HOME/project/course_data/variants"
RESOURCEDIR="${HOME}/project/course_data/resources"

gatk FilterMutectCalls \
-R "$REFDIR"/ref_genome.fa \
-V "$VARIANTDIR"/somatic.vcf.gz \
-O "$VARIANTDIR"/somatic.filtered.vcf.gz

bcftools view -f PASS "$VARIANTDIR"/somatic.filtered.vcf.gz -Oz > "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz
bcftools index --tbi "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz

zcat "$VARIANTDIR"/somatic.filtered.vcf.gz | grep -v "^#" | cut -f 7 | sort | uniq -c | sort -nr | head -n 10
zcat "$VARIANTDIR"/somatic.vcf.gz | grep -v "^#" | wc -l