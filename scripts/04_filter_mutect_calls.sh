#!/usr/bin/env bash

REFDIR=~/project/course_data/reference
VARIANTDIR=~/project/course_data/variants

gatk FilterMutectCalls \
-R "$REFDIR"/ref_genome.fa \
-V "$VARIANTDIR"/somatic.vcf.gz \
-O "$VARIANTDIR"/somatic.filtered.vcf.gz

bcftools view -f PASS "$VARIANTDIR"/somatic.filtered.vcf.gz -Oz \
> "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz

bcftools index --tbi "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz
