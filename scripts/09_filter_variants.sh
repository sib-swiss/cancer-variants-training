#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
REFDIR=/config/data/reference
RESOURCEDIR=/config/data/resources
VARIANTDIR=/config/data/variants

gatk FilterMutectCalls \
-R "$REFDIR"/ref_genome.fa \
-V "$VARIANTDIR"/somatic.vcf.gz \
-O "$VARIANTDIR"/somatic.filtered.vcf.gz