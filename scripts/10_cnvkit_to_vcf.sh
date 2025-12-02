#!/usr/bin/env bash

VARIANTDIR="${HOME}/project/course_data/variants/cnvkit"

cnvkit.py export vcf "$VARIANTDIR"/tumor.call.cns -o "$VARIANTDIR"/tumor.call.vcf

bgzip "$VARIANTDIR"/tumor.call.vcf 
tabix -p vcf "$VARIANTDIR"/tumor.call.vcf.gz
