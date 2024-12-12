#!/usr/bin/env bash

VARIANTDIR=~/project/course_data/variants

# Run VEP ~< 5min
vep -i "$VARIANTDIR"/cancer_variants.vcf \
    -cache \
    --dir /data/.vep \
    --species homo_sapiens \
    --offline \
    --force_overwrite \
    --assembly GRCh38 \
    --format vcf \
    --symbol \
    --sift b \
    --polyphen b \
    --pick \
    --domains \
    --output_file "$VARIANTDIR"/cancer_variants_output.txt

