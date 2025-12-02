#!/usr/bin/env bash

# Define directories
VARIANTDIR="${HOME}/project/course_data/variants"

# Run VEP ~< 5min
vep -i "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz \
    --cache \
    --dir "$HOME"/.vep \
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