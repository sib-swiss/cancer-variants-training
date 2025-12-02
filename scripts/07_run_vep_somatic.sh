#!/usr/bin/env bash

# Define directories
VARIANTDIR="${HOME}/project/course_data/variants"

vep -i "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz \
--cache \
--dir "$HOME"/.vep \
--species homo_sapiens \
--assembly GRCh38 --offline --format vcf \
--symbol --force_overwrite -e \
--output_file "$VARIANTDIR"/somatic.filtered.PASS2_annotated.txt
