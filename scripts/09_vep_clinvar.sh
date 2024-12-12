#!/usr/bin/env bash

# Download ClinVar VCF

VARIANTDIR=~/project/course_data/variants
RESOURCEDIR=/config/project/course_data/resources

# Run VEP with cancer-specific annotations
vep -i "$VARIANTDIR"/somatic.filtered.PASS.vcf.gz \
    --cache \
    --dir /data/.vep \
    --assembly GRCh38 \
    --format vcf \
    --force_overwrite \
    --fork 2 \
    --everything \
    --custom "$RESOURCEDIR"/clinvar/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --output_file "$VARIANTDIR"/somatic.filtered.PASS2_annotated_clinvar.txt

# Filter results
filter_vep \
-i "$VARIANTDIR"/somatic.filtered.PASS2_annotated_clinvar.txt \
-o "$VARIANTDIR"/somatic.filtered.PASS2_annotated_clinvar_filtered.txt \
--force_overwrite \
--filter "(IMPACT is HIGH or IMPACT is MODERATE) or \
(SIFT match deleterious or PolyPhen match probably_damaging) or \
(ClinVar_CLNSIG match pathogenic)"

# Generate summary (using simple grep/awk commands)
echo "Mutation Type Summary:" > summary.txt
grep -v "#" "$VARIANTDIR"/somatic.filtered.PASS2_annotated_clinvar_filtered.txt | cut -f7 | sort | uniq -c >> summary.txt
