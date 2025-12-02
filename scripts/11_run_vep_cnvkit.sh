#!/usr/bin/env bash

VARIANTDIR="${HOME}/project/course_data/variants/cnvkit"
RESOURCEDIR="${HOME}/project/course_data/VEP_dbs"

vep -i "$VARIANTDIR"/tumor.call.vcf.gz \
    --cache \
    --dir $HOME/.vep \
    --assembly GRCh38 \
    --format vcf \
    --fork 2 \
    --max_sv_size 100000000 \
    --everything \
    --regulatory \
    --af_gnomade \
    --force_overwrite \
    --plugin SpliceRegion \
    --plugin REVEL,file="$RESOURCEDIR"/revel/new_tabbed_revel_grch38.tsv.gz \
    --plugin AlphaMissense,file="$RESOURCEDIR"/alphamissense/AlphaMissense_hg38.tsv.gz \
    --custom "$RESOURCEDIR"/clinvar/clinvar.vcf.gz,ClinVar,vcf,overlap,0,CLNSIG,CLNREVSTAT,CLNDN \
    --output_file "$VARIANTDIR"/tumor.call_annotated_all.txt

# Filter the data CNV
filter_vep \
    -i "${VARIANTDIR}/tumor.call_annotated_all.txt" \
    -o "${VARIANTDIR}/tumor.call_annotated_all_IMPACT_filtered.txt" \
    --force_overwrite \
    --filter "(IMPACT is HIGH or IMPACT is MODERATE) or (ClinVar_CLNSIG match pathogenic)"

# filtering by Symbol
filter_vep \
    -i "${VARIANTDIR}/tumor.call_annotated_all.txt" \
    -o "${VARIANTDIR}/chr6ch17_mtumor_annotated_GENE_filtered.txt" \
    --filter "SYMBOL exists" \
    --force_overwrite