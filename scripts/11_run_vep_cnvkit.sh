
VARIANTDIR=~/project/course_data/variants
RESOURCEDIR=/config/project/course_data/resources

vep -i "$VARIANTDIR"/cnvkit/tumor.call.vcf.gz \
    --cache \
    --dir /data/.vep \
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
    --output_file "$VARIANTDIR"/cnvkit/tumor.call_annotated_all.txt

filter_vep \
    -i "$VARIANTDIR"/cnvkit/tumor.call_annotated_all.txt \
    -o "$VARIANTDIR"/cnvkit/tumor.call_annotated_all_filtered.txt \
    --force_overwrite \
    --filter "(IMPACT is HIGH or IMPACT is MODERATE) and \
        (REVEL >= 0.75 or \
            am_class = 'likely_pathogenic' or \
            SIFT = 'deleterious' and PolyPhen = 'probably_damaging')"