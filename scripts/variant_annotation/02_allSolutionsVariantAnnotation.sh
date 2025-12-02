#!/usr/bin/env bash

# Base directories
BASE_DIR="${HOME}/project/course_data"
RES_DIR="${BASE_DIR}/resources"
DBS_DIR="${BASE_DIR}/VEP_dbs"
VAR_DIR="${BASE_DIR}/variants"
CNV_DIR="${VAR_DIR}/cnvkit"
VEP_CACHE="${HOME}/.vep"

# Input Files
MUTECT_VCF="${VAR_DIR}/somatic.filtered.PASS.vcf.gz"
CLINVAR="${DBS_DIR}/clinvar/clinvar.vcf.gz"
REVEL="${DBS_DIR}/revel/new_tabbed_revel_grch38.tsv.gz"
ALPHA="${DBS_DIR}/alphamissense/AlphaMissense_hg38.tsv.gz"

# Create output directories
mkdir -p "${VAR_DIR}"
mkdir -p "${CNV_DIR}"

# Activate environment
mamba activate ngs-tools

# The VCF manual file
echo -e '##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr7\t55242465\trs121913529\tA\tT\t.\tPASS\t.
chr17\t7577121\trs28934578\tC\tT\t.\tPASS\t.
chr13\t32936646\trs28897743\tC\tT\t.\tPASS\t.
chr17\t7674894\trs28934574\tG\tA\t.\tPASS\t.
chr13\t32914438\trs80357090\tT\tC\t.\tPASS\t.
chr17\t41245466\trs28897686\tG\tA\t.\tPASS\t.
chr3\t178936091\trs63751289\tG\tA\t.\tPASS\t.
chr10\t89624230\trs61751507\tG\tA\t.\tPASS\t.
chr11\t108098576\trs386833395\tG\tA\t.\tPASS\t.
chr4\t55141055\trs1800744\tA\tG\t.\tPASS\t.
chr1\t11190510\trs6603781\tG\tA\t.\tPASS\t.
chr7\t140453136\trs121913530\tA\tT\t.\tPASS\t.
chr17\t43094464\trs28897672\tA\tC\t.\tPASS\t.
chr5\t112175770\trs121913343\tC\tT\t.\tPASS\t.' >"$VAR_DIR"/cancer_variants.vcf

# 4. Verify the file was created correctly:
# Check the content
cat "$VAR_DIR"/cancer_variants.vcf

# Make sure you have 13 variants
grep -v '^#' "$VAR_DIR"/cancer_variants.vcf | wc -l # it should show 14

# 1. BASIC VEP CALLS
# Initial call (Example)
vep -i "${VAR_DIR}"/cancer_variants.vcf \
  --cache "${VEP_CACHE}" \
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
  --output_file "${VAR_DIR}/output.txt"

# With everything
vep -i "${MUTECT_VCF}" \
  --cache "${VEP_CACHE}" \
  --species homo_sapiens \
  --assembly GRCh38 \
  --offline \
  --format vcf \
  --symbol \
  --force_overwrite \
  --everything \
  --output_file "${VAR_DIR}/chr6ch17_mutect2_annotated.txt"

# With --fork 2
vep -i "${MUTECT_VCF}" \
  --cache "${VEP_CACHE}" \
  --species homo_sapiens \
  --assembly GRCh38 \
  --offline \
  --format vcf \
  --symbol \
  --force_overwrite \
  --everything \
  --fork 2 \
  --output_file "${VAR_DIR}/chr6ch17_mutect2_annotated.txt"

# First filtering
filter_vep \
  -i "${VAR_DIR}/chr6ch17_mutect2_annotated.txt" \
  -o "${VAR_DIR}/chr6ch17_mutect2_annotated_filtered.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) or (SIFT match deleterious or PolyPhen match probably_damaging)"

# 2. VEP WITH CUSTOM ANNOTATION (ClinVar)
# Run VEP with cancer-specific annotations
vep -i "${MUTECT_VCF}" \
  --cache "${VEP_CACHE}" \
  --assembly GRCh38 \
  --format vcf \
  --force_overwrite \
  --fork 2 \
  --everything \
  --custom "${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" \
  --output_file "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar.txt"

# Filter results
filter_vep \
  -i "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar.txt" \
  -o "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar_filtered.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) or \
    (SIFT match deleterious or PolyPhen match probably_damaging) or \
    (ClinVar_CLNSIG match pathogenic)"

# Generate summary
echo "Mutation Type Summary:" >"${VAR_DIR}/summary.txt"
grep -v "#" "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar_filtered.txt" | cut -f7 | sort | uniq -c >>"${VAR_DIR}/summary.txt"
# 3. COMPLEX VEP (Plugins)
# use all course plugins, but there are more!
vep -i "${MUTECT_VCF}" \
  --cache "${VEP_CACHE}" \
  --assembly GRCh38 \
  --format vcf \
  --fork 2 \
  --everything \
  --af_gnomade \
  --force_overwrite \
  --plugin SpliceRegion \
  --plugin REVEL,file="${REVEL}" \
  --plugin AlphaMissense,file="${ALPHA}" \
  --custom "${CLINVAR},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN" \
  --output_file "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar_and_plugins.txt"

# Filter the data -> Stringent
filter_vep \
  -i "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar_and_plugins.txt" \
  -o "${VAR_DIR}/stringent_filtered.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) and \
        (REVEL >= 0.75 or \
         am_class = 'likely_pathogenic' or \
         SIFT = 'deleterious' and PolyPhen = 'probably_damaging')"

# Filter the data -> Moderate
filter_vep \
  -i "${VAR_DIR}/chr6ch17_mutect2_annotated_with_clinVar_and_plugins.txt" \
  -o "${VAR_DIR}/moderate_filtered.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) and \
           (REVEL >= 0.5 or \
           am_class = 'likely_pathogenic' or \
           SIFT = 'deleterious' or \
           PolyPhen = 'probably_damaging')"

# 4. CNV ANALYSIS (added as supplementary)
# # Exercise CNV setup
# tar -xvzf "${RES_DIR}/CNVkit_data/cnvkit.tar.gz" -C "${CNV_DIR}" --strip-components=1

# Convert CNS to VCF
cnvkit.py export vcf "${CNV_DIR}/tumor.call.cns" -o "${CNV_DIR}/tumor.call.vcf"

bgzip -f "${CNV_DIR}/tumor.call.vcf"

tabix -p vcf -f "${CNV_DIR}/tumor.call.vcf.gz"

# Call VEP with full calls (overlap)
vep -i "${CNV_DIR}/tumor.call.vcf.gz" \
  --cache "${VEP_CACHE}" \
  --assembly GRCh38 \
  --format vcf \
  --fork 2 \
  --max_sv_size 100000000 \
  --everything \
  --regulatory \
  --af_gnomade \
  --force_overwrite \
  --plugin SpliceRegion \
  --plugin REVEL,file="${REVEL}" \
  --plugin AlphaMissense,file="${ALPHA}" \
  --custom "${CLINVAR},ClinVar,vcf,overlap,0,CLNSIG,CLNREVSTAT,CLNDN" \
  --output_file "${CNV_DIR}/tumor.call_annotated_all_overlap.txt"

# Call VEP per gene (exact)
vep -i "$CNV_DIR"/tumor.call.vcf.gz \
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
  --plugin REVEL,file="${REVEL}" \
  --plugin AlphaMissense,file="${ALPHA}" \
  --custom "${CLINVAR}",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
  --output_file "${CNV_DIR}/tumor.call_annotated_all_exact.txt"

## Exact
# Filter the data CNV
filter_vep \
  -i "${CNV_DIR}/tumor.call_annotated_all_exact.txt" \
  -o "${CNV_DIR}/tumor.call_annotated_all_IMPACT_filtered_exact.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) or (ClinVar_CLNSIG match pathogenic)"

# filtering by Symbol
filter_vep \
  -i "${CNV_DIR}/tumor.call_annotated_all_exact.txt" \
  -o "${CNV_DIR}/tumor.call_annotated_GENE_filtered_exact.txt" \
  --filter "SYMBOL exists" \
  --force_overwrite

## Overalp
# Filter the data CNV
filter_vep \
  -i "${CNV_DIR}/tumor.call_annotated_all_overlap.txt" \
  -o "${CNV_DIR}/tumor.call_annotated_all_IMPACT_filtered__overlap.txt" \
  --force_overwrite \
  --filter "(IMPACT is HIGH or IMPACT is MODERATE) or (ClinVar_CLNSIG match pathogenic)"

# filtering by Symbol
filter_vep \
  -i "${CNV_DIR}/tumor.call_annotated_all_overlap.txt" \
  -o "${CNV_DIR}/tumor.call_annotated_GENE_filtered__overlap.txt" \
  --filter "SYMBOL exists" \
  --force_overwrite
