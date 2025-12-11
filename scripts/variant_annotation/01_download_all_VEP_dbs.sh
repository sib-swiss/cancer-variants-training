#!/usr/bin/env bash

# Downloads ClinVar, AlphaMissense, REVEL for VEP v115

set -e

# Setup directories
BASE_DIR="${HOME}/project/course_data/VEP_dbs"
mkdir -p "${BASE_DIR}"/{clinvar,alphamissense,revel}

echo "=== Downloading ClinVar ==="
cd "${BASE_DIR}/clinvar"
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -c https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

echo "=== Downloading AlphaMissense ==="
cd "${BASE_DIR}/alphamissense"
wget -c https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz

echo "=== Downloading REVEL ==="
echo "This takes around >15 min⏰"
cd "${BASE_DIR}/revel"
wget -c https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip
unzip -o revel-v1.3_all_chromosomes.zip
cat revel_with_transcript_ids | tr "," "\t" > tabbed_revel.tsv
sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel.tsv
bgzip -f new_tabbed_revel.tsv

# Prepare for GRCh38
zcat new_tabbed_revel.tsv.gz | head -n1 > h
zgrep -h -v ^#chr new_tabbed_revel.tsv.gz | awk '$3 != "." ' | sort -k1,1 -k3,3n - | cat h - | bgzip -c > new_tabbed_revel_grch38.tsv.gz
tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_grch38.tsv.gz

# Cleanup
rm -f revel_with_transcript_ids tabbed_revel.tsv h revel-v1.3_all_chromosomes.zip

echo ""
echo "=== Done! ==="
echo "ClinVar:       ${BASE_DIR}/clinvar/clinvar.vcf.gz"
echo "AlphaMissense: ${BASE_DIR}/alphamissense/AlphaMissense_hg38.tsv.gz"
echo "REVEL:         ${BASE_DIR}/revel/new_tabbed_revel_grch38.tsv.gz"