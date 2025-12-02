#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
READDIR="$HOME/project/course_data/reads"
RESOURCEDIR="${HOME}/project/course_data/resources"

for sample in tumor normal
do
    gatk BaseRecalibrator \
    --reference  "$REFDIR"/ref_genome.fa \
    --input "$ALIGNDIR"/"$sample".rg.md.bam \
    --known-sites "$RESOURCEDIR"/1000G_phase1.snps.high_confidence.hg38.subset.vcf.gz \
    --known-sites "$RESOURCEDIR"/Mills_and_1000G_gold_standard.indels.hg38.subset.vcf.gz \
    --output "$ALIGNDIR"/"$sample".bqsr.recal.table

    gatk ApplyBQSR \
    --input "$ALIGNDIR"/"$sample".rg.md.bam \
    --bqsr-recal-file "$ALIGNDIR"/"$sample".bqsr.recal.table \
    --output "$ALIGNDIR"/"$sample".recal.bam
done