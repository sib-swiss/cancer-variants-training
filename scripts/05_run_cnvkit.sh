#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
RESOURCEDIR="${HOME}/project/course_data/resources"
VARIANTDIR="${HOME}/project/course_data/variants"

cnvkit.py batch "$ALIGNDIR"/tumor.recal.bam \
--normal "$ALIGNDIR"/normal.recal.bam \
--targets "$REFDIR"/exome_regions.bed.interval_list \
--fasta "$REFDIR"/ref_genome.fa \
--annotate "$RESOURCEDIR"/refFlat.txt \
--output-reference "$VARIANTDIR"/reference.cnn \
--output-dir "$VARIANTDIR"/cnvkit/ \
--processes 4 \
--scatter \
--diagram