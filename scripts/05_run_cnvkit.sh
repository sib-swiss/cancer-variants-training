#!/usr/bin/env bash

ALIGNDIR=~/project/course_data/alignments
REFDIR=~/project/course_data/reference
RESOURCEDIR=~/project/course_data/resources
VARIANTDIR=~/project/course_data/variants

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