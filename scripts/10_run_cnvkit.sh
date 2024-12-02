#!/usr/bin/env bash

ALIGNDIR=~/project/data/alignments
REFDIR=~/project/data/reference
RESOURCEDIR=~/project/data/resources
VARIANTDIR=~/project/data/variants

cnvkit.py batch "$ALIGNDIR"/tumor.rg.md.bam \
--normal "$ALIGNDIR"/normal.rg.md.bam \
--targets "$REFDIR"/exome_regions.bed.interval_list \
--fasta "$REFDIR"/ref_genome.fa \
--annotate "$RESOURCEDIR"/refFlat.txt \
--output-reference "$VARIANTDIR"/reference.cnn \
--output-dir "$VARIANTDIR"/cnvkit/ \
--scatter \
--diagram