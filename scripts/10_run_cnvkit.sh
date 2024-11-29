#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
REFDIR=/config/data/reference
RESOURCEDIR=/config/data/resources
VARIANTDIR=/config/data/variants

cnvkit.py batch "$ALIGNDIR"/tumor.rg.md.bam \
--normal "$ALIGNDIR"/normal.rg.md.bam \
--targets "$REFDIR"/exome_regions.bed.interval_list \
--fasta "$REFDIR"/ref_genome.fa \
--annotate "$RESOURCEDIR"/refFlat.txt.gz \
--output-reference "$VARIANTDIR"/reference.cnn \
--output-dir "$VARIANTDIR"/cnvkit/