#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"
REFDIR="${HOME}/project/course_data/reference"
VARIANTDIR="$HOME/project/course_data/variants"
RESOURCEDIR="${HOME}/project/course_data/resources"

gatk CalculateContamination \
-I "$VARIANTDIR"/tumor.pileups.table \
-matched "$VARIANTDIR"/normal.pileups.table \
-O "$VARIANTDIR"/contamination.table
