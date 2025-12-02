#!/usr/bin/env bash

# Define directories
ALIGNDIR="${HOME}/project/course_data/alignments"

for sample in tumor normal
do
    gatk AddOrReplaceReadGroups \
    --INPUT "$ALIGNDIR"/"$sample".bam \
    --OUTPUT "$ALIGNDIR"/"$sample".rg.bam \
    --RGLB "$sample" \
    --RGPU HWI-ST466.C1TD1ACXX \
    --RGPL ILLUMINA \
    --RGSM "$sample" \
    --RGID HWI-ST466.C1TD1ACXX."$sample"
done