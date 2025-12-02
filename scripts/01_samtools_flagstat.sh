#!/usr/bin/env bash

ALIGNMENT_DIR="${HOME}/project/course_data/alignments"

# Process each sample
for sample in tumor normal; do
    samtools flagstat "${ALIGNMENT_DIR}/$sample.recal.bam" > "${ALIGNMENT_DIR}/$sample.recal.bam.flagstat"
done

echo "Flagstat analysis complete!"