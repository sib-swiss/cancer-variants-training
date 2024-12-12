#!/usr/bin/env bash

cd ~/project/course_data/alignments

for sample in tumor normal
do
    samtools flagstat "$sample".recal.bam > "$sample".recal.bam.flagstat
done