#!/usr/bin/env bash

VARIANTDIR=~/project/course_data/variants

filter_vep \
-i "$VARIANTDIR"/somatic.filtered.PASS2_annotated.txt \
-o "$VARIANTDIR"/somatic.filtered.PASS2_annotated_fitlered.txt \
--force_overwrite \
--filter "(IMPACT is HIGH or IMPACT is MODERATE) or (SIFT match deleterious or PolyPhen match probably_damaging)"
