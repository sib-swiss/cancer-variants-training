
#!/usr/bin/env bash

ALIGNDIR=~/project/data/alignments
RESOURCEDIR=~/project/data/resources

gatk CalculateContamination \
-I "$VARIANTDIR"/tumor.pileups.table \
-matched "$VARIANTDIR"/normal.pileups.table \
-O "$VARIANTDIR"/contamination.table
