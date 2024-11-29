
#!/usr/bin/env bash

ALIGNDIR=/config/data/alignments
RESOURCEDIR=/config/data/resources

gatk CalculateContamination \
-I "$VARIANTDIR"/tumor.pileups.table \
-matched "$VARIANTDIR"/normal.pileups.table \
-O "$VARIANTDIR"/contamination.table
