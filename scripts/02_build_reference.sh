
REFDIR=~/project/data/reference/

mkdir -p "$ALIGNDIR"

bwa index "$REFDIR"/ref_genome.fa
