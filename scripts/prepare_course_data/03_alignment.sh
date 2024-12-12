
REFDIR=~/project/data/reference/
READDIR=~/project/data/reads
ALIGNDIR=~/project/data/alignments

mkdir -p "$ALIGNDIR"

for sample in tumor normal
do
    bwa mem \
    "$REFDIR"/ref_genome.fa \
    "$READDIR"/"$sample"_R1.fastq.gz \
    "$READDIR"/"$sample"_R2.fastq.gz \
    2> "$ALIGNDIR"/$sample.bwa.log \
    | samtools sort \
    | samtools view -bh \
    > "$ALIGNDIR"/"$sample".bam
done