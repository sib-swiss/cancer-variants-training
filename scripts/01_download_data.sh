
mkdir -p data/resources

cd data/resources

# panel of normals

aws s3 \
--no-sign-request --region eu-west-1 \
cp \
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz \
.

aws s3 \
--no-sign-request --region eu-west-1 \
cp \
s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi \
.

cd ../
mkdir reference

cd reference

# exome intervals
wget https://genomedata.org/pmbio-workshop/references/exome/chr6_and_chr17/exome_regions.bed
wget https://genomedata.org/pmbio-workshop/references/exome/chr6_and_chr17/exome_regions.bed.interval_list

# reference genome
wget https://genomedata.org/pmbio-workshop/references/genome/chr6_and_chr17/ref_genome.tar
tar xvf ref_genome.tar
rm ref_genome.tar

# reads
cd ../
mkdir reads
cd reads

wget https://genomedata.org/pmbio-workshop/fastqs/chr6_and_chr17/Exome_Norm.tar
wget https://genomedata.org/pmbio-workshop/fastqs/chr6_and_chr17/Exome_Tumor.tar

tar xvf Exome_Norm.tar
rm Exome_Norm.tar
tar xvf Exome_Tumor.tar
rm Exome_Tumor.tar

mv Exome_Norm/Exome_Norm_R1.fastq.gz normal_R1.fastq.gz
mv Exome_Norm/Exome_Norm_R2.fastq.gz normal_R2.fastq.gz

mv Exome_Tumor/Exome_Tumor_R1.fastq.gz tumor_R1.fastq.gz
mv Exome_Tumor/Exome_Tumor_R2.fastq.gz tumor_R2.fastq.gz

rm -r Exome_Norm
rm -r Exome_Tumor

# subset vcf
cd ../resources
bcftools view -Oz -r chr6,chr17 af-only-gnomad.hg38.vcf.gz > af-only-gnomad.hg38.subset.vcf.gz
bcftools index --tbi af-only-gnomad.hg38.subset.vcf.gz

bcftools view -Oz -r chr6,chr17 1000g_pon.hg38.vcf.gz > 1000g_pon.hg38.subset.vcf.gz
bcftools index --tbi 1000g_pon.hg38.subset.vcf.gz