VARIANTDIR=~/project/course_data/variants

cnvkit.py export vcf "$VARIANTDIR"/cnvkit/tumor.call.cns -o "$VARIANTDIR"/cnvkit/tumor.call.vcf

bgzip "$VARIANTDIR"/cnvkit/tumor.call.vcf 
tabix -p vcf "$VARIANTDIR"/cnvkit/tumor.call.vcf.gz
