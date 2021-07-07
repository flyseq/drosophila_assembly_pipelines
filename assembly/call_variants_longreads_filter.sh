#! /bin/bash
# quality score and repeat region filtering of variant calls

# job parameters
sp="D.melanogaster"               # species/sample name
threads="16"                      # number of threads to use
qv=30                             # minimum site/GT QV threshold

repeats="${sp}.repeats.bed"       # BED intervals of repeats (from RepeatMasker)

# pull out SNPs, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types snps -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${sp}_longread_snps.vcf

# pull out indels, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types indels -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${sp}_longread_indels.vcf

# pull out all sites/intervals passing site quality
bcftools view -i "%QUAL>=${q}" --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.g.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${sp}_longread_allsites.g.vcf