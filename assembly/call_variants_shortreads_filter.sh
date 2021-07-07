#! /bin/bash
# quality score and repeat region filtering of variant calls

# job parameters
sp="D.melanogaster"               # species/sample name
threads="16"                      # number of threads to use
qv=30                             # minimum site/GT QV threshold

# inputs
vcf="${sp}.vcf"                   # VCF containing all sites
repeats="${sp}.repeats.bed"       # BED intervals of repeats (from RepeatMasker)

# outputs
vcf_snps="${sp}_shortread_snps.vcf"           # VCF containing all sites
vcf_indels="${sp}_shortread_indels.vcf"       # VCF containing all sites
vcf_all="${sp}_shortread_allsites.vcf"        # VCF containing all sites

# pull out SNPs passing filters, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types snps -m 2 -M 2 --threads ${threads} ${vcf} \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_snps}

# pull out indels passing filters, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types indels -m 2 -M 2 --threads ${threads} ${vcf} \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_indels}

# pull out all sites passing filters, exclude repeats
bcftools view -i "%QUAL>=${q}" --threads ${threads} ${vcf} \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_all}