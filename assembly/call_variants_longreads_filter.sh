#! /bin/bash
# quality score and repeat region filtering of variant calls

# job parameters
sp="D.melanogaster"               # species/sample name
threads="16"                      # number of threads to use
qv=30                             # minimum site/GT QV threshold

# inputs
repeats="${sp}.repeats.bed"       # BED intervals of repeats (from RepeatMasker)

# outputs
vcf_snps="${sp}_longread_snps.vcf"           # VCF containing all sites
vcf_indels="${sp}_longread_indels.vcf"       # VCF containing all sites
vcf_all="${sp}_longread_allsites.vcf"        # VCF containing all sites

# pull out SNPs, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types snps -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_snps}

# pull out indels, exclude repeats
bcftools view -i "%QUAL>=${q} && MIN(FMT/GQ)>=${qv}" \
    --types indels -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_indels}

# pull out all sites/intervals passing site quality
bcftools view -i "%QUAL>=${q}" --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.g.vcf.gz \
  | bedtools subtract -header -a /dev/stdin -b ${repeats} \
  > ${vcf_all}