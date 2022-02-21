#! /bin/bash
# call variants

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

# inputs
pileup="${sp}.mpileup"                     # read pileup file name

# outputs
vcf="${sp}.vcf"                        # VCF containing all sites

# call variants for all genomic sites
bcftools call -m -Ov --threads ${threads} -f GQ -o ${vcf} ${sp}.mpileup
