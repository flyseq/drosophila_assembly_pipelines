#! /bin/bash
# call variants

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

# inputs
pileup="${sp}.mpileup"                     # read pileup file name

# outputs
vcf_all="${sp}.vcf"                        # VCF containing all sites

# make pileups for each contig
bcftools call -m -Ov --threads ${threads} -f GQ -o ${vcf_all} ${sp}.mpileup