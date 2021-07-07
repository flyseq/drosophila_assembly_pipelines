#! /bin/bash
# generate pileup file for each contig in parallel

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

assm="${sp}.assembly.sm.fasta"             # assembly filename
bam="${sp}.readsToDraft.bam"               # BAM file name
pileup="${sp}.mpileup"                     # read pileup file name

# make pileups for each contig
cat ${assm} \
  | grep ">" \
  | tr -d ">" \
  | parallel -j ${threads} \
        bcftools mpileup -f ${assm} -r {} ${bam} \
        -q 10 -Q 20 -o sr_{}.mpileup

# combine pileups
bcftools concat --threads ${threads} sr_*.mpileup > ${sp}.mpileup \
 && rm sr_contig_*.mpileup