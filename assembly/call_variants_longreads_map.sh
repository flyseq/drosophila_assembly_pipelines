#! /bin/bash
# map Nanopore reads to assembly

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

assm="${sp}.assembly.sm.fasta"             # assembly filename
bam="${sp}.readsToDraft.bam"            # BAM file name
reads="${sp}.passReads.fastq.gz"           # Nanopore reads

minimap2 -ax map-ont -t ${threads} ${assm} ${reads} \
  | samtools sort -@${threads} -O BAM -o ${bam}