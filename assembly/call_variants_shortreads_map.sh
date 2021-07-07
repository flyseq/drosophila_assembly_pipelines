#! /bin/bash
# map Illumina reads to assembly

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

# inputs
assm="${sp}.assembly.sm.fasta"             # assembly filename
read1="${sp}_R1.fastq.gz"                  # Illumina PE reads, forward
read2="${sp}_R2.fastq.gz"                  # Illumina PE reads, reverse

# outputs
bam="${sp}.readsToDraft.bam"               # BAM file name

# map short reads to genome
minimap2 -ax sr -t ${threads} ${assm} ${read1} ${read2} \
  | samtools sort -@${threads} -O BAM -o ${sp}.temp.bam

# remove duplicates
sambamba markdup -r -p -t${threads} ${sp}.temp.bam ${bam} \
 && rm ${sp}.temp.bam