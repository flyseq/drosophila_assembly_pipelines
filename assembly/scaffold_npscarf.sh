#! /bin/bash

# Commands to scaffold a genome with long reads and npScarf

# set up filenames/designators
reads="D.melanogaster.passReads.guppy322.fastq.gz" # Nanopore reads file
assembly="D.melanogaster.clip.fasta" # haplotig purged and clipped genome
sp="D.melanogaster" # species/sample name
threads="16" # number of cores/threads

# map reads to draft genome
minimap2 -ax map-ont --secondary=no -t ${threads} ${assembly} ${reads} \
 | samtools sort -@ ${threads} -m 1G -o reads_to_draft.bam 

# run npScarf
/tools/japsa/bin/jsa.np.npscarf -seq ${assembly} \
    -input reads_to_draft.bam -format bam --support=4 --qual=8