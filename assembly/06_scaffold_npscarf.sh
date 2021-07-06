#! /bin/bash
# commands to scaffold a genome with long reads and npScarf

# job parameters
sp="D.melanogaster"                       # species/sample name
threads="16"                              # number of cores/threads

reads="${sp}.passReads.fastq.gz"          # Nanopore reads file
assembly="${sp}.clip.fasta"               # haplotig purged and clipped genome

# map reads to draft genome
minimap2 -ax map-ont --secondary=no -t ${threads} ${assembly} ${reads} \
 | samtools sort -@ ${threads} -m 1G -o reads_to_draft.bam 

# run npScarf
/tools/japsa/bin/jsa.np.npscarf -seq ${assembly} \
    -input reads_to_draft.bam -format bam --support=4 --qual=8

# clean up
rm reads_to_draft.bam