#! /bin/bash

# This script follows the procedure for running purge_haplotigs as detailed in:
# https://bitbucket.org/mroachawri/purge_haplotigs/src/master/
#
# IMPORTANT: the cutoffs provide to purge_haplotigs cov must be set manually.
#            Please see instructions at the above link for the instructions.

# job parameters
sp="D.melanogaster"                            # species/sample name
threads="16"                                   # number of cores/threads

reads="${sp}.passReads.fastq.gz"               # Nanopore reads file
assembly="${sp}.FlyeAssembly.fasta"            # genome file

# map reads to assembly
minimap2 -ax map-ont -t ${threads} ${assembly} ${reads} --secondary=no \
 | samtools sort -m 1G -o ${sp}.aligned.bam -T ${sp}.tmp.ali

# activate the purge_haplotigs env
conda activate purge_haplotigs

# generate the coverage histogram
purge_haplotigs hist -b ${sp}.aligned.bam -g ${assembly} -t ${threads}

# flag contigs for removal based on coverage
# NOTE: -l -m -h arguments are obtained from the coverage histogram!
purge_haplotigs cov -i ${sp}.aligned.bam -l 1 -m 20 -h 195

# purge haplotigs
purge_haplotigs purge -g ${genome} -c coverage_stats.csv -o ${sp}.purged

# trim overlapping contig ends
purge_haplotigs clip -p ${sp}.purged.fasta -h ${sp}.purged.haplotigs.fasta \
  -o ${sp}.clip

# The purged and clipped output file is ${sp}.clip.fasta 
# (e.g. D.melanogaster.clip.fasta)