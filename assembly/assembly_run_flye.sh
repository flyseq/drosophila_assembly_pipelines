#! /bin/bash
# assemble Nanopore reads with Flye

# job parameters
sp="D.melanogaster"      # sample name/ID
gSize="140m"             # initial guess for genome size
threads="32"             # number of threads to use

# run flye
flye --nano-raw ${sp}.passReads.fastq.gz \
     --genome-size ${gSize} \
     --threads ${threads} \
     --out-dir ${sp}.FlyeAssembly

# copy finished assembly
cp ${sp}.FlyeAssembly/assembly.fasta ${sp}.FlyeAssembly.fasta