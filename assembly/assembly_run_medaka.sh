#! /bin/bash
# polish with one round of Medaka

# job parameters
sp="D.melanogaster"                   # sample name/ID
threads="32"                          # number of threads to use

# activate medaka virtual env
. /medaka/bin/activate

# special setting for RTX 2000 series
# for details, go to https://nanoporetech.github.io/medaka/installation.html
export TF_FORCE_GPU_ALLOW_GROWTH=true 

#run medaka
medaka_consensus -i ${sp}.passReads.fastq.gz -d ${sp}.assembly.racon.fasta \
    -t ${threads} -m r941_min_high -b 40 \
    -o ${sp}.medaka

#copy draft to new file
cp ${sp}.medaka/consensus.fasta ./${sp}.assembly.medaka.fasta