#! /bin/bash

# This script is for computing consensus quality scores with Merqury
# and assumes Merqury is installed in a Conda environment

# job parameters
sp="D.melanogaster"                       # sample name/ID
threads="32"                              # number of threads to use

assm="${sp}.assembly.sm.fasta"            # filename of assembly to qc
read1_trim=${sp}_clean_R1.fastq           # filename of adapter trimmed R1
read2_trim=${sp}_clean_R2.fastq           # filename of adapter trimmed R2
unpaired_trim=${sp}_clean_unpaired.fastq  # filename of adapter trimmed unpaired
outFile="${sp}.merqury.log"               # file to write merqury output to

# activate conda env
conda activate merqury

# get optimal kmer size
# round up
kmer=$( $MERQURY/best_k.sh $(cat ${genome} | grep -v ">" | tr -d '\n' | wc -c) \
         | tail -n1 \
         | awk '{ if($1-int($1) > 0){print int($1) + 1} else {print $1}}' ) 

# build kmer dbs
meryl k=${kmer} count threads=${threads} output read1.meryl \
    ${read1_trim}
meryl k=${kmer} count threads=${threads} output read2.meryl \
    ${read2_trim}
meryl k=${kmer} count threads=${threads} output unpaired.meryl \
    ${unpaired_trim}

# merge kmer dbs
meryl union-sum threads=${threads} output ${sp}.meryl read?.meryl unpaired.meryl

# assess kmer accuracy and recovery, write output to file
merqury ${sp}.meryl ${genome} ${sp}.merqury > ${outFile}

# exit conda env
conda deactivate