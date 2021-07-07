#! /bin/bash

# this script generates a k-mer count histogram with Jellyfish for use
# with GenomeScope

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

reads="${sp}_R1.fastq.gz"          # Illumina PE reads, forward
read2="${sp}_R1.fastq.gz"          # Illumina PE reads, reverse

outFile="${sp}.buscoCoverage.csv"  # filename to write BUSCO output to

# make kmer histogram
jellyfish count -C -m 21 -s 16000000000 -t ${threads} \
    ${read1} ${read2} -o ${sp}.reads.jf
jellyfish histo -t${threads} ${sp}.reads.jf > ${sp}.reads.hist

# sample GenomeScope R command
# genomescope.R ${sp}.reads.hist 21 150 ${sp}.gsoutput