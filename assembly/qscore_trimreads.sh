#! /bin/bash

# Short read adapter trimming in preparation for running Merqury

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

assm="${sp}.assembly.sm.fasta"     # filename of assembly to qc
read1="${sp}_R1.fastq"             # Illumina PE reads, forward
read2="${sp}_R2.fastq"             # Illumina PE reads, reverse

# trim Illumina adapters
/tools/bbmap/bbduk.sh \
    in=${read1} \
    in2=${read2} \
    out=${sp}_clean_R1.fastq \
    out2=${sp}_clean_R2.fastq \
    outs=${sp}_clean_unpaired.fastq \
    ref=/tools/bbmap/resources/nextera.fa.gz,kapa \
    ktrim=r k=23 mink=11 hdist=1 ftm=5 threads=${threads} \
    tpe tbo