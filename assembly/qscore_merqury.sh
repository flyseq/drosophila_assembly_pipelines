#! /bin/bash

cat ~/reads_list.txt |
while read line; do
sp=$(echo ${line} | sed -E 's/([A-Za-z0-9\.]+)\.passReads.+/\1/' | sed 's/.R10//') 
genome="${sp}.assembly.sm.fasta"
threads="75"
read1="${sp}_R1.fastq"
read2="${sp}_R2.fastq"

# copy file
zcat ~/active-data/repeatMasker/fastas/${genome}.gz > ./${genome}
zcat ~/active-data/illuminaReads/${sp}/${read1}.gz > ./${read1}
zcat ~/active-data/illuminaReads/${sp}/${read2}.gz > ./${read2}

# STEP 1: prepare meryl dbs

# 1.1 get optimal kmer size
#round
#kmer=$($MERQURY/best_k.sh $(cat ${genome} | grep -v ">" | tr -d '\n' | wc -c) | tail -n1 | awk '{print int($1+0.5)}')

#round up
kmer=$($MERQURY/best_k.sh $(cat ${genome} | grep -v ">" | tr -d '\n' | wc -c) | tail -n1 | awk '{ if($1-int($1) > 0){print int($1) + 1} else {print $1}}' ) 

# 1.2 build kmer db

# trim adapters
~/tools/bbmap/bbduk.sh \
    in=${read1} \
    in2=${read2} \
    out=${sp}_clean_R1.fastq \
    out2=${sp}_clean_R2.fastq \
    outs=${sp}_clean_unpaired.fastq \
    ref=~/tools/bbmap/resources/nextera.fa.gz,kapa \
    ktrim=r k=23 mink=11 hdist=1 ftm=5 threads=${threads} \
    tpe tbo

# build kmer dbs
meryl k=${kmer} count threads=${threads} output read1.meryl ${sp}_clean_R1.fastq
meryl k=${kmer} count threads=${threads} output read2.meryl ${sp}_clean_R2.fastq
meryl k=${kmer} count threads=${threads} output unpaired.meryl ${sp}_clean_unpaired.fastq

# merge kmer dbs
meryl union-sum threads=${threads} output ${sp}.meryl read?.meryl unpaired.meryl

# STEP 2: assess kmer accuracy and recovery
merqury ${sp}.meryl ${genome} ${sp}.merqury > ${sp}.merqury.log

mkdir -p ~/Dropbox/writing/2020_drosophila_genomes_manuscript/merqury/${sp}
cp ${sp}.merqury* ~/Dropbox/writing/2020_drosophila_genomes_manuscript/merqury/${sp}/

rm -r ./*
done