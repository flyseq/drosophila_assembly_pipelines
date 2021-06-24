#! /bin/bash

mkdir -p ~/active-data/genomesize_kmer
echo $(ls ~/active-data/illuminaReads/) | tr ' ' '\n' > specieslist.txt
threads="60"

cat ~/active-data/genomesize_kmer/specieslist.txt |
while read sp; do
  echo "doing species: ${sp}"
  cp ~/active-data/illuminaReads/${sp}/${sp}_R?.fastq.gz ./
  ls *.fastq.gz | parallel -j2 gunzip {}
  
  # make kmer histogram
  nice jellyfish count -C -m 21 -s 50000000000 -t ${threads} *.fastq -o reads.jf
  nice jellyfish histo -t${threads} reads.jf > reads.histo

  # run GenomeScope
  len=$(($(cat ${sp}_R?.fastq | head -n2 | tail -1 | wc -c) - 1))
  genomescope.R reads.histo 21 ${len} ${sp}.genomescope \
    > ~/active-data/genomesize_kmer/${sp}.genomescope.txt
  mv ${sp}.genomescope ~/active-data/genomesize_kmer/${sp}.genomescope
  mv reads.histo ~/active-data/genomesize_kmer/${sp}.genomescope/

  rm ./*
done