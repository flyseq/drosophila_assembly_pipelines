#! /bin/bash

ls ~/working/updated_genomes \
 | sed 's/.assembly.sm.fasta.gz//' \
 > ~/active-data/genomesize/genomes_list.txt

cat ~/active-data/genomesize/genomes_list.txt |
while read sp; do
cp ~/active-data/buscos_v4/${sp}.buscov4.tar.gz ./
tar zxf ${sp}.buscov4.tar.gz
cp ~/active-data/reads/${sp}.passReads.guppy???.fastq.gz ./
reads=$(ls ${sp}.passReads.guppy???.fastq.gz | tr '\n' ' ' | cut -d' ' -f1)
cat \
  ${sp}.buscov4/run_diptera_odb10/busco_sequences/single_copy_busco_sequences/*.fna \
  > ${sp}.buscocds.fasta

ncomp=$(cat ${sp}.buscocds.fasta | grep ">" | wc -l)

minimap2 -ax map-ont --secondary no -t 78 ${sp}.buscocds.fasta ${reads} \
    | samtools sort -@78 -o reads_to_busco.sam

stats=$(samtools depth -aa reads_to_busco.sam \
         | awk '{OFS=","}{sum+=$3} END { print sum,NR }')

readl=$(zcat ${reads} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c )

mkdir -p ~/active-data/genomesize
echo "${sp},${ncomp},${stats},${readl}" \
  > ~/active-data/genomesize/${sp}.coverage.csv

rm -r ./*
done