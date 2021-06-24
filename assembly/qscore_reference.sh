#! /bin/bash

# prepare genes to be run with Nanopore's Pomoxis
# code adapted from 
# https://github.com/nanoporetech/pomoxis/blob/master/scripts/assess_assembly
# custom chunking to assess at only genes
threads="75"
outDir="$HOME/active-data/qc2"

cat ${outDir}/genome_list.txt | \
while read assm refseq annot; do
  sp=$(echo ${assm} | sed 's/.assembly.sm.fasta//')
  repeat=$(echo ${refseq} | sed 's/_genomic.fna.gz/_rm.out.gz/')
  
  #download reference genome
  wget -q ${refseq} -O- | gunzip > ${sp}.ref.fasta

  #download reference annotations
  wget -q ${annot} -O- | gunzip > ${sp}.ref.gff

  #download reference repeats
  wget -q ${repeat} -O- | gunzip | awk '{OFS="\t"}{print $5,$6-1,$7}' \
    | tail +4 > ${sp}.repeats.bed

  #copy over nanopore assm
  zcat ~/active-data/repeatMasker/fastas/${assm}.gz > ${assm}
  
  # get longest transcripts from GFF
  conda activate agat
  agat_convert_sp_gxf2gxf.pl \
    --gff ${sp}.ref.gff --merge_loci -o ${sp}.ref.merged.gff
  agat_sp_keep_longest_isoform.pl -gff ${sp}.ref.merged.gff \
    -o ${sp}.longestIsoform.gff
  conda deactivate

  # get chrom sizes
  conda activate medaka
  faidx ${sp}.ref.fasta -i chromsizes | sort -k1,1 -k2,2n > ${sp}.chromSizes
  awk '{print $1"\t"0"\t"$2}' ${sp}.chromSizes | sort -k1,1 -k2,2n \
   > ${sp}.chromSizes.bed
  # get intergenic
  cat ${sp}.longestIsoform.gff \
   | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' \
   | awk '($3 != "region"){print $0}' > in_sorted.gff
  bedtools complement -i in_sorted.gff -g ${sp}.chromSizes \
   | awk '{print $0"\tNA"}'> ${sp}.intergenic.bed
  # get exons
  cat in_sorted.gff \
   | awk '{OFS="\t"} $1 ~ /^#/ {next} {if ($3 == "exon") {split($9, a, ";"); print $1, $4-1, $5, a[3]}}' \
   | sed -E 's/Dbxref=([A-Za-z=]+:[^,]+).+/\1/' > ${sp}.exons.bed
  # get introns
  bedtools complement \
    -i <(cat ${sp}.exons.bed ${sp}.intergenic.bed | sort -k1,1 -k2,2n) \
    -g ${sp}.chromSizes > ${sp}.introns.bed
  # get CDS
  cat in_sorted.gff \
   | awk '{OFS="\t"} $1 ~ /^#/ {next} {if ($3 == "CDS") {split($9, a, ";"); print $1, $4-1, $5, a[3]}}' \
   | sed -E 's/Dbxref=([A-Za-z=]+:[^,]+).+/\1/' > ${sp}.exons.bed
  # run Pomoxis on introns, intergenic
  rm ${sp}.ref.fasta.??i
  for region in "introns" "intergenic" "repeats"; do
    assess_assembly -r ${sp}.ref.fasta -i ${assm} -t ${threads} -c 100000 \
      -l 10000 -b ${sp}.${region}.bed > ${sp}.${region}.qv
  done
  # run Pomoxis on whole genome
  assess_assembly -r ${sp}.ref.fasta -i ${assm} -t ${threads} -c 100000 \
      -l 10000 > ${sp}.wholegenome.qv
  # run Pomoxis on exons
  assess_assembly -r ${sp}.ref.fasta -i ${assm} -t ${threads} -c 100000 \
      -l 1 -b ${sp}.exons.bed > ${sp}.exons.qv
  
  # get indel BED intervals
  cat assm_indel_ge1.txt | grep -vE "^type" \
   | awk '{OFS="\t"}{print $3,$4,$5,$1,$2,$9}' > ${sp}.indels.bed

  # get CDS
  cat in_sorted.gff \
   | awk '{OFS="\t"} $1 ~ /^#/ {next} {if ($3 == "CDS") {split($9, a, ";"); print $1, $4-1, $5, a[3]}}' \
   | sed -E 's/Dbxref=([A-Za-z=]+:[^,]+).+/\1/' > ${sp}.cds.bed
  echo -e "chr\tstart\tstop\tgene\ttype\tlen\tstrand" > ${sp}_indel_gene_list.tsv
  bedtools intersect -wa -wb -a ${sp}.cds.bed -b ${sp}.indels.bed \
   | awk '{OFS="\t"}{print $5,$6,$7,$4,$8,$9,$10}' \
   ${sp}_indel_gene_list.tsv

  # get all unique gene IDs
  cat ${sp}.cds.bed | awk '{print $NF}' | sort | uniq > ${sp}_gene_list.tsv

  #calculate proportion of genes without indel
  count=$(cat ${sp}_indel_gene_list.tsv | tail +2 | awk '{print $4}' | sort | uniq | wc -l)
  tot=$(cat ${sp}_gene_list.tsv | wc -l)
  echo "scale=2; ${count}/${tot}" | bc > ${sp}_gene_indelprop.txt
  
  for i in $(ls -d */); do 
    sp=$(echo ${i} | tr -d '/')
    cd ${sp}
    count=$(cat ${sp}_indel_gene_list.tsv | tail +2 | awk '{print $4}' | sort | uniq | wc -l)
    tot=$(cat ${sp}_gene_list.tsv | wc -l)
    echo "scale=2; ${count}/${tot}" | bc > ${sp}_gene_indelprop.txt
    cd ..
  done

  mkdir -p ${outDir}/${sp}
  cp ${sp}.*.qv ${outDir}/${sp}
  cp ${sp}_indel_gene_list.tsv ${outDir}/${sp}
  cp ${sp}_gene_list.tsv ${outDir}/${sp}
  cp ${sp}_gene_indelprop.txt ${outDir}/${sp}

  conda deactivate

  rm -r ./*
done

# collect all scores in csv
( echo "species,region,type,mean,q10,q50,q90"
for sp in $(ls -d */ | tr -d '/'); do
  for region in "exons" "introns" "intergenic" "repeats" "wholegenome"; do
    for type in "err_ont" "err_bal" "iden" "del" "ins"; do
      stats=$(cat ${sp}/${sp}.${region}.qv | grep -A6 "#  Q Scores" \
       | grep "${type}" | awk '{print $2,$3,$4,$5}' | tr ' ' ',')
      echo "${sp},${region},${type},${stats}"
    done
  done 
done ) > ref_qscore_summary.csv