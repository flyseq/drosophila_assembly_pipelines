#! /bin/bash

# Quality assessment with Nanopore's Pomoxis tool
# this assumes you have Pomoxis installed in a Conda environment

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

assm="${sp}.assembly.sm.fasta"     # filename of assembly to qc
refgff="${sp}.longestIsoform.gff"  # filename of reference GFF, longest isoforms
refassm="${sp}.reference.fasta"    # reference assembly
refrep="${sp}.repeats.bed"         # BED file of repeats (from refSeq *.out.gz)

conda activate medaka

# get chrom sizes
faidx ${refassm} -i chromsizes | sort -k1,1 -k2,2n > ${sp}.chromSizes
awk '{print $1"\t"0"\t"$2}' ${sp}.chromSizes | sort -k1,1 -k2,2n \
 > ${sp}.chromSizes.bed

# get BED intervals of intergenic regions
cat ${sp}.longestIsoform.gff \
 | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' \
 | awk '($3 != "region"){print $0}' > in_sorted.gff
bedtools complement -i in_sorted.gff -g ${sp}.chromSizes \
 | awk '{print $0"\tNA"}'> ${sp}.intergenic.bed

# get BED intervals of exons
cat in_sorted.gff \
 | awk '{OFS="\t"} $1 ~ /^#/ {next} {if ($3 == "exon") {split($9, a, ";"); print $1, $4-1, $5, a[3]}}' \
 | sed -E 's/Dbxref=([A-Za-z=]+:[^,]+).+/\1/' > ${sp}.exons.bed

# get BED intervals of introns
bedtools complement \
  -i <(cat ${sp}.exons.bed ${sp}.intergenic.bed | sort -k1,1 -k2,2n) \
  -g ${sp}.chromSizes > ${sp}.introns.bed

# get BED intervals of CDS
cat in_sorted.gff \
 | awk '{OFS="\t"} $1 ~ /^#/ {next} {if ($3 == "CDS") {split($9, a, ";"); print $1, $4-1, $5, a[3]}}' \
 | sed -E 's/Dbxref=([A-Za-z=]+:[^,]+).+/\1/' > ${sp}.exons.bed

# run Pomoxis for introns, intergenic regions, and repeats
rm ${sp}.reference.fasta.??i
for region in "introns" "intergenic" "repeats"; do
  assess_assembly -r ${refassm} -i ${assm} -t ${threads} -c 100000 \
    -l 0 -b ${sp}.${region}.bed > ${sp}.${region}.qv
done

# run Pomoxis on whole genome
assess_assembly -r ${refassm} -i ${assm} -t ${threads} -c 100000 \
    -l 0 > ${sp}.wholegenome.qv

# run Pomoxis on exons
assess_assembly -r ${refassm} -i ${assm} -t ${threads} -c 100000 \
    -l 1 -b ${sp}.exons.bed > ${sp}.exons.qv

# get indel BED intervals
cat assm_indel_ge1.txt | grep -vE "^type" \
 | awk '{OFS="\t"}{print $3,$4,$5,$1,$2,$9}' > ${sp}.indels.bed

# get overlap of CDS intervals with indels
echo -e "chr\tstart\tstop\tgene\ttype\tlen\tstrand" > ${sp}_indel_gene_list.tsv
bedtools intersect -wa -wb -a ${sp}.cds.bed -b ${sp}.indels.bed \
 | awk '{OFS="\t"}{print $5,$6,$7,$4,$8,$9,$10}' \
 > ${sp}_indel_gene_list.tsv

# get all unique gene IDs
cat ${sp}.cds.bed | awk '{print $NF}' | sort | uniq > ${sp}_gene_list.tsv

#calculate proportion of genes with indel
count=$(cat ${sp}_indel_gene_list.tsv | tail +2 | awk '{print $4}' | sort | uniq | wc -l)
tot=$(cat ${sp}_gene_list.tsv | wc -l)
echo "scale=4; ${count}/${tot}" | bc > ${sp}_gene_indelprop.txt

# collect all scores in csv
( echo "species,region,type,mean,q10,q50,q90"
  for region in "exons" "introns" "intergenic" "repeats" "wholegenome"; do
    for type in "err_ont" "err_bal" "iden" "del" "ins"; do
      stats=$(cat ${sp}.${region}.qv | grep -A6 "#  Q Scores" \
       | grep "${type}" | awk '{print $2,$3,$4,$5}' | tr ' ' ',')
      echo "${sp},${region},${type},${stats}"
    done
  done 
) > ${sp}.referenceQV.csv

conda deactivate