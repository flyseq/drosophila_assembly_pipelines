#! /bin/bash

# this script estimates genome size from depth of coverage over single-copy
# BUSCOs

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use
buscoRun="${sp}.buscov4"           # path to folder containing BUSCOv4 run

reads="${sp}.passReads.fastq.gz"   # nanopore reads (can use Illumina instead)
assm="${sp}.assembly.sm.fasta"     # filename of assembly BUSCO was run with
outFile="${sp}.buscoCoverage.csv"  # filename to write BUSCO output to

# get BED interval of complete and single-copy BUSCO coordinates
cat ${buscoRun}/run*/full_table.tsv | grep -v "#" \
 | awk '{OFS="\t"}($2=="Complete"){print $3,$4-1,$5,$1}' > ${sp}.buscos.bed

# get number of complete and single-copy BUSCOs
ncomp=$( cat ${sp}.buscos.bed | wc -l )

# map reads to genome
minimap2 -ax map-ont --secondary no -t ${threads} ${assm} ${reads} \
  | samtools sort -@${threads} -o ${sp}_reads_to_busco.sam

# compute coverage statistics
# total depth over all busco sites
sum=$(samtools depth -b ${sp}.buscos.bed ${sp}_reads_to_busco.sam \
        | awk '{sum+=$3} END { print sum","NR }')
# total number of bp in BUSCOs
len=$(cat ${sp}.buscos.bed \
        | awk '{sum+=($3-$2)} END { print sum }')
# total number of bp in reads
readl=$(zcat ${reads} \
          | paste - - - - \
          | cut -f 2 \
          | tr -d '\n' \
          | wc -c )

# write info to file
# species name, num complete BUSCOS, total depth, total BUSCO bp, total read bp
# genome size is computed as: (total read bp)/(total depth/total BUSCO bp)
echo "${sp},${ncomp},${sum},${len},${readl}" \
  > ${outFile}
