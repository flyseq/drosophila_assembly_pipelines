#! /bin/bash

# This script removes all tagged contaminant contigs from the Pilon-polished
# assembly.

# job parameters
sp="D.melanogaster"                # sample name/ID

assm="${sp}.assembly.pilon.fasta"     # Pilon assembly filename
remove="${sp}_remove_contigs.txt"     # File with contigs to exclude
outassm="${sp}.assembly.clean.fasta"  # Filename to write cleaned assembly to

# get complement of contigs tagged for removal
cat $assm \
  | grep ">" \
  | tr -d ">" \
  | grep -Fxvf D.melanogaster_remove_contigs.txt - \
  > ${sp}_keep_contigs.txt

# get new fasta minus contaminant contigs
# rename contigs to simplify
# then clean up files
seqtk subseq -l80 ${assm} ${sp}_keep_contigs.txt \
  | awk '/^>/{print ">contig_" ++i; next}{print}' \
  > ${outassm} \
 && rm D.melanogaster_remove_contigs.txt \
 && rm ${sp}_keep_contigs.txt