#! /bin/bash

# Sample BLASTn command

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

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
seqtk subseq -l80 ${sp}.pilon.fasta ${sp}_nobusco_contigs.txt \
  | awk '/^>/{print ">contig_" ++i; next}{print}' \
  > ${outassm} \
 && rm ${sp}_nobusco_contigs.txt