#! /bin/bash

# Sample BLASTn command

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

assm="${sp}.pilon.nobusco.fasta"  # Pilon assembly minus BUSCO contigs filename

# a sample BLAST command
# you may want to increase -max_target_seqs for more results
blastn -task megablast -query ${assm} -db nt \
    -outfmt '6 qseqid sskingdoms sblastnames sscinames scomnames stitle evalue bitscore' \
    -max_target_seqs 1 -max_hsps 1 -num_threads 70 -evalue 1e-25 \
    -out ${assm}.out