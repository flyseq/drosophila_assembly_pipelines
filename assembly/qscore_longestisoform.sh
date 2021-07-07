#! /bin/bash

# Get longest isoforms for use with Pomoxis

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

# reference annotations 
refgff="GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff" 

conda activate agat

# merge overlapping transcripts of the same type
agat_convert_sp_gxf2gxf.pl \
  --gff ${refgff} --merge_loci -o ${sp}.ref.merged.gff

# keep only longest isoforms
agat_sp_keep_longest_isoform.pl -gff ${sp}.ref.merged.gff \
  -o ${sp}.longestIsoform.gff

conda deactivate