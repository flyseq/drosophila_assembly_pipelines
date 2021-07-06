#! /bin/bash

# this script excludes any contigs associated with a BUSCO from being considered
# for removal during the decontamination step

# job parameters
sp="D.melanogaster"                # sample name/ID
threads="32"                       # number of threads to use

assm="${sp}.assembly.pilon.fasta"  # Pilon-polished assembly filename

# get BUSCO annotations
# this assumes BUSCOv3 is installed in a Conda environment
# IMPORTANT: path to lineage database is specific to your installation:
#        -l /path/to/busco/database
conda activate buscov3
run_BUSCO.py -i ${assm} -o ${sp}.pilon.busco \
    -l diptera_odb9 --species fly \
    -m geno -c ${threads}
conda deactivate

# get list of contig names where BUSCOs are present
cat ${sp}.pilon.busco/full_table_${sp}*.tsv \
  | grep -v "^#" \
  | grep -v "Missing" \
  | awk '{print $3}' \
  | sort \
  | uniq \
  > ${sp}_keep_contigs.txt

# get contigs without BUSCOs (complement)
cat $assm \
  | grep ">" \
  | tr -d ">" \
  | grep -Fxvf ${sp}_keep_contigs.txt - \
  > ${sp}_nobusco_contigs.txt

# get fasta excluding contigs with BUSCOs
seqtk subseq -l80 ${sp}.pilon.fasta ${sp}_nobusco_contigs.txt \
  > ${sp}.pilon.nobusco.fasta \
 && rm ${sp}_nobusco_contigs.txt