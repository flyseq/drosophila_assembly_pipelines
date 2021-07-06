#! /bin/bash
# basecalling Nanopore fast5s

# job parameters
sp="D.melanogaster"      # sample name/ID
nanoRaw="D.melanogaster" # path of parent folder containing fast5s
threads="32"             # number of threads to use

# run guppy with HAC model
guppy_basecaller -i ${nanoRaw} -s ${sp}.basecalled --recursive \
  -c dna_r9.4.1_450bps_hac.cfg --device "cuda:0" \
  --trim_strategy dna --qscore_filtering --calib_detect

#gather reads passing default quality filter and compress
cat ./${sp}.basecalled/pass/*.fastq \
  | pigz -p ${threads} \
  > ${sp}.passReads.fastq.gz