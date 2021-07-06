#! /bin/bash
# polish Flye assembly with two rounds of Racon

# job parameters
sp="D.melanogaster"      # sample name/ID
threads="32"             # number of threads to use

# copy Flye assembly to temp file
cp ${sp}.FlyeAssembly.fasta temp_draft.fa

# polish twice with Racon
for i in {1..2}; do
    minimap2 -x map-ont -t 32 temp_draft.fa ${sp}.passReads.fastq.gz \
      > temp_reads_to_draft.paf
    racon -t ${threads} -c 4 -m 8 -x -6 -g -8 -w 500 ${sp}.passReads.fastq.gz \
    	temp_reads_to_draft.paf temp_draft.fa \
      > temp_draft_new.fa
    mv temp_draft_new.fa temp_draft.fa
done

#clean up files
mv temp_draft.fa D.melanogaster.assembly.racon.fasta
rm temp_reads_to_draft.paf