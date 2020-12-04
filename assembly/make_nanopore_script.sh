# This script creates a set of scripts to perform the Nanopore-based portion of
# the genome assembly pipeline.
#
# To use it, place it into your working directory alongside a folder containing
# the Nanopore fast5 files. The folder should be named with the species/strain
# name, that name should match the sp="D.melanogaster" name below. 
#
# The final output from this pipeline is the Medaka-polished assembly e.g.,
# D.melanogaster.assembly.medaka.fasta

# set these parameters for your assembly
sp="D.melanogaster"
genomeSize="140m"

#set the number of CPU threads to use (should probably set this manually)
threads=$(nproc)
# if you have a lot of ultra-long reads, flye sometimes runs out of memory. Set
# flyeThreads < threads if this is an issue.
flyeThreads="78"

#
# Make basecaller script
# 
cat > ./01_run_guppy.sh <<EOF
#run guppy in HAC mode
guppy_basecaller -i ${sp} \\
  -s ${sp}.basecalled --recursive \\
  -c dna_r9.4.1_450bps_hac.cfg \\
  --device "cuda:0" \\
  --trim_strategy dna \\
  --qscore_filtering --calib_detect

#gather reads passing default quality filter and gzip
cat ./${sp}.basecalled/pass/*.fastq | pigz -p ${threads} > ${sp}.passReads.fastq.gz
EOF

#
# Make assembly script
#
cat > ./02_run_flye.sh <<EOF
#run flye
flye --nano-raw ${sp}.passReads.fastq.gz \\
     --genome-size ${genomeSize} \\
     --threads ${flyeThreads} \\
     --out-dir ${sp}.FlyeAssembly
EOF

#
# Make racon polishing script
#
cat > ./03_run_racon.sh <<EOF
#copy Flye assembly to temp file
cp ${sp}.FlyeAssembly/assembly.fasta ./temp_draft.fa

#polish with Racon twice
for i in {1..2}; do
    minimap2 -x map-ont -t ${threads} temp_draft.fa ${sp}.passReads.fastq.gz > \\
        temp_reads_to_draft.paf
    racon -t ${threads} -c 4 -m 8 -x -6 -g -8 -w 500 ${sp}.passReads.fastq.gz \\
    	temp_reads_to_draft.paf temp_draft.fa > temp_draft_new.fa
    mv temp_draft_new.fa temp_draft.fa

#clean up files
mv temp_draft.fa ${sp}.assembly.racon.fasta
rm temp_reads_to_draft.paf
EOF

#
# Make Medaka polishing script
#
cat > ./04_run_medaka.sh <<EOF
#activate medaka virtual env
. /medaka/bin/activate

#special setting for RTX 2000 series
export TF_FORCE_GPU_ALLOW_GROWTH=true

#run medaka
medaka_consensus -i ${sp}.passReads.fastq.gz -d ${sp}.assembly.racon.fasta \\
    -t ${threads} -m r941_min_high -b 40 \\
    -o ${sp}.medaka

#copy draft to new file
cp ${sp}.medaka/consensus.fasta ./${sp}.assembly.medaka.fasta
EOF
