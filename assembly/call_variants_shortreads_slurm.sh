#! /bin/bash
#
#SBATCH --job-name=het_calc_master
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov
#SBATCH --ntasks=1
#SBATCH --array=1-91
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# load stuff
. /home/users/bkim331/.bashrc
ml system ncurses

# read in sp info
read -r reads mod <<< $( sed -n "${SLURM_ARRAY_TASK_ID}p" data_list.tsv )

#set up files
cwd=$(pwd)
sp=$(echo $reads | sed -E 's/\.passReads\.guppy...\.fastq\.gz//' | sed -E 's/.R10//')
ref="${sp}.assembly.sm.fasta"
mkdir -p ${cwd}/${sp} \
 && cd ${cwd}/${sp}

#Qscore threshold
q="30"

# download
rclone copy stanfordbox:100x100/assemblies/repeat_masked/${sp}.repeatMasker.tar.gz ./
rclone copy stanfordbox:100x100/assembly_data/illuminaForAssembly/${sp}/ ./

# check if split; cat if
if [ -f ${sp}_R1.fastq.gz.000 ]; then
    cat ${sp}_R1.fastq.gz.0?? > ${sp}_R1.fastq.gz
fi
if [ -f ${sp}_R2.fastq.gz.000 ]; then
    cat ${sp}_R2.fastq.gz.0?? > ${sp}_R2.fastq.gz
fi

#extract info
tar zxvfO ${sp}.repeatMasker.tar.gz \
    ${sp}.repeatMasker/${sp}.assembly.cleaned.fasta.out | \
    awk '{OFS="\t"}(NR>=4){print $5,$6-1,$7}' | bedtools sort | 
    bedtools merge -i - > repeats.bed
tar zxvfO ${sp}.repeatMasker.tar.gz \
    ${sp}.repeatMasker/${sp}.assembly.cleaned.fasta.masked > ${ref}

####################################
#                                  #
#   do short read variant calling  #
#                                  #
####################################

threads=16
cat << EOF > sr1_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=sr1_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=3G

. /home/users/bkim331/.bashrc
ml system
ml ncurses

cd ${cwd}/${sp}

# map short reads to genome
minimap2 -ax sr -t ${threads} ${ref} ${sp}_R1.fastq.gz ${sp}_R2.fastq.gz \
    > sr_reads_to_draft.sam
sambamba view -S -t ${threads} sr_reads_to_draft.sam -f bam | \
    sambamba sort -t ${threads} -o sr_reads_to_draft_sorted.bam /dev/stdin

# remove duplicates
sambamba markdup -r -p -t${threads} sr_reads_to_draft_sorted.bam \
    sr_reads_to_draft_dedup.bam
EOF

sr1ID=$(sbatch sr1_${sp}.sh | sed 's/Submitted batch job //')

threads="16"
cat << EOF > sr2_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=sr2_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=3G
#SBATCH --dependency=afterok:${sr1ID}

. /home/users/bkim331/.bashrc
ml system ncurses parallel

cd ${cwd}/${sp}

# make pileups for each contig (faster)
cat ${ref} | grep ">" | tr -d ">" | parallel -j ${threads} \
    bcftools mpileup -f ${ref} -r {} sr_reads_to_draft_dedup.bam \
    -q 10 -Q 20 -o sr_{}.mpileup
bcftools concat --threads ${threads} sr_*.mpileup > ${sp}.mpileup \
 && rm sr_contig_*.mpileup
EOF

sr2ID=$(sbatch sr2_${sp}.sh | sed 's/Submitted batch job //')

threads="16"
cat << EOF > sr3_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=sr3_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=3G
#SBATCH --dependency=afterok:${sr2ID}

. /home/users/bkim331/.bashrc
ml system
ml ncurses

cd ${cwd}/${sp}

# call all sites
bcftools call -m -Ov --threads ${threads} -o sr_calls.g.vcf ${sp}.mpileup
EOF

sr3ID=$(sbatch sr3_${sp}.sh | sed 's/Submitted batch job //')

threads="16"
cat << EOF > sr4_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=sr4_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=3G
#SBATCH --dependency=afterok:${sr3ID}

. /home/users/bkim331/.bashrc
ml system
ml ncurses

cd ${cwd}/${sp}

# filter gVCF for QV
bcftools sort sr_calls.g.vcf | bcftools view -i '%QUAL>=30 && DP>=10' > sr_filt.g.vcf

# remove gVCF sites in repetitive regions
cat sr_filt.g.vcf | grep -E "^#" > sr_filt_norep.g.vcf
bedtools subtract -a sr_filt.g.vcf -b repeats.bed >> sr_filt_norep.g.vcf

# filter for SNPs with QV
bcftools view -i '%QUAL>=30' --types snps -m 2 -M 2 --threads ${threads} \
    sr_filt.g.vcf | bcftools view -i 'MIN(FMT/GQ)>=30' > sr_filt_snps.vcf

# remove SNPs in repetitive regions
cat sr_filt_snps.vcf | grep -E "^#" > sr_norepeat_snps.vcf
bedtools subtract -a sr_filt_snps.vcf -b repeats.bed >> sr_norepeat_snps.vcf

# filter for indels with QV
bcftools view -i '%QUAL>=30' --types indels -m 2 -M 2 --threads ${threads} \
    sr_filt.g.vcf | bcftools view -i 'MIN(FMT/GQ)>=30' > sr_filt_indels.vcf

# remove indels in repetitive regions
cat sr_filt_indels.vcf | grep -E "^#" > sr_norepeat_indels.vcf
bedtools subtract -a sr_filt_indels.vcf -b repeats.bed >> sr_norepeat_indels.vcf
EOF

sr4ID=$(sbatch sr4_${sp}.sh | sed 's/Submitted batch job //')