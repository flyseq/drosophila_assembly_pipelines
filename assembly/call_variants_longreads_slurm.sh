#! /bin/bash
#
#SBATCH --job-name=deepvariant
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --array=1-105
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# load stuff
. /home/users/bkim331/.bashrc
ml system ncurses

# read in sp info
read -r reads mod <<< $( sed -n "${SLURM_ARRAY_TASK_ID}p" data_list.tsv )

#set up files
sp=$(echo $reads | sed -E 's/\.passReads\.guppy...\.fastq\.gz//' | sed -E 's/.R10//')
ref="${sp}.assembly.sm.fasta"
mkdir -p /scratch/users/bkim331/deepvariant/${sp} \
 && cd /scratch/users/bkim331/deepvariant/${sp}
cwd=$(pwd)

#Qscore threshold
q="30"

# download
rclone copy stanfordbox:100x100/assemblies/repeat_masked/${sp}.repeatMasker.tar.gz ./
rclone copy petrovlab:100x100/data_assembly/nanoporeReads/${reads} ./

#extract info
tar zxvfO ${sp}.repeatMasker.tar.gz \
    ${sp}.repeatMasker/${sp}.assembly.cleaned.fasta.out | \
    awk '{OFS="\t"}(NR>=4){print $5,$6-1,$7}' | bedtools sort | 
    bedtools merge -i - > repeats.bed
tar zxvfO ${sp}.repeatMasker.tar.gz \
    ${sp}.repeatMasker/${sp}.assembly.cleaned.fasta.masked > ${ref}

####################################
#                                  #
#   do long read variant calling   #
#                                  #
####################################

threads="16"
cat << EOF > lr1_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=lr1_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=normal,hns,dpetrov,owners
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=3G

. /home/users/bkim331/.bashrc
ml system
ml ncurses

cd ${cwd}

# map long reads to genome
minimap2 -ax map-ont -t ${threads} ${ref} ${reads} > lr_reads_to_draft.sam
sambamba view -S -t ${threads} lr_reads_to_draft.sam -f bam | \
    sambamba sort -t ${threads} -o lr_reads_to_draft_sorted.bam /dev/stdin
EOF

lr1ID=$(sbatch lr1_${sp}.sh | sed 's/Submitted batch job //')

threads="16"
cat << EOF > lr2_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=lr2_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=owners,dpetrov,hns,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=2G
#SBATCH --dependency=afterok:${lr1ID}

. /home/users/bkim331/.bashrc
ml system ncurses

cd ${cwd}

simg="/home/groups/dpetrov/bernard/singularity_images/deepvariant.simg"

# call variants with DeepVariant
[ -d deepvariant ] && rm -r deepvariant
singularity exec \${simg} \
    run_pepper_margin_deepvariant call_variant -b lr_reads_to_draft_sorted.bam \
    -f ${sp}.assembly.sm.fasta -o deepvariant -t ${threads} --gvcf --ont
EOF

lr2ID=$(sbatch lr2_${sp}.sh | sed 's/Submitted batch job //')

threads="8"
cat << EOF > lr3_${sp}.sh
#! /bin/bash
#
#SBATCH --job-name=lr3_${sp}
#SBATCH --time=24:00:00
#SBATCH --partition=owners,dpetrov,hns,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${threads}
#SBATCH --mem-per-cpu=2G
#SBATCH --dependency=afterok:${lr2ID}

. /home/users/bkim331/.bashrc
ml system ncurses

cd ${cwd}

# filter for SNPs with QV
bcftools view -i '%QUAL>=${q}' --types snps -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
    | bcftools view -i 'MIN(FMT/GQ)>=30'> lr_filtered_snps.vcf

# remove SNPs in repetitive regions
cat lr_filtered_snps.vcf | grep -E "^#" > lr_norepeat_snps.vcf
bedtools subtract -a lr_filtered_snps.vcf -b repeats.bed >> lr_norepeat_snps.vcf

# filter for indels with QV
bcftools view -i '%QUAL>=${q}' --types indels -m 2 -M 2 --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.vcf.gz \
    | bcftools view -i 'MIN(FMT/GQ)>=30' > lr_filtered_indels.vcf

# remove indels in repetitive regions
cat lr_filtered_indels.vcf | grep -E "^#" > lr_norepeat_indels.vcf
bedtools subtract -a lr_filtered_indels.vcf -b repeats.bed >> lr_norepeat_indels.vcf

# get all sites passing QV
bcftools view -i '%QUAL>=${q}' --threads ${threads} \
    deepvariant/PEPPER_MARGIN_DEEPVARIANT_FINAL_OUTPUT.g.vcf.gz \
    > lr_all_filt.vcf

# remove SNPs in repetitive regions
cat lr_all_filt.vcf | grep -E "^#" > lr_all_filt_norep.vcf
bedtools subtract -a lr_all_filt.vcf -b repeats.bed >> lr_all_filt_norep.vcf
EOF

lr3ID=$(sbatch lr3_${sp}.sh | sed 's/Submitted batch job //')