#! /bin/bash
# call variants with mapped Nanopore reads

# job parameters
sp="D.melanogaster"                        # species/sample name
threads="16"                               # number of threads to use

# inputs
assm="${sp}.assembly.sm.fasta"             # assembly filename
bam="${sp}.readsToDraft.bam"               # BAM file name

# path to the PEPPER-Margin-DeepVariant Singularity image
simg="pepper_deepvariant_r0.4.sif"         

# call variants
singularity exec ${simg} \
    run_pepper_margin_deepvariant call_variant -b ${bam} \
    -f ${assm} -o deepvariant -t ${threads} --gvcf --ont