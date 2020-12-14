# Three rounds of Pilon polishing

#! /bin/bash

#edit these variables to match your assembly parameters
input="D.melanogaster.assembly.medaka.fasta"
output="D.melanogaster.assembly.pilon.fasta"

threads="78"
memory="200G"

read1="${species}_R1.fastq"
read2="${species}_R2.fastq"

cp ${input} temp_draft.fa
for i in {1..3}; do
	# map short reads to draft
    minimap2 -ax sr -t ${threads} temp_draft.fa ${read1} ${read2} | samtools sort -o reads_to_draft.bam --threads ${threads}
    samtools index reads_to_draft.bam
    # run pilon
    java -Xmx${memory} -jar /tools/pilon.jar --genome temp_draft.fasta --frags reads_to_draft.bam --outdir pilon --threads ${threads}
    # clean up
    cp ./pilon/pilon.fasta temp_draft.fasta
    rm  -r ./pilon
done

#final cleanup of temp files
mv temp_draft.fasta ${output}
rm reads_to_draft.bam
rm reads_to_draft.bam.bai