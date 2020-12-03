sp="L.clarofinis"
genomeSize="500m"
threads="78"
flyeThreads="78"
#guppy version
gv="guppy351"

cwd=$(pwd)

cat > ./${sp}/run_nanopore_workflow.sh << EOM
#singularity shell --nv -B ${cwd}/${sp}:/scratch/ -B /media/bernardkim/active-data/:/active-data/ ~/singularity_images/guppy.3.2.4_raconGPU_medakaGPU.simg

cd /scratch/
guppy_basecaller -i ${sp} \\
    -s ${sp}.basecalled --recursive \\
    -c dna_r9.4.1_450bps_hac.cfg \\
    --device "cuda:0" \\
    --trim_strategy dna \\
    --qscore_filtering --calib_detect

# for HAC base caller use
#    --config dna_r9.4.1_450bps_hac.cfg \\
# for modified base caller modify the config and enable fast5 output
#    --config dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg \\
#    --fast5_out \\
# for R10 use
#    --config dna_r10_450bps_hac.cfg
# removed for Guppy 3.2.4
#    --gpu_runners_per_device 4

#compile fastq into a single gzipped file
cat ./${sp}.basecalled/pass/*.fastq | pigz -p 70 > ${sp}.passReads.${gv}.fastq.gz
cat ./${sp}.basecalled/fail/*.fastq ./${sp}.basecalled/pass/*.fastq | pigz -p 70 > \\
    ${sp}.allReads.${gv}.fastq.gz

#compute an md5 hash for each fastq.gz
md5sum ${sp}.passReads.${gv}.fastq.gz | awk {'print \$1'} > ${sp}.passReads.${gv}.fastq.gz.md5
md5sum ${sp}.allReads.${gv}.fastq.gz | awk {'print \$1'} > ${sp}.allReads.${gv}.fastq.gz.md5

#compute read lengths for pass reads
zcat ${sp}.passReads.${gv}.fastq.gz | awk '{if(NR%4==2) print length(\$1)}' > ${sp}.readLengths.txt

#make transfer directory
mkdir -p transfer_box_reads

#split files if bigger than 15GB and move to transfer dir
for readType in "passReads" "allReads"; do
    if [ \$readType == "passReads" ]; then
        cp ${sp}.\${readType}.${gv}.fastq.gz /active-data/reads/
    fi
    fs=\$(stat --printf="%s" ${sp}.\${readType}.${gv}.fastq.gz)
    if [ \$fs -gt 15000000000 ]; then
        split -d -a3 -b 10000m ${sp}.\${readType}.${gv}.fastq.gz ${sp}.\${readType}.${gv}.fastq.gz.
	mv ${sp}.\${readType}.${gv}.fastq.gz.* ./transfer_box_reads
	rm ${sp}.\${readType}.${gv}.fastq.gz
    else
	mv ${sp}.\${readType}.${gv}.fastq.gz* ./transfer_box_reads
    fi
done

#get ready for assembly and polishing
cp /active-data/reads/${sp}.passReads.${gv}.fastq.gz ./
gunzip ${sp}.passReads.${gv}.fastq.gz

#assemble with Flye
flye --nano-raw ${sp}.passReads.${gv}.fastq \\
     --genome-size ${genomeSize} \\
     --threads ${flyeThreads} \\
     --out-dir ${sp}.FlyeAssembly

#save Flye output
mkdir -p transfer_box_flye
cp ./${sp}.FlyeAssembly/assembly.fasta ./transfer_box_flye/${sp}.FlyeAssembly.fasta
tar cf - ${sp}.FlyeAssembly | pigz -p ${threads} > ./transfer_box_flye/${sp}.FlyeAssembly.tar.gz
cp ./transfer_box_flye/${sp}.FlyeAssembly.fasta /active-data/FlyeAssemblies/${sp}.FlyeAssembly.fasta

#polish with racon
cp /active-data/FlyeAssemblies/${sp}.FlyeAssembly.fasta ./temp_draft.fa
for i in {1..2}; do
    minimap2 -x map-ont -t ${threads} temp_draft.fa ${sp}.passReads.${gv}.fastq > \\
        temp_reads_to_draft.paf
    racon -t ${threads} -c 4 -m 8 -x -6 -g -8 -w 500 ${sp}.passReads.${gv}.fastq \\
    	temp_reads_to_draft.paf temp_draft.fa > temp_draft_new.fa
    mv temp_draft_new.fa temp_draft.fa
done

#clean up and copy to active data drive
rm temp_reads_to_draft.paf
mv temp_draft.fa ${sp}.assembly.racon.fasta
cp ${sp}.assembly.racon.fasta /active-data/racon/

#polish with medaka
. /medaka/bin/activate
export TF_FORCE_GPU_ALLOW_GROWTH=true
medaka_consensus -i ${sp}.passReads.${gv}.fastq -d ${sp}.assembly.racon.fasta \\
    -t ${threads} -m r941_min_high -b 40 \\
    -o ${sp}.medaka

# for R9 revD
# -m r941_min_high
# for R10 
# -m r10_min_high

#save drafts
cp ./${sp}.medaka/consensus.fasta /active-data/medaka/${sp}.assembly.medaka.fasta
cp ./${sp}.medaka/consensus.fasta ./${sp}.assembly.medaka.fasta
EOM

echo "singularity shell --nv -B ${cwd}/${sp}:/scratch/ -B /media/bernardkim/active-data/:/active-data/ ~/singularity_images/guppy.3.2.4_raconGPU_medakaGPU.simg"