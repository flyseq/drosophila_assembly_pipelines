# Genome assembly workflow
The genome assembly pipeline.

## Installing BUSCO
We used [BUSCO v3](https://gitlab.com/ezlab/busco/-/tree/3.0.2) to assess assembly completeness through each step. BUSCO v4 was released during the course of this work and should be used, but we maintained v3 through the pipeline for internal consistency.

BUSCO was installed in a Conda environment:
```bash
conda create --name buscov3
git clone --branch "3.0.2" https://gitlab.com/ezlab/busco.git
cd busco
python setup.py install --user
```

and otherwise run according to the developer's instructions.

## Nanopore-based assembly
Edit parameters at top of `make_nanopore_script.sh` and run to spawn a set of smaller job scripts. These should be run sequentially. Changes should be made to threads requested if varying between tasks (e.g. submitting to different nodes on a cluster). The draft sequence is generated following [ONT's recommendations](https://nanoporetech.github.io/medaka/draft_origin.html#how-should-i-create-my-draft-sequence). 
1. `01_run_guppy.sh`: Run Guppy basecaller in high-accuracy mode
1. `02_run_flye.sh`: Run Flye, generate initial draft assembly
1. `03_run_racon.sh`: Polish twice with Racon
1. `04_run_medaka.sh`: Polish once with Medaka

## Haplotig identification and removal (not run in a container)
Duplicate contigs in the assembly (representing alternative haplotypes, or haplotigs) were identified and removed with the [Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) pipeline. 

## Pilon polishing
Three rounds of Pilon polishing were performed. Edit parameters at the top of `polish_pilon.sh` and run.

## Decontamination
The NCBI BLAST database was downloaded locally and is too big to provide in a container. A Docker Image with BLAST+ applications pre-installed can be obtained by: 
```bash
docker pull ncbi/blast
```

## Repeat masking
A repeat masking image is not provided because the RepBase RepeatMasker library cannot be freely provided. The instructions for running RepeatMasker are provided [at this link](http://www.repeatmasker.org/RMDownload.html).