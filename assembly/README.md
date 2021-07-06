# Genome assembly workflow
The genome assembly pipeline is provided, step by step, in individual scripts. To ensure reproducibility and consistency of compute environments, each step was run either in a container or a Conda environment.

Note that older versions of many of these programs have been specified here for the sake of consistency across our assemblies. In many cases there have been significant updates to these programs. Please use the most up-to-date versions for best results.

## Conda environment setups

## 1. BUSCO
We used [BUSCO v3](https://gitlab.com/ezlab/busco/-/tree/3.0.2) to assess assembly completeness through each step of the assembly. [BUSCO v4](https://gitlab.com/ezlab/busco/-/tree/4.1.4) was released during the course of this work and was used to assess the completeness of the final assemblies.

### BUSCO v3:
```bash
conda create --name buscov3
git clone --branch "3.0.2" https://gitlab.com/ezlab/busco.git
cd busco
python setup.py install --user
```

### BUSCO v4:
```bash
conda create --name buscov4 -c bioconda -c conda-forge busco=4.1.4
```

## 2. Purge_haplotigs

```bash
conda create --name purge_haplotigs -c bioconda -c conda-forge \
    purge_haplotigs=1.1.1
```

## Nanopore-based assembly
Edit parameters at top of `make_nanopore_script.sh` and run to spawn a set of smaller job scripts. These should be run sequentially. Changes should be made to threads requested if varying between tasks (e.g. submitting to different nodes on a cluster). The draft sequence is generated following [ONT's recommendations](https://nanoporetech.github.io/medaka/draft_origin.html#how-should-i-create-my-draft-sequence). 
1. `01_run_guppy.sh`: Run Guppy basecaller in high-accuracy mode
1. `02_run_flye.sh`: Run Flye, generate initial draft assembly
1. `03_run_racon.sh`: Polish twice with Racon
1. `04_run_medaka.sh`: Polish once with Medaka

## Haplotig identification and removal (not run in a container)
Duplicate contigs in the assembly (representing alternative haplotypes, or haplotigs) were identified and removed with the [Purge Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) pipeline. 

Purge_haplotigs was installed in a Conda environment:
```bash
conda create --name purge_haplotigs -c bioconda -c conda-forge purge_haplotigs
```

There are manual steps to the workflow in `purge_haplotigs.sh`. Please see the purge_haplotigs repository for instructions.

## Scaffolding
If haplotig purging was performed, we attempted to re-scaffold the assembly using long reads. The Dockerfile includes an installation of [npScarf](https://github.com/mdcao/npScarf). 

## Pilon polishing
Three rounds of Pilon polishing were performed. Edit parameters at the top of `polish_pilon.sh` and run.

## Decontamination
The NCBI BLAST database was downloaded locally and is too big to provide in a container. A Docker Image with BLAST+ applications pre-installed can be obtained by: 
```bash
docker pull ncbi/blast
```

## Repeat masking
A repeat masking image is not provided because the RepBase RepeatMasker library cannot be freely provided. The instructions for running RepeatMasker are provided [at this link](http://www.repeatmasker.org/RMDownload.html).

## Variant calling
Short read variant calling tools are available in the Docker image. [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) already provides Docker and Singularity images.

## Quality assessment
[Merqury](https://github.com/marbl/merqury) and [Pomoxis](https://github.com/nanoporetech/pomoxis) were installed in Conda environments.

Merqury:
```bash
conda create --name merqury -c bioconda merqury
```

Pomoxis:
```bash
conda create --name pomoxis -c bioconda pomoxis
```

## Genome size estimation
Although Jellyfish is provided in the Docker, we have not included [GenomeScope](https://github.com/schatzlab/genomescope) as it can simply be run as an R script. A [web interface is also available](http://qb.cshl.edu/genomescope/).