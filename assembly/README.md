# Genome assembly and quality assessment workflows
The genome assembly pipeline is provided, step by step, in individual scripts. To ensure reproducibility and consistency of compute environments, each step was run either in a container or a Conda environment.

Note that older versions of many of these programs have been specified here for the sake of consistency across our assemblies. In many cases there have been significant updates to these programs. Please use the most up-to-date versions for best results.

## Setup
Setup instructions for conda enviroments and containers.

### 1. BUSCO

[BUSCO v3](https://gitlab.com/ezlab/busco/-/tree/3.0.2) was used to assess assembly completeness through each step of the assembly. 

```bash
conda create --name buscov3 -c bioconda -c conda-forge busco=3.0.2
```

[BUSCO v4](https://gitlab.com/ezlab/busco/-/tree/4.1.4) was released during the course of this work and was used to assess the completeness of the final assemblies.

```bash
conda create --name buscov4 -c bioconda -c conda-forge busco=4.0.6
```

### 2. Purge_haplotigs

```bash
conda create --name purge_haplotigs -c bioconda -c conda-forge \
    purge_haplotigs=1.1.1
```

### 3. RepeatMasker
Although we set up RepeatMasker locally, it is now provided as part of the [Dfam-TETools container](https://github.com/Dfam-consortium/TETools). A Singularity image can be built with the following command:

```bash
singularity build tetools.simg docker://dfam/tetools:latest
```

### 4. NCBI BLAST+

```bash
conda create --name blast -c bioconda blast=2.10.1
```

To download a local copy of the Nucleotide (NT) database,

```bash
perl update_blastdb.pl --decompress nt
```

Then make sure the `BLASTDB` environmental variable is set to the directory containing the NT database files.
```bash
export BLASTDB="/path/to/blastdb"
```

### 5. Merqury

```bash
conda create --name merqury -c bioconda merqury
```

### 6. Pomoxis
Due to compatibility issues between the [Pomoxis](https://github.com/nanoporetech/pomoxis) and [AGAT](https://github.com/NBISweden/AGAT) Conda installations, we recommend installing them into separate Conda environments.

Pomoxis:
```bash
conda create --name pomoxis -c bioconda pomoxis
```

AGAT (this takes a while):
```bash
conda create --name agat -c bioconda agat
```

### 7. PEPPER-Margin-DeepVariant
[PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) already provides a Docker image. A local copy of the Singularity image is generated with the command

```bash
singularity pull docker://kishwars/pepper_deepvariant:r0.4
```

## Genome assembly workflow

### Assembly with Nanopore reads
Edit parameters at top of each script for each assembly, then run scripts sequentially. Changes should be made to threads requested if varying between tasks (e.g. submitting to different nodes on a cluster). The long-read based draft sequence is generated following [ONT's recommendations](https://nanoporetech.github.io/medaka/draft_origin.html#how-should-i-create-my-draft-sequence).
* `assembly_run_guppy.sh`: Run [Guppy basecaller](https://community.nanoporetech.com/downloads) (requires login) in high-accuracy mode
* `assembly_run_flye.sh`: Run [Flye](https://github.com/fenderglass/Flye), generate initial draft assembly
* `assembly_run_racon.sh`: Polish twice with [Racon](https://github.com/isovic/racon)
* `assembly_run_medaka.sh`: Polish once with [Medaka](https://github.com/nanoporetech/medaka)

### Haplotig identification and removal (run in Conda)
If the BUSCO duplication rate exceeds 1% at this step, duplicate contigs in the assembly (representing alternative haplotypes, or haplotigs) are identified and removed with [purge_haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/). 

Haplotig purging is performed on the Flye assembly, not the Medaka-polished assembly. The BUSCO assessment is done after Medaka polishing rather than immediately on the Flye assembly since the complete BUSCO %s can be underestimated from the Flye sequences, as they are less accurate.

* `purge_haplotigs.sh`: Haplotig purging, following the workflow provided by the developer. 

Coverage cutoffs must be set manually during the analysis. Please see the purge_haplotigs repository [for specific instructions](https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial).

### Scaffolding
If haplotig purging was performed, we attempted to re-scaffold the assembly using long reads. The Dockerfile includes an installation of [npScarf](https://github.com/mdcao/npScarf).

* `scaffold.sh`: Scaffolding, bridge contigs only if supported by at least 4 long reads

### Pilon polishing
One round of [Pilon](https://github.com/broadinstitute/pilon) polishing was performed. 

* `polish_pilon.sh`

### Decontamination
[NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) tools were used to align contigs against the nucleotide database. Any contigs with extensive matches to fungal, bacterial, or viral sequences were tagged for removal. A local copy of the NT database is recommended. We do not provide this as an image due to the large size of the NT database.

* `decont_exclude_busco_contigs.sh`: Identify contigs with single-copy genes present, exclude them from the BLASTn search
* `decont_blastn.sh`: Run `blastn` to align non-BUSCO contigs against NCBI NT database
* With manual curation of BLAST results, identify contigs that look like microbial contamination. A list of contigs to remove, with one contig name per line, should be provided as `<sample ID>_remove_contigs.txt`.
* `decont_remove_contamin.sh`: Remove contaminant sequences from the assembly.

### Repeat masking
The RepBase RepeatMasker library itself cannot be provided thus we cannot provide a self-contained repeatMasker image. Please use the [Dfam-TETools Docker Container](https://github.com/Dfam-consortium/TETools#using-the-container) and follow their instructions for [running RepeatMasker with the RepBase library](https://github.com/Dfam-consortium/TETools#using-repbase-repeatmasker-edition).

Sample command:
```bash
RepeatMasker --species 7214 -xsmall -pa ${threads} ${assemblyFile}
```

## Genome QC workflows
Scripts used for various quality assessments.

### Genome size estimation with long reads and BUSCOs
BUSCO should be run on the genome of interest before running the provided script.

* `genomesize_busco.sh`

### Genome size estimation from short reads
Although Jellyfish is included in the Docker setup, we have not included [GenomeScope](https://github.com/schatzlab/genomescope) as it is simply run as an R script. Only short read data are needed for the following scripts.

* `genomesize_jellyfish.sh`: generate a k-mer count histogram from short reads, use `<sample>.reads.hist` for genomeScope analysis
* [Instructions for running GenomeScope](https://github.com/schatzlab/genomescope#running-genomescope-on-the-command-line)

A [web interface is also available for GenomeScope](http://qb.cshl.edu/genomescope/).

### Variant with short reads
Short read variant calling tools are available in the Docker image. 

* `call_variants_shortreads_map.sh`: Map Illumina reads to a reference
* `call_variants_shortreads_pileup.sh`: Generate read pileup with [BCFtools](https://github.com/samtools/bcftools)
* `call_variants_shortreads_call.sh`: Call variants with BCFtools
* `call_variants_longreads_filter.sh`: Filter with BCFtools, site/genotype quality score of at least 30, exclude repetitive regions
  + `<sample>_shortread_snps.vcf`: VCF file containing SNPs 
  + `<sample>_shortread_indels.vcf`: VCF file containing indels
  + `<sample>_shortread_allsites.vcf`: VCF file containing all sites passing quality filter

### Variant calling with long reads
[PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) already provides Docker and Singularity images.

* `call_variants_longreads_map.sh`: Map Nanopore reads to reference 
* `call_variants_longreads_call.sh`: Call variants with PEPPER-Margin-DeepVariant
* `call_variants_longreads_filter.sh`: Filter variants with [BCFtools](https://github.com/samtools/bcftools). Filter on site/genotype quality score of at least 30, exclude repetitive regions
  + `<sample>_longread_snps.vcf`: VCF file containing SNPs 
  + `<sample>_longread_indels.vcf`: VCF file containing indels
  + `<sample>_longread_allsites.g.vcf`: gVCF file containing sites/intervals of sites passing quality filter (*note*: this is a gVCF)

### Quality assessment with Merqury
[Merqury](https://github.com/marbl/merqury) was used to perform k-mer based quality assessment. A genome assembly and short reads are required to run the following scripts.

* `qscore_trimreads.sh`: adapter trimming short reads in preparation for Merqury
* `qscore_merqury.sh`: script to run Merqury quality assessment

### Quality assessment with Pomoxis
A high-quality reference genome is required for reference-based quality assessments with Pomoxis. The following code is set up to work with NBCI RefSeq Drosophila genomes.

* `qscore_longestisoform.sh`: prepare RefSeq GFF for use with Pomoxis
* `qscore_reference.sh`: compute quality scores for different genomic regions. Creates the following useful output files
  + `${sp}_indel_gene_list.tsv`: list of all coding sequence indels
  + `${sp}_gene_list.tsv`: list of all protein-coding genes
  + `${sp}.referenceQV.csv`: table of quality scores