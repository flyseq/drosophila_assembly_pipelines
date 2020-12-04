# Genome assembly workflow
The genome assembly pipeline.

## Installing BUSCO
placeholder
-we used buscov3
-small changes will need to be made to accommodate buscov4 but should easily be doable

## Nanopore-based assembly
Edit parameters at top of `make_nanopore_script.sh` and run to spawn a set of
smaller job scripts. These should be run sequentially. Changes should be made to
threads requested if varying between tasks (e.g. submitting to different nodes
on a cluster). The draft sequence is generated following [ONT's
recommendations](https://nanoporetech.github.io/medaka/draft_origin.html#how-should-i-create-my-draft-sequence). 
1. `01_run_guppy.sh`: Run Guppy basecaller in high-accuracy mode
1. `02_run_flye.sh`: Run Flye, generate initial draft assembly
1. `03_run_racon.sh`: Polish twice with Racon
1. `04_run_medaka.sh`: Polish once with Medaka

## Haplotig identification and removal (not run in a container)
Duplicate contigs in the assembly (representing alternative haplotypes, or
haplotigs) were identified and removed with the [Purge
Haplotigs](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/)
pipeline. 



placeholder

## Pilon polishing
placeholder

## Decontamination
placeholder

## 