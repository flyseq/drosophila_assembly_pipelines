# BUSCO-based analyses
Workflows for BUSCO-based analyses.

## Installing BUSCO
Although the BUSCO developers provide Docker images, we ran all BUSCO analyses
in a Conda evironment. To set up the environment:

```bash
conda create -n busco -c bioconda -c conda-forge busco=4.0.6 python=3
```

## Running BUSCO on individual genomes (not done in a container)
BUSCO analysis is run for each assembly and output saved in a tarball:
```bash
assm_name="D.melanogaster.assembly.sm.fasta"
sp=${assm_name%%.assembly*}

busco -i ${assm_name} -o ${sp}.busco -c $(nproc) -l diptera_odb10 \
    --augustus_species fly \
 && tar zcvf -o ${sp}.busco.tar.gz ${sp}.busco
```

All BUSCO output tarballs are saved to a single directory. Subsequent analysis
needs the path this directory to work.

## Workflow: building a protein tree from BUSCOs
`busco_phylo.sh`
1. Extract protein sequences of single-copy BUSCOs
1. Align with MAFFT
1. Infer gene trees with RAxML-NG
1. Infer species tree with ASTRAL-MP
1. Extract protein sequences of single-copy BUSCOs
1. Align with MAFFT
1. Infer gene trees with RAxML-NG
1. Infer species tree with ASTRAL-MP

`data/busco_gene_trees_ml.tree`: Newick-formatted individual gene trees (RAxML)

`data/busco_species_astral.tree`: Newick-formatted species tree (ASTRAL)

## Workflow: BUSCO synteny analysis (not done in a container)
1. `compile_busco_markers.sh`: Extract position information for each complete BUSCO
1. `busco_synteny_figure.py`: Take BUSCO position information `complete_busco_locations.csv` and parse
connections into node/edge tables; parse D. melanogaster BUSCO information using
provided `dmel_complete_buscos.csv`
1. `busco_synteny_figure.R`: Apply ForceAtlas2 layout algorithm and plot

`data/busco_graph_edges.csv` and `data/busco_graph_nodes.csv`: node/edge tables
