# Additional data
This directory contains the raw data for figures presented in the manuscript. Additional explanations are provided for data columns with ambiguous meanings.

### Figure 1
* `busco_summary.csv`: BUSCOv4 scores for each assembly
* `contiglens*.csv`: Files containing contig lengths from which N50, NG50, etc. can be computed. 
    + `contig_rank`: Size rank of contigs in each genome. Largest contigs have lowest rank. Ties are resolved by random choice.
    + `contig_prop`: Proportion of the genome assembly contained in each contig.
    + `contig_cumsum`: Cumulative length of the genome assembly represented by each contig when sorted according to `contig_rank`.
    + `contig_cumsum_prop`: Cumulative proportion of the genome assembly represented by each contig when sorted according to `contig_rank`.
* `genome_size.csv`: Genome size estimates for each assembly
    + `nbusco`: Number of complete single-copy BUSCOs used for genome size estimation
    + `lr_sumcoverage`: Read depth per site summed across all sites in BUSCO genes
    + `lr_coveredbases`: Number of sites in BUSCO genes with read depth of at least 1
    + `lr_totalbases`: Total number of sites in BUSCO genes
    + `lr_nanoporelen`: Total bases in Nanopore reads
    + `lr_genomesize`: Genome size estimated from read depth in BUSCOs
    + `contam`: Proportion of the original assembly tagged as contamination and removed
    + `lr_genomesize_adj`: `lr_genomesize` \* `(1-contam)`
    + `sr_genomesize`: Genome size estimated with JellyFish + GenomeScope from short read k-mer histogram
    + `assm_size`: Length of the assembled sequences

### Figure 2
Estimates of heterozygosity in the samples used for assembly.
* `het_summary.csv`: Summary table of variant calls for long reads (`lr`) and short reads (`sr`).
    + `num_callable_*`: Number of callable (variant+invariant) sites passing quality filtering
    + `snp_het_*`: Number of SNPs called as heterozygote
    + `snp_hom_*`: Number of SNPs called as homozygous non-reference
    + `indel_het_*`: Number of indels called as heterozygote
    + `indel_hom_*`: Number of indels called as homozygous non-reference

* `het_vs_contiguity.csv`: Comparisons of measured diversity with genome contiguity.
    + `srhet`: SNP diversity estimated from short reads
    + `lrhet`: SNP diversity estimated from long reads
    + `50kb`: Proportion of the Nanopore bases contained in reads >50kb

### Figure 3
Consensus quality scores.

* `qv_reference_summary.csv`: Quality scores calculated by comparisons to NCBI RefSeq genomes. Genomes were divided into 100kb windows before quality assessment. Mean, 10th, 50th, and 90th quantile scores are reported.
    + `err_ont`: Overall QV, computed using the length of the alignment
    + `err_bal`: Overall QV, computed using the reference span
    + `iden`: QV(SNVs)
    + `del`: QV(deletions)
    + `ins`: QV(insertions)

* `qv_summary.csv`: Overall mean quality scores per genome. QV estimates that could not be computed because of the lack of short reads or a reference assembly are represented by `NA`.
    + `shortreadsource`: Did short reads come from the same sample or a different sample?
    + `merqury`: Merqury-estimated QV
    + `lr_qv`: QV estimated from long-read variant calls
    + `sr_qv`: QV estimated from short-read variant calls
    + `pomoxis`: QV estimated with a reference-based comparison
    + `refsource`: Is the reference genome the same strain as our assembly?

### Figure 4
Files for BUSCO synteny network analysis. Only single-copy and complete BUSCOs were used for these analyses.

* `busco_graph_nodes.csv`: List of nodes (BUSCOs)
* `busco_graph_edges.csv`: Edge weights for the graph
* `complete_busco_locations.csv`: The locations of each BUSCO annotation in each assembly.
    + `ID`: BUSCO ID
* `dmel_complete_buscos.csv`: Locations of BUSCO annotations in the D.melanogaster reference genome

### Figure 5
Phylogenetic tree files and the RepeatMasker output summary. Please note that the RepeatMasker table may be inaccurate for species without well-characterized repeats present in the Dfam 3.1 and RepBase RepeatMasker libraries.

* `busco_genes_ml.tree`: Gene trees inferred with RAxML
* `busco_genes_astral.tree`: Species tree inferred from gene trees with ASTRAL
* `rm_summary.csv`: Summary table of RepeatMasker output. 

### Figure 6
Contig length information for assemblies from lower coverage read sets downsampled from the full *D. jambulina* Nanopore reads.

* `contiglens_coverage.csv`: Files containing contig lengths from which N50, NG50, etc. can be computed. Columns have the same meanings as `contiglens*.csv` files for **Figure 1**.