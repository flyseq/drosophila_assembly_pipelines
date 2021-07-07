# Figure data
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

### Figure 3

### Figure 4

### Figure 5

### Figure 6
* `contiglens_coverage.csv`: Files containing contig lengths from which N50, NG50, etc. can be computed. Columns have the same meanings as `contiglens*.csv` files for *Figure 1*.