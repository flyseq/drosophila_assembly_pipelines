setwd("~/Dropbox/writing/2020_drosophila_genomes_manuscript/figures/busco_synteny_figure")

#install.packages("devtools")
#devtools::install_github("analyxcompany/ForceAtlas2")

library(igraph)
library(tidyverse)
library(ggraph)
library(ForceAtlas2)

# ForceAtlas2 R package info can be found at the following links
#https://github.com/analyxcompany/ForceAtlas2
#https://rdrr.io/github/analyxcompany/ForceAtlas2/man/layout.forceatlas2.html

# read in edges
edges <- read.csv("edges.csv")
names(edges) <- c("from","to","type","weights")
edges <- edges[,!names(edges) == "type"]

# read in nodes
nodes <- read.csv("nodes_with_attrs.csv")
names(nodes) <- c("ID","label","chr")
nodes <- nodes[,!names(nodes) == "label"]

# perform layout
g <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
layout <- layout.forceatlas2(g, iterations=3000, plotstep=10,
                             tolerance=1, gravity=1)

# plot
q <- ggraph(g, layout=layout) + 
  geom_edge_density() +
  geom_edge_link(alpha=0.01) +
  geom_node_point(aes(color=chr), alpha=0.7) +
  theme_classic(base_size=14) +
  scale_color_manual(element_blank(),values=c("#F1BB7B", "#F06266", "#739F89", 
                                       "#D57237", "#E2D236", "#7294D5", "black")) +
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.ticks=element_blank(), legend.position="top")

ggsave("busco_synteny.pdf", q, width=6, height=6.5, units="in")