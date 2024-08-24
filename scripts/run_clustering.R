library(tidyverse)
library(ggrepel)
library(proxy)
library(uwot)
source("/home/rheinnec/schwab/repos/ookinete_clustering/scripts/functions.R")

clustering_data <- read_tsv("/home/rheinnec/schwab/repos/ookinete_clustering/data/clustering_data.tsv")

final_figure_plot(clustering_data, 
                  cutoff=0.18,  
                  outdir="/g/schwab/Marco/ook_plots/cfinal")
