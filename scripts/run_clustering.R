library(tidyverse)
library(ggrepel)
library(proxy)
library(uwot)
source("/home/rheinnec/schwab/repos/ookinete_clustering/scripts/functions.R")

clustering_data <- read_tsv("/home/rheinnec/schwab/repos/ookinete_clustering/data/clustering_data.tsv")

size_tot=5

final_figure_plot(clustering_data, 
                  cutoff=0.18,  
                  outdir=paste0("/g/schwab/Marco/ook_plots/new_cfinal_size", size_tot, "_repel_kreiskreis"),
                  size_pts=size_tot,
                  size_txt=size_tot)




