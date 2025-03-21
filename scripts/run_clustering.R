library(tidyverse)
library(ggrepel)
library(proxy)
library(uwot)
source("/g/schwab/Marco/repos/ookinete_clustering/scripts/functions.R")

clustering_data <- read_tsv("/g/schwab/Marco/repos/ookinete_clustering/data/clustering_data.tsv") 


full_data <- read_tsv("/g/schwab/Marco/repos/ookinete_clustering/data/full_data.tsv")


size_tot=4


outdir <- "/g/schwab/Marco/ook_plots/cfinal_210325"

final_figure_plot(clustering_data, 
                  cutoff=0.18,  
                  outdir=paste0(outdir),
                  size_pts=size_tot,
                  size_txt=size_tot)



final_staging <- read_tsv(file.path("/g/schwab/Marco/repos/ookinete_clustering/data/final_staging.tsv"))



cell_id_mapping <- read_csv(file.path(outdir, "cell_id_mapping.csv")) %>%
  bind_rows(
    tibble(cell_id=c("z1", "z5", "z16"),
           seriation_id=c(26,27,28))
  )

norm_between_one_zero <- full_data %>%
  select(-cell_id) %>%
  as.matrix() %>% apply(., 2, normalize_column) 


feature_sorter <- tibble(
  
  feature=c(
    
    'has_nuctail', 
    'vol_sc', 
    'vol_nuc', 
    'vol_cell', 
    'Apicoplast Area',   
    'area_nuc', 
     'Apicoplast Volume',
     'sphericity_nuc',
    'sphericity_cell',
    'area_cell', 
    'Apicoplast branching',   
    'vol_nucleo', 
    'has_micronuc',
    'retort_length',
    'n_cryst',
    'has_micronemes', 
    'is_full',
    "seriation_id"),
  
  
  neds_name= c(
    
    "nuclear tail",
    "synaptonemal complexes",
    "nuclear volume",
    "cellular volume",
    'apicoplast area',
    "nuclear area",
    'apicoplast volume',
    "nuclear sphericity",
    "cellular sphericity",
    "cellular area",
    'apicoplast branching',
    "nucleolar volume",
    "nuclear appendix",
    "retort length",
    "crystalloids",
    "micronemes",
    "nuclear relocalization",
    "seriation id"
    
  ),
  
  neds_name_table= c(
    
    "nuclear tail",
    "synaptonemal complexes volume",
    "nuclear volume",
    "cellular volume",
    'apicoplast area',
    "nuclear area",
    'apicoplast volume',
    "nuclear sphericity",
    "cellular sphericity",
    "cellular area",
    'apicoplast branching',
    "nucleolar volume",
    "nuclear appendix",
    "retort length",
    "crystalloids",
    "micronemes",
    "nuclear relocalization",
    "cell id"
    
  )
  
)




hm_raw_f_data <- norm_between_one_zero %>% 
  bind_cols(
    full_data %>% select(cell_id)
  ) %>% 
  left_join(cell_id_mapping) %>%
  left_join(final_staging) %>%
  mutate(
    #stage=
  ) %>%
  # mutate(cell_id=rownames(hm_raw_f_data)) %>%
  pivot_longer(cols = colnames(norm_between_one_zero), names_to = "feature", values_to = "val") %>%
  
  mutate(
    is_complete=ifelse(seriation_id >25, "T", "F")
  )
  
  # mutate(#cell_id=factor(str_replace_all(cell_id, "z", ""), levels=str_replace_all(cell_id_order$cell_id, "z", "")),
  #   cell_id=factor(seriation_id, levels=cell_id_order$seriation_id),
  #   feature=factor(feature, levels=rev(feature_vec)))

hm_raw_features <-
  hm_raw_f_data %>%
  ggplot(.)+
    facet_grid(~stage,  scales = "free", space="free")+
  geom_tile(
    aes(y=factor(feature, levels = rev(feature_sorter$feature)), 
        x=factor(seriation_id, levels = as.character(c(1:28))), fill=val)
  )+
  
  scale_fill_gradientn(
    name="norm. val.",
    limits=c(-max(hm_raw_f_data$val),max(hm_raw_f_data$val)), colors=c("white","firebrick"))+
  scale_y_discrete(breaks=feature_sorter$feature, 
                   labels=feature_sorter$neds_name)+
  theme_bw()+
  #    theme(legend.position = "bottom",
  #          plot.margin = unit(c(0,1,0,0), "cm"))+
  xlab("cell id")+
    ylab("feature")




pdf(file.path(outdir, "supp_full_feature_heatmap.pdf"), width=12.5, height=4.5)
print(hm_raw_features)
dev.off()



make_heatmap_figure(select(full_data %>% filter(!cell_id %in% c("z1", "z16", "z5")), all_of(feature_sorter$feature)))+
  scale_x_discrete(breaks=feature_sorter$feature, 
                   labels=feature_sorter$neds_name)+
  scale_y_discrete(breaks=feature_sorter$feature, 
                   labels=feature_sorter$neds_name)+
  #scale_fil
  theme(
    axis.text.x = element_text(angle=315, hjust=0)
  )



final_ov_table_pre <- full_data %>%
  select(-has_sc) %>%
  left_join(cell_id_mapping)


new_names <- lapply(names(final_ov_table_pre), function(ON){
  print(ON)
  new <- feature_sorter %>%
    filter(feature==ON) 
  
  print(new)
  
  if(nrow(new)==0){
    return(ON)
  } else {
    return(pull(new, neds_name))
  }
  
}) %>%
  unlist()


colnames(final_ov_table_pre) <- new_names


final_ov_table <- final_ov_table_pre %>%
  select(-cell_id) %>%
  .[c(ncol(.), 2:(ncol(.)-1))]


write_tsv(final_ov_table,
          file=file.path(outdir, "final_overview.tsv"))




