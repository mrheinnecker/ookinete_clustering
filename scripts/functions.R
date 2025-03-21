
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}



make_heatmap_figure <- function(num_data){
  
  cormat <- round(cor(num_data, method = "spearman"),2)
  
  upper_tri <- get_upper_tri(cormat)
  
  melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)
  # Heatmap
  
  heatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    geom_text(aes(label=value))+
    scale_fill_gradient2(low = "blue", high = "firebrick", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Spearman\nCorrelation") +
    theme_minimal()+ 
    coord_fixed()
  
  return(heatmap)
  
}




normalize_column <- function(col) {
  (col - min(col, na.rm=T)) / (max(col, na.rm=T) - min(col, na.rm=T))
}


dendroplot <- function(my_dat, dend, cutoff, size_txt){
  
  stage_labels <- tibble(
    lab=c("Stage I", "Stage II", "Stage III", "Stage IV", "Stage V", "Stage VI", "Stage VII"),
    x=c(2.5, 6, 9, 12.5, 16, 19.5, 23.5),
    col=c("black", "black",rep("white", 5))
  )
  
  dendro <- ggdendro::ggdendrogram(dend, rotate = F, theme_dendro = T)+
    geom_hline(yintercept=cutoff, linetype=2)+
    geom_tile(
      data=my_dat,
      height=0.25,
      aes(x=x, y=-0.125, fill=as.character(cluster)),
      show.legend = F
    )+
    scale_fill_manual(values= gcol_clust(),
                      breaks=c(1:7))+
    geom_text(
      data=stage_labels,
      inherit.aes = F,
      show.legend = F,
      size=size_txt,
      aes(x=x,
          y=-0.06,
          color=lab,
          label=lab)
    )+
    scale_color_manual(breaks = stage_labels$lab, values=stage_labels$col)+
    # ggnewscale::new_scale_fill()+
    # geom_tile(
    #   data=my_dat,
    #   height=0.5,
    #   aes(x=x, y=-0.75, fill=as.character(sec_col)),
    #   show.legend = T
    # )+
    # scale_fill_manual(values = c("gray80", "gray20", "gray40", "gray60"))+
    
    #ifelse(seriation_id %in% c(1:4, 6,8,10), "b", "w")  
    ggnewscale::new_scale_color()+
    geom_text(
      data=my_dat %>% 
        mutate(txt_col=ifelse(
          label %in% c(1:4, 6,8,10),
          #label %in% c("z13", "z14", "z15", "z17"), 
          "b", "w")),
      aes(x=x, 
          y=-0.14, 
          label=str_replace_all(label, "z", ""), 
          color=txt_col
      ),
      show.legend = F,
      size=size_txt,
      #color="white",
      vjust=1,
      hjust=0.5,
      #angle=270
    )+
    scale_color_manual(breaks = c("b", "w"), values = c("black", "white"))+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          legend.position = "none")
  return(dendro)
}


gcol_clust <- function(){
  
  return(
    #     c("#f8fd98", 
    #       "#bde79d",
    #       "#86cda3",
    #       "#5db0a8",
    #       "#47909f",
    #       "#326c8a",
    #       "#2d4c6b"   
    # )
    
    c("#326c8a",
      "#47909f",
      "#2d4c6b",  
      "#5db0a8",
      "#f8fd98",
      "#86cda3",
      "#bde79d"
      
      
      
    )
    
  )
  
}

final_figure_plot <- function(rel_feature_data, cutoff=1, outdir="/home/rheinnec", size_txt=2, size_pts=5){
  
  METHOD <- "complete"
  
  num_data <- rel_feature_data %>%
    select(-cell_id, -sec_id, -time) %>%
    mutate_all(
      .funs = as.numeric
    )
  
  norm_between_one_zero <- num_data  %>% as.matrix() %>% apply(., 2, normalize_column) 
  
  rownames(norm_between_one_zero) <- factor(rel_feature_data$cell_id)
  
  scaled_to_zero <- scale(num_data)
  
  norm_data <- scaled_to_zero
  rownames(norm_data) <- factor(rel_feature_data$cell_id)
  
  
  # euclidean_dist <- dist(norm_data)
  #  gower_dist <- vegan::vegdist(norm_data, method = "gower")
  
  num_matrix <- num_data %>% as.matrix()
  rownames(num_matrix) <- factor(rel_feature_data$cell_id)
  
  
  
  
  #gower_dist <- vegan::vegdist(num_matrix, method = "gower")
  gower_dist <- proxy::simil(num_matrix, method = "gower")
  
  S <- as.matrix(gower_dist)
  L <- diag(S) - S
  E <- eigen(L)
  vec.df <- data.frame(index = 1:nrow(L), 
                       value = E$vectors[,(nrow(L)-1)])
  
  #vec.df <- tibble(value=sort(E$vectors[,(nrow(L)-1)]))
  
  sorted.idx <- order(vec.df$value, decreasing = TRUE)
  
  vec.df$value <- vec.df$value[sorted.idx]
  
  
  vec.df$time <- rel_feature_data$time[sorted.idx]  
  
  
  
  
  
  #ggplot(vec.df, aes(x=-value, y=1))+geom_point()+theme_bw()
  #geom_text_repel(aes(label=index))
  
  
  seriation_id_mapping <- rel_feature_data %>% mutate(sec_id=factor(sec_id, levels=sorted.idx)) %>% 
    arrange(sec_id) %>% rownames_to_column("seriation_id") %>%
    select(cell_id, seriation_id)
  
  
  ### umap
  
  # 
  # set.seed(20240822)
  # 
  # u <- umap(as.dist(gower_dist), n_neighbors = 5)
  # colnames(u) <- c("X1", "X2")
  # 
  # u <- as_tibble(u) %>%
  #   mutate(cell_id=rel_feature_data$cell_id) %>%
  #   left_join(seriation_id_mapping)
  # 
  # #p <- 
  # umap_plot <- 
  #   ggplot(u, aes(y = X1, x= X2)) + 
  #   geom_point(size=5, shape=1, aes(color=)) #+ 
  #    geom_text_repel(aes(label = seriation_id), 
  #       show_guide = FALSE, 
  # #color="white"
  # size=3,
  # direction = "both",
  # max.overlaps = 100
  # ) + theme_bw()
  
  
  # 
  # pdf("output/umap.pdf")
  # print(umap_plot)
  # dev.off()
  
  
  
  #df <- num_data
  ## JK code
  #df$time <- as.factor(df$time)
  
  #features <- df[, 1:8] %>% as.data.frame()
  #rownames(features) <- rel_feature_data$cell_id
  
  # S <- proxy::simil(num_matrix, method = "gower")
  # 
  # S <- as.matrix(S)
  # L <- diag(S) - S
  # E <- eigen(L)
  # vec.df <- data.frame(index = 1:nrow(L), 
  #                      value = E$vectors[,(nrow(L)-1)])
  # sorted.idx <- order(vec.df$value, decreasing = TRUE)
  # 
  # rel_feature_data %>% mutate(sec_id=factor(sec_id, levels=sorted.idx)) %>% 
  #   arrange(sec_id) %>% rownames_to_column("seriation_ids")
  
  
  
  # vec.df$value <- vec.df$value[sorted.idx]
  # #p <- 
  # ggplot(vec.df, aes(x = index, y = value, colour = df$time[sorted.idx])) + 
  #   geom_point() + 
  #   geom_text_repel(aes(label = rel_feature_data$cell_id[sorted.idx]), 
  #                   show_guide = FALSE, 
  #                   max.overlaps = 10) + 
  #   guides(color = guide_legend(title = "Time (h)"))
  # 
  # vec.df %>%
  #   arrange(value)
  
  ##
  
  hc <- hclust(1-gower_dist, method=METHOD)
  
  # Convert the hierarchical clustering result to a dendrogram object
  dend_raw <- as.dendrogram(hc)
  
  
  #order.dendrogram(dend_raw)
  
  
  dend_old <- reorder(dend_raw, wts=c(
    10000000*c(1, 1, 1, 1),
    1000000*c(10000000,0,0,1,0,0,0),
    1*c(1, 100000, 100, 0),
    100000*c(1,1,1),
    200000*c(10000000,0,1,1),
    2*c(10000000,10000000,1)
  ))
  #plot(dend)
  
  order.dendrogram(dend_raw)
  
  
  
  rel_feature_data %>%
    left_join(seriation_id_mapping) %>%
    
    mutate(sec_id=factor(sec_id, levels=order.dendrogram(dend_raw))) %>%
    arrange(sec_id) %>%
    pull(seriation_id) %>%
    as.numeric()
  
  
  
  dend <- reorder(dend_raw, wts=c(
    0, # ändert alles 
    0, ## z18z3-z6
    100000, ## z2-z4 
    100, # z18z3-z6
    
    
    100000000, ## ändert alles
    0, #z8-z27
    0, #z7z10-z9z11
    
    
    0, #z7z10-z9z11
    0, #z7z10-z9z11
    10000, # z21z20 - z12z24 & z12-z24
    
    
    0, # z13-z15z14z17 (wenn man auf 10 ändert)
    0, # z15-z14z17 (wenn man auf 10 ändert)
    1, # z15-z14z17 (wenn man auf 0 ändert)
    0, # z15-z14z17 (wenn man auf 10 ändert)
    10, ## öndert viel
    
    # 19
    1000000, #z12z14-z13z11
    
    
    # 20 21
    0,
    100000, #z12z14-z13z11
    #z12z14-z13z11
    
    ## 22 23
    10000, # stage2 -stage3
    1000, ## öndert alles
    
    
    # 24 25
    2000, ## z12z14-z13z11
    20, ## z28z25-z26
    
    
    200000, ## z28z25-z26
    200000, ## z2z4-z8z27
    0   ## 25 - 28
  ))
  
  
  plot(dend)
  
  cell_id_order <- rel_feature_data %>%
    #rownames_to_column("sec_id") %>%
    mutate(sec_id=factor(sec_id, levels=order.dendrogram(dend))) %>%
    arrange(sec_id) %>%
    left_join(seriation_id_mapping, by="cell_id")
  #cutoff <- 0.2
  clusters <- cutree(hc, h = cutoff)
  
  # Convert cluster IDs to factor for coloring
  clusters_factor <- as.factor(clusters)
  
  
  my_dat <- tibble(cluster=clusters, 
                   label=rel_feature_data$cell_id,
                   order=factor(c(1:length(clusters)), level=order.dendrogram(dend))) %>%
    arrange(order) %>%
    mutate(x=c(1:length(clusters)))
  
  ### seriation plot
  
  n_neigh <- 2
  
  mn_time_data <- 
    lapply(seq(1,nrow(vec.df)), function(INDEX){
      
      val <- vec.df %>%
        filter(
          index %in% seq(INDEX-n_neigh, INDEX+n_neigh)
        ) %>%
        pull(time) %>%
        as.character() %>%
        as.numeric() %>%
        mean()
      
      tibble(
        index=INDEX,
        mn_time=val
      ) %>%
        return()
      
    }) %>%
    bind_rows() %>%
    left_join(vec.df,
              by="index") %>%
    left_join(seriation_id_mapping %>% mutate(seriation_id=as.numeric(seriation_id)), 
              by=c("index"="seriation_id")) %>%
    left_join(
      my_dat,
      by=c("cell_id"="label")
    )
  
  
  
  seriation_plot <- ggplot(mn_time_data,
                           aes(x=-value,
                           ))+
    geom_line(aes(y=mn_time),
              linetype=2)+
    geom_point(#y=13.75, 
      size=size_pts,
      aes(color=as.character(cluster),
          y=time)
    )+
    
    scale_color_manual(values=gcol_clust())+
    ggnewscale::new_scale_color()+
    
    # geom_text(aes(label=index,
    #               y=time,
    #               color=ifelse(index %in% c(1:4, 6,8,10), "b", "w")),
    #           size=size_txt,
    #           show.legend=F)+
    geom_text_repel(aes(label=index,
                  y=time
                  #color=ifelse(index %in% c(1:4, 6,8,10), "b", "w")
                  ),
                  max.overlaps = 500,
              size=size_txt,
              show.legend=F)+
    scale_color_manual(breaks = c("b", "w"),
                       values = c("black", "white"))+
    geom_point(#y=13.75,
               shape=1,
               size=size_pts,
               color="black",
               aes(y=time)
               #aes(color=as.character(cluster)

    )+
    
    
    scale_y_continuous(#limits = c(8,20), 
      breaks = c(8,12,16,20),
      expand = c(0.1,0.1))+
    scale_x_continuous(breaks = seq(-0.3, 0.3, 0.1),
                       labels = round((-1)*seq(-0.3, 0.3, 0.1),1))+
    theme_bw()+
    theme(legend.position = "none")+
    ylab("time (h)")+
    xlab("spectral seriation value")
  
  
  
  
  
  ####
  ## umap
  
  
  set.seed(20240822)
  
  u <- umap(as.dist(gower_dist), n_neighbors = 5)
  colnames(u) <- c("X1", "X2")
  
  u <- as_tibble(u) %>%
    mutate(cell_id=rel_feature_data$cell_id) %>%
    left_join(seriation_id_mapping) %>%
    left_join(my_dat %>% select(cell_id=label, cluster))
  
  #p <- 
  umap_plot <- 
    ggplot(u, aes(y = X1, x= X2)) + 
    geom_point(size=size_pts, #shape=1, 
               aes(color=as.character(cluster)),
               show.legend = F)+
    geom_point(#y=13.75,
      shape=1,
      size=size_pts,
      color="black",
      #aes(y=time)
      #aes(color=as.character(cluster)
      
    )+
    scale_color_manual(values = gcol_clust())
  
  
  ####
  
  
  
  heatmap <- make_heatmap_figure(num_data)+
    theme(axis.title = element_blank(),
          axis.text.x = element_text(angle = 315, hjust = 0))
  
  scatter_data <- rel_feature_data %>%
    left_join(my_dat %>% select(cell_id=label, cluster)) %>%
    select(-time) %>%
    mutate_at(.vars = "cluster", .funs=as.character)
  
  
  scatter_plot_data <- scatter_data %>%
    pivot_longer(
      cols=names(num_data)
    )
  
  jitterplot <-
    ggplot(scatter_plot_data,
           aes(x=cluster,
               y=as.numeric(value),
               color=cluster))+
    facet_wrap(~name, scales = "free_y", ncol=1, strip.position = "right")+
    geom_jitter()+
    geom_boxplot(alpha=0, outlier.alpha = 0)+
    theme_bw()
  
  
  
  # Plot the dendrogram with clusters colored
  
  dendro_data <- my_dat %>%
    left_join(seriation_id_mapping, by=c("label"="cell_id")) %>%
    mutate(label=seriation_id)
  
  dendro <- dendroplot(dendro_data, 
                       dend, cutoff, size_txt)
  
  
  #plot(dend)
  
  
  feature_vec <- c("has_nuctail", "vol_cell", "has_micronuc", "has_nucleo", "retort_length", "n_cryst", "has_micronemes", "is_full")
  #   
  hm_raw_f_data <- norm_between_one_zero %>% as.data.frame() %>% 
    rownames_to_column("cell_id") %>%
    left_join(seriation_id_mapping, by="cell_id") %>%
    select(-cell_id) %>%
    # mutate(cell_id=rownames(hm_raw_f_data)) %>%
    pivot_longer(cols = colnames(num_data), names_to = "feature", values_to = "val") %>%
    mutate(#cell_id=factor(str_replace_all(cell_id, "z", ""), levels=str_replace_all(cell_id_order$cell_id, "z", "")),
      cell_id=factor(seriation_id, levels=cell_id_order$seriation_id),
      feature=factor(feature, levels=rev(feature_vec)))
  
  hm_raw_features <-
    hm_raw_f_data %>%
    ggplot(.)+
    geom_tile(
      aes(y=feature, x=cell_id, fill=val)
    )+
    
    scale_fill_gradientn(limits=c(-max(hm_raw_f_data$val),max(hm_raw_f_data$val)), colors=c("blue","white","firebrick"))+
    #scale_y_discrete(breaks=rownames(norm_data), labels=str_replace_all(rownames(norm_data),"z", ""))+
    theme_bw()+
    #    theme(legend.position = "bottom",
    #          plot.margin = unit(c(0,1,0,0), "cm"))+
    xlab("cell id")
  
  
  
  
  
  
  #pca_result <- prcomp(norm_data, scale. = F)
  
  
  pca_result <- prcomp(num_data, scale. = T)
  
  
  total_variance <- sum(pca_result$sdev^2)
  
  # Calculate the percentage of variance covered by each principal component
  variance_percentage <- (pca_result$sdev^2 / total_variance) * 100
  
  
  
  plot_data <- bind_cols(
    rel_feature_data %>% select(cell_id),
    as_tibble(pca_result$x)
  ) %>%
    left_join(my_dat %>% select(cell_id=label, cluster)) %>%
    # left_join(
    #   tibble(cell_id=names(kmc$cluster), km_cluster=kmc$cluster)
    # ) %>%
    mutate_at(.vars = c("cluster"#, 
                        #"km_cluster"
    ), .funs=as.character) %>%
    left_join(
      seriation_id_mapping,
      by="cell_id"
    )
  
  
  cluster_plot <- 
    ggplot(plot_data,
           aes(y=PC1,
               x=PC2,
               #size=PC3,
               color=PC3))+
    geom_point(size=5
    )+
    geom_point(size=5,
               shape=1, 
               color="black")+
    
    scale_color_gradientn(limits=c(
      -max(abs(plot_data$PC3)),max(abs(plot_data$PC3))
    ), 
    colors=c("blue", "white","firebrick"),
    name=paste0("PC3 (", round(variance_percentage[3]), "%)"))+
    
    ggnewscale::new_scale_color()+
    geom_text(
      show.legend = F,
      #color="black",
      size=2,
      aes(label=seriation_id,#str_replace_all(cell_id, "z", ""), #color=cluster
          color=ifelse(
            between(PC3, -0.45, 0.45),
            "1",
            "2"
          )
          #color=cluster
      ),
      # direction = "both",
      #  max.overlaps = 1000
    )+
    scale_color_manual(values= c("black", "white"),
                       breaks=c("1", "2"))+
    ylab(paste0("PC1 (", round(variance_percentage[1]), "%)"))+
    xlab(paste0("PC2 (", round(variance_percentage[2]), "%)"))+theme_bw()+
    theme(
      #plot.margin = unit(c(0,0,0,0), "cm")
    )+coord_fixed()
  
  # pdf("pca.pdf", width=5, height = 5)
  # cluster_plot
  # dev.off()
  
  pcs_data <- 
    pca_result$rotation %>%
    as_tibble() %>%
    mutate(feature=rownames(pca_result$rotation)) %>%
    mutate(
      
      feature=factor(feature, levels=rev(feature_vec))) %>%
    pivot_longer(colnames(pca_result$rotation), names_to = "pc", values_to = "val") #%>%
  
  
  nminusone_features <- length(feature_vec)-1
  
  var_data <- tibble(
    pc=paste0("PC", c(1:8)),
    var=variance_percentage
  ) %>%
    mutate(
      var_end=(var+(100/nminusone_features))/(100/nminusone_features)
    )
  
  
  #      filter(pc %in% c("PC1", "PC2", "PC3", "PC4")) %>%
  contrib_plot <- 
    ggplot(
      pcs_data,
      aes(y=feature, x=pc), 
      #color=pc, group=pc)
    )+
    
    geom_point(aes(size=abs(val), 
                   color="variance\ncoverage (%)"
    ))+    
    
    
    geom_segment(data = var_data,
                 inherit.aes = F,
                 linewidth=10,
                 #color="black",
                 alpha=1,
                 aes(x=pc, 
                     y=1, 
                     color="variance\ncoverage (%)",
                     xend=pc, 
                     yend=var_end))+
    
    scale_color_manual(values = c("gray80"), name="")+
    
    geom_text(
      data = var_data,
      inherit.aes = F,
      #linewidth=10,
      color="gray80",
      alpha=1,
      vjust=0,
      aes(x=pc, 
          y=var_end,
          #color=ifelse(between())
          label=paste0(round(var)),#"%")
          #xend=pc, 
          #yend=var_end
      )
      
    )+
    #+
    ggnewscale::new_scale_color()+ 
    geom_point(aes(size=abs(val), color=val))+
    scale_color_gradientn(name="contribution",
                          limits=c(-max(abs(pca_result$rotation)),
                                   max(abs(pca_result$rotation))),
                          breaks=c(-max(abs(pca_result$rotation))+0.1,
                                   0,
                                   max(abs(pca_result$rotation))-0.1),
                          labels=c("negative","0", "positive"),
                          #limits=c(-1,1),
                          #limits=c(-max(norm_data),max(norm_data)),
                          colors=c("blue", "white","firebrick"))+
    
    scale_size(guide = 'none')+
    
    
    #geom_line()+
    theme_bw()+
    xlab("principal component")
  
  
  #+
  # theme(axis.title.y = element_blank(),
  #       legend.position = "bottom",
  #                   axis.ticks.y = element_blank(),
  #                   axis.text.y = element_blank(),
  #                   plot.margin = unit(c(0,0,0,0), "cm")
  #   )#+
  #ylab("contribution")+
  #coord_flip()
  
  
  
  
  
  ###### combine all plots    
  
  ## dendro + heatmap
  
  neds_labs <- tibble(
    
    feature=feature_vec,
    neds_lab=c(
      "nuclear tail",
      "cellular volume",
      "nuclear appendix",
      "nucleolus",
      "retort length",
      "crystalloids",
      "micronemes",
      "nuclear relocalization"
    )
    
  )
  
  p_empty <- ggplot(tibble(lab=""))+
    geom_text(aes(x=1, y=1, label=lab))+
    theme_void() 
  
  hm_dendro <- 
    cowplot::plot_grid(
      
      dendro+
        scale_x_discrete(limits=c("1","25"))+
        theme_void(),
      p_empty,
      hm_raw_features+
        scale_y_discrete(
          position = "right",
          breaks=neds_labs$feature,
          labels=neds_labs$neds_lab
        )+
        theme(legend.position = "none",
              axis.text.y.right = element_text(hjust=0.5)),
      ncol=1,
      align="v",
      axis = "tblr",
      rel_heights=c(2,0.5,2)
      
    )
  
  hm_dendro_in <- 
    cowplot::plot_grid(
      
      dendro+
        scale_x_discrete(limits=c("1","25"))+
        theme_void()+
        theme(
          #text = element_text(family = "Arial")
          ),
      p_empty,
      hm_raw_features+
        scale_y_discrete(
          position = "right",
          breaks=neds_labs$feature,
          labels=neds_labs$neds_lab
        )+
        theme(legend.position = "none",
              #text = element_text(family = "Arial"),
              axis.text.y.right  = element_text(hjust=0.5)),
      ncol=1,
      align="v",
      axis = "tblr",
      rel_heights=c(0.5,0.18,0.7)
      
    )
  
  # comb_plots_out <- 
  #   cowplot::plot_grid(
  #     
  #     umap_plot+
  #       geom_text_repel(aes(label=seriation_id),
  #                 size=size_txt,
  #                 max.overlaps = 500)+
  #       # geom_text(aes(label=seriation_id),
  #       #           size=size_txt)+
  #       theme_bw(),
  #     hm_dendro,
  #     ncol=2,
  #     rel_widths = c(1,2.6)
  #     
  #   ) 
  
  comb_plots_in <- 
    cowplot::plot_grid(
      
      umap_plot+
        ggnewscale::new_scale_color()+
        # geom_text(aes(label=seriation_id,
        #               color=ifelse(seriation_id %in% c(1:4, 6,8,10), "b", "w")),
        #           size=size_txt,
        #           show.legend=F)+
        geom_text_repel(aes(label=seriation_id,
                      #color=ifelse(seriation_id %in% c(1:4, 6,8,10), "b", "w")
                      ),
                  size=size_txt,
                  show.legend=F,
                  max.overlaps = 500)+
        scale_color_manual(breaks = c("b", "w"),
                           values = c("black", "white"))+
        scale_y_continuous(breaks=c(-8:5))+
        scale_x_continuous(breaks = c(-4:5))+
        theme_bw()+
        theme(panel.grid.minor = element_blank()),#+
        #coord_fixed(),
      hm_dendro_in,
      ncol=2,
      rel_widths = c(1,2.6),
      axis = "tl",
      align="hv"
      
    ) 
  
  ser_in_main_fig <-
    
    cowplot::plot_grid(
      seriation_plot+theme(panel.grid.minor = element_blank()),
      comb_plots_in,
      ncol=1,
      rel_heights = c(2,2)
      #axis = "l",
      #align="v"
    )
  
  
  ## new arrangement... seriation inside
  
  
  ser_in_comb <-
    cowplot::plot_grid(
      
      seriation_plot+
        scale_y_continuous(
          position = "right",
          breaks = c(8,12,16,20),
          expand = c(0.2,0.2)
        )+
        theme(
          #text = element_text(family = "Arial")
          #axis.title.x = element_blank(),
                           ),
      hm_dendro_in,
      ncol=1,
      rel_heights = c(1,3)
      #axis = "tl",
      #align="hv"
      
    ) %>%
    cowplot::plot_grid(
      
      umap_plot+
        ggnewscale::new_scale_color()+
        # geom_text(aes(label=seriation_id,
        #               color=ifelse(seriation_id %in% c(1:4, 6,8,10), "b", "w")),
        #           size=size_txt,
        #           show.legend=F)+
        geom_text_repel(aes(label=seriation_id),
                         max.overlaps=1000,
                            #color=ifelse(seriation_id %in% c(1:4, 6,8,10), "b", "w")
        
        size=size_txt,
        show.legend=F,
        #max.overlaps = 500
        )+
        scale_color_manual(breaks = c("b", "w"),
                           values = c("black", "white"))+
        scale_y_continuous(breaks=seq(-8,5,2))+
        scale_x_continuous(breaks = seq(-4, 5, 2))+
        theme_bw()+
        theme(panel.grid.minor = element_blank()),
              #text = element_text(family = "Arial")),#+
        #coord_fixed(),
      p_empty,
      .,
      ncol=3,
      rel_widths = c(1,0.1,5)
      
    )
  
  
  
  dir.create(outdir)
  
  write_csv(seriation_id_mapping,
            file=file.path(outdir, "cell_id_mapping.csv"))
  
  
  pdf(file.path(outdir, "combined_with_seriation_inside.pdf"), width=12.5, height=6.5)
  print(ser_in_comb)
  dev.off()
  
  pdf(file.path(outdir, "combined_with_seriation.pdf"), width=18, height=11.46497)
  print(ser_in_main_fig)
  dev.off()
  
  pdf(file.path(outdir, "combined_in.pdf"), width=18, height=8)
  print(comb_plots_in)
  dev.off()
  
  pdf(file.path(outdir, "dendrogram.pdf"), width=10, height=7)
  print(dendro)
  dev.off()
  
  pdf(file.path(outdir, "pca.pdf"), width=10, height=7)
  print(cluster_plot)
  dev.off()
  
  pdf(file.path(outdir, "contributors.pdf"), width=7, height=7)
  print(contrib_plot)
  dev.off()
  
  pdf(file.path(outdir, "heatmap.pdf"), width=10, height=7)
  print(hm_raw_features)
  dev.off()
  
  pdf(file.path(outdir, "corr_heatmap.pdf"), width=10, height=7)
  print(heatmap)
  dev.off()
  
  
  pdf(file.path(outdir, "jitterplot.pdf"), width=10, height=15)
  print(jitterplot)
  dev.off()
  
  
  WIDTH <- 25
  HEIGHT <- 15
  
  
  ggsave(file.path(outdir, "combined_with_seriation_new.svg"), plot=ser_in_comb, width = 12.5*0.7, height= 6.5*0.7, unit="in")
  
  ggsave(file.path(outdir, "pca.svg"), plot=cluster_plot, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "dendrogram.svg"), plot=dendro, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "contributors.svg"), plot=contrib_plot, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "heatmap.svg"), plot=hm_raw_features, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "combined_with_seriation.svg"), plot=ser_in_main_fig, width = WIDTH, height= HEIGHT, unit="cm")
  
  
  ggsave(file.path(outdir, "combined_in.svg"), plot=comb_plots_in, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "corr_heatmap.svg"), plot=heatmap, width = WIDTH, height= HEIGHT, unit="cm")
  
  ggsave(file.path(outdir, "jitterplot.svg"), plot=jitterplot, width = WIDTH, height= HEIGHT, unit="cm")
  
  return(ser_in_comb)
}