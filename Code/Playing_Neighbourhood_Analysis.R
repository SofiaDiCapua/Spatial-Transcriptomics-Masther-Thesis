library(Seurat) 
library(dplyr)
library(magrittr)
library(ggplot2)
library(raster)
#library(scclusteval)
library(ComplexHeatmap)
library(dbscan)
library(ggplotify)
# For finding the best K
library(cluster)
library(purrr)
library(tidyr)
library(tibble)
library(circlize)

set.seed(1234)

setwd("C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO") # In Rmd this line has to be executed in the terminal

fl <- c("DU53")
path <- paste0("Neighbourhood/", fl)

# Search range for K
k_min <- 2
k_max <- 15

# Majority vote helper
get_mode_k <- function(x) {
  ux <- sort(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

# Efficient helper: run kmeans once per K and reuse results
compute_kmeans_grid <- function(x, k_vec, nstart = 25, iter.max = 100) {
  res <- vector("list", length(k_vec))
  names(res) <- as.character(k_vec)
  for (kk in k_vec) {
    km <- kmeans(x, centers = kk, nstart = nstart, iter.max = iter.max)
    res[[as.character(kk)]] <- km
  }
  res
}

# Simple internal Kneedle for elbow detection
# Works well for decreasing WSS curves
find_kneedle <- function(x, y) {
  x_n <- (x - min(x)) / (max(x) - min(x))
  y_n <- (y - min(y)) / (max(y) - min(y))
  
  # For a decreasing curve, elbow = max distance from straight line
  # line from (0,1) to (1,0)  -> use y' = 1 - x
  d <- (1 - x_n) - y_n
  
  x[which.max(d)]
}

# Neighbourhood with KMEANS
for(i in fl){
  pvp <- readRDS(paste0("C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/Annotations/", i, "/RDS/", i, ".rds")) 
  dir.create(file.path("Neighbourhood", i), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path("Neighbourhood", i, "RDS"), showWarnings = FALSE, recursive = TRUE)
  
  # GETTING NEIGHBOURHOOD MATRIX ----------------------------------------------
  mat<- pvp@images$fov$centroids@coords[,1:2] 
  mat<- as.matrix(mat)
  rownames(mat)<- pvp@images$fov$centroids@cells
  mat<- mat[Cells(pvp), ]
  
  # Find the closest cells for each cell within 25 um radius.
  # Documentation has shown similar outputs from both fixed radius and knn algorithms
  eps <- 25
  nn <- dbscan::frNN(x= mat, eps = eps)
  # Create the Neighbourhood count matrix.
  nn_df <- purrr::map2_dfr(
    nn$id, seq_along(nn$id), function(neighbours, cell_index) {
      if (length(neighbours) > 0) {
        data.frame(values = neighbours, ind = cell_index)
      } else {
        data.frame(values = NA_integer_,ind = cell_index)
      }})
  # get the annotations for cells that are within the 25 um of each cell.
  cluster_ids<- pvp$singler_labels %>% unname()
  nn_df$cluster_id<- cluster_ids[nn_df$values]
  nn_df$cluster_id<- factor(nn_df$cluster_id)
  nn_count<- nn_df %>%
    group_by(ind) %>%
    count(cluster_id, .drop = FALSE)
  # pivot it to a wide format and create a cell x cluster_id matrix
  nn_count <- nn_count %>% tidyr::pivot_wider(names_from = cluster_id,
                                              values_from = n,values_fill = 0)
  # Keep real cell IDs (one per row)
  cell_ids <- rownames(mat)[nn_count$ind]
  # convert to a matrix
  nn_mat<- nn_count[,-1] %>% as.matrix() 
  # We should scale this matrix as it is new formed data, there's nothing from 
  # the scaling done in QC in this matrix. We do it below when necessary
  
  rownames(nn_mat) <- cell_ids
  
  # GETTING BEST K -------------------------------------------------------------
  k_vec <- k_min:k_max
  
  # Precompute kmeans once per K to reuse across methods
  km_grid <- compute_kmeans_grid(scale(nn_mat), k_vec, nstart = 25, iter.max = 100)
  
  # Elbow method
  # Computing within-cluster sum of squares
  wss <- sapply(k_vec, function(k) km_grid[[as.character(k)]]$tot.withinss)
  
  # Kneedle algorithm to get the best K
  best_k_elbow <- find_kneedle(k_vec, wss)
  
  elbow_df <- data.frame(k = k_vec, wss = wss)
  
  p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
    geom_line() + geom_point() +
    {if (!is.na(best_k_elbow)) geom_vline(xintercept = best_k_elbow, linetype = 2, color = "red")} +
    {if (!is.na(best_k_elbow)) geom_point(data = elbow_df[elbow_df$k == best_k_elbow, ],
                                          aes(x = k, y = wss), color = "red", size = 3)} +
    labs(
      title = paste0(i, " - Elbow method with Kneedle"),
      subtitle = paste0("Best K = ", best_k_elbow),
      x = "Number of clusters (K)",
      y = "Total within-cluster sum of squares"
    ) +
    theme_bw()
  
  ggsave(filename = paste0("Neighbourhood/", i, "/Elbow_Kneedle.png"), plot = p_elbow, width = 7, height = 5)
  
  # Silhouette method - Too much memory needed
  #dmat <- dist(nn_mat) 
  #sil_width <- sapply(k_vec, function(k) {
  #  cl <- km_grid[[as.character(k)]]$cluster
  #  ss <- cluster::silhouette(cl, dmat)
  #  mean(ss[, "sil_width"])
  #})
  #best_k_sil <- k_vec[which.max(sil_width)]
  #sil_df <- data.frame(k = k_vec, silhouette = sil_width)
  
  #p_sil <- ggplot(sil_df, aes(x = k, y = silhouette)) +
  #  geom_line() +
  #  geom_point() +
  #  geom_vline(xintercept = best_k_sil, linetype = 2, color = "red") +
  #  geom_point(data = sil_df[sil_df$k == best_k_sil, ],
  #             aes(x = k, y = silhouette), color = "red", size = 3) +
  #  labs(
  #    title = paste0(i, " - Silhouette method"),
  #    subtitle = paste0("Best K = ", best_k_sil),
  #    x = "Number of clusters (K)",
  #    y = "Average silhouette width"
  #  ) + theme_bw()
  
  #ggsave(filename = paste0("Neighbourhood/", i, "_Silhouette.png"), plot = p_sil, width = 7, height = 5)
  
  
  # Gap statistic - Too much memory needed (869.1 GB)
  # It repeatedly runs kmeans on the real data + many bootstrap reference datasets,
  # which can easily cause:
  # "cannot allocate vector of size XXX Gb"
  #
  # If you want to use it again, better run it on a subset of cells:
  # nn_mat_subset <- nn_mat[sample(nrow(nn_mat), 3000), ]
  
  # gap_stat <- cluster::clusGap(
  #   nn_mat,
  #   FUN = function(x, k) kmeans(x, centers = k, nstart = 25, iter.max = 100),
  #   K.max = k_max,
  #   B = 50
  # )
  #
  # best_k_gap <- cluster::maxSE(
  #   f = gap_stat$Tab[, "gap"],
  #   SE.f = gap_stat$Tab[, "SE.sim"],
  #   method = "firstSEmax"
  # )
  #
  # p_gap <- factoextra::fviz_gap_stat(gap_stat) +
  #   geom_vline(xintercept = best_k_gap, linetype = 2, color = "red") +
  #   ggtitle(paste0(i, " - Gap statistic")) +
  #   labs(subtitle = paste0("Best K = ", best_k_gap))
  #
  # ggsave(
  #   filename = paste0("Neighbourhood/", i, "_GapStatistic.png"),
  #   plot = p_gap,
  #   width = 7,
  #   height = 5
  # )
  
  # Since Gap statistic is commented, use only Elbow + Silhouette
  #k_votes <- c(best_k_elbow, best_k_sil)
  #k_votes <- k_votes[!is.na(k_votes)]
  #nk <- get_mode_k(k_votes)
  nk <- best_k_elbow
  
  # Cleaner print
  cat("\n")
  cat("###############################################\n")
  cat("### BEST K SELECTION :", i, "\n")
  cat("###############################################\n")
  cat("Elbow (Kneedle)   -> K =", best_k_elbow, "\n")
  cat("-----------------------------------------------\n")
  cat("FINAL SELECTED K  -> K =", nk, "\n")
  cat("###############################################\n\n")
  
  # CLUSTERING -----------------------------------------------------------------
  # k-means clustering
  # The k slected from the previous methods should be the one used (max voting strategy)
  k_means_res <- km_grid[[as.character(nk)]] # Reusing the previously computed one
  
  k_means_id <- data.frame(
    cell_id = rownames(nn_mat),
    kmeans_cluster = k_means_res$cluster
  )
  
  rownames(k_means_id) <- k_means_id$cell_id
  
  pvp@meta.data$cell_id <- rownames(pvp@meta.data)
  pvp@meta.data[["kmeans_cluster"]] <- k_means_id[rownames(pvp@meta.data), "kmeans_cluster"]
  

  h <- (round(max(pvp@images$fov@boundaries$centroids@coords[,1]))/1000)*2 # height
  w <- (round(max(pvp@images$fov@boundaries$centroids@coords[,2]))/1000)*2 # width
  
  # Plotting communities in the tissue
  # Clusters = communities as they represent a set of cellular types
  p2<- ImageDimPlot(pvp, fov = "fov", cols = "polychrome", axes = TRUE, 
                    group.by = "kmeans_cluster", border.color = NA,
                    size =0.5, alpha =1, border.size = 0)+
    theme(legend.text = element_text(size = 14),
      legend.key.size = unit(0.8, "cm"),legend.spacing.x = unit(0.3, "cm"))
  p1<- ImageDimPlot(pvp, fov = "fov", cols = "polychrome", axes = TRUE, 
                    group.by = "singler_labels", border.color = NA,
                    size =0.5, alpha =1, border.size = 0)+
    theme(legend.text = element_text(size = 14),
          legend.key.size = unit(0.8, "cm"),legend.spacing.x = unit(0.3, "cm"))
  p <- p1 | p2
  
  # Plotting hetmap 
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.rect(x = x, y = y, width = width *0.99, 
                    height = height *0.99,
                    gp = grid::gpar(col = "grey", 
                                    fill = fill, lty = 1, lwd = 0.5))
  }
  
  col_fun=circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  
  
  mat3<- scale(table(pvp$kmeans_cluster, pvp$singler_labels))
  hm<-Heatmap(t(scale(mat3)), 
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              rect_gp = grid::gpar(type = "none"),
              cell_fun = cell_fun)
  hm <- as.ggplot(hm)
  ggsave(filename = paste0("Neighbourhood/",i,"/Spatial_Neighbourhood_K", nk, ".png"),plot = p,height = h,width =w*2, limitsize = F)
  ggsave(filename = paste0("Neighbourhood/",i,"/Heat_K", nk, ".png"),plot = hm,height = 10,width =10, limitsize = F)
  
  # Replicating images e and f from figure 5 -----------------------------------
  # Pretty colors like the figure
  cols_rbp <- c("Poor" = "#D84FB3", "Baseline" = "#7F7F7F", "Rich" = "#3C84C6")
  cols_epi <- c("Baseline" = "black", "Epithelial-rich" = "#3C84C6", "Epithelial-poor" = "#D84FB3")
  cols_fib <- c("Baseline" = "black", "Fibroblasts-rich" = "#3C84C6", "Fibroblasts-poor" = "#D84FB3")
  
  # COMMUNITY X CELL-TYPE ENRICHEMNT -------------------------------------------
  # The right code I must use 
  nn_obj <- CreateSeuratObject(counts = t(nn_mat), min.features = 1)
  # Now nCount is how many neighbouring cells each cell has
  # and nFeature is how many cell types each cell has as neighbour 
  nn_obj <- SCTransform(nn_obj, vst.flavor = "v2")
  
  # Adding metadata info (clusters)
  k_means_df <- data.frame(
    cell_id = rownames(nn_mat),
    kmeans_cluster = k_means_res$cluster
  )
  rownames(k_means_df) <- k_means_df$cell_id
  pvp@meta.data$cell_id <- rownames(pvp@meta.data)
  pvp@meta.data[["kmeans_cluster"]] <- k_means_df[rownames(pvp@meta.data), "kmeans_cluster"]
  nn_obj <- AddMetaData(nn_obj, metadata = k_means_df)
  
  Idents(nn_obj) <- "kmeans_cluster"
  
  # Like this we are treating each cell type as a "gene" and the ask:
  # Which neighbouring cell types are enriched in each k-means community ?
  # FindAllMarkers() is designed to compare numerical feature values across 
  # groups of cells. It needs a matrix shaped like: features x cells
  markers <- FindAllMarkers(
    nn_obj, slot = "data", # Used the normalized data, no raw counts
    features = rownames(nn_obj), # Features are the cell types 
    test.use = "wilcox", logfc.threshold = 0, min.pct = 0, min.diff.pct = 0)
  
  # Convert marker results to matrix:
  # rows = cell types
  # columns = communities
  enrich_df <- markers %>%
    dplyr::select(kmeans_cluster = cluster, cell_type = gene, avg_log2FC ) %>%
    pivot_wider(id_cols = cell_type, names_from = kmeans_cluster,
      values_from = avg_log2FC, values_fill = 0)
  
  mat3_t <- enrich_df %>% column_to_rownames("cell_type") %>% as.matrix()
  
  # Threshold: logfold > 0 means enriched compared with background (logfold=0)
  thr <- 0
  
  # Compact ternary heatmap matrix: Poor / Baseline / Rich
  mat_disc <- ifelse(mat3_t > thr, 1, ifelse(mat3_t < -thr, -1, 0))
  
  df <- as.data.frame(mat_disc)
  df$CellType <- rownames(mat_disc)
  
  df_long <- pivot_longer(df, cols = -CellType, names_to = "Community", values_to = "value")
  
  df_long$group <- factor(df_long$value,
                          levels = c(0, 1, -1),
                          labels = c("Baseline", "Rich", "Poor"))
  
  df_long$CellType <- factor(df_long$CellType, levels = rev(rownames(mat_disc)))
  
  p <- ggplot(df_long, aes(Community, CellType, fill = group)) +
    geom_tile(color = "white", linewidth = 1) +
    scale_fill_manual(values = cols_rbp, name = NULL) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(face = "bold", color = "black"),
          axis.ticks = element_blank(),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.text = element_text(face = "bold"))
  ggsave(paste0("Neighbourhood/", i, "/Community_RichPoor_Heatmap.png"), p, width = 8, height = 5, dpi = 300)
  
  # Helper: assign rich/poor/baseline to each cell from its community
  state_from_comm <- function(lineage, prefix) {
    v <- mat3_t[lineage, ]
    s <- ifelse(v > thr, paste0(prefix, "-rich"),
                ifelse(v < -thr, paste0(prefix, "-poor"), "Baseline"))
    out <- s[as.character(pvp$kmeans_cluster)]
    names(out) <- Cells(pvp)
    
    factor(out,
           levels = c("Baseline", paste0(prefix, "-rich"), paste0(prefix, "-poor")))
  }
  
  pvp[["epithelial_state"]] <- state_from_comm("Epithelial cells", "Epithelial")
  pvp[["fibroblast_state"]] <- state_from_comm("Fibroblasts", "Fibroblasts")
  
  # Spatial plots like (f) and (g)
  p_epi <- ImageDimPlot(
    pvp, fov = "fov", group.by = "epithelial_state",
    cols = cols_epi, axes = FALSE, dark.background = TRUE, border.color = NA,
    size =0.5, alpha =1, border.size = 0
  ) + theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.title = element_blank(), legend.position = "top",
    legend.text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))
  
  p_fib <- ImageDimPlot(
    pvp, fov = "fov", group.by = "fibroblast_state",
    cols = cols_fib, axes = FALSE, dark.background = TRUE, border.color = NA,
    size =0.5, alpha =1, border.size = 0
  ) + theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.title = element_blank(), legend.position = "top",
    legend.text = element_text(size = 14), legend.key.size = unit(0.8, "cm"))
  
  ggsave(paste0("Neighbourhood/", i, "/RichPoor_Communities.png"), p_epi | p_fib, height = 8, width = 12, limitsize = FALSE)
  
  saveRDS(object = pvp,file = paste0("Neighbourhood/",i,"/RDS/",i, ".rds"))
}

