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

############################################################
# Auxiliary functions
############################################################
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

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(...), "\n", sep = "")
}

############################################################
# Setting working directory and data access
############################################################
# To be able to execute it both in cluster and locally :
# ---- Paths: default (local) + override (cluster) ----
args <- commandArgs(trailingOnly = TRUE)

# Defaults for your laptop (edit these 2 lines once)
raw_root_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/RAW/XENIUM"
out_root_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/QC"

# If you run on the cluster, pass two args: raw_root out_root
raw_root <- if (length(args) >= 1) args[1] else raw_root_default
out_root <- if (length(args) >= 2) args[2] else out_root_default

raw_root <- normalizePath(raw_root, mustWork = TRUE)
out_root <- normalizePath(out_root, mustWork = FALSE)

p <- file.path

# Safety: don't ever write into raw
if (identical(raw_root, out_root)) stop("out_root cannot be the same as raw_root")
raw_root_sep <- paste0(raw_root, .Platform$file.sep)
out_root_sep <- paste0(out_root, .Platform$file.sep)
if (startsWith(out_root_sep, raw_root_sep)) stop("out_root cannot be inside raw_root")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# Sample list = folder names inside raw_root
file_list <- list.dirs(raw_root, full.names = FALSE, recursive = FALSE)

log_msg("RAW root directory: %s", raw_root)
log_msg("Output directory: %s", out_root)
log_msg("Samples detected: %s", paste(file_list, collapse = ", "))

options(future.globals.maxSize = 8000 * 1024^5) # For big size samples this parameter has to be modified or sent to the cluster
future.seed=TRUE
counts_list <- vector("list", length(file_list))
names(counts_list) <- file_list

############################################################
# Output folder + accumulator for downstream Cohen's D analysis
############################################################
# communities.txt (written at the end) holds the per-(sample, Kclusters)
# cell-type ENRICHMENT matrix (avg_log2FC from FindAllMarkers). CohensD.R
# reads this directly to label communities rich/poor.
communities_dir <- file.path(out_root, "communities")
dir.create(communities_dir, recursive = TRUE, showWarnings = FALSE)

# One enrichment table per sample is appended here, then bound and saved.
communities_all <- list()

# Neighbourhood with KMEANS
for(i in file_list){
  if (i == "DU68"){ next } # DU68 has empty files
  path_in  <- p(raw_root, i)
  path_out <- p(out_root, i)
  
  log_msg("--------------------------------------------------")
  log_msg("Processing sample: %s", i)
  log_msg("Input path: %s", path_in)
  log_msg("Output path: %s", path_out)
  
  # ---- Resume support ----
  # If this sample's enrichment file already exists, it was completed in a
  # previous (possibly killed) run. Skip recomputing it so the job can be
  # executed in batches. Delete the file to force reprocessing.
  enrich_file <- file.path(communities_dir, i, paste0(i, "_community_enrichment.txt"))
  if (file.exists(enrich_file)) {
    log_msg("Sample %s already processed (%s exists). Skipping.", i, enrich_file)
    next
  }
  
  # Creates the folders needed to save the images and data
  dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, i), recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, "RDS"), recursive = TRUE, showWarnings = FALSE)
  dir.create(p(communities_dir, i), recursive = TRUE, showWarnings = FALSE)
  
  pvp <- readRDS(file.path(path_in, "RDS", paste0(i, ".rds")))
  
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
  
  ggsave(filename = paste0(path_out, "/Elbow_Kneedle.png"), plot = p_elbow, width = 7, height = 5)
  
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
    kmeans_cluster = k_means_res$cluster,
    singler_labels = pvp$singler_labels[match(rownames(nn_mat), Cells(pvp))]
  )
  
  rownames(k_means_id) <- k_means_id$cell_id
  
  pvp@meta.data$cell_id <- rownames(pvp@meta.data)
  pvp@meta.data[["kmeans_cluster"]] <- k_means_id[rownames(pvp@meta.data), "kmeans_cluster"]
  
  # Save per-cell community assignment for downstream Cohen's D analysis.
  # CohensD.R reads this file at communities/<sample>/<sample>_Kclusters.txt
  # to count cell types inside each rich community.
  write.table(
    k_means_id,
    file = file.path(communities_dir, i, paste0(i, "_Kclusters.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
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
  ggsave(filename = paste0(path_out,"/Spatial_Neighbourhood_K", nk, ".png"),plot = p,height = h,width =w*2, limitsize = F)
  ggsave(filename = paste0(path_out, "/", i, "/Heat_K", nk, ".png"),plot = hm,height = 10,width =10, limitsize = F)
  
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
  
  # Save sample-level community ENRICHMENT table for downstream Cohen's D.
  #   rows of mat3_t  : cell types
  #   columns of mat3_t: kmeans communities
  # We transpose so each row is one (sample, Kclusters) community and each
  # column is a cell type holding its avg_log2FC enrichment.
  community_sample <- as.data.frame(t(mat3_t), check.names = FALSE)
  community_sample$Kclusters <- rownames(community_sample)
  community_sample$sample    <- i
  community_sample <- community_sample %>%
    dplyr::select(sample, Kclusters, dplyr::everything())
  
  communities_all[[i]] <- community_sample
  
  # Persist this sample's enrichment to disk IMMEDIATELY. This way, if the job
  # is later killed (time/memory limit), every sample finished so far is
  # already saved and the run can be resumed in batches without redoing them.
  write.table(
    community_sample,
    file = file.path(communities_dir, i, paste0(i, "_community_enrichment.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
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
  ggsave(paste0(path_out, "/Community_RichPoor_Heatmap.png"), p, width = 8, height = 5, dpi = 300)
  
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
  
  ggsave(paste0(path_out, "/RichPoor_Communities.png"), p_epi | p_fib, height = 8, width = 12, limitsize = FALSE)
  
  saveRDS(object = pvp,file = paste0(path_out,"/RDS/",i, ".rds"))
}

############################################################
# Save combined community enrichment table for all samples
############################################################
# IMPORTANT: we rebuild communities.txt by reading every per-sample
# enrichment file present ON DISK, not from the in-memory `communities_all`
# list. This makes the consolidation robust to batched / interrupted runs:
# whatever samples have finished so far (in this run or previous ones) are
# included, and re-running the script always regenerates a complete
# communities.txt from all enrichment files available.
#
# Each per-sample file: communities/<sample>/<sample>_community_enrichment.txt
# One row per (sample, Kclusters); one numeric column per cell type holding
# its avg_log2FC enrichment. Different samples may contain different cell
# types, so bind_rows fills the gaps with NA; we set those to 0 (baseline,
# i.e. not enriched) so communities.txt is a clean numeric matrix that
# CohensD.R can threshold directly (>0 rich, <0 poor, 0 baseline).

enrich_files <- list.files(
  communities_dir,
  pattern    = "_community_enrichment\\.txt$",
  recursive  = TRUE,
  full.names = TRUE
)

if (length(enrich_files) == 0) {
  warning("No per-sample enrichment files found in ", communities_dir,
          ". communities.txt was not written.")
} else {
  communities_all_disk <- lapply(enrich_files, function(f) {
    read.table(f, sep = "\t", header = TRUE,
               check.names = FALSE, stringsAsFactors = FALSE)
  })
  
  communities_all_df <- dplyr::bind_rows(communities_all_disk)
  communities_all_df[is.na(communities_all_df)] <- 0
  
  write.table(
    communities_all_df,
    file = file.path(communities_dir, "communities.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  log_msg("Consolidated %d per-sample enrichment files into: %s  (%d rows x %d cols)",
          length(enrich_files),
          file.path(communities_dir, "communities.txt"),
          nrow(communities_all_df), ncol(communities_all_df))
}