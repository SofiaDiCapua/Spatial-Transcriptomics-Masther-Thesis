#!/usr/bin/env Rscript

library(Seurat)
library(future)
library(ggplot2)
library(dbscan)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(reshape2)
library(Matrix)

dashline <- function() paste(rep("-", 50), collapse = "")

xlsx_available <- tryCatch({ library(xlsx); TRUE }, error = function(e) FALSE)
if (!xlsx_available) library(readxl)

options(future.globals.maxSize = 8000 * 1024^5)
future.seed <- TRUE
set.seed(1234)

log_msg <- function(fmt, ...) {
  message(sprintf(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", fmt), ...))
}

args <- commandArgs(trailingOnly = TRUE)

running_on_cluster <- length(args) >= 3

if (running_on_cluster) {
  
  neighbors_dir <- args[1]
  metadata_file <- args[2]
  raw_root      <- args[3]
  
  log_msg("Running in CLUSTER mode")
  
} else {
  
  log_msg("Running in LOCAL mode")
  
  neighbors_dir <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/COHENSD"
  
  metadata_file <- paste0(
    "C:/Users/sofia/Documents/Health_Data_Science/",
    "TFM/Code/DUTRENEO/Neighbourhood/",
    "treatment_response_table.xlsx"
  )
  
  raw_root <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/Neighbourhood"
}

dir.create(neighbors_dir, recursive = TRUE, showWarnings = FALSE)

neighbors_dir <- normalizePath(neighbors_dir, mustWork = FALSE)
metadata_file <- normalizePath(metadata_file, mustWork = TRUE)
raw_root      <- normalizePath(raw_root, mustWork = TRUE)

log_msg("Output directory : %s", neighbors_dir)
log_msg("Metadata file    : %s", metadata_file)
log_msg("Raw root         : %s", raw_root)

cohens_dir <- neighbors_dir

dir.create(cohens_dir, recursive = TRUE, showWarnings = FALSE)

# Keep only sample folders such as DU53, DU54, etc. This avoids treating an RDS folder as a sample.
file_list <- list.dirs(raw_root, full.names = FALSE, recursive = FALSE)
file_list <- file_list[grepl("^DU[0-9]+$", file_list)]

if (length(file_list) == 0) {
  stop("No sample folders found in raw_root. raw_root must contain folders like DU53/RDS/DU53.rds")
}

log_msg("Samples found: %s", paste(file_list, collapse = ", "))

if (xlsx_available) {
  metadata <- read.xlsx(metadata_file, sheetIndex = 1, header = TRUE)
} else {
  metadata <- as.data.frame(read_excel(metadata_file, sheet = 1))
}

colnames(metadata) <- tolower(colnames(metadata))

required_metadata_cols <- c("sample", "type", "treatment", "response")
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata))
if (length(missing_metadata_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_metadata_cols, collapse = ", "))
}

metadata$sample <- as.character(metadata$sample)
metadata$response <- as.character(metadata$response)
metadata$treatment <- as.character(metadata$treatment)

celltypes <- c(
  "Neurons", "B-cells", "Endothelial cells", "Plasma cells",
  "Granulocytes", "Epithelial cells", "Fibroblasts", "Macrophages",
  "Macrophages M1", "Macrophages M2", "Tregs", "CD4+ T-cells",
  "Smooth muscle", "CD8+ T-cells", "NK cells", "Monocytes"
)

# Use tibble() instead of data.frame() to avoid logical/character bind_rows conflicts.
communities <- tibble()

for (i in file_list) {
  
  if (i == "DU68") {
    log_msg("Skipping %s.", i)
    next
  }
  
  log_msg(dashline())
  log_msg("Processing sample: %s", i)
  
  sample_id <- i
  sample_dir <- file.path(neighbors_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  rds_path  <- file.path(raw_root, sample_id, "RDS", paste0(sample_id, ".rds"))
  
  if (!file.exists(rds_path)) {
    log_msg("WARNING: RDS not found for %s. Skipping.", sample_id)
    next
  }
  
  enrich_file <- file.path(
    sample_dir,
    paste0(sample_id, "_community_enrichment.txt")
  )
  
  # Load previous per-sample results only if they exist and have rows.
  if (file.exists(enrich_file)) {
    log_msg("Sample %s already processed. Loading existing enrichment.", sample_id)
    markers.df <- read.table(enrich_file, sep = "\t", header = TRUE, check.names = FALSE)
    if (nrow(markers.df) > 0) {
      markers.df$sample <- as.character(markers.df$sample)
      markers.df$kmeans_cluster <- as.character(markers.df$kmeans_cluster)
      communities <- bind_rows(communities, as_tibble(markers.df))
      next
    } else {
      log_msg("Existing enrichment for %s has 0 rows. Recomputing.", sample_id)
      file.remove(enrich_file)
    }
  }
  
  pvp <- readRDS(rds_path)
  
  mat <- pvp@images$fov$centroids@coords[, 1:2]
  rownames(mat) <- rownames(pvp@meta.data)
  
  pvp@meta.data$cell_id <- rownames(pvp@meta.data)
  
  needed_meta <- c("cell_id", "singler_labels", "kmeans_cluster")
  missing_obj_cols <- setdiff(needed_meta, colnames(pvp@meta.data))
  if (length(missing_obj_cols) > 0) {
    log_msg("WARNING: Missing columns in %s: %s. Skipping.", sample_id, paste(missing_obj_cols, collapse = ", "))
    next
  }
  
  mdata <- pvp@meta.data[, needed_meta]
  colnames(mdata)[colnames(mdata) == "singler_labels"] <- "singleR_final"
  
  mdata$cell_id <- as.character(mdata$cell_id)
  mdata$singleR_final <- as.character(mdata$singleR_final)
  mdata$kmeans_cluster <- as.character(mdata$kmeans_cluster)
  
  mdata <- mdata[!mdata$singleR_final %in% c("Unassigned", "Erythrocytes"), ]
  
  stopifnot(all(mdata$cell_id %in% rownames(mat)))
  mat <- mat[mdata$cell_id, , drop = FALSE]
  stopifnot(all(mdata$cell_id == rownames(mat)))
  
  # --------------------------------------------------------------------------
  # Build neighborhood composition matrix: rows = cells, columns = neighbour cell types
  # --------------------------------------------------------------------------
  R <- 25
  nn <- frNN(x = mat, eps = R)
  
  stopifnot(all(mdata$cell_id == names(nn$id)))
  
  nn_df <- stack(nn$id)
  colnames(nn_df) <- c("values", "ind")
  
  cell_ids <- mdata$cell_id
  
  # frNN returns integer neighbour positions, not barcodes. Convert positions -> cell_ids -> labels.
  nn_df$neighbor_celltype <-
    setNames(as.character(mdata$singleR_final), cell_ids)[
      cell_ids[as.integer(nn_df$values)]
    ]
  
  nn_df$neighbor_cluster <-
    setNames(as.character(mdata$kmeans_cluster), cell_ids)[
      cell_ids[as.integer(nn_df$values)]
    ]
  
  # For the target heatmap, count neighbouring cell types, not neighbouring k-means clusters.
  nn_count <- nn_df %>%
    filter(!is.na(neighbor_celltype)) %>%
    group_by(ind) %>%
    count(neighbor_celltype, .drop = FALSE) %>%
    pivot_wider(
      names_from = neighbor_celltype,
      values_from = n,
      values_fill = 0
    )
  
  nn_mat <- as.matrix(nn_count[, -1, drop = FALSE])
  rownames(nn_mat) <- as.character(nn_count$ind)
  colnames(nn_mat) <- as.character(colnames(nn_count)[-1])
  storage.mode(nn_mat) <- "double"
  nn_mat[is.na(nn_mat)] <- 0
  
  if (nrow(nn_mat) == 0 || ncol(nn_mat) < 2) {
    log_msg("WARNING: sample %s has insufficient neighbour cell-type features. Skipping.", sample_id)
    next
  }
  
  k_means_df <- data.frame(
    cell_id         = mdata$cell_id,
    kmeans_cluster = as.character(mdata$kmeans_cluster),
    singleR_final  = as.character(mdata$singleR_final),
    row.names      = mdata$cell_id,
    stringsAsFactors = FALSE
  )
  
  sample_dir <- file.path(neighbors_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.table(
    k_means_df,
    file = file.path(sample_dir, paste0(sample_id, "_Kclusters.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Seurat expects features x cells. Here features = neighbour cell types.
  nn_mat_t <- t(nn_mat)
  rownames(nn_mat_t) <- make.unique(as.character(rownames(nn_mat_t)))
  colnames(nn_mat_t) <- make.unique(as.character(colnames(nn_mat_t)))
  nn_mat_t <- as(nn_mat_t, "dgCMatrix")
  
  nn_obj <- CreateSeuratObject(counts = nn_mat_t, min.features = 1)
  nn_obj <- SCTransform(nn_obj, vst.flavor = "v2")
  nn_obj <- AddMetaData(nn_obj, metadata = k_means_df)
  
  Idents(nn_obj) <- "kmeans_cluster"
  
  markers <- FindAllMarkers(
    nn_obj,
    slot            = "data",
    features        = rownames(nn_obj),
    test.use        = "wilcox",
    logfc.threshold = 0,
    min.pct         = 0,
    min.diff.pct    = 0
  )
  
  if (nrow(markers) == 0) {
    log_msg("WARNING: FindAllMarkers returned 0 rows for %s. Skipping.", sample_id)
    next
  }
  
  # Keep one row per k-means community and one column per neighbour cell type.
  # Do NOT transpose here; transposing caused V1-V6 columns.
  markers.df <- markers[, c("gene", "cluster", "avg_log2FC")] %>%
    mutate(
      gene = as.character(gene),
      cluster = as.character(cluster)
    ) %>%
    pivot_wider(
      names_from = gene,
      values_from = avg_log2FC,
      values_fill = 0
    )
  
  colnames(markers.df)[colnames(markers.df) == "cluster"] <- "kmeans_cluster"
  markers.df$kmeans_cluster <- as.character(markers.df$kmeans_cluster)
  markers.df$sample <- as.character(sample_id)
  
  sample_dir <- file.path(neighbors_dir, sample_id)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.table(
    markers.df,
    file = enrich_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  communities <- bind_rows(communities, as_tibble(markers.df))
  
  rm(pvp, nn, nn_df, nn_mat, nn_obj, markers, markers.df, k_means_df, mdata, mat)
  gc()
  
  log_msg("Sample %s done.", sample_id)
}

if (nrow(communities) == 0) {
  stop("No communities were generated. Check RDS paths, sample folders, and FindAllMarkers output.")
}

communities$sample <- as.character(communities$sample)
communities$kmeans_cluster <- as.character(communities$kmeans_cluster)

communities <- merge(
  as.data.frame(communities),
  metadata[, c("sample", "type", "treatment", "response")],
  by = "sample"
)

write.table(
  communities,
  file = file.path(neighbors_dir, "communities.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_msg("communities.txt saved.")

# ------------------------------------------------------------------------------
# Cohen's D: cell-type enrichment by community
# ------------------------------------------------------------------------------

log_msg("Building xenium_metadata_all from RDS objects...")

xenium_metadata_all <- tibble()

for (i in file_list) {
  
  if (i == "DU68") next
  
  sample_id <- i
  rds_path <- file.path(raw_root, sample_id, "RDS", paste0(sample_id, ".rds"))
  
  if (!file.exists(rds_path)) {
    log_msg("WARNING: RDS not found for %s. Skipping in metadata build.", sample_id)
    next
  }
  
  pvp <- readRDS(rds_path)
  pvp@meta.data$cell_id <- rownames(pvp@meta.data)
  
  tmp_meta <- pvp@meta.data[, c("cell_id", "singler_labels", "kmeans_cluster")]
  colnames(tmp_meta)[colnames(tmp_meta) == "singler_labels"] <- "singleR_final"
  
  tmp_meta$sample <- sample_id
  tmp_meta$cell_id <- as.character(tmp_meta$cell_id)
  tmp_meta$singleR_final <- as.character(tmp_meta$singleR_final)
  tmp_meta$kmeans_cluster <- as.character(tmp_meta$kmeans_cluster)
  
  tmp_meta <- tmp_meta[!tmp_meta$singleR_final %in% c("Unassigned", "Erythrocytes"), ]
  
  xenium_metadata_all <- bind_rows(xenium_metadata_all, as_tibble(tmp_meta))
  
  rm(pvp, tmp_meta)
  gc()
}

if (nrow(xenium_metadata_all) == 0) {
  stop("xenium_metadata_all has 0 rows. Check RDS files and filtering.")
}

xenium_metadata_all <- merge(
  as.data.frame(xenium_metadata_all),
  metadata[, c("sample", "treatment", "response")],
  by = "sample"
)

xenium_metadata_all <- xenium_metadata_all[
  xenium_metadata_all$treatment == "DUTRENEO",
]

xenium_metadata_all$response[xenium_metadata_all$response == "PARTIAL"] <- "COMPLETE"

colnames(xenium_metadata_all)[
  colnames(xenium_metadata_all) == "kmeans_cluster"
] <- "Kclusters"

communities <- read.table(
  file.path(neighbors_dir, "communities.txt"),
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

if (nrow(communities) == 0) {
  stop("communities.txt has 0 rows. Re-run the sample processing section.")
}

colnames(communities)[
  colnames(communities) == "kmeans_cluster"
] <- "Kclusters"

ct_cols <- celltypes[celltypes %in% colnames(communities)]

if (length(ct_cols) == 0) {
  stop(
    "No cell-type columns found in communities.txt. Current columns are: ",
    paste(colnames(communities), collapse = ", ")
  )
}

for (ct in ct_cols) {
  cname <- paste0("C_", ct)
  communities[[cname]] <- rep("mix", nrow(communities))
  communities[[cname]][communities[[ct]] > 0] <- "rich"
  communities[[cname]][communities[[ct]] < 0] <- "poor"
}

xenium_metadata_all <- merge(
  xenium_metadata_all,
  communities[, c(
    "sample",
    "Kclusters",
    grep("^C_", colnames(communities), value = TRUE)
  )],
  by = c("sample", "Kclusters")
)

if (nrow(xenium_metadata_all) == 0) {
  stop("After merging metadata with communities, there are 0 rows. Check sample/Kclusters matching.")
}

# Relative abundance of cell types inside each rich community.
df_relative <- tibble()

community_cols <- grep("^C_", colnames(communities), value = TRUE)

for (c in community_cols) {
  
  community_name <- sub("^C_", "", c)
  
  tmp <- xenium_metadata_all %>%
    select(sample, singleR_final, response, all_of(c)) %>%
    rename(status = all_of(c)) %>%
    filter(status == "rich") %>%
    filter(singleR_final != community_name) %>%
    group_by(sample, response) %>%
    mutate(total_count = n()) %>%
    group_by(sample, response, singleR_final, total_count) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      proportion = count / total_count,
      community = community_name
    )
  
  if (nrow(tmp) > 0) {
    df_relative <- bind_rows(df_relative, tmp)
  }
}

if (nrow(df_relative) == 0) {
  stop("df_relative is empty: no cells found in rich communities. Check C_* columns in communities.")
}

# Complete missing sample/community/celltype combinations with zeros while preserving response.
sample_response <- metadata[, c("sample", "response")]
sample_response$response[sample_response$response == "PARTIAL"] <- "COMPLETE"
sample_response <- sample_response[!duplicated(sample_response$sample), ]

all_samples <- unique(df_relative$sample)
all_communities <- unique(df_relative$community)
all_celltypes <- celltypes

df_relative <- df_relative %>%
  select(sample, response, community, singleR_final, count, proportion) %>%
  complete(
    sample = all_samples,
    community = all_communities,
    singleR_final = all_celltypes,
    fill = list(count = 0, proportion = 0)
  ) %>%
  select(-response) %>%
  left_join(sample_response, by = "sample") %>%
  mutate(response = ifelse(response == "PARTIAL", "COMPLETE", response))

cohen_d_fun <- function(x, g) {
  x_R  <- x[g == "COMPLETE"]
  x_NR <- x[g == "NO"]
  
  if (length(x_R) < 2 || length(x_NR) < 2) return(NA_real_)
  
  sd_R  <- sd(x_R)
  sd_NR <- sd(x_NR)
  
  pooled_sd <- sqrt(
    ((length(x_R) - 1) * sd_R^2 + (length(x_NR) - 1) * sd_NR^2) /
      (length(x_R) + length(x_NR) - 2)
  )
  
  if (is.na(pooled_sd) || pooled_sd == 0) return(0)
  
  (mean(x_R) - mean(x_NR)) / pooled_sd
}

CD <- df_relative %>%
  group_by(community, singleR_final) %>%
  summarise(
    cohenD = cohen_d_fun(proportion, response),
    pvalue = tryCatch(
      wilcox.test(proportion ~ response)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

CD$cohenD[is.na(CD$cohenD)] <- 0
CD$cohenD[CD$cohenD > -1 & CD$cohenD < 1] <- 0
CD$cohenD[CD$cohenD > 2.4] <- 2.4
CD$cohenD[CD$cohenD < -2.4] <- -2.4

write.table(
  CD,
  file = file.path(cohens_dir, "cohensD_values.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

heatmap_matrix <- acast(
  CD,
  singleR_final ~ community,
  value.var = "cohenD",
  fill = 0
)

heatmap_matrix <- heatmap_matrix[
  intersect(celltypes, rownames(heatmap_matrix)),
  intersect(celltypes, colnames(heatmap_matrix)),
  drop = FALSE
]

colnames(heatmap_matrix) <- paste0(colnames(heatmap_matrix), "-rich")

col_fun <- colorRamp2(
  c(-2.4, -1, 0, 1, 2.4),
  c("orange", "orange", "white", "darkolivegreen3", "darkolivegreen")
)

ht <- Heatmap(
  heatmap_matrix,
  name = "CR/PR\n\nCohen's D\nof relative\nabundance\n\nNR",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  na_col = "white",
  rect_gp = gpar(col = "white", lwd = 0.5),
  row_title = "Cell Type",
  column_title = NULL,
  row_names_side = "right",
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2.4, -1, 0, 1, 2.4),
    labels = c("-2.4", "-1", "0", "1", "2.4")
  )
)

pdf(
  file.path(cohens_dir, "CD_relative_abundance_heatmap.pdf"),
  width = 12,
  height = 8
)
draw(ht, heatmap_legend_side = "right")
dev.off()

png(
  file.path(cohens_dir, "CD_relative_abundance_heatmap.png"),
  width = 3600,
  height = 2400,
  res = 300
)
draw(ht, heatmap_legend_side = "right")
dev.off()

log_msg("Done! Heatmap saved to: %s", cohens_dir)
