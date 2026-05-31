# How you run it in both places
## Laptop (no args, uses defaults)_
##    Rscript qc.R
## Cluster (pass paths, no code changes):
##    Rscript Final_QC_Xenium.R /mnt/beegfs/scapua/DUTRENEO /mnt/beegfs/scapua/DUTRENEO_QC
## To run it in R inside the cluster:
## source("Final_Annotation_Xenium.R")

############################################################
# Loading libraries
############################################################
library(Seurat)
library(ggplot2)
library(progressr)
library(ggplotify)
library(viridis)
library(patchwork)
library(BiocParallel)
library(scrapper)
library(cowplot)
library(knitr)
library(celldex)
library(ensembldb)
library(SingleR)
library(clustree)
library(ggtext)
library(tidyverse)
library(pheatmap)
library(UCell)
library(dplyr)
library(tidyr)
library(purrr)
############################################################
# Helper functions
############################################################
log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(...), "\n", sep = "")
}

pretty_print_markers <- function(markers_list) {
  for (name in names(markers_list)) {
    cat("\n", paste0("=== ", name, " ==="), "\n")
    cat(paste(markers_list[[name]], collapse = ", "), "\n")
  }
}

# PanglaoDB green genes (canonical genes per cell type)
panglao_canonical <- list(
  "DC"      = c("ITGAX","ZBTB46", "CD86","LAMP3","CD83","CD1A","SLAMF7","CX3CR1","AXL","DAB2"),
  "Monocytes"      = c("CD14","APOBEC3A", "CCR2","CSF3R","FCGR3A","IFITM3"),
  "Fibroblasts"      = c("VIM","PDGFRB", "LUM","COL6A2","VTN","COL1A2","COL1A1","SERPINH1","POSTN","ASPN"),
  "Plasma cells"      = c("MZB1","IGHG1", "IGKC","JCHAIN","SPAG4","TGM5","SIK1"),
  "NK cells"      = c("NKG7","KLRF1", "KLRD1","GNLY","NCR1","GZMA","HOPX","ITGAM","TGFB1","GZMB"),
  "Smooth muscle"      = c("ACTA2","MYL9", "RGS5","HHIP","TAGLN","MYH11","GJA4","CNN1","DES","MYLK"),
  "Endothelial cells"  = c("CD93","PECAM1","VWF","KDR","EMCN","EGFL7","FLT1","ID3","CDH5","CLDN5"),
  "Epithelial cells"  = c("KRT14","EPCAM","CDH1","CD24","MUC1","KRT3","ANPEP","SCGB1A1","TP63"),
  "CD4+ T-cells"       = c("CD4","GZMK","IL10"),  
  "CD8+ T-cells"       = c("CD8A","ZNF683","GNG4","PDCP1","TOX2","SATB1","CCR9"),
  "Tregs"              = c("FOXP3","IL2RA","CTLA4","IKZF2","TIGIT","TNFRSF18","CCR8","ENTPD1", "LAG3","BATF"), # From Chatgpt
  "Granulocytes"       = c("CSF3R","FCGR3B","CXCR2","S100A8","S100A9","S100A12","MPO","ELANE","LCN2","CTSG"), # From Chatgpt
  "Macrophages"        = c("CD68","FCGR1","NAAA","TYROBP","CCL12","LYZ2","AIF1","SEPP1"),
  "B-cells"            = c("MS4A1","CD79A","PXK","CD19","IGHD","CD74","BANK1","IGHM")
)

# ---- Function ----
check_top_markers_vs_reference <- function(top_markers_tbl,
                                           reference_markers,
                                           cluster_col = "cluster",
                                           gene_col = "gene") {
  # It intersects the canonical markers per cell type with the found markers
  # Its aim is to simply checking if the markers make sense
  tbl <- top_markers_tbl %>%
    dplyr::mutate(
      .cluster = as.character(.data[[cluster_col]]),
      .gene = toupper(as.character(.data[[gene_col]]))
    ) %>%
    dplyr::filter(.cluster %in% names(reference_markers))
  
  tbl %>%
    dplyr::group_by(.cluster) %>%
    dplyr::summarise(
      top_genes = list(unique(.gene)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      canonical_genes = purrr::map(.cluster, ~ toupper(reference_markers[[as.character(.x)]])),
      overlap_genes   = purrr::map2(top_genes, canonical_genes, intersect),
      missing_genes   = purrr::map2(canonical_genes, overlap_genes, setdiff),
      n_overlap       = purrr::map_int(overlap_genes, length),
      overlap_pct     = round(100 * n_overlap / purrr::map_int(top_genes, length), 1),
      overlap_genes   = purrr::map_chr(overlap_genes, ~ if (length(.x) == 0) "-" else paste(.x, collapse = ", ")),
      missing_genes   = purrr::map_chr(missing_genes, ~ if (length(.x) == 0) "-" else paste(.x, collapse = ", "))
    ) %>%
    dplyr::rename(cell_type = .cluster) %>%
    dplyr::select(cell_type, n_overlap, overlap_pct, overlap_genes, missing_genes)
}

# Build a top-N markers-per-cluster table from either a FindAllMarkers result
# (cluster column present) or a pairwise FindMarkers result (needs idents).
top_n_markers <- function(obj, idents = NULL, n = 10,
                          logfc.threshold = 0.25, min.pct = 0.25,
                          p_adj_cutoff = 0.05) {
  if (is.null(idents) || length(idents) > 2) {
    # Multi-class (within-subset) DE
    m <- FindAllMarkers(obj, only.pos = TRUE,
                        min.pct = min.pct, logfc.threshold = logfc.threshold)
    m %>%
      dplyr::filter(p_val_adj < p_adj_cutoff, avg_log2FC > 0) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(avg_log2FC, n = n, with_ties = FALSE) %>%
      dplyr::ungroup()
  } else {
    # Pairwise DE; split by sign and re-label
    m <- FindMarkers(obj, ident.1 = idents[1], ident.2 = idents[2],
                     only.pos = FALSE, min.pct = min.pct,
                     logfc.threshold = logfc.threshold)
    m$gene <- rownames(m)
    pos <- m %>% dplyr::filter(p_val_adj < p_adj_cutoff, avg_log2FC > 0) %>%
      dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
      dplyr::mutate(cluster = idents[1]) %>% dplyr::slice_head(n = n)
    neg <- m %>% dplyr::filter(p_val_adj < p_adj_cutoff, avg_log2FC < 0) %>%
      dplyr::arrange(avg_log2FC) %>%
      dplyr::mutate(cluster = idents[2], avg_log2FC = -avg_log2FC) %>%
      dplyr::slice_head(n = n)
    dplyr::bind_rows(pos, neg) %>%
      dplyr::select(cluster, gene, avg_log2FC, p_val_adj, dplyr::everything())
  }
}

# Identify "unconcrete" groups from a plot_df of (signature, singler_label, sum_score).
# For each candidate ct: if the top-scoring label for ct's signature is NOT ct,
# the group is {ct} ∪ {all labels scoring strictly above ct}. Skips NA and
# "Unassigned". Deduplicates groups by sorted label set.
find_unconcrete_groups <- function(plot_df, candidate_cts) {
  groups <- Filter(Negate(is.null), lapply(candidate_cts, function(ct) {
    sub <- plot_df %>%
      dplyr::filter(signature == ct, !is.na(sum_score)) %>%
      dplyr::mutate(singler_label = as.character(singler_label)) %>%
      dplyr::filter(!is.na(singler_label), singler_label != "Unassigned") %>%
      dplyr::arrange(dplyr::desc(sum_score))
    if (nrow(sub) == 0) return(NULL)
    if (sub$singler_label[1] == ct) return(NULL)            # ct wins -> concrete
    ct_score <- sub$sum_score[sub$singler_label == ct]
    if (length(ct_score) == 0) return(NULL)
    above <- setdiff(sub$singler_label[sub$sum_score > ct_score[1]], ct)
    if (length(above) == 0) return(NULL)
    c(ct, above)
  }))
  if (length(groups) == 0) return(list())
  groups[!duplicated(
    vapply(groups, function(g) paste(sort(g), collapse = "__"), character(1)))]
}

# Score a Seurat object with per-cluster signatures, then save one bar plot
# per signature showing mean score across SingleR labels.
score_and_plot_signatures <- function(obj, top_markers, out_dir,
                                      score_prefix = "Score_",
                                      plot_title_suffix = "",
                                      width = 12, height = 6,
                                      nbin = 24, ctrl = 10) {
  sigs <- top_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(genes = list(unique(gene)), .groups = "drop")
  sigs <- setNames(sigs$genes, as.character(sigs$cluster))
  sigs <- sigs[lengths(sigs) > 0]
  if (length(sigs) == 0) return(obj)
  
  raw_prefix <- paste0(score_prefix, "_raw_")
  # AddModuleScore can fail on small/sparse subsets when nbin > distinct
  # average-expression values. Try the requested nbin first, then fall back.
  ams_attempt <- function(nbin_try) {
    tryCatch(
      AddModuleScore(obj, features = sigs, ctrl = ctrl,
                     nbin = nbin_try, name = raw_prefix),
      error = function(e) NULL
    )
  }
  scored <- ams_attempt(nbin)
  if (is.null(scored)) {
    for (fallback in c(10, 5, 3)) {
      if (fallback >= nbin) next
      log_msg("    AddModuleScore failed at nbin=%d; retrying with nbin=%d",
              nbin, fallback)
      scored <- ams_attempt(fallback)
      if (!is.null(scored)) { nbin <- fallback; break }
    }
  }
  if (is.null(scored)) {
    log_msg("    AddModuleScore failed at all nbin values; skipping scoring.")
    return(list(obj = obj, plot_df = NULL))
  }
  obj <- scored
  raw_cols <- grep(paste0("^", raw_prefix), colnames(obj@meta.data), value = TRUE)
  colnames(obj@meta.data)[match(raw_cols, colnames(obj@meta.data))] <-
    paste0(score_prefix, names(sigs))
  
  plot_df <- obj@meta.data %>%
    dplyr::mutate(singler_label = singler_labels) %>%
    dplyr::select(singler_label, dplyr::starts_with(score_prefix)) %>%
    tidyr::pivot_longer(dplyr::starts_with(score_prefix),
                        names_to = "signature", values_to = "score",
                        names_prefix = score_prefix) %>%
    dplyr::filter(signature %in% names(sigs)) %>%
    dplyr::group_by(signature, singler_label) %>%
    dplyr::summarise(sum_score = mean(score, na.rm = TRUE), .groups = "drop")
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  for (ct in names(sigs)) {
    ct_fs <- gsub("[^A-Za-z0-9_-]+", "_", ct)
    p <- ggplot(plot_df %>% dplyr::filter(signature == ct),
                aes(x = fct_reorder(singler_label, sum_score), y = sum_score)) +
      geom_col() + coord_flip() + theme_bw() +
      labs(x = "SingleR label", y = "Mean of AddModuleScore",
           title = paste0(ct, " signature", plot_title_suffix))
    ggsave(file.path(out_dir, paste0(ct_fs, "_Signature.png")),
           p, width = width, height = height)
  }
  list(obj = obj, plot_df = plot_df)
}

############################################################
# Setting working directory and data access
############################################################
# To be able to execute it both in cluster and locally :
# ---- Paths: default (local) + override (cluster) ----
args <- commandArgs(trailingOnly = TRUE)

# Defaults for your laptop (edit these 2 lines once)
raw_root_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/QC"
out_root_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/Annotations"

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
counts_list <- vector("list", length(file_list))
names(counts_list) <- file_list
# DU68 does has empty files
# file_list <- c("DU54","DU71")
for(fl in file_list){
  
  path_in  <- p(raw_root, fl)
  path_out <- p(out_root, fl)
  
  log_msg("--------------------------------------------------")
  log_msg("Processing sample: %s", fl)
  log_msg("Input path: %s", path_in)
  log_msg("Output path: %s", path_out)
  
  # Creates the folders needed to save the images and data
  dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, "Signature"), recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, "Signature_UNCONCRETE"), recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, "RDS"), recursive = TRUE, showWarnings = FALSE)
  
  ############################################################
  # Loading the data
  ############################################################
  pvp <- readRDS(paste0(path_in,"/RDS/",fl,".rds"))
  cat("[", fl, "] After readRDS - Sum: ", sum(pvp$nCount_Xenium), " Length(entries) : ",length(pvp$nCount_Xenium),"\n")
  
  # CODE OF PART 1 (OPTIONAL):
  # Apropiate clustering --> Clustree
  
  # CODE OF PART 2: 
  
  ############################################################
  # SingleR
  ############################################################
  # UPLOADING THE REFERENCE DATASET 
  
  ref.data.bpe <- celldex::BlueprintEncodeData(ensembl=FALSE) 
  print(ref.data.bpe)
  # Keep only selected main labels
  keep <- ref.data.bpe$label.main %in% c("Epithelial cells","Fibroblasts",
                                         "Smooth muscle", "Pericytes",
                                         "Monocytes", "CD4+ T-cells","CD8+ T-cells",
                                         "B-cells", "NK cells", "Macrophages",
                                         "Endothelial cells", "Neutrophils",
                                         "Eosinophils", "DC")
  
  ref.sub <- ref.data.bpe[, keep]
  print(unique(ref.sub$label.main))
  # Extract expression data (normalized data is usually in the "data" slot)
  test_data <- GetAssayData(pvp, assay = "Xenium", layer = "data")
  
  # MODIFYING IT TO OUR TASTE
  ref.sub$label.fine[ref.sub$label.fine == "Pericytes"] <- "Smooth muscle"
  ref.sub$label.fine[ref.sub$label.fine == "mv Endothelial cells"] <- "Endothelial cells"
  ref.sub$label.fine[ref.sub$label.fine %in% c("CD4+ Tcm", "CD4+ Tem")] <- "CD4+ T-cells"
  ref.sub$label.fine[ref.sub$label.fine %in% c("CD8+ Tcm", "CD8+ Tem")] <- "CD8+ T-cells"
  ref.sub$label.fine[ref.sub$label.fine %in% c("Tregs")] <- "Tregs"
  ref.sub$label.fine[ref.sub$label.fine %in% c("Neutrophils", "Eosinophils")] <- "Granulocytes"
  ref.sub$label.fine[ref.sub$label.fine %in% c("Macrophages M1", "Macrophages M2")] <- "Macrophages"
  ref.sub$label.fine[ref.sub$label.fine %in% c("naive B-cells", "Class-switched memory B-cells", "Memory B-cells")] <- "B-cells"
  
  # EXECUTING SINGLER
  # set parallelization
  bp <- if (.Platform$OS.type == "windows") {
    SnowParam(workers = 4) # cluster
  } else {
    MulticoreParam(workers = 4) # Windows
  }
  
  # perform label transfer at the single cell-level,
  # using pseudo-bulk Chromium profiles as reference
  res <- SingleR(
    test=test_data, ref=ref.sub, 
    labels=ref.sub$label.fine, 
    aggr.ref=TRUE, BPPARAM=bp)
  
  # Saving the results in the Seurat object
  pvp$singler_labels <- res$labels
  
  # MAKING NICE PLOTS OF SINGLER RESULTS
  low_confi <- pruneScores(res, min.diff.med = 0.1)  
  res$labels[low_confi] <- "Unassigned"
  tab_labels <- sort(table(res$labels), decreasing = TRUE) 
  log_msg("SingleR results: %s", paste(names(tab_labels), tab_labels, sep = "=", collapse = ", "))
  
  
  pvp$singler_labels <- res$labels[match(rownames(pvp@meta.data), rownames(res))]
  
  p1 <- ImageDimPlot(
    object = pvp,
    fov = "fov",
    boundaries = "centroids",
    group.by = "singler_labels",
    cols = "polychrome",
    axes = TRUE
  ) 
  ggsave(filename = paste0(path_out,"/AfterPruneScore.png"), plot = p1, width = 12, height = 12, dpi = 300)
  
  umap_SR <- DimPlot(pvp, group.by = "singler_labels", cols = "polychrome", label = TRUE)
  
  ggsave(filename = paste0(path_out,"/UMAP_SingleR.png"),
         umap_SR,
         width = 10,
         height = 6)
  
  
  ############################################################
  # Get markers for each cell type
  ############################################################
  
  # GETTING ALL MARKERS
  
  # 1. Use those labels as identities
  Idents(pvp) <- "singler_labels" # So that the clusters in FindAllMarkers are the cell types
  
  # 2. Use RNA assay for marker finding
  DefaultAssay(pvp) <- "SCT"
  
  markers_subset <- FindAllMarkers(
    pvp, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 3. Top 10 genes per cell type
  # Keep only statistically significant positive markers
  # - p_val_adj < 0.05: significant after multiple-testing correction
  # - abs(avg_log2FC) > 0: gene is more expressed in that cell type than in the others
  markers_sig <- markers_subset %>%
    dplyr::filter(
      p_val_adj < 0.05, # p_val corrected for multiple testing (e.g. Bonferroni)
      avg_log2FC > 1 # Highly expressed, lower would mean "Unassigned"
    )
  
  # For each cell type, keep the top 10 genes with the highest positive log fold-change
  top10_genesXcellType <- markers_sig %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()
  
  # INSPECTING RESULT
  log_msg("Top 10 genes per cell type:\n%s", paste(capture.output(print(top10_genesXcellType)), collapse = "\n"))
  
  # Saving the intersection with the canonical genes per cell type from PanglaoDB
  panglao_check <- check_top_markers_vs_reference(
    top_markers_tbl = top10_genesXcellType,
    reference_markers = panglao_canonical)
  write.csv(panglao_check, file.path(path_out, "panglao_check_results.csv"), row.names = FALSE)
  
  # Making dot plot
  # To solve the error non-unique values when setting 'row.names' 
  # The error was due to NaNs cause by it not locking at the sketched data (the one that was clusterized)
  # When looking at the results of the clustering the full data didn't have clusters assigns to all the entries
  DefaultAssay(pvp) <- "SCT" 
  
  top5 <- markers_subset %>% 
    group_by(cluster) %>% 
    dplyr::filter((avg_log2FC) > 1) %>% 
    slice_max(order_by = (avg_log2FC), n = 5, with_ties = FALSE) %>% 
    ungroup()
  
  genes <- unique(top5$gene)
  
  DP <- DotPlot(object = pvp, features = genes,group.by = "singler_labels") + 
    RotatedAxis() + 
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(filename = paste0(path_out,"/dotplot_wide.png"), plot = DP, width = 30, height = 12, dpi = 300)
  
  
  # GETTING SCORE FOR THE TOP 10 MARKERS PER CELL TYPE
  gene_sets <- top10_genesXcellType %>%
    group_by(cluster) %>% # Cluster here refers to the cell type
    summarise(genes = list(unique(gene)), .groups = "drop")
  
  signatures <- setNames(gene_sets$genes, gene_sets$cluster)
  
  # Computing Score per cell type
  pvp <- AddModuleScore(pvp, features = signatures, ctrl = 10, name = "Top10Score_")
  
  score_cols <- grep("^Top10Score_", colnames(pvp@meta.data), value = TRUE)
  celltypes <- names(signatures)
  colnames(pvp@meta.data)[match(score_cols, colnames(pvp@meta.data))] <- paste0("Score_", celltypes)
  
  # MAKING NICE PLOTS OF THE RESULTS
  plot_df <- pvp@meta.data %>%
    mutate(singler_label = singler_labels) %>%
    dplyr::select(singler_label, starts_with("Score_")) %>%
    pivot_longer(
      cols = starts_with("Score_"),
      names_to = "signature",
      values_to = "score",
      names_prefix = "Score_"
    ) %>%
    dplyr::filter(signature %in% names(signatures)) %>%
    group_by(signature, singler_label) %>%
    summarise(sum_score = mean(score, na.rm = TRUE), .groups = "drop")
  
  for (cell_type in names(signatures) ){
    plot <- ggplot(plot_df %>% dplyr::filter(signature == cell_type),
                   aes(x = fct_reorder(singler_label, sum_score), y = sum_score)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      labs(
        x = "SingleR label",
        y = "Mean of AddModuleScore",
        title = paste0(cell_type, " signature score across SingleR labels")
      )
    ggsave(filename = paste0(path_out,"/Signature/", cell_type,"_Signature.png"),
           plot, width = 24, height = 6)
    log_msg("Differential Expression for %s", cell_type)
  }
  
  ############################################################
  # Per-group analysis of unconcrete cell types
  # Rule (applied here and again at the recheck step below):
  #   ct is unconcrete iff the top-scoring label for ct's signature is NOT ct.
  #   group = {ct} ∪ {all labels scoring above ct}.
  # Size dispatch:
  #   - 1 label above ct -> pair  -> FindMarkers
  #   - 2+ labels above  -> trio+ -> FindAllMarkers within subset
  ############################################################
  mismatched_groups <- find_unconcrete_groups(plot_df, names(signatures))
  
  log_msg("Unconcrete groups detected: %d", length(mismatched_groups))
  for (g in mismatched_groups) log_msg("  group: [%s]", paste(g, collapse = ", "))
  
  if (length(mismatched_groups) == 0) {
    log_msg("No unconcrete cell types detected for sample %s", fl)
  } else for (group_labels in mismatched_groups) {
    group_key    <- paste(sort(group_labels), collapse = "_vs_")
    group_key_fs <- gsub("[^A-Za-z0-9_-]+", "_", group_key)
    group_dir    <- file.path(path_out, "Signature_UNCONCRETE", group_key_fs)
    log_msg("Processing unconcrete group: %s", group_key)
    
    pvp_group <- subset(pvp, subset = singler_labels %in% group_labels)
    label_counts <- table(pvp_group$singler_labels)
    if (length(label_counts) < 2 || any(label_counts < 3)) {
      log_msg("  Skipping (insufficient cells): %s",
              paste(names(label_counts), label_counts, sep = "=", collapse = ", "))
      next
    }
    Idents(pvp_group)    <- "singler_labels"
    DefaultAssay(pvp_group) <- "SCT"
    
    top10 <- top_n_markers(pvp_group,
                           idents = if (length(group_labels) == 2) group_labels else NULL)
    if (nrow(top10) == 0) { log_msg("  No significant markers; skipping."); next }
    
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(top10,
              file.path(group_dir, paste0("top10_markers_", group_key_fs, ".csv")),
              row.names = FALSE)
    write.csv(check_top_markers_vs_reference(top10, panglao_canonical),
              file.path(group_dir, paste0("panglao_check_", group_key_fs, ".csv")),
              row.names = FALSE)
    
    sps <- score_and_plot_signatures(
      obj              = pvp_group,
      top_markers      = top10,
      out_dir          = group_dir,
      score_prefix     = paste0("GroupScore_", group_key_fs, "_"),
      plot_title_suffix = paste0(" (group: ", group_key, ")"),
      nbin             = 10
    )
    pvp_group     <- sps$obj
    plot_df_group <- sps$plot_df
    if (is.null(plot_df_group)) {
      log_msg("  Skipping recheck for %s: no scoring data.", group_key)
      next
    }
    
    # ---- Recheck: apply the SAME rule using the NEW per-group scores ----
    # For each ct in the group, build {ct} ∪ {labels above ct} from plot_df_group.
    # Each resulting recheck-group becomes one violin plot.
    recheck_groups <- find_unconcrete_groups(plot_df_group, group_labels)
    
    if (length(recheck_groups) == 0) {
      log_msg("  All cell types in group %s are now concrete.", group_key)
    } else for (rgrp in recheck_groups) {
      rgrp_key    <- paste(sort(rgrp), collapse = "_vs_")
      rgrp_key_fs <- gsub("[^A-Za-z0-9_-]+", "_", rgrp_key)
      log_msg("  Violin step for recheck group: [%s]", paste(rgrp, collapse = ", "))
      
      # Build per-label marker lists (top 20, |avg_log2FC| > 0.5, p_val_adj < 0.05)
      if (length(rgrp) == 2) {
        # PAIR -> single FindMarkers, split by sign
        ident_1 <- rgrp[1]; ident_2 <- rgrp[2]
        mdf <- FindMarkers(pvp_group, ident.1 = ident_1, ident.2 = ident_2,
                           logfc.threshold = 0.25)
        mdf$gene <- rownames(mdf)
        m1 <- mdf %>%
          dplyr::filter(avg_log2FC >  0.5, p_val_adj < 0.05) %>%
          dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
          head(20) %>% dplyr::pull(gene)
        m2 <- mdf %>%
          dplyr::filter(avg_log2FC < -0.5, p_val_adj < 0.05) %>%
          dplyr::arrange(avg_log2FC) %>%
          head(20) %>% dplyr::pull(gene)
        sigs_violin <- setNames(list(m1, m2), c(ident_1, ident_2))
      } else {
        # TRIO+ -> FindAllMarkers within recheck subset, top 20 per label
        pvp_rgrp <- subset(pvp_group, subset = singler_labels %in% rgrp)
        rcounts <- table(pvp_rgrp$singler_labels)
        if (length(rcounts) < 2 || any(rcounts < 3)) {
          log_msg("    Skipping (insufficient cells): %s",
                  paste(names(rcounts), rcounts, sep = "=", collapse = ", "))
          next
        }
        Idents(pvp_rgrp)       <- "singler_labels"
        DefaultAssay(pvp_rgrp) <- "SCT"
        fam <- FindAllMarkers(pvp_rgrp, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
        fam_sig <- fam %>%
          dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
          dplyr::group_by(cluster) %>%
          dplyr::arrange(dplyr::desc(avg_log2FC), .by_group = TRUE) %>%
          dplyr::slice_head(n = 20) %>%
          dplyr::ungroup()
        sigs_violin <- split(fam_sig$gene, as.character(fam_sig$cluster))
        # Keep only labels in rgrp and with non-empty signatures
        sigs_violin <- sigs_violin[intersect(rgrp, names(sigs_violin))]
      }
      
      empty <- vapply(sigs_violin, function(x) length(x) == 0, logical(1))
      if (any(empty) || length(sigs_violin) < 2) {
        log_msg("    Skipping violin: incomplete signatures (%s).",
                paste(names(sigs_violin),
                      vapply(sigs_violin, length, integer(1)),
                      sep = "=", collapse = ", "))
        next
      }
      
      pretty_print_markers(sigs_violin)
      pvp_group <- AddModuleScore_UCell(pvp_group, features = sigs_violin, name = NULL)
      
      vln <- VlnPlot(pvp_group, features = names(sigs_violin),
                     group.by = "singler_labels", pt.size = 0)
      ggsave(file.path(group_dir, paste0(rgrp_key_fs, "_Violin_plot.png")),
             vln,
             width  = max(12, 4 * length(sigs_violin)),
             height = 8)
    }
  }
  
  
  # Saving pvp data
  saveRDS(pvp, file = paste0(path_out,"/RDS/",fl,".rds"))
  log_msg("RDS file saved.")
}