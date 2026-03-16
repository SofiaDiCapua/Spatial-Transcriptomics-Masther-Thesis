# How you run it in both places
## Laptop (no args, uses defaults)_
##    Rscript qc.R
## Cluster (pass paths, no code changes):
##    Rscript Final_QC_Xenium.R /mnt/beegfs/scapua/DUTRENEO /mnt/beegfs/scapua/DUTRENEO_QC

############################################################
# Loading libraries
############################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(progressr)
library(ggplotify)
library(viridis)
library(patchwork)
#library(SummarizedExperiment)
#library(scuttle)
# For the distance exclusion
library(RANN)
library(scales)

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

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      sprintf(...), "\n", sep = "")
}

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
  dir.create(p(path_out, "Clustering"), recursive = TRUE, showWarnings = FALSE)
  dir.create(p(path_out, "RDS"), recursive = TRUE, showWarnings = FALSE)
  ############################################################
  # Loading the data
  ############################################################
  xenium.obj <- LoadXenium(path_in, fov = "fov", segmentations = "cell")
  ## Use the log so that the counts distribution is not skewed and can be plotted prettier
  xenium.obj$log_count <- log1p(xenium.obj$nCount_Xenium)
  xenium.obj$log_features <- log1p(xenium.obj$nFeature_Xenium)
  
  ############################################################
  # Quality Control
  ############################################################
  
  ## Visualize the number of counts per regions
  Img_nCount <- ImageFeaturePlot(object = xenium.obj,
                                 features = "log_count",
                                 fov = "fov",
                                 axes = T,
                                 size = 0.8,
                                 alpha = 1,
                                 border.size =NA) +
    scale_fill_viridis(option = "D", direction = 1)
  Img_nFeatures <- ImageFeaturePlot(object = xenium.obj,
                                    features = "log_features",
                                    fov = "fov",
                                    axes = T,
                                    size = 0.8,
                                    alpha = 1,
                                    border.size =NA) +
    scale_fill_viridis(option = "D", direction = 1)
  ggsave(p(path_out,"Intensity_x_region.png"),plot = wrap_plots(Img_nCount, Img_nFeatures),height = 7,width = 7)
  log_msg("Intensity plot saved.")
  ## Lower threshold for counts
  vln_c <- VlnPlot(xenium.obj,features = "nCount_Xenium",assay = "Xenium",group.by = "orig.ident",pt.size = 0) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0,hjust = 0.5)) +
    xlab("") +
    scale_x_discrete(breaks=c("SeuratProject"),labels=c(fl))
  
  # Violin Feature Plot
  vln_f <- VlnPlot(xenium.obj,features = "nFeature_Xenium",assay = "Xenium",group.by = "orig.ident",pt.size = 0) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 0,hjust = 0.5)) +
    xlab("") +
    scale_x_discrete(breaks=c("SeuratProject"),labels=c(fl))
  
  wrap_plots(vln_c, vln_f)
  ggsave(p(path_out,"Violin_plot.png"),plot = wrap_plots(vln_c, vln_f),height = 7,width = 7)
  log_msg("Violin plot saved.")
  # Excluding cells with 0 counts or too few counts
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 5)
  
  ## Plot the distribution of negative probes per cell and cropp (> 10% is good threshold)
  ### Compute negative counts
  # Total negative probes per cell
  xenium.obj$neg_total <-
    xenium.obj$nCount_BlankCodeword +
    xenium.obj$nCount_ControlCodeword +
    xenium.obj$nCount_ControlProbe
  # Percentage relative to total counts
  xenium.obj$neg_pct <-
    xenium.obj$neg_total / xenium.obj$nCount_Xenium * 100
  ### Distribution
  neg_plot <- ggplot(xenium.obj@meta.data,
                     aes(x = neg_total)) +
    geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
    theme_minimal() +
    labs(x = "Negative probe counts per cell",
         y = "Number of cells") + scale_x_log10()
  neg_percent <- ggplot(xenium.obj@meta.data,
                        aes(x = neg_pct)) +
    geom_histogram(bins = 100, fill = "darkred", alpha = 0.7) +
    geom_vline(xintercept = 10, linetype = "dashed", size = 1) +
    theme_minimal() +
    labs(x = "Negative probes (%) per cell",
         y = "Number of cells")
  neg_combined <- wrap_plots(neg_plot, neg_percent, ncol = 2) +
    plot_annotation(
      title = paste0("Negative Probe Quality Control - ", fl),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  ggsave(p(path_out,"Negative_probes.png"),plot = neg_combined,height = 6,width = 12)
  log_msg("Negative probe QC plot saved.")
  # Excluding
  xenium.obj <- subset(xenium.obj, subset = neg_pct <= 10)
  
  ## Remove detached regions (isolated cells) which would be the ones that have no neighbors in 50~70 microns (cells outside the tissue)
  # Computing which cells are isolated
  coords <- GetTissueCoordinates(xenium.obj[["fov"]]) %>% as.data.frame()
  rownames(coords) <- coords$cell
  # nearest-neighbor distance per sample
  xenium.obj$nn_dist_um <- NA_real_
  for (s in unique(xenium.obj$orig.ident)) {
    cells <- intersect(colnames(xenium.obj)[xenium.obj$orig.ident == s], rownames(coords))
    if (length(cells) < 3) next
    d <- nn2(as.matrix(coords[cells, c("x","y")]), k = 2)$nn.dists[,2] # nn2 = knn
    xenium.obj$nn_dist_um[cells] <- d
  }
  xenium.obj$is_isolated_50um <- xenium.obj$nn_dist_um > 70
  # Summary BEFORE excluding
  xenium.obj@meta.data %>%
    group_by(orig.ident) %>%
    summarise(
      total = n(),
      excluded = sum(is_isolated_50um, na.rm = TRUE),
      excluded_pct = percent(excluded / total, accuracy = 0.1)
    ) %>%
    arrange(desc(excluded)) %>%
    print(n = Inf)
  # Visualizing
  plot_df <- coords %>%
    transmute(
      cell = cell, x = x, y = y,
      sample = xenium.obj$orig.ident[cell],
      isolated = xenium.obj$is_isolated_50um[cell]
    )
  
  iso_plot <- ggplot() +
    geom_point(data = subset(plot_df, !isolated), aes(x = x, y = y),
               size = 0.15, alpha = 0.02) +
    geom_point(data = subset(plot_df, isolated), aes(x = x, y = y),
               size = 1.0, alpha = 1) +
    coord_equal() + theme_void() +
    facet_wrap(~sample) +
    labs(title = "Isolated cells (no neighbor within 70 µm)")
  ggsave(p(path_out,"Isolated_cells.png"),plot = iso_plot,height = 7,width = 7)
  log_msg("Isolated cells plot saved.")
  # Excluding
  xenium.obj <- subset(xenium.obj, subset = !is_isolated_50um | is.na(is_isolated_50um))
  
  ## Plot distributions of count in each sample (color de histogram by sample)
  counts_list[[fl]] <- tibble(sample = fl, nCount_Xenium = xenium.obj$nCount_Xenium)
  
  ############################################################
  # Preprocessing
  ############################################################
  # Annotate the dimensions of the sample to create with correct proportions the plots
  h <- (round(max(xenium.obj@images$fov@boundaries$centroids@coords[,1]))/1000)*2 # height
  w <- (round(max(xenium.obj@images$fov@boundaries$centroids@coords[,2]))/1000)*2 # width
  # Sketch Clustering #
  xenium.obj <- NormalizeData(xenium.obj)
  # CLUSTERING AND DIM REDUCT
  xenium.obj <- FindVariableFeatures(xenium.obj)
  xenium.obj <- ScaleData(xenium.obj)
  # we select 50,0000 cells and create a new 'sketch' assay
  xenium.obj <- SketchData(
    assay = "Xenium",
    object = xenium.obj,
    ncells = round(dim(xenium.obj@assays$Xenium$data)[2]/6), # Between 1/4 y 1/6 of the total amount of cells
    method = "LeverageScore",
    sketched.assay = "sketch"
  )
  
  # perform clustering workflow
  xenium.obj <- FindVariableFeatures(xenium.obj)
  xenium.obj <- ScaleData(xenium.obj)
  xenium.obj <- RunPCA(xenium.obj, assay = "sketch", reduction.name = "pca.sketch")
  xenium.obj <- FindNeighbors(xenium.obj, assay = "sketch", reduction = "pca.sketch", dims = 1:20)
  xenium.obj <- FindClusters(xenium.obj, cluster.name = "seurat_cluster.sketched", resolution = 2.6)
  xenium.obj <- RunUMAP(xenium.obj, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:20)
  xenium.obj <- ProjectData(
    object = xenium.obj,
    assay = "Xenium",
    full.reduction = "full.pca.sketch",
    sketched.assay = "sketch",
    sketched.reduction = "pca.sketch",
    umap.model = "umap.sketch",
    dims = 1:50,
    refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
  )
  DefaultBoundary(xenium.obj[["fov"]]) <- "segmentations" # This way the cells are outlined
  
  Idents(xenium.obj) <- "seurat_cluster.projected"
  IDP <- ImageDimPlot(xenium.obj,
                      group.by = "seurat_cluster.projected",
                      fov = "fov",
                      border.color = NA,
                      axes = T,
                      size = 0.5,
                      alpha = 1,
                      border.size =0)
  ggsave(filename = p(path_out, "Clustering", "Image_Dim_Plot_General.png"), plot = IDP, height = h, width = w)
  log_msg("ImageDimPlot saved.")
  DP <- DimPlot(object = xenium.obj,group.by = "seurat_cluster.projected",pt.size = 0.5,reduction = "full.umap.sketch")
  ggsave(filename = p(path_out, "Clustering", "Dim_Plot_General.png"), plot = DP, height = h, width = w)
  log_msg("DimPlot saved.")
  # SCT v5 Calculation AND DIM REDUCT
  xenium.obj <- SCTransform(xenium.obj,vst.flavor="v2",return.only.var.genes = F,assay = "Xenium")
  
  Img_nCount <- ImageFeaturePlot(object = xenium.obj,
                                 features = "nCount_Xenium",
                                 fov = "fov",
                                 axes = T,
                                 size = 0.5,
                                 alpha = 1,
                                 border.size =0) +
    scale_fill_viridis(option = "D", direction = 1)
  
  ggsave(filename = p(path_out,"Image_nCount_Plot_AfterQC.png"),plot = Img_nCount, height = h, width = w)
  
  # nFeatures Plot
  Img_nFeature <- ImageFeaturePlot(object = xenium.obj,
                                   features = "nFeature_Xenium",
                                   fov = "fov",
                                   axes = T,
                                   size = 0.5,
                                   alpha = 1,
                                   border.size =0) +
    scale_fill_viridis(option = "D", direction = 1)
  ggsave(filename = p(path_out,"Image_nFeatures_Plot_AfterQC.png"),plot = Img_nFeature, height = h, width = w)
  
  # Saving processed data
  saveRDS(xenium.obj, file = p(path_out, "RDS", paste0(fl, ".rds")))
  log_msg("RDS file saved.")
}

dist_count <- ggplot(bind_rows(counts_list), aes(x = nCount_Xenium, fill = sample, color = sample)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, alpha = 0.25, position = "identity") +
  geom_density(alpha = 0.8) +
  scale_x_log10() +
  geom_vline(xintercept = 4000, linetype = "dashed") +
  theme_minimal()
ggsave(p(out_root, "Distribution_counts.png"), plot = dist_count, height = 7, width = 7)