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
library(BiocParallel)
library(scrapper)
library(cowplot)
library(knitr)
library(dplyr)
library(celldex)
library(ensembldb)
library(SingleR)
library(clustree)
library(ggtext)
library(tidyverse)
library(pheatmap) 
library(clustree)
library(SingleR)
library(UCell)

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
  pvp <- readRDS(path_in)
  cat("[", fl, "] After readRDS - Sum: ", sum(pvp$nCount_Xenium), " Length(entries) : ",length(pvp$nCount_Xenium),"\n")
  
  ############################################################
  # Clustree
  ############################################################
  
  resolutions <- c(0.1,0.5, 0.8,1,1.5,2,2.5,3,3.5,4,4.5,5)
  for (res in resolutions){
    DefaultAssay(pvp) <- "sketch"
    pvp <- FindClusters(object = pvp,cluster.name = paste0("ReducedRES.",res), resolution = res, algorithm = 4)
  }
  
  cols <- grep("^sketchedRES\\.", colnames(pvp@meta.data), value = TRUE)
  keep <- complete.cases(pvp@meta.data[, cols, drop = FALSE])
  pvp_clean <- subset(pvp, cells = rownames(pvp@meta.data)[keep])
  
  p <- clustree::clustree(pvp_clean, prefix = "sketchedRES.")
  ggsave(filename = paste0(path,"/Reduced_Clustree_", fl,".png"), plot = p, width = 12, height = 20)

  # CODE OF PART 2:
  # SingleR
  
  # Get concrete clusters if clustering and singleR coincide
  
  # Log2Fold para confirmar (parejas de cell types que se confundan)
  
  # FindAllMarkers to further confirm ya en general (trios o más)
  
  # Si el top de addmodulescore no es el cell type por el cual lo etsamos haciendo volverlo a hacer todo solo para esos cell types rarillos
  
  # Si no se pueden concretar pues unasigned y a chuparla
  
}


