library(dplyr)
library(tidyr)
library(tibble)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(xlsx)

set.seed(1234)

############################################################
# Paths
############################################################

args <- commandArgs(trailingOnly = TRUE)

neighbors_dir_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/QC"
metadata_file_default <- "C:/Users/sofia/Documents/Health_Data_Science/TFM/Code/DUTRENEO/treatment_response_table.xlsx"

# Cluster:
# Rscript CohensD.R /mnt/beegfs/scapua/DUTRENEO_Neighbourhood /mnt/beegfs/scapua/treatment_response_table.xlsx
neighbors_dir <- if (length(args) >= 1) args[1] else neighbors_dir_default
metadata_file <- if (length(args) >= 2) args[2] else metadata_file_default

neighbors_dir <- normalizePath(neighbors_dir, mustWork = TRUE)
metadata_file <- normalizePath(metadata_file, mustWork = TRUE)

communities_dir <- file.path(neighbors_dir, "communities")
cohens_dir      <- file.path(communities_dir, "cohensD_plots")

dir.create(cohens_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# Read sample treatment / response metadata
############################################################

metadata <- read.xlsx(
  metadata_file,
  sheetIndex = 1,
  header     = TRUE
)

# Excel columns come in UPPERCASE (SAMPLE, RESPONSE, TYPE, TREATMENT, ...).
# Normalize to lowercase so the rest of the pipeline (which expects 'sample',
# 'treatment', 'response') works identically to the original 5_communities.R.
colnames(metadata) <- tolower(colnames(metadata))

required_metadata_cols <- c("sample", "treatment", "response")
missing_metadata_cols <- setdiff(required_metadata_cols, colnames(metadata))

if (length(missing_metadata_cols) > 0) {
  stop(
    "Metadata file is missing columns: ",
    paste(missing_metadata_cols, collapse = ", ")
  )
}

metadata <- metadata %>% filter(treatment == "DUTRENEO")
# Change to "STANDARD" to compute the heatmap for the other treatment.

# Do NOT collapse PARTIAL into COMPLETE in `metadata`. In the original
# 5_communities.R the PARTIAL -> COMPLETE replacement is applied only to the
# per-cell table, so PARTIAL samples are excluded from Cohen's D
# (cohens_d() looks only at "COMPLETE" vs "NO").

############################################################
# Read community enrichment table
############################################################

communities_file <- file.path(communities_dir, "communities.txt")

# If the consolidated communities.txt is missing (e.g. the neighbourhood job
# was killed mid-loop before consolidation), rebuild it from the per-sample
# enrichment files written during the run. Columns are aligned by name and
# missing cell types filled with 0.
if (!file.exists(communities_file)) {
  
  enrich_files <- list.files(
    communities_dir,
    pattern    = "_community_enrichment\\.txt$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  if (length(enrich_files) == 0) {
    stop(
      "Missing communities.txt at: ", communities_file,
      "\nand no *_community_enrichment.txt files found to rebuild it.",
      "\nRun Final_Neighbourhood_Analysis.R (or merge_communities.R) first."
    )
  }
  
  message("communities.txt not found; rebuilding from ",
          length(enrich_files), " per-sample enrichment files.")
  
  tables <- lapply(enrich_files, function(f) {
    df <- read.table(f, sep = "\t", header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE)
    if (nrow(df) == 0 ||
        !all(c("sample", "Kclusters") %in% colnames(df))) return(NULL)
    df %>% mutate(sample = as.character(sample),
                  Kclusters = as.character(Kclusters))
  })
  tables <- tables[!vapply(tables, is.null, logical(1))]
  
  if (length(tables) == 0) {
    stop("All per-sample enrichment files were empty or invalid.")
  }
  
  rebuilt <- dplyr::bind_rows(tables)
  rebuilt[is.na(rebuilt)] <- 0
  
  write.table(rebuilt, file = communities_file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("Wrote consolidated communities.txt: ",
          nrow(rebuilt), " rows x ", ncol(rebuilt), " cols.")
}

communities <- read.table(
  communities_file,
  sep         = "\t",
  header      = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# communities.txt is now the ENRICHMENT MATRIX written by
# Final_Neighbourhood_Analysis.R:
#   one row per (sample, Kclusters), one numeric column per cell type holding
#   the avg_log2FC enrichment (from Seurat FindAllMarkers).
if (!all(c("sample", "Kclusters") %in% colnames(communities))) {
  stop(
    "communities.txt must contain columns 'sample' and 'Kclusters' plus one ",
    "numeric column per cell type (the avg_log2FC enrichment matrix written ",
    "by Final_Neighbourhood_Analysis.R)."
  )
}

communities <- communities %>%
  mutate(
    sample    = as.character(sample),
    Kclusters = as.character(Kclusters)
  )

# Any remaining NA (cell type absent in a sample) -> 0 (baseline, not enriched)
celltype_cols <- setdiff(colnames(communities), c("sample", "Kclusters"))
communities[celltype_cols][is.na(communities[celltype_cols])] <- 0

############################################################
# Add rich / poor / mix labels to communities
############################################################
# Same thresholding as Final_Neighbourhood_Analysis.R (thr = 0):
#   avg_log2FC > 0 -> "rich"   (enriched vs background)
#   avg_log2FC < 0 -> "poor"   (depleted vs background)
#   avg_log2FC = 0 -> "mix"    (baseline)

for (ct in celltype_cols) {
  cname <- paste0("C_", ct)
  communities[[cname]] <- "mix"
  communities[[cname]][communities[[ct]] > 0] <- "rich"
  communities[[cname]][communities[[ct]] < 0] <- "poor"
}

############################################################
# Read per-cell community assignments
############################################################
# Files written by Final_Neighbourhood_Analysis.R contain:
#   cell_id, kmeans_cluster, singler_labels
# The original 5_communities.R expects the cell-type column to be called
# `singleR_final`, so we rename on the fly.

samples_to_use <- intersect(metadata$sample, unique(communities$sample))

if (length(samples_to_use) == 0) {
  stop("No overlapping samples between metadata and communities.txt")
}

xenium_metadata_all <- data.frame()

for (s in samples_to_use) {
  
  k_file <- file.path(communities_dir, s, paste0(s, "_Kclusters.txt"))
  
  if (!file.exists(k_file)) {
    warning("Missing Kcluster file for sample: ", s)
    next
  }
  
  tmp <- read.table(
    k_file,
    sep              = "\t",
    header           = TRUE,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  
  # Accept either `singleR_final` (original naming) or `singler_labels`
  # (Final_Neighbourhood_Analysis.R naming).
  if (!"singleR_final" %in% colnames(tmp) &&
      "singler_labels" %in% colnames(tmp)) {
    tmp <- dplyr::rename(tmp, singleR_final = singler_labels)
  }
  
  required_k_cols <- c("cell_id", "kmeans_cluster", "singleR_final")
  missing_k_cols  <- setdiff(required_k_cols, colnames(tmp))
  
  if (length(missing_k_cols) > 0) {
    stop(
      "File ", k_file, " is missing columns: ",
      paste(missing_k_cols, collapse = ", "),
      "\nExpected: cell_id, kmeans_cluster, and one of singleR_final / singler_labels."
    )
  }
  
  tmp <- tmp %>%
    mutate(
      sample    = s,
      Kclusters = as.character(kmeans_cluster)
    )
  
  xenium_metadata_all <- bind_rows(xenium_metadata_all, tmp)
}

if (nrow(xenium_metadata_all) == 0) {
  stop("xenium_metadata_all is empty after reading per-sample files.")
}

############################################################
# Merge metadata and filter cells
############################################################
# Sanity check: `metadata` must be unique per `sample` so the join with the
# per-cell table is many-to-one (and not a cartesian product).

dup_meta <- metadata %>% count(sample) %>% filter(n > 1)
if (nrow(dup_meta) > 0) {
  warning(
    "metadata has duplicated rows for samples: ",
    paste(dup_meta$sample, collapse = ", "),
    ". Keeping the first row per sample."
  )
  metadata <- metadata %>% distinct(sample, .keep_all = TRUE)
}

xenium_metadata_all <- xenium_metadata_all %>%
  left_join(metadata, by = "sample", relationship = "many-to-one") %>%
  filter(
    treatment == "DUTRENEO",
    !singleR_final %in% c("Unassigned", "Erythrocytes")
  )

############################################################
# Merge rich / poor / mix community labels into cells
############################################################
# The original 5_communities.R does a `merge` here and assumes
# (sample, Kclusters) is UNIQUE in `communities`. If it isn't (e.g.
# communities.txt was built differently and has duplicated keys), the join
# would produce a cartesian explosion (millions/billions of rows). We
# diagnose and deduplicate explicitly.

c_cols <- grep("^C_", colnames(communities), value = TRUE)

communities_for_join <- communities %>%
  select(sample, Kclusters, all_of(c_cols)) %>%
  mutate(
    sample    = as.character(sample),
    Kclusters = as.character(Kclusters)
  )

n_before <- nrow(communities_for_join)
dup_keys <- communities_for_join %>%
  count(sample, Kclusters) %>%
  filter(n > 1)

if (nrow(dup_keys) > 0) {
  warning(
    "communities.txt has ", nrow(dup_keys),
    " (sample, Kclusters) combinations with multiple rows. ",
    "Keeping the first row per combination. ",
    "First duplicates:\n",
    paste(
      utils::capture.output(print(utils::head(dup_keys, 10))),
      collapse = "\n"
    )
  )
  communities_for_join <- communities_for_join %>%
    distinct(sample, Kclusters, .keep_all = TRUE)
}

cat(
  "communities.txt: ", n_before, " rows -> ",
  nrow(communities_for_join), " unique (sample, Kclusters) rows used for join.\n",
  sep = ""
)

xenium_metadata_all <- xenium_metadata_all %>%
  mutate(
    sample    = as.character(sample),
    Kclusters = as.character(Kclusters)
  ) %>%
  left_join(
    communities_for_join,
    by = c("sample", "Kclusters"),
    relationship = "many-to-one"
  )

# Mirror the original script: collapse PARTIAL into COMPLETE only in the
# per-cell table. The metadata-side `response` is left untouched, so PARTIAL
# samples remain excluded from the Cohen's D calculation.
xenium_metadata_all$response[xenium_metadata_all$response == "PARTIAL"] <- "COMPLETE"

############################################################
# Relative abundance per rich community
############################################################
# IMPORTANT: do NOT group by `response` here. Match the original script so
# totals are per sample (not per sample x response). `response` is attached
# afterwards via merge with metadata.

df_relative <- data.frame()

rich_cols <- grep("^C_", colnames(communities), value = TRUE)

for (c in rich_cols) {
  
  community_celltype <- sub("^C_", "", c)
  
  tmp <- xenium_metadata_all %>%
    select(sample, singleR_final, all_of(c)) %>%
    filter(.data[[c]] == "rich") %>%
    filter(singleR_final != community_celltype) %>%
    group_by(sample) %>%
    mutate(total_count = n()) %>%
    group_by(sample, singleR_final, total_count) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(
      proportion = count / total_count,
      community  = community_celltype
    )
  
  df_relative <- bind_rows(df_relative, tmp)
}

if (nrow(df_relative) == 0) {
  stop("df_relative is empty. Check whether any communities were labelled as rich.")
}

# Attach response from metadata (PARTIAL stays PARTIAL here, so those samples
# are naturally excluded from the Cohen's D computation below).
df_relative <- merge(
  df_relative,
  metadata[, c("sample", "response")],
  by = "sample"
)

cat("Response groups present in df_relative:\n")
print(table(df_relative$response))

############################################################
# Cohen's D helper
############################################################

cohens_d <- function(values, groups) {
  
  x_complete <- values[groups == "COMPLETE"]
  x_no       <- values[groups == "NO"]
  
  n_complete <- length(x_complete)
  n_no       <- length(x_no)
  
  if (n_complete < 2 || n_no < 2) {
    return(NA_real_)
  }
  
  sd_complete <- sd(x_complete)
  sd_no       <- sd(x_no)
  
  pooled_sd <- sqrt(
    ((n_complete - 1) * sd_complete^2 +
       (n_no - 1) * sd_no^2) /
      (n_complete + n_no - 2)
  )
  
  if (is.na(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  
  (mean(x_complete) - mean(x_no)) / pooled_sd
}

############################################################
# Cohen's D and Wilcoxon test
############################################################
# Compute per (community, singleR_final) group. We pass the proportion and
# response vectors explicitly to small helpers so there is no reliance on
# cur_data() column ordering.

wilcox_p <- function(values, groups) {
  keep <- groups %in% c("COMPLETE", "NO")
  v <- values[keep]
  g <- groups[keep]
  
  if (length(unique(g)) < 2) return(NA_real_)
  if (any(table(g) < 1))      return(NA_real_)
  
  tryCatch(
    stats::wilcox.test(v ~ g)$p.value,
    error = function(e) NA_real_
  )
}

CD <- df_relative %>%
  group_by(community, singleR_final) %>%
  summarise(
    cohenD = cohens_d(proportion, response),
    pvalue = wilcox_p(proportion, response),
    .groups = "drop"
  )

CD$cohenD[CD$cohenD > -1 & CD$cohenD < 1] <- 0

write.table(
  CD,
  file      = file.path(cohens_dir, "CohensD_results.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

############################################################
# Heatmap matrix
############################################################

heatmap_matrix <- reshape2::acast(
  CD,
  singleR_final ~ community,
  value.var = "cohenD"
)

heatmap_matrix[is.na(heatmap_matrix)] <- 0

col_fun <- circlize::colorRamp2(
  c(-2.4, -0.5, 0, 0.5, 2.4),
  c("orange", "white", "white", "white", "darkolivegreen")
)

ht <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name            = "Cohen's D",
  col             = col_fun,
  cluster_rows    = TRUE,
  cluster_columns = TRUE,
  na_col          = "white"
)

pdf(
  file.path(cohens_dir, "CD_HD_heatmap_with_partials.pdf"),
  width  = 7,
  height = 6
)

ComplexHeatmap::draw(ht, heatmap_legend_side = "right")

dev.off()