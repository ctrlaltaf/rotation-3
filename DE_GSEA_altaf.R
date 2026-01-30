# based off DE_GSEA_pipe_EXAMPLE.R file


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# 2) Set Bioconductor version compatible with your R
BiocManager::install(version = "3.22", ask = FALSE)

# 3) Install required Bioconductor packages
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db"
), ask = FALSE)

# 4) Install CRAN packages used in your script
install.packages(c(
  "msigdbr",
  "openxlsx",
  "dplyr",
  "tibble"
), repos = "https://cloud.r-project.org")

stopifnot(exists("untreated_all"))

source("/Users/altafbarelvi/Desktop/DE_GSEA_functions.R")

out_root <- "/Users/altafbarelvi/Code/rotation-3/outputs/gsea/"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

out_de_dir   <- file.path(out_root, "DE")
out_gsea_dir <- file.path(out_root, "GSEA")
dir.create(out_de_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_gsea_dir, showWarnings = FALSE, recursive = TRUE)


de_one_vs_all_path <- run_de_one_vs_all_findallmarkers(
  seurat_obj = untreated_all,
  group_by = "seurat_clusters",     # change to any metadata column (e.g., "Dose")
  assay = "RNA",                    # set if multi-assay; otherwise NULL
  slot = "data",                    # typical for log-normalized ("data")
  test_use = "wilcox",           # Consider other tests for DE and which we should use here
  min_pct = 0.1,
  logfc_threshold = 0.25,
  only_pos = FALSE,
  out_de_xlsx = "Clusters_OneVsAll_DE.xlsx",
  out_dir = out_de_dir
)
