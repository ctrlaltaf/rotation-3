#!/usr/bun/env Rscript

# load in required libraries

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

# setup dirs 

out_dir   <- "outputs"
rds_dir   <- file.path(out_dir, "rds")
plots_dir <- file.path(out_dir, "plots")
tabs_dir  <- file.path(out_dir, "tables")
logs_dir  <- file.path(out_dir, "logs")

dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- "/Volumes/spencer_group/Rotation Students/Altaf/scRNAseq_raw_data"

#load in the h5 file and create seurat object

untreated_h5_file <- file.path(data_dir, "Untreated_FullDepth/filtered_feature_bc_matrix.h5")

stopifnot(file.exists(untreated_h5_file))

untreated.data <- Read10X_h5(untreated_h5_file)
untreated <- CreateSeuratObject(counts = untreated.data, project = "unt", assay = "RNA", min.cells =  3, min.features =  200 )

# QC by mitrochondial counts
untreated[["percent.mt"]] <- PercentageFeatureSet(untreated, pattern = "^MT-")

head(untreated@meta.data, 5)

mt_vln_qc <- VlnPlot(untreated, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

ggsave(filename = "outputs/plots/qc_mitochondial_violin.png",
       plot = mt_vln_qc, 
       width = 10,
       height = 10,
       dpi = 300
       )

mt_count_rna_plt <- FeatureScatter(untreated, feature1 = "nCount_RNA", feature2 = "percent.mt")
count_rna_feature_rna_plt <- FeatureScatter(untreated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
combined_scatter_plt_qc <- mt_count_rna_plt + count_rna_feature_rna_plt
ggsave(filename = "outputs/plots/qc_scatter_count_mt_feature.png",
       plot = combined_scatter_plt_qc,
       width = 10,
       height = 10,
       dpi = 300
       )

pre_filter_cell_count <- ncol(untreated)

untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pre_filter_cell_count - ncol(untreated)
