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

untreated_h5_file <- file.path(data_dir, "Untreated_FullDepth/filtered_feature_bc_matrix.h5")
untreated.data <- Read10X_h5(untreated_h5_file)
untreated <- CreateSeuratObject(counts = untreated.data, project = "untreated", 
                                assay = "RNA", min.cells =  3, min.features =  200 )

etop_2_5_h5_file <- file.path(data_dir, "2_5uM_Etop_FullDepth/filtered_feature_bc_matrix.h5")
etop_2_5.data <- Read10X_h5(etop_2_5_h5_file)
etop_2_5 <- CreateSeuratObject(counts = etop_2_5.data, project = "etop_2_5", 
                                assay = "RNA", min.cells =  3, min.features =  200 )

etop_10_h5_file <- file.path(data_dir, "10uM_Etop_FullDepth/filtered_feature_bc_matrix.h5")
etop_10.data <- Read10X_h5(etop_10_h5_file)
etop_10 <- CreateSeuratObject(counts = etop_10.data, project = "etop_10", 
                               assay = "RNA", min.cells =  3, min.features =  200 )

etop_25_h5_file <- file.path(data_dir, "25uM_Etop_FullDepth/filtered_feature_bc_matrix.h5")
etop_25.data <- Read10X_h5(etop_25_h5_file)
etop_25 <- CreateSeuratObject(counts = etop_25.data, project = "etop_25", 
                              assay = "RNA", min.cells =  3, min.features =  200 )

untreated$condition <- "untreated"
etop_2_5$condition  <- "etop_2.5uM"
etop_10$condition   <- "etop_10uM"
etop_25$condition   <- "etop_25uM"


combined <- merge(
  x = untreated,
  y = list(etop_2_5, etop_10, etop_25),
  add.cell.ids = c("UNT", "E2.5", "E10", "E25"),
  project = "etoposide_full_depth"
)

combined
table(combined$condition)

#mitochondrial qc
combined[["percent.mt"]] <- PercentageFeatureSet(combined, patter = "^MT-")
vln_rna_mt_plt <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = file.path(plots_dir, "combined/vln_rna_mt.png"),
       plot = vln_rna_mt_plt,
       height = 20,
       width = 20,
       dpi = 300)

count_rna_mt_plt <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
count_rna_feature_plt <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
count_rna_mt_plt + count_rna_feature_plt 

ggsave(filename = file.path(plots_dir, "combined/combined_feature_rna_mt.png"),
       plot = count_rna_mt_plt + count_rna_feature_plt ,
       height = 15,
       width = 15,
       dpi = 300)

# filter out
combined_pre_count <- ncol(combined)
combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
combined_pre_count - ncol(combined)

# normalizing data
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

#feature selection
combined <- FindVariableFeatures(combined, selection.method = "vst", nFeatures = 2000)
top10 <- head(VariableFeatures(combined), 10)

std_var_features_plt <- VariableFeaturePlot(combined)
labeled_feature_plt <- LabelPoints(plot = std_var_features_plt, points = top10, repel = TRUE)
std_var_features_plt + labeled_feature_plt