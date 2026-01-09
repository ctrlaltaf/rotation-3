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


