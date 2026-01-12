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

rm(etop_10, etop_10.data, etop_2_5, etop_2_5.data, etop_25, etop_25.data, untreated, untreated.data)
gc()

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
combined <- subset(combined, subset = nFeature_RNA > 3500 & percent.mt < 15)
# originally filtered by nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
combined_pre_count - ncol(combined)

# normalizing data
combined <- NormalizeData(combined)
# above normalization is the default params

#feature selection
combined <- FindVariableFeatures(combined)
top10 <- head(VariableFeatures(combined), 10)

std_var_features_plt <- VariableFeaturePlot(combined)
labeled_feature_plt <- LabelPoints(plot = std_var_features_plt, points = top10, repel = TRUE)
std_var_features_plt + labeled_feature_plt

ggsave(filename = file.path(plots_dir, "combined/std_var_features_plt.png"),
       plot = std_var_features_plt + labeled_feature_plt ,
       height = 10,
       width = 15,
       dpi = 300)

gc()
# scale data
combined <- ScaleData(combined)
combined <- ScaleData(combined, features = VariableFeatures(combined)) #just scale a set of the features not all
mem.maxVSize()

#PCA
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
pca_viz_plt <- VizDimLoadings(combined, dims = 1:3, reduction = "pca")
ggsave(filename = "outputs/plots/combined/pca_viz_plt.png",
       plot = pca_viz_plt,
       height = 10,
       width = 10,
       dpi = 300)
combined_pca <- DimPlot(combined, reduction = "pca") + NoLegend()
ggsave(filename = "outputs/plots/combined/pca_dim.png",
       plot = combined_pca,
       height = 10,
       width = 10,
       dpi = 300)

png("outputs/plots/combined/pca_dim_heatmap.png", width = 20, height = 10, units = "in", res = 300)
DimHeatmap(combined, dims = 1:3, cells = 500, balanced = TRUE)
dev.off()

#find dimensionality of dataset
elbow_plot <- ElbowPlot(combined, ndims = 50)
ggsave(filename = "outputs/plots/combined/elbow.png",
       plot = elbow_plot,
       height = 10,
       width = 10,
       dpi = 300,
       bg = "white")

#calculation from ref script
pct <- combined[["pca"]]@stdev / sum(combined[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

co3 <- which(cumu > 80)[1]

co3


# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs


#cluster cells

combined <- FindNeighbors(combined, dims = 1:30) # set originally at 10.
combined <- FindClusters(combined, resolution = 0.8) #originally 0.5
head(Idents(combined), 5)

combined <- RunUMAP(combined, dims = 1:30)
combined <- RunTSNE(combined)

cluster_colors <- c(
  "untreated" = "#F37B6F",
  "etop_2.5uM" = "#7BB438",
  "etop_10uM" = "#30C1BF",
  "etop_25uM" = "#CA79FB"
)


cluster_condition_plot <- DimPlot(combined, reduction = "umap", group.by = "condition", pt.size = 0.5, label.size = 9, cols = cluster_colors) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
# we group by condition to get the dose cluster

cluster_default_plot <- DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9) +
  theme(
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
cluster_final <- cluster_condition_plot + cluster_default_plot
ggsave(filename = "outputs/plots/combined/cluster_condition.png",
       plot = cluster_condition_plot,
       height = 10,
       width = 10,
       dpi = 300)
ggsave(filename = "outputs/plots/combined/cluster_default.png",
       plot = cluster_default_plot,
       height = 8,
       width = 8,
       dpi = 300)
ggsave(filename = "outputs/plots/combined/cluster_final.png",
       plot = cluster_final,
       height = 10,
       width = 20,
       dpi = 300)

