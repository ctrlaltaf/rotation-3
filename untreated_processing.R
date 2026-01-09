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

#pre filtering step
pre_filter_cell_count <- ncol(untreated)

untreated <- subset(untreated, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pre_filter_cell_count - ncol(untreated)

#normalizing 
untreated <- NormalizeData(untreated, normalization.method = "LogNormalize", scale.factor = 10000)

# feature selection
nrow(untreated) # number of feautures/genes
untreated <- FindVariableFeatures(untreated, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(untreated), 10)
top10_feature_plt <- VariableFeaturePlot(untreated)
top_10_feature_label <- LabelPoints(plot = top10_feature_plt, points = top10, repel = TRUE)
combined_top_feature_plt <- top10_feature_plt + top_10_feature_label
ggsave(filename = "outputs/plots/top10_features.png",
       plot = combined_top_feature_plt,
       width = 20,
       height = 10,
       dpi = 300)


# data scaling

all.genes <- rownames(untreated)
untreated <- ScaleData(untreated, features = all.genes)

# PCA
untreated <- RunPCA(untreated, features = VariableFeatures(object = untreated))
print(untreated[["pca"]], dims = 1:5, nfeatures = 5)
test_dim_loadings <- VizDimLoadings(untreated, dims = 1:2, reduction = "pca")
# for PC 1 and PC 2, which are. the genes that contribute most to these components. higher magniture means more influence
pca_dim_plot <- DimPlot(untreated, reduction = "pca") + NoLegend()
# i guess this is a sanity check in this case to see if the cells are similar (which they are since same type)
pca_dim_heatmap <- DimHeatmap(untreated, dims = 1:3, cells = 500, balanced = TRUE)

ggsave(filename = "outputs/plots/test_dim_loadings.png",
      plot = test_dim_loadings,
      width = 20,
      height = 10,
      dpi = 300)

ggsave(filename = "outputs/plots/pca_dim_plot.png",
       plot = pca_dim_plot,
       width = 20,
       height = 10,
       dpi = 300)

png("outputs/plots/pca_dim_heatmap.png", width = 20, height = 10, units = "in", res = 300)
print(DimHeatmap(untreated, dims = 1:3, cells = 500, balanced = TRUE))
dev.off()

#dimensionality check
untreated_elbow_plot <- ElbowPlot(untreated)
ggsave(filename = "outputs/plots/untreated_elbow_plot.png",
       plot = untreated_elbow_plot,
       width = 20,
      height = 10,
      dpi = 300,
      bg = "white"
)

# clustering
untreated <- FindNeighbors(untreated, dims = 1:8)
untreated <- FindClusters(untreated, resolution = 0.5)
head(Idents(untreated), 5)


#umap

untreated <- RunUMAP(untreated, dims = 1:8)

untreated_umap <- DimPlot(untreated, reduction = "umap")
ggsave(filename = "outputs/plots/untreated_umap.png",
       plot = untreated_umap,
       width = 10,
       height = 10,
       dpi = 300)
saveRDS(untreated, file = "outputs/rds/untreated.rds")

# find cluster biomarkers in cluster 2 compared to all other cluster
cluster2.markers <- FindMarkers(untreated, ident.1 = 2)
head(cluster2.markers, n = 5)

untreated.markers <- FindAllMarkers(untreated, only.pos = TRUE)
untreated.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

umap_gene_plot <- FeaturePlot(untreated, features = c("SSR1", "TSNAX", "C7orf50"))
ggsave(filename = "outputs/plots/umap_gene_plot.png",
       plot = umap_gene_plot,
       height = 10,
       width = 10,
       dpi = 300)
