#!/usr/bun/env Rscript

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(patchwork)
library(readxl)
library(ggrepel)


# look at the cluster that we are interested in
merged_seurat_filtered <- readRDS("/Users/altafbarelvi/Code/rotation-3/merged_seurat_filtered.rds")
# saveRDS(merged_seurat_filtered, file = "/Volumes/spencer_group/Rotation Students/Altaf/merged_seurat_filtered.rds")

colnames(merged_seurat_filtered@meta.data)

# checking to see if we have the same umaps
original_cluster_plot <- DimPlot(merged_seurat_filtered, reduction = 'umap', pt.size = .01, shuffle = TRUE, label = TRUE, label.size = 9) + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), axis.line = element_blank()) + 
  xlim(-8, 8) + ylim(-8, 8) +
  labs(title = NULL, x = NULL, y = NULL)
original_cluster_plot
ggsave(filename = "outputs/plots/subset/original_cluster.png",
       plot = original_cluster_plot,
       height = 10,
       width = 10,
       dpi = 300
       )

original_cluster_phase_plot <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = "Phase", pt.size = .01, shuffle = TRUE, label = FALSE, label.size = 9) + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), axis.line = element_blank()) + 
  xlim(-8, 8) + ylim(-8, 8) +
  labs(title = NULL, x = NULL, y = NULL)
original_cluster_phase_plot
ggsave(filename = "outputs/plots/subset/original_cluster_phase_plot.png",
       plot = original_cluster_phase_plot,
       height = 10,
       width = 10,
       dpi = 300
)

combined_original_cluster <- original_cluster_plot + original_cluster_phase_plot
ggsave(filename = "outputs/plots/subset/combined_original_cluster.png",
       plot = combined_original_cluster,
       width = 20,
       height = 10,
       dpi = 300)

#subset for only the untreated cells
unique(merged_seurat_filtered$Dose)
untreated_all <- subset(merged_seurat_filtered, subset = Dose == "UT")

#get barcodes for original ut 13 cells
orig_ut13_cells <- WhichCells(
  merged_seurat_filtered,
  expression = Dose == "UT" & seurat_clusters == 13
)

#get barcodes for original ut 15 cells
orig_ut15_cells <- WhichCells(
  merged_seurat_filtered,
  expression = Dose == "UT" & seurat_clusters == 15
)


#QC
vln_rna_mt_plt <- VlnPlot(untreated_all, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
vln_rna_mt_plt

FeatureScatter(untreated_all, feature1 = "nCount_RNA", feature2 = "mitoPercent") +
  ggtitle("nCount vs mitoPercent (UT only)")

FeatureScatter(untreated_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("nCount vs nFeature (UT only)")

FeaturePlot(untreated_all, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)


# analyze sillouette grid results from cluster_stability_altaf.R

res <- read.csv("outputs/cluster_stability/silhouette_grid_results.csv") %>%
  mutate(
    pcs = as.integer(pcs),
    resolution = as.numeric(resolution),
    mean_silhouette = as.numeric(mean_silhouette)
  )

# 5) Plot heatmap (easy to see best region + drops)
plot_df <- res %>%
  mutate(
    pcs = factor(pcs, levels = sort(unique(pcs))),
    resolution = factor(resolution, levels = sort(unique(resolution)))
  )

p_heat <- ggplot(plot_df, aes(x = resolution, y = pcs, fill = mean_silhouette)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", mean_silhouette)), size = 3) +
  labs(
    title = "Mean silhouette across (PCs x resolution)",
    x = "Resolution",
    y = "# PCs (dims)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(filename = "outputs/plots/subset/mean_silhouette.png",
       plot = p_heat,
       width = 10,
       height = 8,
       dpi = 300)

# process
untreated_all <- FindVariableFeatures(untreated_all)
untreated_all <- ScaleData(untreated_all, features = VariableFeatures(untreated_all))
untreated_all <- RunPCA(untreated_all, features = VariableFeatures(untreated_all))

ElbowPlot(untreated_all, ndims = 50)

# sweep through pc and res params (ONLY DO THIS IF YOU NEED)

source("/Users/altafbarelvi/Desktop/DE_GSEA_functions.R")

pc_params  <- c(15, 20, 25, 30)
res_params <- c(0.2, 0.4, 0.6)

# Output roots
out_root <- "/Users/altafbarelvi/Code/rotation-3/outputs/gsea/"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

out_de_dir   <- file.path(out_root, "DE")
out_gsea_dir <- file.path(out_root, "GSEA")
dir.create(out_de_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Plot output
plot_dir <- "/Users/altafbarelvi/Code/rotation-3/outputs/plots/subset/param_sweep"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# Summary table
sweep_results <- data.frame(
  pcs = integer(),
  resolution = numeric(),
  ut13_overlap = integer(),
  n_clusters = integer(),
  de_xlsx = character(),
  stringsAsFactors = FALSE
)

# One long PDF of all UMAP outputs
pdf(file.path(plot_dir, "UMAP_parameter_sweep_all.pdf"), width = 18, height = 12)

for (pc in pc_params) {
  for (res in res_params) {
    
    message("Running: PCs = ", pc, " | Resolution = ", res)
    obj <- untreated_all
    
    # Cluster/UMAP
    obj <- FindNeighbors(obj, dims = 1:pc, verbose = FALSE)
    obj <- FindClusters(obj, resolution = res, verbose = FALSE)
    obj <- RunUMAP(obj, dims = 1:pc, verbose = FALSE)
    
    # UT13 overlap + cluster count
    ut13_in_new <- intersect(orig_ut13_cells, colnames(obj))
    ut13_count <- length(ut13_in_new)
    n_clust <- length(unique(Idents(obj)))
    
    # UT15 overlap + cluster count
    ut15_in_new <- intersect(orig_ut15_cells, colnames(obj))
    ut15_count <- length(ut15_in_new)
    n_clust <- length(unique(Idents(obj)))
    
    p_clusters <- DimPlot(
      obj, reduction = "umap",
      group.by = "seurat_clusters",
      label = TRUE, label.size = 6,
      pt.size = 0.4
    ) + ggtitle(paste0("PC=", pc, " | res=", res))
    
    p_phase <- DimPlot(
      obj, reduction = "umap",
      group.by = "Phase",
      pt.size = 0.3,
      shuffle = TRUE
    ) + ggtitle("Cell Cycle Phase")
    
    p_highlight_13 <- DimPlot(
      obj, reduction = "umap",
      cells.highlight = ut13_in_new,
      cols = "lightgrey",
      cols.highlight = "red",
      pt.size = 0.4
    ) + ggtitle("Original UT13 highlighted")
    
    p_highlight_15 <- DimPlot(
      obj, reduction = "umap",
      cells.highlight = ut15_in_new,
      cols = "lightgrey",
      cols.highlight = "red",
      pt.size = 0.4
    ) + ggtitle("Original UT15 highlighted")

    
    combined_untreated_umap <- p_clusters + p_phase + p_highlight_13 + p_highlight_15
      plot_annotation(
        title = paste0(
          "PCs=", pc,
          " | res=", res,
          " | UT13 overlap=", ut13_count,
          " | #Clusters=", n_clust
        )
      )

    print(combined_untreated_umap)
    ggsave(
      filename = file.path(plot_dir, paste0("umap_PC", pc, "_res", res, ".png")),
      plot = combined_untreated_umap,
      width = 16, height = 12, dpi = 300
    )
    

    cfg_tag <- paste0("pc", pc, "_res", res)
    cfg_de_dir <- file.path(out_de_dir, cfg_tag)
    dir.create(cfg_de_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Run DE one-vs-all for this config

    de_path <- run_de_one_vs_all_findallmarkers(
      seurat_obj = obj,
      group_by = "seurat_clusters",
      assay = "RNA",
      slot = "data",
      test_use = "wilcox",
      min_pct = 0.1,
      logfc_threshold = 0.25,
      only_pos = FALSE,
      out_de_xlsx = "Clusters_OneVsAll_DE.xlsx",
      out_dir = cfg_de_dir
    )
    
    # Store summary row
    sweep_results <- rbind(
      sweep_results,
      data.frame(
        pcs = pc,
        resolution = res,
        ut13_overlap = ut13_count,
        n_clusters = n_clust,
        de_xlsx = de_path,
        stringsAsFactors = FALSE
      )
    )
  }
}

dev.off()

# Save summary table
write.csv(sweep_results, file.path(plot_dir, "parameter_sweep_summary.csv"), row.names = FALSE)

sweep_results



# highlight where the untreated cluster outside of the main untreated clusters are
unique(merged_seurat_filtered$Dose)
original_cluster_plot
top_island_clusters <- c(12, 9, 4, 7, 3, 2, 10, 1, 8, 5)

Idents(merged_seurat_filtered) <- "seurat_clusters"
untreated_cells <- WhichCells(merged_seurat_filtered, expression = Dose == "UT")
top_island_cells <- WhichCells(merged_seurat_filtered, idents = top_island_clusters)
untreated_top_island_cells <- intersect(untreated_cells, top_island_cells)

length(untreated_top_island_cells)
head(untreated_top_island_cells)

untreated_outside_main_island <- DimPlot(
  merged_seurat_filtered,
  reduction = "umap",
  cells.highlight = untreated_top_island_cells,
  cols = c("lightgrey", "#F37B6F"),
  pt.size = 0.01,
  shuffle = TRUE
) +
  ggtitle("Untreated cells that are outside of main untreated island") +
  theme(
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.line  = element_blank(),
  )

ggsave(filename = "outputs/plots/subset/untreated_outside_main_island.png",
       plot = untreated_outside_main_island,
       width = 10,
       height = 10,
       dpi = 300)

# use umap dim = 25 and res = 0.4

untreated_all <- FindNeighbors(untreated_all, dims = 1:25, verbose = FALSE)
untreated_all <- FindClusters(untreated_all, resolution = 0.4, verbose = FALSE)
untreated_all <- RunUMAP(untreated_all, dims = 1:25, verbose = FALSE)

# UT13 overlap + cluster count
ut13_in_new <- intersect(orig_ut13_cells, colnames(untreated_all))
ut13_count <- length(ut13_in_new)
n_clust <- length(unique(Idents(untreated_all)))

# UT15 overlap + cluster count
ut15_in_new <- intersect(orig_ut15_cells, colnames(untreated_all))
ut15_count <- length(ut15_in_new)
n_clust <- length(unique(Idents(untreated_all)))

p_clusters <- DimPlot(
  untreated_all, reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE, label.size = 6,
  pt.size = 0.4
) + ggtitle(paste0("PC=", 25, " | res=", 0.4))

p_phase <- DimPlot(
  untreated_all, reduction = "umap",
  group.by = "Phase",
  pt.size = 0.3,
  shuffle = TRUE
) + ggtitle("Cell Cycle Phase")

p_highlight_13 <- DimPlot(
  untreated_all, reduction = "umap",
  cells.highlight = ut13_in_new,
  cols = "lightgrey",
  cols.highlight = "red",
  pt.size = 0.4
) + ggtitle("Original UT13 highlighted")

p_highlight_15 <- DimPlot(
  untreated_all, reduction = "umap",
  cells.highlight = ut15_in_new,
  cols = "lightgrey",
  cols.highlight = "red",
  pt.size = 0.4
) + ggtitle("Original UT15 highlighted")


combined_untreated_umap <- p_clusters + p_phase + p_highlight_13 + p_highlight_15
plot_annotation(
  title = paste0(
    "PCs=", 25,
    " | res=", 0.4,
    " | UT13 overlap=", ut13_count,
    " | #Clusters=", n_clust
  )
)
ggsave(filename = "outputs/plots/subset/combined_chosen_untreated_umap.png",
       plot = combined_untreated_umap,
       height = 10, 
       width = 14, 
       dpi = 300)

# do GSEA

de_output_file <- "/Users/altafbarelvi/Code/rotation-3/outputs/gsea/de/pc25_res0.4/Clusters_OneVsAll_DE.xlsx"
gsea_output_file <- "/Users/altafbarelvi/Code/rotation-3/outputs/gsea/GSEA/Clusters_OneVsAll_GSEA.xlsx"
file.exists(gsea_output_file)

# make volcano plots for each cluster of interest

p_volcano_5 <- plot_volcano_from_de_workbook(
  de_xlsx_path = de_output_file,
  cluster_id   = 5,
  sheet_prefix = "Cluster_",
  padj_thresh  = 0.05,
  lfc_thresh   = 1.5,
  label_top_n  = 10,
  title_prefix = "Volcano"
)
de_top_genes_c_5 <- c("IGFL2-AS1","AC004231.1", "AL139042.1","KCNK15", "FPR1", 
                      "S100P", "PTGFR", "AC092810.3", "ADH1B", "DCN", "CAMK2B", "CLCA2",
                      "LINC01503", "PLAT", "MT2A", "IRF6", "S100A2", "SERINC2", "FAT2", 
                      "IGFBP7")

ggsave(filename = "outputs/plots/subset/volcano_5.png",
       plot = p_volcano_5,
       width = 10,
       height = 10,
       dpi = 300)



p_volcano_6 <- plot_volcano_from_de_workbook(
  de_xlsx_path = de_output_file,
  cluster_id   = 6,
  sheet_prefix = "Cluster_",
  padj_thresh  = 0.05,
  lfc_thresh   = 1.5,
  label_top_n  = 10,
  title_prefix = "Volcano"
)
de_top_genes_c_6 <- c("PYCARD", "LINC01127", "KRT14", "PARVB", "SNHG18", "IL1R2",
                      "TENM1", "TMEM216", "SOX5", "SYT14", "IGFBP7", "BST2", "FUCA2",
                      "ESRP1", "SMOC1", "IGFBP2", "HLA-B", "ITGA10", "LINC02289", "PRSS21")

ggsave(filename = "outputs/plots/subset/volcano_6.png",
       plot = p_volcano_6,
       width = 10,
       height = 10,
       dpi = 300)


p_volcano_7 <- plot_volcano_from_de_workbook(
  de_xlsx_path = de_output_file,
  cluster_id   = 7,
  sheet_prefix = "Cluster_",
  padj_thresh  = 0.05,
  lfc_thresh   = 1.5,
  label_top_n  = 10,
  title_prefix = "Volcano"
)

de_top_genes_c_7 <- c("WFDC2", "AC090204.1", "NUDT11", "EMCN", "RYR1", "PARVB", 
                      "RRAGD", "LCP1", "C10orf67", "CHI3L1", "BST2", "CCND2", "PRSS21", "TACC2", 
                      "LIMCH1", "COL8A1", "CCBE1", "CENPF", "FBLN1", "TBL1X")

ggsave(filename = "outputs/plots/subset/volcano_7.png",
       plot = p_volcano_7,
       width = 10,
       height = 10,
       dpi = 300)



# feature plot the top genes

de_top_feature_c_5 <- FeaturePlot(untreated_all, features = de_top_genes_c_5)
de_top_feature_c_6 <- FeaturePlot(untreated_all, features = de_top_genes_c_6)
de_top_feature_c_7 <- FeaturePlot(untreated_all, features = de_top_genes_c_7)

ggsave(filename = "outputs/plots/subset/de_top_feature_c_5.png",
       plot = de_top_feature_c_5,
       height = 20,
       width = 20,
       dpi = 300
       )
ggsave(filename = "outputs/plots/subset/de_top_feature_c_6.png",
       plot = de_top_feature_c_6,
       height = 20,
       width = 20,
       dpi = 300
)
ggsave(filename = "outputs/plots/subset/de_top_feature_c_7.png",
       plot = de_top_feature_c_7,
       height = 20,
       width = 20,
       dpi = 300
)

gsea_plots_c_5 <- plot_gsea_dotplot_from_workbook(
  gsea_xlsx_path = gsea_output_file,
  cluster_id     = 5,
  sheet_prefix   = "Cluster_",
  which          = c("GO_BP", "Hallmark"),
  top_n          = 20,
  title_prefix   = "GSEA dot plot"
)

print(gsea_plots_c_5$GO_BP)
ggsave(filename = "outputs/plots/subset/gsea_c_5.png",
       plot = gsea_plots_c_5$GO_BP,
       width = 10,
       height = 10,
       dpi = 300)


gsea_plots_c_6 <- plot_gsea_dotplot_from_workbook(
  gsea_xlsx_path = gsea_output_file,
  cluster_id     = 6,
  sheet_prefix   = "Cluster_",
  which          = c("GO_BP", "Hallmark"),
  top_n          = 20,
  title_prefix   = "GSEA dot plot"
)

print(gsea_plots_c_6$GO_BP)
ggsave(filename = "outputs/plots/subset/gsea_c_6.png",
       plot = gsea_plots_c_6$GO_BP,
       width = 10,
       height = 10,
       dpi = 300)

gsea_plots_c_7 <- plot_gsea_dotplot_from_workbook(
  gsea_xlsx_path = gsea_output_file,
  cluster_id     = 7,
  sheet_prefix   = "Cluster_",
  which          = c("GO_BP", "Hallmark"),
  top_n          = 20,
  title_prefix   = "GSEA dot plot"
)

print(gsea_plots_c_7$GO_BP)
ggsave(filename = "outputs/plots/subset/gsea_c_7.png",
       plot = gsea_plots_c_7$GO_BP,
       width = 10,
       height = 10,
       dpi = 300)

