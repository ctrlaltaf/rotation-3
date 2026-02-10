#!/usr/bun/env Rscript

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(patchwork)
library(readxl)
library(ggrepel)
library(janitor)


# look at the cluster that we are interested in
merged_seurat_filtered <- readRDS("/Users/altafbarelvi/Code/rotation-3/merged_seurat_filtered.rds")
#saveRDS(merged_seurat_filtered, file = "/Volumes/spencer_group/Rotation Students/Altaf/merged_seurat_filtered.rds")
#saveRDS(untreated_all, file = "/Volumes/spencer_group/Rotation Students/Altaf/untreated_all.rds")

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
  de_de_output_file = de_output_file,
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
  de_de_output_file = de_output_file,
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
  de_de_output_file = de_output_file,
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

de_top_feature_c_5 <- FeaturePlot(untreated_all, features = de_top_genes_c_5, min.cutoff = "q2", max.cutoff =  "q98")
de_top_feature_c_6 <- FeaturePlot(untreated_all, features = de_top_genes_c_6, min.cutoff = "q2", max.cutoff =  "q98")
de_top_feature_c_7 <- FeaturePlot(untreated_all, features = de_top_genes_c_7, min.cutoff = "q2", max.cutoff =  "q98")

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
  gsea_de_output_file = gsea_output_file,
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
  gsea_de_output_file = gsea_output_file,
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
  gsea_de_output_file = gsea_output_file,
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




# investigate p21

gene_interest <-"CDKN1A"
p21_df <- FetchData(
  object = untreated_all,
  vars   = c(gene_interest, "seurat_clusters"),
  layer   = "data"
)

p21_df %>%
  group_by(seurat_clusters) %>%
  summarise(
    n = n(),
    mean_expr   = mean(.data[[gene_interest]]),
    median_expr = median(.data[[gene_interest]]),
    pct_expr    = mean(.data[[gene_interest]] > 0) * 100
  ) %>%
  arrange(as.integer(as.character(seurat_clusters)))

expr_by_cluster <- split(p21_df[[gene_interest]], p21_df$seurat_clusters)

expr_by_cluster[["3"]][1:10]


avg <- AverageExpression(
  object   = p21_df,
  assays   = DefaultAssay(p21_df),
  features = gene_interest,
  slot     = "data",
  group.by = "seurat_clusters"
)

# look at the top and bottom genes and their feature plot


top_n <- 100                 # how many "top" genes per cluster
bottom_n <- 100              # how many "bottom" genes per cluster
padj_cutoff <- 0.05         # used only if an adjusted p-value column exists
require_significant <- TRUE # if FALSE, will ignore padj cutoff
# -----------------------

# Helper: find first matching column from a set of candidates
pick_col <- function(df, candidates) {
  hits <- intersect(names(df), candidates)
  if (length(hits) == 0) return(NA_character_)
  hits[[1]]
}

# Common column name variants (after clean_names())
gene_candidates <- c("gene", "genes", "feature", "features", "symbol", "gene_symbol")
logfc_candidates <- c("avg_log2fc", "avg_logfc", "log2fc", "logfc", "avg_log_2fc", "avg_log_2_fc")
padj_candidates <- c("p_val_adj", "p_adj", "padj", "p_val_adj_x", "p_val_adj_y", "p_adjust", "pvalue_adj", "p_val_adj_1")

# Read sheet names (each should correspond to a cluster)
sheets <- readxl::excel_sheets(de_output_file)
if (length(sheets) == 0) stop("No sheets found in workbook: ", de_output_file)

message("Found ", length(sheets), " sheet(s): ", paste(sheets, collapse = ", "))

# Containers
cluster_outputs <- list()
all_top <- list()
all_bottom <- list()

for (sh in sheets) {
  message("\n--- Processing sheet: ", sh, " ---")
  
  df_raw <- tryCatch(
    readxl::read_xlsx(de_output_file, sheet = sh),
    error = function(e) {
      message("FAILED reading sheet: ", sh)
      message("Reason: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(df_raw)) next
  if (nrow(df_raw) == 0) {
    message("  (empty sheet; skipping)")
    next
  }
  
  df_raw <- readxl::read_xlsx(de_output_file, sheet = sh)
  if (nrow(df_raw) == 0) {
    message("  (empty sheet; skipping)")
    next
  }
  
  # Clean names and keep a copy of original names if you want to inspect
  df <- janitor::clean_names(df_raw)
  
  # Some Seurat exports store gene names as rownames; try to recover if needed
  # If there is no obvious gene column, we’ll use first column if it looks like gene symbols.
  gene_col <- pick_col(df, gene_candidates)
  
  if (is.na(gene_col)) {
    # fallback: if first column is character-ish, treat it as gene
    first_col <- names(df)[1]
    if (is.character(df[[first_col]]) || is.factor(df[[first_col]])) {
      gene_col <- first_col
      message("  Using first column as gene column: ", gene_col)
    } else {
      stop("Could not identify gene column in sheet '", sh, "'. Columns: ",
           paste(names(df), collapse = ", "))
    }
  }
  
  logfc_col <- pick_col(df, logfc_candidates)
  if (is.na(logfc_col)) {
    stop("Could not identify logFC/log2FC column in sheet '", sh, "'. Columns: ",
         paste(names(df), collapse = ", "))
  }
  
  padj_col <- pick_col(df, padj_candidates)
  
  # Basic cleaning: drop rows with missing gene or logFC
  df2 <- df %>%
    mutate(
      gene = as.character(.data[[gene_col]]),
      logfc = suppressWarnings(as.numeric(.data[[logfc_col]]))
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(logfc))
  
  # Apply significance filter if possible/desired
  if (!is.na(padj_col)) {
    df2 <- df2 %>%
      mutate(p_adj = suppressWarnings(as.numeric(.data[[padj_col]])))
    if (require_significant) {
      df2 <- df2 %>% filter(!is.na(p_adj), p_adj <= padj_cutoff)
      message("  Using adjusted p-value filter: ", padj_col, " <= ", padj_cutoff)
    } else {
      message("  Adjusted p-value column found (", padj_col, "), but filter disabled.")
    }
  } else {
    message("  No adjusted p-value column found; ranking purely by logFC.")
  }
  
  if (nrow(df2) == 0) {
    message("  No rows left after filtering; skipping output for this sheet.")
    next
  }
  
  # Top upregulated and bottom downregulated
  top_tbl <- df2 %>%
    arrange(desc(logfc)) %>%
    slice_head(n = top_n) %>%
    mutate(cluster = sh, direction = "top_up")
  
  bottom_tbl <- df2 %>%
    arrange(logfc) %>%
    slice_head(n = bottom_n) %>%
    mutate(cluster = sh, direction = "bottom_down")
  
  # Save combined per-cluster table (keep original columns too)
  # We'll bind top and bottom, but also keep gene/logfc/p_adj if present.
  keep_cols <- c("cluster", "direction", "gene", "logfc")
  if ("p_adj" %in% names(top_tbl)) keep_cols <- c(keep_cols, "p_adj")
  
  per_cluster <- bind_rows(top_tbl, bottom_tbl) %>%
    dplyr::select(any_of(keep_cols), everything())
  
  cluster_outputs[[sh]] <- per_cluster
  
  all_top[[sh]] <- top_tbl
  all_bottom[[sh]] <- bottom_tbl
  
  message("  Kept ", nrow(per_cluster), " rows (top ", min(top_n, nrow(top_tbl)),
          " + bottom ", min(bottom_n, nrow(bottom_tbl)), ")")
}


top_df <- bind_rows(all_top) %>% dplyr::select(cluster, gene, logfc, any_of("p_adj"), everything())
bottom_df <- bind_rows(all_bottom) %>% dplyr::select(cluster, gene, logfc, any_of("p_adj"), everything())


clusters <- c("Cluster_5", "Cluster_6", "Cluster_7")


for (c in clusters) {
  
  top_cluster_subset <- top_df %>%
    filter(cluster == c)
  
  bottom_cluster_subset <- bottom_df %>%
    filter(cluster == c)
  
  top_feature <- FeaturePlot(untreated_all, features = top_cluster_subset$gene,  min.cutoff = "q2", max.cutoff =  "q98")
  top_output_name <- paste0("outputs/plots/subset/", c, "_top_gene_features.pdf")
  
  ggsave(
    filename = top_output_name,
    plot     = top_feature,
    width    = 20,
    height   = 100,
    units    = "in",
    useDingbats = FALSE,
    limitsize = FALSE
  )
  
  bottom_feature <- FeaturePlot(untreated_all, features = bottom_cluster_subset$gene, min.cutoff = "q2", max.cutoff =  "q98")
  bottom_output_name <- paste0("outputs/plots/subset/", c, "_bottom_gene_features.pdf")
  
  ggsave(
    filename = bottom_output_name,
    plot     = bottom_feature,
    width    = 20,
    height   = 100,
    units    = "in",
    useDingbats = FALSE,
    limitsize = FALSE
  )
}

# combine potential marker genes for each cluster

c_5_markers <- c("IGFL2-AS1","DKK1", "S100P", "DCN","CFB", "FAT2", "CLCA2", "CDKN1C", "S100A4", "MUC1", "CACNA2D3")

c_6_markers <- c("PYCARD", "KRT14", "SOX5","TGM1", "ESRP1")

c_7_markers <- c("CHI3L1", "FEZ1", "MME", "TACC2")

c_6_7_markers <- c("PARVB", "BST2", "GLDC")

all_markers <- c(
  c_5_markers,
  c_6_markers,
  c_7_markers,
  c_6_7_markers
)

potential_markers <- FeaturePlot(untreated_all, features = all_markers, min.cutoff = "q2", max.cutoff =  "q98")
ggsave(
  filename = "outputs/plots/subset/potential_markers_feature_plt.pdf",
  plot     = potential_markers,
  width    = 30,
  height   = 30,
  units    = "in",
  useDingbats = FALSE,
  limitsize = FALSE
)

print(all_markers)

p53_family_features <- FeaturePlot(untreated_all, features = c("TP53", "TP73", "TP63"))

ggsave(filename = "outputs/plots/subset/p53_family_feature.png",
       plot = p53_family_features,
       width = 14,
       height = 10,
       dpi = 300)

gcn2geneslist <- c("CHAC1", "SLC3A2", "HSPA5", "DDIT3", "XBP1", "ATF3", "ATF4", "ATF6", "CEBPA", "CEBPB", "CEBPD",
                   "CEBPG", "CEBPZ", "NFE2L2", "JUN", "FOS", "TRIB3", "ASNS", "SESN2", "SLC7A5", "PSAT1",
                   "PHGDH", "VLDLR", "GCN1", "EIF2AK4", "EIF2S1", "PPP1R15A")

gcn2_list_features <- FeaturePlot(untreated_all, features = gcn2geneslist,  min.cutoff = "q2", max.cutoff =  "q98")

ggsave(filename = "outputs/plots/subset/gcn2_list_feature.pdf",
       plot = gcn2_list_features,
       width    = 30,
       height   = 40,
       units    = "in",
       limitsize = FALSE,
       useDingbats = FALSE
       )

# pairwise clustering analysis 

ut_subset_clusters <- c("5", "6", "7")
ut_main_clusters <- c("0", "1")

for (i in ut_main_clusters){
  for (j in ut_subset_clusters) {
    output_de_name <- paste("cluster_", i,"-",j,"_pair.xlsx", sep = "")
    de_pairwise_path <- run_de_pairwise_findmarkers(
      seurat_obj = untreated_all,
      group_by = "seurat_clusters",           # metadata column that defines groups
      group_1 = j,            # label in obj@meta.data$Treatment
      group_2 = i,                  # label in obj@meta.data$Treatment
      assay = "RNA",
      slot = "data",
      test_use = "wilcox",
      min_pct = 0.1,
      logfc_threshold = 0.25,
      only_pos = FALSE,
      out_de_xlsx = output_de_name,
      out_dir = "outputs/gsea/pairwise_de"
    )
    
    output_gsea_name <- paste("cluster_", i,"-",j,"_pair.xlsx", sep = "")
    gsea_out <- gsea_from_de_workbook_fgsea(
      de_xlsx_path = paste("outputs/gsea/pairwise_de/", output_de_name, sep = ""),
      out_gsea_xlsx = output_gsea_name,
      out_dir = "outputs/gsea/pairwise_gsea",
      species = "Homo sapiens",
      seed = 1
    )
    
    gsea_plots <- plot_gsea_dotplot_from_workbook(
      gsea_xlsx_path = paste("outputs/gsea/pairwise_gsea/", output_gsea_name, sep = ""),
      sheet = paste(j, "_vs_", i, sep = ""),
      which = c("GO_BP", "Hallmark"),
      top_n_each = 15,
      title_prefix = paste("Cluster", j, "vs", i)
    )
    
    go_bp_output <- paste("outputs/plots/subset/gsea/go_bp_cluster_", i, "-", j,".png", sep = "")
    hallmark_output <- paste("outputs/plots/subset/gsea/hallmark_cluster_", i, "-", j,".png", sep = "")
    
    ggsave(go_bp_output, gsea_plots$GO_BP, width = 20, height = 8, dpi = 300)
    ggsave(hallmark_output, gsea_plots$Hallmark, width = 20, height = 8, dpi = 300)
    
  }
}

# merge cluster 0 and 1 together

untreated_01_merged <- untreated_all

old <- as.integer(as.character(untreated_01_merged$seurat_clusters))

new <- ifelse(old %in% c(0, 1), 0L, old - 1L)

# store as factor with clean levels 0..K
new <- factor(new, levels = sort(unique(new)))

# save into metadata (so it’s persistent in the new object)
untreated_01_merged$seurat_clusters <- new

# set active identities to these clusters
Idents(untreated_01_merged) <- "seurat_clusters"

# sanity check
table(old = untreated_all$seurat_clusters, new = untreated_01_merged$seurat_clusters)

DimPlot(untreated_01_merged, reduction = 'umap', pt.size = .01, shuffle = TRUE, label = TRUE, label.size = 9) + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(), axis.line = element_blank()) + 
  xlim(-8, 8) + ylim(-8, 8) +
  labs(title = NULL, x = NULL, y = NULL)


ut_subset_clusters <- c("4", "5", "6")
ut_main_clusters <- c("0")

for (i in ut_main_clusters){
  for (j in ut_subset_clusters) {
    output_name <- paste("cluster_", i,"-",j,"_pair_merged_main_ut.xlsx", sep = "")
    de_pairwise_path <- run_de_pairwise_findmarkers(
      seurat_obj = untreated_01_merged,
      group_by = "seurat_clusters",           # metadata column that defines groups
      group_1 = j,            # label in obj@meta.data$Treatment
      group_2 = i,                  # label in obj@meta.data$Treatment
      assay = "RNA",
      slot = "data",
      test_use = "wilcox",
      min_pct = 0.1,
      logfc_threshold = 0.25,
      only_pos = FALSE,
      out_de_xlsx = output_name,
      out_dir = "outputs/gsea/pairwise_de/"
    )
    
    gsea_output_name <- paste("cluster_", i, "-", j, "_pair_merged_main_GSEA.xlsx", sep = "")
    gsea_out <- gsea_from_de_workbook_fgsea(
      de_xlsx_path = paste("outputs/gsea/pairwise_de/", output_name, sep = ""),
      out_gsea_xlsx = gsea_output_name,
      out_dir = "outputs/gsea/pairwise_gsea",
      species = "Homo sapiens",
      seed = 1
    )
    
    gsea_plots <- plot_gsea_dotplot_from_workbook(
      gsea_xlsx_path = paste("outputs/gsea/pairwise_gsea/", gsea_output_name, sep = ""),
      sheet = paste(j, "_vs_", i, sep = ""),
      which = c("GO_BP", "Hallmark"),
      top_n_each = 15,
      title_prefix = paste("Cluster", j, "vs", i)
    )
    
    go_bp_output <- paste("outputs/plots/subset/gsea/go_bp_cluster_", i, "-", j,"merged_main",".png", sep = "")
    hallmark_output <- paste("outputs/plots/subset/gsea/hallmark_cluster_", i, "-", j,"merged_main",".png", sep = "")
    
    ggsave(go_bp_output, gsea_plots$GO_BP, width = 20, height = 8, dpi = 300)
    ggsave(hallmark_output, gsea_plots$Hallmark, width = 20, height = 8, dpi = 300)
    
  }
}



