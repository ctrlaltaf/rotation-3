############################################################
## clustree cluster stability analysis (RNA assay)
## Resolutions: 0.5 → 1.0
############################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(clustree)
  library(ggplot2)
  library(patchwork)
})

# palbo <- readRDS('~/Desktop/scRNAseq/palbo_analysis/palbo_RERUN/datasets/palbo_base_allcells_postQC-DimRed.rds')
#untreated_all <- readRDS("/Volumes/spencer_group/Rotation Students/Altaf/untreated_all.rds")

colnames(untreated_all@meta.data)

stopifnot(exists("untreated_all"))
obj <- untreated_all
#obj <- etop_UT
############################################################
## save_plot helper (provided)
############################################################
save_plot <- function(plot, filename, path = ".", height = 6, width = 8,
                      resolution = 300, format = "pdf") {
  file_path <- file.path(path, paste0(filename, ".", format))
  
  if (format == "png") {
    png(file = file_path, height = height, width = width,
        units = "in", res = resolution)
    print(plot)
    dev.off()
    
  } else if (format == "pdf") {
    pdf(file = file_path,
        height = height,
        width = width,
        useDingbats = FALSE,
        compress = FALSE)
    print(plot)
    dev.off()
    message("Vectorized PDF saved for Illustrator editing: ", file_path)
    
  } else {
    stop("Unsupported format. Use 'png' or 'pdf'.")
  }
  
  message("Plot saved as ", file_path)
}

############################################################
## Output directory
############################################################
# fig_dir <- "~/Desktop/scRNAseq/untreated_all_analysis/untreated_all_RERUN/untreated_all_figures_REMAKE/cluster_stability"
fig_dir <- "/Users/altafbarelvi/Code/rotation-3/outputs/cluster_stability"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## Parameters
############################################################
DefaultAssay(obj) <- "RNA"
dims_use <- 1:30
res_grid <- seq(0.5, 1.0, by = 0.1)
prefix   <- "RNA_snn_res."

############################################################
## Build graph once (RNA-based)
############################################################
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)

############################################################
## Run clustering across resolutions
############################################################
for (res in res_grid) {
  message("Running FindClusters at resolution = ", res)
  obj <- FindClusters(obj, resolution = res, verbose = FALSE)
}

############################################################
## Confirm clustering columns exist
############################################################
clust_cols <- grep(
  paste0("^", gsub("\\.", "\\\\.", prefix)),
  colnames(obj@meta.data),
  value = TRUE
)
print(clust_cols)

############################################################
## clustree visualization
############################################################
p_tree <- clustree(obj@meta.data, prefix = prefix) +
  ggtitle("Cluster stability across resolutions (RNA, 0.5 → 1.0)")
p_tree
save_plot(
  plot = p_tree,
  filename = "clustree_RNA_res_0p5_to_1p0",
  path = fig_dir,
  width = 12,
  height = 8,
  format = "pdf"
)

############################################################
## Number of clusters vs resolution
############################################################
clust_cols_in_range <- clust_cols[
  as.numeric(sub(prefix, "", clust_cols)) %in% res_grid
]

nclust_df <- data.frame(
  resolution = as.numeric(sub(prefix, "", clust_cols_in_range)),
  n_clusters = sapply(
    clust_cols_in_range,
    function(cc) length(unique(obj@meta.data[[cc]]))
  )
)
nclust_df <- nclust_df[order(nclust_df$resolution), ]

p_nclust <- ggplot(nclust_df, aes(x = resolution, y = n_clusters)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Number of clusters vs resolution (RNA)",
    x = "Resolution",
    y = "Number of clusters"
  )
p_nclust

save_plot(
  plot = p_nclust,
  filename = "nClusters_vs_resolution_RNA_0p5_to_1p0",
  path = fig_dir,
  width = 7,
  height = 5,
  format = "pdf"
)

############################################################
## Optional: UMAP preview at each resolution
############################################################
# Assumes UMAP already exists; uncomment if needed:
# obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

umap_plots <- lapply(res_grid, function(res) {
  colname <- paste0(prefix, res)
  stopifnot(colname %in% colnames(obj@meta.data))
  
  DimPlot(
    obj,
    group.by = colname,
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(paste0("RNA clustering @ resolution ", res))
})

p_umap <- wrap_plots(umap_plots, ncol = 2)
p_umap
save_plot(
  plot = p_umap,
  filename = "UMAP_RNA_res_0p5_to_1p0",
  path = fig_dir,
  width = 14,
  height = 18,
  format = "pdf"
)


###############################################################################
# Goal:
#   For a Seurat object `obj`, compute mean silhouette scores across a grid of:
#     - number of PCs used for FindNeighbors (dims = 1:np)
#     - FindClusters resolution values
#   Then plot a heatmap 
#
# Notes (practical):
#   - Silhouette requires a distance matrix; for large datasets, a full NxN
#     distance matrix is too big. This script computes silhouette on a
#     random subsample of cells (default 5000) using Euclidean distances in
#     PCA space. This is typical for stability/quality sweeps.
#
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(cluster)  # silhouette()
  library(viridis)  # viridis color scales
})

# -----------------------------
# 0) Basic checks
# -----------------------------
stopifnot(exists("obj"))
stopifnot(inherits(obj, "Seurat"))

set.seed(1)

# -----------------------------
# 1) User parameters
# -----------------------------
# Grid to evaluate
pc_grid  <- c(10, 15, 20, 25, 30, 35, 40)
res_grid <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)

# Silhouette computation settings
pca_reduction <- "pca"
max_cells_for_silhouette <- 5000  # subsample cap for silhouette distance matrix
silhouette_metric <- "euclidean"  # Euclidean in PCA space

# Optional: run clustering on all cells or only on the silhouette subsample
# - TRUE  (recommended): cluster full object, then compute silhouette on subset
# - FALSE: cluster subset only (faster, but cluster labels may differ)
cluster_full_object <- TRUE

# Performance tip: reuse the same neighbor graph name across sweeps
# (FindNeighbors will overwrite this graph each time)
graph_name <- "sweep_snn"

# Plot options
plot_title <- "Resolution parameter sweeping across PCA dims and FindClusters\nSilhouette grid search"
digits_to_label <- 2

# Output (optional)
out_dir <- "/Users/altafbarelvi/Code/rotation-3/outputs/cluster_stability"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
heatmap_pdf <- file.path(out_dir, "silhouette_grid_heatmap.pdf")
results_csv <- file.path(out_dir, "silhouette_grid_results.csv")

# -----------------------------
# 2) Ensure PCA exists (minimal)
# -----------------------------
if (!(pca_reduction %in% Reductions(obj))) {
  message("PCA reduction not found. Running minimal preprocessing + PCA...")
  if (!("RNA" %in% Assays(obj))) stop("Expected an RNA assay. Please set DefaultAssay(obj) appropriately.")
  DefaultAssay(obj) <- "RNA"
  
  # Minimal pipeline: if you already have these steps done, you can remove this block
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
}

# Determine max PCs available
pca_embed_all <- Embeddings(obj, reduction = pca_reduction)
max_pcs_available <- ncol(pca_embed_all)
if (max(pc_grid) > max_pcs_available) {
  stop(sprintf(
    "pc_grid requests up to %d PCs, but %s only has %d PCs. Re-run RunPCA(npcs=...) or adjust pc_grid.",
    max(pc_grid), pca_reduction, max_pcs_available
  ))
}

# -----------------------------
# 3) Choose silhouette cells (subsample)
# -----------------------------
all_cells <- colnames(obj)
n_all <- length(all_cells)

if (n_all <= max_cells_for_silhouette) {
  sil_cells <- all_cells
} else {
  sil_cells <- sample(all_cells, size = max_cells_for_silhouette, replace = FALSE)
}
message(sprintf("Silhouette will be computed on %d / %d cells.", length(sil_cells), n_all))

# Pre-extract PCA embeddings for silhouette cells (we’ll slice per # PCs later)
pca_embed_sil <- pca_embed_all[sil_cells, , drop = FALSE]

# -----------------------------
# 4) Helper: compute mean silhouette for one (pcs, resolution)
# -----------------------------
compute_mean_silhouette <- function(seu, pcs, res) {
  dims_use <- 1:pcs
  
  # ---- clustering
  if (cluster_full_object) {
    seu_run <- seu
    
    seu_run <- FindNeighbors(
      seu_run,
      reduction = pca_reduction,
      dims = dims_use,
      graph.name = graph_name,
      verbose = FALSE
    )
    
    # Use a unique cluster column name for this sweep
    clust_col <- sprintf("sweep_res%.2f_pcs%d", res, pcs)
    seu_run <- FindClusters(
      seu_run,
      graph.name = graph_name,
      resolution = res,
      algorithm = 1,
      verbose = FALSE
    )
    
    # FindClusters writes to Idents and typically to "seurat_clusters" unless you set it.
    # We will grab Idents immediately and store it as a metadata column for clarity.
    seu_run[[clust_col]] <- Idents(seu_run)
    
    # cluster labels for silhouette cells
    cl <- seu_run@meta.data[sil_cells, clust_col, drop = TRUE]
  } else {
    # Cluster only the silhouette subset (faster; may differ from full-object clustering)
    seu_sub <- subset(seu, cells = sil_cells)
    
    seu_sub <- FindNeighbors(
      seu_sub,
      reduction = pca_reduction,
      dims = dims_use,
      graph.name = graph_name,
      verbose = FALSE
    )
    
    seu_sub <- FindClusters(
      seu_sub,
      graph.name = graph_name,
      resolution = res,
      algorithm = 1,
      verbose = FALSE
    )
    
    cl <- Idents(seu_sub)
  }
  
  # Need at least 2 clusters to compute silhouette
  cl <- as.integer(as.factor(cl))
  if (length(unique(cl)) < 2) {
    return(NA_real_)
  }
  
  # ---- distances in PCA space (silhouette subset only)
  X <- pca_embed_sil[, dims_use, drop = FALSE]
  d <- dist(X, method = silhouette_metric)
  
  sil <- cluster::silhouette(cl, d)
  mean(sil[, "sil_width"], na.rm = TRUE)
}

# -----------------------------
# 5) Run the grid search
# -----------------------------
grid_df <- expand.grid(
  pcs = pc_grid,
  resolution = res_grid
) %>%
  arrange(pcs, resolution)

message(sprintf("Running %d combinations (pcs x resolution)...", nrow(grid_df)))

# Use a simple loop (more transparent + easy to interrupt)
mean_sil <- numeric(nrow(grid_df))

for (i in seq_len(nrow(grid_df))) {
  pcs_i <- grid_df$pcs[i]
  res_i <- grid_df$resolution[i]
  
  message(sprintf("  [%d/%d] pcs=%d, res=%.2f", i, nrow(grid_df), pcs_i, res_i))
  mean_sil[i] <- compute_mean_silhouette(obj, pcs = pcs_i, res = res_i)
}

results <- grid_df %>%
  mutate(mean_silhouette = mean_sil)

# Save numeric results
write.csv(results, results_csv, row.names = FALSE)
message("Saved results: ", results_csv)

# -----------------------------
# 6) Plot heatmap like your example
# -----------------------------
# Match your figure orientation: y = number of PCs, x = resolution
# Ensure factors preserve numeric ordering
plot_df <- results %>%
  mutate(
    pcs = factor(pcs, levels = sort(unique(pcs))),
    resolution = factor(resolution, levels = sort(unique(resolution)))
  )

p <- ggplot(plot_df, aes(x = resolution, y = pcs, fill = mean_silhouette)) +
  geom_tile(color = NA) +
  geom_text(
    aes(label = ifelse(is.na(mean_silhouette), "NA",
                       format(round(mean_silhouette, digits_to_label), nsmall = digits_to_label))),
    size = 3
  ) +
  scale_fill_viridis(
    option = "viridis",
    direction = -1,
    na.value = "grey90",
    name = "Mean\nsilhouette"
  ) +
  labs(
    title = plot_title,
    x = "Clustering Resolution",
    y = "Number of PCs (dims)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p)



###############################################################################
# End
###############################################################################


############################################################
## Save table + object
############################################################
write.csv(
  nclust_df,
  file = file.path(fig_dir, "nClusters_vs_resolution_RNA_0p5_to_1p0.csv"),
  row.names = FALSE
)

saveRDS(
  obj,
  file = file.path(fig_dir, "untreated_all_with_RNA_res_0p5_to_1p0.RDS")
)

message("Cluster stability analysis complete.")