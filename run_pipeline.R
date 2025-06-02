#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(Seurat)
  library(tidyverse)
  library(data.table)
  library(clusterProfiler)
})

files <- list.files("R", full.names = TRUE, pattern = "\\.R$")
sapply(files, source)

# Parse config file from command line
option_list <- list(
  make_option("--config", type = "character", help = "Path to config.yaml", metavar = "file")
)

opt <- parse_args(OptionParser(option_list = option_list))
config <- yaml::read_yaml(opt$config)

# Load expression data (10X filtered_feature_bc_matrix directory)
message("Loading 10X gene expression data...")
data <- Read10X(data.dir = config$expression_data_dir)
seurat_obj <- CreateSeuratObject(counts = data)

# Basic QC filtering (optional)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nCount_RNA < 25000)

# Normalize and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# PCA and clustering
seurat_obj <- RunPCA(seurat_obj)

# Calculate variance explained
pca_stdev <- seurat_obj@reductions$pca@stdev
pca_var <- pca_stdev^2
pca_var_explained <- pca_var / sum(pca_var) * 100

# Convert to data frame for easy plotting
pca_df <- data.frame(
  PC = 1:length(pca_var_explained),
  VarianceExplained = pca_var_explained
)

# Plot elbow (variance explained)
library(ggplot2)

p <- ggplot(pca_df[1:20,], aes(x = PC, y = VarianceExplained)) +
  geom_point(size = 2, color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red") +
  labs(title = "PCA Variance Explained",
       x = "Principal Component",
       y = "Percent Variance Explained") +
  theme_minimal()

# Save plot
ggsave("pca_variance_explained.png", p, width = 7, height = 5)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5) # Louvain clustering algorithm
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Immune Cell Clusters")
ggsave("umap_clusters.png", p)

# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fwrite(markers, file = "cluster_markers.csv")

# Get all cluster IDs
clusters <- unique(markers$cluster)

# Create folder for enrichment results
dir.create("enrichment_results", showWarnings = FALSE)

# Loop over clusters
for (clust in clusters) {
  message(glue::glue("Running enrichment for cluster {clust}"))

  top_genes <- markers %>% filter(cluster == clust) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 100) %>% pull(gene)


  # Only proceed if there are enough genes
  if (length(top_genes) > 10) {
    enrich_result <- enrichGO(
      gene = top_genes,
      OrgDb = "org.Hs.eg.db",
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH"
    )

    # Save enrichment result
    write.csv(as.data.frame(enrich_result),
              file = sprintf("enrichment_results/cluster_%s_go_enrichment.csv", clust))
  }
}
message("Pipeline completed successfully.")
