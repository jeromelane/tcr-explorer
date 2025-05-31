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
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# Plot UMAP
p <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle("Immune Cell Clusters")
ggsave("umap_clusters.png", p)

# Find markers for each cluster
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fwrite(markers, file = "cluster_markers.csv")

# Select top genes for enrichment (e.g. cluster 0)
top_genes <- markers %>% filter(cluster == 0) %>% pull(gene)

# Run enrichment analysis using clusterProfiler (KEGG or GO example)
message("Running pathway enrichment...")
enrich_result <- enrichGO(gene = top_genes,
                          OrgDb = "org.Hs.eg.db",
                          keyType = "SYMBOL",
                          ont = "BP", pAdjustMethod = "BH")

# Save enrichment results
write.csv(as.data.frame(enrich_result), file = "go_enrichment.csv")

message("Pipeline completed successfully.")
