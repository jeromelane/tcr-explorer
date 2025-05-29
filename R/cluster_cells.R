cluster_cells <- function(seurat_obj, config) {
  message("Clustering cells...")
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  seurat_obj <- Seurat::RunPCA(seurat_obj)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- Seurat::FindClusters(seurat_obj)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)
  return(seurat_obj)
}
