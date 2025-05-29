load_expression_data <- function(config) {
  message("Loading expression data...")
  data <- Seurat::Read10X(data.dir = config$expression_data_dir)
  seurat_obj <- Seurat::CreateSeuratObject(data)
  return(seurat_obj)
}
