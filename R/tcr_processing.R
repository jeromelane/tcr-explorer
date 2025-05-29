annotate_tcrs <- function(seurat_obj, config) {
  message("Annotating TCR data...")
  tcr_data <- read.csv(config$tcr_data_file)
  seurat_obj <- scRepertoire::combineTCR(tcr_data, samples = "sample1", cells = "T-AB") %>%
    scRepertoire::combineExpression(seurat_obj)
  return(seurat_obj)
}
