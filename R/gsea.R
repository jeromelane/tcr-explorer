run_gsea <- function(seurat_obj, config) {
  message("Running GSEA...")
  markers <- Seurat::FindAllMarkers(seurat_obj)
  m_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7")  # immunologic signatures
  gsea_result <- clusterProfiler::enricher(
    gene = markers$gene,
    TERM2GENE = m_df[, c("gs_name", "gene_symbol")]
  )
  print(gsea_result)
}
