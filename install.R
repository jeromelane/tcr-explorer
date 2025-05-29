# install.R

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install CRAN packages
cran_packages <- c(
  "Seurat",
  "scRepertoire",
  "yaml",
  "optparse",
  "ggplot2",
  "msigdbr"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
bioc_packages <- c("clusterProfiler")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

message("All required packages have been installed.")
