#!/usr/bin/env Rscript

# Fail on warnings to catch install issues early
options(warn=2)

# Install package managers first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

BiocManager::install()


# Install CRAN packages
cran_packages <- c(
  "Seurat",
  "yaml",
  "optparse",
  "ggplot2",
  "tidyverse",
  "data.table"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
bioc_packages <- c("clusterProfiler", "org.Hs.eg.db")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

message("All required packages have been installed.")
