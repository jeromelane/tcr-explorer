#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tcr.explorer)
  library(yaml)
})

option_list <- list(
  make_option("--config", type="character", help="Path to config.yaml", metavar="file")
)

opt <- parse_args(OptionParser(option_list=option_list))
config <- yaml::read_yaml(opt$config)

sce <- load_expression_data(config)
sce <- cluster_cells(sce, config)
sce <- annotate_tcrs(sce, config)
run_gsea(sce, config)
run_cross_reactivity(sce, config)
