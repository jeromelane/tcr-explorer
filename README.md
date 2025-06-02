# tcr-explorer

`tcr_explorer` is an R-based pipeline for TCR repertoire analysis of 10X Genomics single-cell RNA-seq data.

---

## 🧪 Scientific Objective

- Identify distinct **immune cell populations** from single-cell transcriptomic data.
- Detect **gene expression markers** that characterize each cluster.
- Perform **pathway enrichment analysis** to gain functional insight into the biological processes active in specific immune cell subsets.
- Provide data that can contribute to studies in:
  - Immuno-oncology
  - Autoimmunity
  - Vaccine response
  - Tumor microenvironment profiling
  - Immunotherapy biomarker discovery

---

## 📊 Pipeline Overview

### 1️⃣ Input

- A 10X Genomics filtered_feature_bc_matrix directory.
- A `config.yaml` file specifying data location.

### 2️⃣ Quality Control (QC)

- Filter cells based on:
  - Minimum expressed genes (`nFeature_RNA > 200`)
  - Maximum total UMI count (`nCount_RNA < 25000`)
- These filters help remove:
  - Empty droplets
  - Poor quality cells
  - Multiplets or doublets

### 3️⃣ Normalization & Scaling

- Normalize gene expression per cell.
- Identify highly variable genes.
- Scale data to remove unwanted sources of variation.

### 4️⃣ Dimensionality Reduction & Clustering

- Perform PCA for initial dimensionality reduction.
- Build nearest-neighbor graph.
- Cluster cells using Louvain community detection.
- Visualize clusters using UMAP embedding.

### 5️⃣ Visualization

- Export high-resolution UMAP plot (`umap_clusters.png`), showing labeled immune cell clusters.

### 6️⃣ Marker Gene Identification

- Detect cluster-specific marker genes using `FindAllMarkers()` in Seurat.
- Save marker list to `cluster_markers.csv`.

### 7️⃣ Functional Enrichment (Pathway Analysis)

- Select top marker genes from Cluster 0.
- Perform Gene Ontology (GO) enrichment using `clusterProfiler`.
- Export pathway enrichment results to `enrichment_results` folder.

## Installation

### Install R

Make sure R (version ≥ 4.0) is installed on your system.

For Ubuntu example:

```bash
sudo apt update
sudo apt install r-base
```

For other Linux distributions or OS, please refer to the official CRAN installation instructions:
👉 https://cran.r-project.org/

### Install system dependencies

Some R packages require system libraries to compile:

```bash
sudo apt update
sudo apt install libboost-all-dev liblzma-dev libbz2-dev libpcre2-dev libcurl4-openssl-dev libxml2-dev libssl-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libgit2-dev
```

### Install R packages

All R dependencies can be installed automatically via the provided script.

From an R terminal:

```R
source("install.R")
```

From a bash terminal:

```bash
Rscript install.R
```

## Run script

### Download example test data

You can download the test sample with:

```bash
wget https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/vdj_v1_hs_nsclc_multi_5gex_t_b/vdj_v1_hs_nsclc_multi_5gex_t_b_count_filtered_feature_bc_matrix.tar.gz
```

### Execution of the pipeline

Currently the program uses a yaml config file as input to run.
expression_data_dir key should provide the path to expression data such as vdj_v1_hs_nsclc_multi_5gex_t_b_count_filtered_feature_bc_matrix.tar.gz. An example of config file is given in config folder.

From bash terminal:

```bash
Rscript run_pipeline.R --config /your/path/to/config.yaml
```

Outputs:

- umap_clusters.png (cell clusters)
- cluster_markers.csv (differential gene markers)
- cluster_[id]_go_enrichment.csv (pathway enrichment for the cluster [id]) in enrichment_results folder