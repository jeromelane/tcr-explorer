# tcr-explorer

`tcr_explorer` is an R-based pipeline for TCR repertoire analysis.

---

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
