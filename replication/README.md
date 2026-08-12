# Replication package — *msPCA: An R Package for Sparse PCA with Multiple Components*

This folder contains the code and data needed to reproduce the results in the paper and in the [benchmarking article](https://jeanpauphilet.github.io/msPCA/articles/benchmarking.html). Scripts are self-contained: each sources only base R and the packages listed below.

It is excluded from the package build (`.Rbuildignore`), so it ships with the GitHub repository but not with the CRAN tarball. The `msPCA` package itself is the source tree one level up.

## Folder structure

```
replication/
├── worked_example/
│   └── worked_example.R      # Section 4 — basic workflow on mtcars
├── case_study/
│   ├── case_study_SnP.R      # Section 5 — S&P 500 factor analysis
│   ├── snp_varyingk_results.csv         # Pre-computed sweep over sparsity levels
│   ├── snp_k10_orthogonality_loadings.txt  # Loadings at k = 10, orthogonality constraint
│   └── snp_k10_correlation_loadings.txt    # Loadings at k = 10, zero-correlation constraint
├── benchmarking/
│   ├── notebook_mtcars.R     # Section 6 — Table 2 (mtcars, p = 11)
│   ├── notebook_pitprops.R   # Section 6 — Table 3 (Pitprops, p = 13)
│   ├── notebook_breast.R     # Section 6 — Table 4 (breast cancer, p = 500)
│   ├── notebook_riboflavin.R # Section 6 — Table 5 (riboflavin, p = 4088)
│   ├── benchmarking_results_mtcars.csv
│   ├── benchmarking_results_pitprops.csv
│   ├── benchmarking_results_breast.csv
│   └── benchmarking_results_riboflavin.csv
├── run_all.R                 # Master script — runs all scripts in order
├── environment_details.R     # Captures R version and package versions
├── package_versions.csv      # Output of environment_details.R
└── r_platform_info.csv       # Output of environment_details.R
```

### One file is not included

`case_study/SnP_returns_cleaned.csv` (daily log-returns, 423 stocks, 2010–2019, about 36 MB) is too large to version-control and is not in this folder. Two consequences:

- `case_study/case_study_SnP.R` will not run end to end as shipped. Rebuild the returns file from the [source Kaggle dataset](https://www.kaggle.com/datasets/yash16jr/s-and-p500-daily-update-dataset) (CC0 1.0) first; section (B) of `../data-raw/snp500.R` documents the derivation.
- Everything the S&P 500 vignette needs is already in the package: the deflated correlation matrix ships as `data(snp500)`, and the pre-computed sparsity sweep as `inst/vignette-data/snp_varyingk_results.csv`.

## Requirements

macOS on Apple M2. Required packages and the versions used for the results reported in the benchmarking article:

| Package | Version |
|---------|---------|
| msPCA | 0.5.0 |
| elasticnet | 1.3 |
| PMA | 1.2.4 |
| sparsepca | 0.1.2 |
| mixOmics | 6.36.0 |
| nsprcomp | 0.5.1.2 |
| amanpg | 0.1.0 |
| reticulate | 1.40.0 |
| scikit-learn (Python 3.11) | 1.6.1 |

Additional packages used in specific scripts: `RSpectra`, `readr`, `dplyr`, `MASS`, `hdi`, `datasets` (base R).

`sessionInfo.txt`, `package_versions.csv` and `r_platform_info.csv` in this folder are the authoritative record of the environment used; the table above is a summary.

## Reproducing results

### Run everything at once

The easiest way to reproduce all results is to use the master script `run_all.R`.
From within the `replication/` directory, run from the R console:

```r
setwd("/path/to/replication")   # set to the location of this folder
source("run_all.R")
```

Or from the command line:

```bash
cd /path/to/replication
Rscript run_all.R
```

`run_all.R` sets BLAS/OpenMP thread counts to 1 to maximize cross-platform reproducibility. Using 2 or 4 threads may improve performance, but can introduce small numerical differences in non-convex optimization results.

This sources all scripts in the correct order:
1. `environment_details.R` — captures R version and package versions
2. `worked_example/worked_example.R` — Section 4 worked example
3. `benchmarking/notebook_mtcars.R` — Section 5 benchmarking (mtcars)
4. `benchmarking/notebook_pitprops.R` — Section 5 benchmarking (Pitprops)
5. `benchmarking/notebook_breast.R` — Section 5 benchmarking (breast cancer)
6. `benchmarking/notebook_riboflavin.R` — Section 5 benchmarking (riboflavin)
7. `case_study/case_study_SnP.R` — Section 6 case study

Note: the riboflavin and breast cancer benchmarking notebooks are computationally intensive and may take several minutes to complete.

### Run scripts individually

Each script can also be sourced individually from within the `replication/` directory. Each benchmarking notebook writes its results table to the corresponding CSV file in `benchmarking/`. The case study script writes figures and tables to `case_study/`.

To verify the R environment, run `environment_details.R` first; it writes `package_versions.csv` and `r_platform_info.csv`.
