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
│   ├── snp_nsprcomp_results.csv         # Cached nsprcomp comparison (needs raw returns to recompute)
│   ├── snp_k10_orthogonality_loadings.txt  # Loadings at k = 10, orthogonality constraint
│   └── snp_k10_correlation_loadings.txt    # Loadings at k = 10, zero-correlation constraint
├── benchmarking/
│   ├── bench_utils.R         # Runtime + peak-memory harness (see below)
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

`case_study/SnP_returns_cleaned.csv` (daily log-returns, 423 stocks, 2010–2019, about 36 MB) is too large to version-control and is not in this folder. Consequences:

- `case_study/case_study_SnP.R` still runs end to end as shipped. Every `msPCA` result is computed from the market-deflated correlation matrix distributed with the package as `data(snp500)`, which is all those methods need.
- Only the `nsprcomp` reference curve requires the raw returns, since it takes a data matrix rather than a correlation matrix. The script tests for the returns file: if present it rebuilds the deflated data matrix `XR`, runs `nsprcomp`, and refreshes `case_study/snp_nsprcomp_results.csv`; if absent it reads that cached file instead; if neither exists it warns and omits `nsprcomp` from the tables and figures.
- To recompute the `nsprcomp` numbers from scratch, rebuild the returns file from the [source Kaggle dataset](https://www.kaggle.com/datasets/yash16jr/s-and-p500-daily-update-dataset) (CC0 1.0) into `case_study/`; section (B) of `../data-raw/snp500.R` documents the derivation.
- The vignette is likewise self-contained: the pre-computed sparsity sweep ships as `inst/vignette-data/snp_varyingk_results.csv`.

## Requirements

macOS on Apple M2. Required packages and the versions used for the results reported in the benchmarking article:

| Package | Version |
|---------|---------|
| msPCA | 0.5.1 |
| elasticnet | 1.3 |
| PMA | 1.2.4 |
| sparsepca | 0.1.2 |
| mixOmics | 6.36.0 |
| nsprcomp | 0.5.1.2 |
| amanpg | 0.1.0 |
| reticulate | 1.40.0 |
| scikit-learn (Python 3.11) | 1.6.1 |

Additional packages used in specific scripts: `ggplot2`, `dplyr`, `tidyr`, `tibble`, `MASS`, `hdi`, `datasets` (base R). `RSpectra` is needed only for the `nsprcomp` branch of the case study, i.e. when the raw returns file is present. `callr` is required by all four benchmarking notebooks.

`sessionInfo.txt`, `package_versions.csv` and `r_platform_info.csv` in this folder are the authoritative record of the environment used; the table above is a summary.

### How runtime and memory are measured

`benchmarking/bench_utils.R` runs each method in a fresh R subprocess (`callr::r()`) and records both timing and peak resident set size, read from `getrusage(RUSAGE_SELF).ru_maxrss` via a small compiled probe. R-level memory metrics (`gc()`, `bench::mark(mem_alloc)`) are deliberately not used: `msPCA` works in RcppEigen buffers and `scikit-learn` in the embedded CPython heap, neither of which R's garbage collector sees, so an R-level metric would under-report those two and fully charge the pure-R packages.

**One repetition per process.** `ru_maxrss` is a monotone high-water mark, so several methods in one session — or even several repetitions of one method — contaminate each other. Every repetition therefore gets its own fresh R process. The algorithm is invoked exactly as many times as it would be otherwise, so this costs only the extra R startups, and in exchange each repetition yields an independent peak-RSS sample.

That independence is necessary, because peak RSS is not deterministic even for deterministic code: it depends on when the garbage collector happens to run. Across two identical runs of Pitprops, with identical seeds and bit-identical FVE and orthogonality, `nscumcomp`'s working set came out as 72.4 MB and then 211.2 MB, and `mixOmics`'s as 0.6 MB and then 13.5 MB. Memory is therefore reported as a median over repetitions with `working_set_min_mb` / `working_set_max_mb` alongside; check the range before leaning on the median.

Tuning sweeps run in the parent and are excluded from all reported figures.

There is deliberately **no warm-up call** before the baseline reading, so `mem_delta_mb` includes whatever each package allocates on its first call — roughly 4 MB of namespace, compiled-code and S4 method-table loading. Warming up to remove that cost was tried and rejected: any allocation before the baseline leaves R's garbage collector with higher heap growth targets, so the measured run then collects less often and reaches a higher peak, by more than the baseline rose. On Pitprops this inflated `mixOmics` from 0.6 to 11.6 MB and `nscumcomp` from 72.4 to 275.8 MB, while `msPCA` went *down*, 0.8 to 0.1 MB — because `msPCA` allocates in Eigen rather than on R's heap. That biases the comparison in `msPCA`'s favour, which is the one thing this harness must not do. Carrying ~4 MB of loading cost in every row is the lesser evil, and it is small next to the figures that carry the comparison. See the note at the top of `bench_utils.R`.

### Repetitions and seeds

Each method is run `REPS` times with a **different seed per repetition** (43, 44, 45, …). The harness calls `set.seed()` in the child immediately before each repetition, so benchmarked functions must not set their own seed; `scikit-learn`, which does not read R's RNG, receives the seed explicitly as `random_state`. `REPS` is 5 for mtcars and Pitprops and 3 for breast cancer and riboflavin, where individual runs take seconds to minutes.

Runtime, memory, FVE and the orthogonality violation are all reported as medians over the repetitions, each with its range alongside. Report the ranges: neither wall-clock time nor peak RSS is reproducible to the precision a bare median implies. Runtime varies with machine conditions and with seed-dependent iteration counts — `msPCA (Sigma)` on riboflavin has come out at 58.6 s and at 85.7 s in different sittings — and peak RSS varies with garbage-collector timing. FVE and orthogonality are medians with the spread reported alongside, so that solution stability under reseeding is visible: deterministic methods (`elasticnet`, `PMA`, `sparsepca`, `mixOmics`, `amanpg`, `prcomp`) collapse to a zero-width FVE range, while `nsprcomp` and `nscumcomp` use random initialisation and do not.

### Columns in each results CSV

| Column | Meaning |
|---|---|
| `runtime_s` | Median elapsed time over `REPS` repetitions |
| `runtime_min_s`, `runtime_max_s` | Range of elapsed time across repetitions; wall-clock timing varies with machine conditions as well as with seed-dependent iteration counts |
| `working_set_mb` | `input_mb + mem_delta_mb`, median over repetitions — the headline memory figure |
| `working_set_min_mb`, `working_set_max_mb` | Range of the above across repetitions; a wide range means the median should not be leaned on |
| `input_mb` | Size of the matrix handed to the method (Σ is p × p, X is n × p) |
| `mem_delta_mb` | Peak RSS minus post-setup baseline: the solver's own working memory, excluding first-call package loading |
| `peak_rss_mb` | Absolute high-water RSS, including the R interpreter and dependencies |
| `baseline_rss_mb` | RSS after packages, `setup()` and inputs, before the method is called |
| `gc_max_mb` | R-level `gc()` maximum; a cross-check only, not comparable across methods |
| `nnz_pc1` … | Nonzeros per component, from the representative repetition |
| `fve`, `fve_min`, `fve_max` | Median fraction of variance explained, and its range over repetitions |
| `orth`, `orth_max` | Median and worst-case orthogonality violation over repetitions |
| `n_fail` | Repetitions in which the method failed (`nscumcomp` can fail on co-linear axes) |

The representative repetition is the one whose FVE is closest to the median; `nnz_pc*` is read from it, so those counts always come from a single real solution rather than an average.

Compare methods on `working_set_mb`. Do not compare `peak_rss_mb` across languages — the `sklearn` row carries a whole CPython heap in its baseline, which is a fact about `reticulate` rather than about the algorithm.

Requires macOS or Linux. On Windows the memory columns are returned as `NA`; runtimes are unaffected.

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

`run_all.R` runs the contents of this folder only. Between them these scripts regenerate the benchmarking tables, the paper figures, and the two sets consumed elsewhere in the package: `../vignettes/figures/snp_*.png` and `../inst/vignette-data/snp_varyingk_results.csv`. The case-study script writes both destinations on every run, so those need no copying by hand.

Two things it deliberately leaves alone:

- **The synthetic benchmark.** `notebooks/notebook_synthetic.R` produces the figures shown in the top-level README and on the pkgdown home page (`man/figures/synthetic_*.png`). It lives outside `replication/` and is run separately, from the repository root. It writes to `notebooks/`, so the two PNGs must then be copied into `man/figures/` — otherwise the README keeps displaying the previous run's figures.
- **The rendered site.** After a re-run, update the tables and prose in `vignettes/` to match the new numbers, then re-render `docs/` with `pkgdown::build_site()`.

Note: the riboflavin and breast cancer benchmarking notebooks are computationally intensive and may take several minutes to complete.

### Run scripts individually

Each script can also be sourced individually from within the `replication/` directory. Each benchmarking notebook writes its results table to the corresponding CSV file in `benchmarking/`.

The case study script writes to `case_study/` **and into the package sources**: the four figures are saved as PDFs here for the paper and as 150-dpi PNGs in `../vignettes/figures/` for the vignette and pkgdown site, and the sparsity sweep is copied to `../inst/vignette-data/snp_varyingk_results.csv`, which the vignette reads via `system.file()`. Both destinations are written on every run, so re-running the case study refreshes the website too. This matters because the two previously drifted apart silently — the vignette spent several days describing factors that a re-run had already replaced. Because the script assumes it is sourced from `replication/`, run it from there rather than from `case_study/`.

To verify the R environment, run `environment_details.R` first; it writes `package_versions.csv` and `r_platform_info.csv`.
