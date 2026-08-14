## ============================================================
## run_all.R
## Master replication script for "msPCA: Multiple Sparse PCA in R"
##
## Run this file from within the code/ directory:
##   setwd("/path/to/code")
##   source("run_all.R")
##
## Scripts are sourced in the following order:
##   1. Environment details (package versions, session info)
##   2. Worked example (Section 4)
##   3. Benchmarking notebooks (Section 5)
##   4. Case study: S&P 500 (Section 6)
##   5. Synthetic benchmark (README / website figures)
##
## Running this file regenerates every derived file that the paper, the
## vignettes and the pkgdown site display:
##
##   benchmarking/benchmarking_results_*.csv   benchmarking article tables
##   case_study/snp_*.pdf                      paper figures
##   ../vignettes/figures/snp_*.png            case-study vignette figures
##   ../inst/vignette-data/snp_varyingk_results.csv
##                                             read by the case-study vignette
##   ../man/figures/synthetic_*.png            README and pkgdown home page
##   ../notebooks/msPCA_synthetic_results.csv  synthetic benchmark results
##
## The three destinations outside replication/ matter: writing only to the
## replication folder is how the vignette and README came to display figures
## older than the results behind them. The scripts now write to both places,
## so nothing needs to be copied by hand.
##
## What this does NOT do is rebuild the rendered website in docs/. After a
## re-run, the tables and prose in vignettes/ have to be updated to match the
## new numbers, and then pkgdown::build_site() re-renders docs/.
## ============================================================

## Reproducibility settings: force single-threaded BLAS/OpenMP behavior
## to reduce cross-platform numerical differences.
Sys.setenv(
	OPENBLAS_NUM_THREADS = "1",
	OMP_NUM_THREADS = "1",
	MKL_NUM_THREADS = "1",
	BLIS_NUM_THREADS = "1",
	VECLIB_MAXIMUM_THREADS = "1"
)

## 1. Environment details
cat("\n============================================================\n")
cat(" [1/8] environment_details.R\n")
cat("============================================================\n\n")
source("environment_details.R")

## 2. Worked example
cat("\n============================================================\n")
cat(" [2/8] worked_example/worked_example.R\n")
cat("============================================================\n\n")
source("worked_example/worked_example.R")

## 3. Benchmarking
cat("\n============================================================\n")
cat(" [3/8] benchmarking/notebook_mtcars.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_mtcars.R")

cat("\n============================================================\n")
cat(" [4/8] benchmarking/notebook_pitprops.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_pitprops.R")

cat("\n============================================================\n")
cat(" [5/8] benchmarking/notebook_breast.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_breast.R")

cat("\n============================================================\n")
cat(" [6/8] benchmarking/notebook_riboflavin.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_riboflavin.R")

## 4. Case study
cat("\n============================================================\n")
cat(" [7/8] case_study/case_study_SnP.R\n")
cat("============================================================\n\n")
source("case_study/case_study_SnP.R")

## 5. Synthetic benchmark (README / pkgdown home-page figures)
##
## This one lives in ../notebooks/ rather than here, because its outputs are
## package front-page assets rather than paper results. It is included so that
## a single run of this script refreshes every generated file the README, the
## vignettes and the website display -- man/figures/synthetic_*.png among them,
## which previously went two months stale against the results behind them.
##
## Note: unlike the scripts above, it installs any missing dependencies
## (mvtnorm, readr, dplyr, ggplot2, RColorBrewer) rather than failing.
## It also locates the package root itself, so sourcing it from here is safe.
cat("\n============================================================\n")
cat(" [8/8] ../notebooks/notebook_synthetic.R\n")
cat("============================================================\n\n")
source("../notebooks/notebook_synthetic.R")

cat("\n============================================================\n")
cat(" All scripts completed.\n")
cat("============================================================\n\n")