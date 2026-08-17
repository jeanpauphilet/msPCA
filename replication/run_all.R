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
##
## This script runs the contents of replication/ only. It regenerates:
##
##   benchmarking/benchmarking_results_*.csv   benchmarking article tables
##   case_study/snp_*.pdf                      paper figures
##   ../vignettes/figures/snp_*.png            case-study vignette figures
##   ../inst/vignette-data/snp_varyingk_results.csv
##                                             read by the case-study vignette
##
## The last two live outside this folder because the case-study vignette and
## the website read them from there; case_study_SnP.R writes both destinations
## on every run, so nothing needs copying by hand.
##
## NOT covered here: the synthetic benchmark behind the README and pkgdown
## home-page figures (man/figures/synthetic_*.png). That is
## notebooks/notebook_synthetic.R, which sits outside replication/ and is run
## separately from the repository root.
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

## Python side (the scikit-learn comparison). Sourced first because
## reticulate::py_require() can only pin package and interpreter versions
## before Python is initialised, and environment_details.R below initialises
## it. See python_setup.R for what this does and how to override it.
source("python_setup.R")

## 1. Environment details
cat("\n============================================================\n")
cat(" [1/7] environment_details.R\n")
cat("============================================================\n\n")
source("environment_details.R")

## 2. Worked example
cat("\n============================================================\n")
cat(" [2/7] worked_example/worked_example.R\n")
cat("============================================================\n\n")
source("worked_example/worked_example.R")

## 3. Benchmarking
cat("\n============================================================\n")
cat(" [3/7] benchmarking/notebook_mtcars.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_mtcars.R")

cat("\n============================================================\n")
cat(" [4/7] benchmarking/notebook_pitprops.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_pitprops.R")

cat("\n============================================================\n")
cat(" [5/7] benchmarking/notebook_breast.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_breast.R")

cat("\n============================================================\n")
cat(" [6/7] benchmarking/notebook_riboflavin.R\n")
cat("============================================================\n\n")
source("benchmarking/notebook_riboflavin.R")

## 4. Case study
cat("\n============================================================\n")
cat(" [7/7] case_study/case_study_SnP.R\n")
cat("============================================================\n\n")
# source("case_study/case_study_SnP.R")

cat("\n============================================================\n")
cat(" All scripts completed.\n")
cat("============================================================\n\n")