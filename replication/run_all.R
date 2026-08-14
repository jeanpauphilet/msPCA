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