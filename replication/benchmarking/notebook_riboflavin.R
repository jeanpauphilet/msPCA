############################################################
## Benchmarking notebook: msPCA vs competing packages on Riboflavin
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: hdi::riboflavin (p = 4,088 genes, n = 71 bacterial cultures)
## Task:    r = 2 sparse PCs with k = 20 nonzeros each
##
## Methods compared: msPCA (Sigma and X input), elasticnet, mixOmics,
##                   nsprcomp, dense PCA
## For mixOmics::spca(), keepX = rep(k, r) matches the cardinality
## budget exactly. For nsprcomp::nsprcomp(), k = rep(k, r) gives the
## per-component nonzero budget.
##
## NOTE: p >> n (4,088 >> 71). The correlation matrix S is rank-
## deficient (rank <= 70). Methods that operate on S directly
## (msPCA type = "Sigma") may produce fewer meaningful components;
## results should be interpreted accordingly.
##
## Each method is timed AND memory-profiled by bench_utils.R, which runs
## it in a fresh R subprocess and records peak resident set size from
## getrusage(). See bench_utils.R for why R-level metrics (gc(),
## bench::mark) would bias the comparison in msPCA's favour.
##
## NOTE ON MEMORY: this is the notebook where the memory comparison is
## decisive. The correlation matrix S is 4088^2 x 8 bytes = about
## 134 MB; the data matrix X is 71 x 4088 x 8 = about 2.3 MB, roughly
## 57x smaller. msPCA's type = "X" path applies Sigma implicitly as
## X^T(X v)/(n-1) and never allocates the p x p matrix, so the saving
## shows up in BOTH input_mb (S is not an input) and mem_delta_mb (no
## Eigen copy of S is made). Compare the working_set_mb column, which
## sums the two without double-counting.
############################################################

library("msPCA")
library("elasticnet")
library("mixOmics")
library("nsprcomp")

source("benchmarking/bench_utils.R")

## Helper: unit-normalise columns
unit_norm <- function(L) {
  for (j in seq_len(ncol(L))) {
    nm <- sqrt(sum(L[, j]^2))
    if (nm > 0) L[, j] <- L[, j] / nm
  }
  L
}

## Load Riboflavin production data.
## hdi::riboflavin is a list with $x (71 x 4088 log-expression matrix)
## and $y (71-element log-riboflavin-production response vector).
data("riboflavin", package = "hdi")
X_ribo <- riboflavin$x               # n = 71, p = 4088
X_ribo <- scale(X_ribo)              # centre and scale

## Build correlation matrix (rank-deficient: rank <= n-1 = 70).
## This is ~134 MB and is needed in the parent to evaluate FVE and
## orthogonality for every method, and as the input to msPCA (Sigma).
cat("Computing correlation matrix (p =", ncol(X_ribo), ")... ")
S <- cor(riboflavin$x)
cat("done. S is", round(as.numeric(object.size(S)) / 1024^2, 1), "MB;",
    "X is", round(as.numeric(object.size(X_ribo)) / 1024^2, 1), "MB.\n")

r <- 2L; k <- 20L

## Repetitions per method. Each repetition runs in its OWN process, so every
## one yields an independent peak-RSS sample; runtime, memory and FVE are all
## reported as medians over repetitions. Five is the minimum that gives the
## memory median any resistance to outliers -- peak RSS varies substantially
## run to run for methods that allocate heavily on R's heap (see bench_utils.R).
REPS <- 5L

MAXITER <- 100L


############################################################
## Measured runs -- one fresh subprocess per method.
## Each fun() returns a plain p x r loadings matrix.
############################################################

## --- msPCA, Sigma input (dense p x p covariance operator) ---
cat("Running msPCA (Sigma)...\n")
b_mspca_S <- bench_method(
  fun = function(S, r, k, maxIter) {
    mspca(S, type = "Sigma", r = r, ks = rep(k, r), verbose = FALSE,
          maxIter = maxIter, feasibilityConstraintType = 0)$x_best
  },
  inputs = list(S = S, r = r, k = k, maxIter = MAXITER),
  packages = "msPCA", reps = REPS, label = "msPCA (Sigma)")

## --- msPCA, X input (matrix-free Gram operator, never forms p x p) ---
cat("Running msPCA (X)...\n")
b_mspca_X <- bench_method(
  fun = function(X, r, k, maxIter) {
    mspca(X, type = "X", r = r, ks = rep(k, r), verbose = FALSE,
          maxIter = maxIter, feasibilityConstraintType = 0)$x_best
  },
  inputs = list(X = X_ribo, r = r, k = k, maxIter = MAXITER),
  packages = "msPCA", reps = REPS, label = "msPCA (X)")

## --- elasticnet (data matrix input; the Gram path is impractical at p = 4088) ---
cat("Running elasticnet...\n")
b_enet <- bench_method(
  fun = function(X, r, k) {
    elasticnet::spca(X, K = r, sparse = "varnum",
                     para = rep(k, r), type = "predictor")$loadings
  },
  inputs = list(X = X_ribo, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (X)")

## --- mixOmics::spca (keepX = rep(k, r) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
cat("Running mixOmics...\n")
b_mixo <- bench_method(
  fun = function(X, r, k) {
    unit_norm(mixOmics::spca(X, ncomp = r, keepX = rep(k, r),
                             center = TRUE)$loadings$X)
  },
  inputs = list(X = X_ribo, r = r, k = k),
  globals = list(unit_norm = unit_norm),
  packages = "mixOmics", reps = REPS, label = "mixOmics")

## --- nsprcomp (k nonzeros per PC; basic deflation) ---
cat("Running nsprcomp...\n")
b_nspr <- bench_method(
  fun = function(X, r, k) {
    nsprcomp::nsprcomp(X, ncomp = r, k = rep(k, r),
                       nneg = FALSE, center = TRUE)$rotation
  },
  inputs = list(X = X_ribo, r = r, k = k),
  packages = "nsprcomp", reps = REPS, label = "nsprcomp")

## --- Dense PCA ---
cat("Running dense PCA...\n")
b_pca <- bench_method(
  fun = function(X, r) prcomp(X, scale. = TRUE)$rotation[, seq_len(r)],
  inputs = list(X = riboflavin$x, r = r),
  packages = character(), reps = REPS, label = "Dense PCA")


############################################################
## Console report and summary table
############################################################

benches <- list(b_mspca_S, b_mspca_X, b_enet, b_mixo, b_nspr, b_pca)

cat("\n")
for (b in benches) report_bench(b, S)

results <- cbind(bench_table(benches), quality_table(benches, S, r))

print(results, row.names = FALSE)
write.csv(results, "benchmarking/benchmarking_results_riboflavin.csv", row.names = FALSE)
