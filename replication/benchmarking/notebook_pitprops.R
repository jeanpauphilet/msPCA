## ============================================================
## Pitprops benchmark: msPCA vs competing packages
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: Pitprops correlation matrix (Jeffers 1967, p = 13)
##   13 physical measurements on 180 timber pit-prop specimens.
##   Classic sparse PCA benchmark; reference results in
##   Zou, Hastie & Tibshirani (2006, Table 3).
##
## Task: r = 6 sparse PCs with k = 4 nonzeros each
##
## Note: mixOmics::spca() and nsprcomp::nscumcomp() require a raw
## data matrix. Since only the correlation matrix is available for
## this dataset, we generate pseudo-data via MASS::mvrnorm() with
## the Pitprops matrix as the population covariance. FVE and
## orthogonality are then evaluated on the original pitprops matrix.
##
## Each method is timed AND memory-profiled by bench_utils.R, which
## runs it in a fresh R subprocess and records peak resident set size
## from getrusage(). See bench_utils.R for why R-level metrics (gc(),
## bench::mark) would bias the comparison in msPCA's favour.
##
## NOTE ON MEMORY: p = 13, so the numerical working set of every
## method is a few kilobytes. msPCA, PMA, mixOmics, amanpg and dense
## PCA accordingly report close to zero. elasticnet, sparsepca and
## nscumcomp do not: they carry a fixed working-memory floor of
## roughly 10-25 MB that is present even at this size, on top of
## which cost scales with p. Note also that the msPCA (X) row here
## consumes MORE memory than msPCA (Sigma): the pseudo-data matrix is
## 500 x 13 while the correlation matrix is only 13 x 13. That is the
## expected ordering whenever n > p, and is worth reporting honestly
## alongside the n << p results in notebook_riboflavin.R.
## ============================================================

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("amanpg")
library("mixOmics")
library("nsprcomp")
library("reticulate")

source("benchmarking/bench_utils.R")

## Python packages (installed via pip: scikit-learn)
py_bin <- Sys.which("python3")
use_python(py_bin, required = TRUE)
sk_decomp    <- import("sklearn.decomposition")

## Helper: unit-normalise columns
unit_norm <- function(L) {
  for (j in seq_len(ncol(L))) {
    nm <- sqrt(sum(L[, j]^2))
    if (nm > 0) L[, j] <- L[, j] / nm
  }
  L
}

## Pitprops correlation matrix from Jeffers (1967), Table 1.
## Variables: topdiam, length, moist, testsg, ovensg, ringtop,
##            ringbut, bowmax, bowdist, whorls, clear, knots, diaknot
pitprops <- matrix(c(
	1.000,  0.954,  0.364,  0.342, -0.129,  0.313,  0.496,  0.424,
	0.592,  0.545,  0.084, -0.019,  0.134,
	0.954,  1.000,  0.297,  0.284, -0.118,  0.291,  0.503,  0.419,
	0.648,  0.569,  0.076, -0.036,  0.144,
	0.364,  0.297,  1.000,  0.882, -0.148,  0.153, -0.029, -0.054,
	0.125, -0.081,  0.162,  0.220,  0.126,
	0.342,  0.284,  0.882,  1.000,  0.220,  0.381,  0.174, -0.059,
	0.137, -0.014,  0.097,  0.169,  0.015,
 -0.129, -0.118, -0.148,  0.220,  1.000,  0.364,  0.296,  0.004,
 -0.039,  0.037, -0.091, -0.145, -0.208,
	0.313,  0.291,  0.153,  0.381,  0.364,  1.000,  0.813,  0.090,
	0.211,  0.274, -0.036,  0.024, -0.329,
	0.496,  0.503, -0.029,  0.174,  0.296,  0.813,  1.000,  0.372,
	0.465,  0.679, -0.113, -0.232, -0.424,
	0.424,  0.419, -0.054, -0.059,  0.004,  0.090,  0.372,  1.000,
	0.482,  0.557,  0.061, -0.357, -0.172,
	0.592,  0.648,  0.125,  0.137, -0.039,  0.211,  0.465,  0.482,
	1.000,  0.526,  0.085, -0.127, -0.175,
	0.545,  0.569, -0.081, -0.014,  0.037,  0.274,  0.679,  0.557,
	0.526,  1.000, -0.319, -0.368, -0.545,
	0.084,  0.076,  0.162,  0.097, -0.091, -0.036, -0.113,  0.061,
	0.085, -0.319,  1.000,  0.029,  0.054,
 -0.019, -0.036,  0.220,  0.169, -0.145,  0.024, -0.232, -0.357,
 -0.127, -0.368,  0.029,  1.000,  0.220,
	0.134,  0.144,  0.126,  0.015, -0.208, -0.329, -0.424, -0.172,
 -0.175, -0.545,  0.054,  0.220,  1.000
), nrow = 13, byrow = TRUE)
varnames <- c("topdiam", "length", "moist", "testsg", "ovensg",
							"ringtop", "ringbut", "bowmax", "bowdist", "whorls",
							"clear", "knots", "diaknot")
rownames(pitprops) <- colnames(pitprops) <- varnames
p <- 13L; r <- 6L; k <- 4L

## Pseudo-data for methods requiring a raw data matrix (mixOmics, nsprcomp).
## Generated from the pitprops correlation matrix; n = 500 gives stable loadings.
set.seed(42)
X_pseudo <- MASS::mvrnorm(n = 500L, mu = rep(0, p), Sigma = pitprops)
colnames(X_pseudo) <- varnames

## Repetitions per method. Each repetition runs in its OWN process, so every
## one yields an independent peak-RSS sample; runtime, memory and FVE are all
## reported as medians over repetitions. Five is the minimum that gives the
## memory median any resistance to outliers -- peak RSS varies substantially
## run to run for methods that allocate heavily on R's heap (see bench_utils.R).
REPS <- 5L

## Tuning sweeps are exploratory and are NOT part of the measured runs.
## They stay in the parent process so their allocations cannot contaminate
## any method's peak-RSS reading. Set to FALSE to skip.
RUN_TUNING_SWEEPS <- TRUE

## Values selected from the sweeps below.
PMA_SUMABSV     <- 1.727
SPARSEPCA_ALPHA <- 0.004
AMANPG_LAMBDA1  <- c(1, 1, 5, 0.4, 0.4, 0)   # tuned for ~4 NNZ per PC
SKLEARN_ALPHA   <- 3.1

## nscumcomp() may fail on datasets with high inter-variable correlations
## ("Co-linear principal axes"); it is retried up to this many times.
NSCUM_MAX_ATTEMPTS <- 20L


## ============================================================
## Tuning sweeps (untimed, unprofiled)
## ============================================================

if (RUN_TUNING_SWEEPS) {
  ## --- PMA (sumabsv tuned for ~4 NNZ per PC on Pitprops) ---
  cat("PMA sumabsv sweep:\n")
  for (sv in seq(1.7, 1.8, by = 0.001)) {
    pma_try <- PMA::SPC(pitprops, sumabsv = sv, K = r, orth = TRUE, trace = FALSE)
    cat("  sumabsv =", sv,
        "| NNZ:", paste(colSums(abs(pma_try$v) > 0), collapse = " "),
        "| FVE:", round(fraction_variance_explained(pitprops, pma_try$v), 4), "\n")
  }

  ## --- sparsepca (alpha tuned for ~4 NNZ per PC) ---
  cat("sparsepca alpha sweep:\n")
  for (a in seq(0.001, 0.005, by = 0.0005)) {
    sp_try <- unit_norm(
      sparsepca::spca(pitprops, k = r, alpha = a, verbose = FALSE, scale = TRUE)$loadings)
    cat("  alpha =", a,
        "| NNZ:", paste(colSums(abs(sp_try) > 0), collapse = " "),
        "| FVE:", round(fraction_variance_explained(pitprops, sp_try), 4), "\n")
  }

  ## --- sklearn SparsePCA (alpha tuned for ~4 NNZ per PC) ---
  cat("sklearn alpha sweep:\n")
  for (a in seq(2.5, 3.5, by = 0.05)) {
    m <- sk_decomp$SparsePCA(n_components = r, alpha = a, random_state = 43L)
    m$fit(X_pseudo)
    L <- unit_norm(t(m$components_))
    cat("  alpha =", a, "| NNZ:", paste(colSums(abs(L) > 0), collapse = " "), "\n")
  }
}


## ============================================================
## Measured runs -- one fresh subprocess per method.
## Each fun() returns a plain p x r loadings matrix.
## ============================================================

#### R METHODS ####

## --- msPCA, Sigma input (dense p x p covariance operator) ---
b_mspca_S <- bench_method(
  fun = function(Sigma, r, k) {
    mspca(Sigma, r = r, ks = rep(k, r),
          verbose = FALSE, feasibilityConstraintType = 0)$x_best
  },
  inputs = list(Sigma = pitprops, r = r, k = k),
  packages = "msPCA", reps = REPS, label = "msPCA (Sigma)")

## --- msPCA, X input (matrix-free Gram operator, never forms p x p) ---
b_mspca_X <- bench_method(
  fun = function(X, r, k) {
    mspca(X, type = "X", r = r, ks = rep(k, r),
          verbose = FALSE, feasibilityConstraintType = 0)$x_best
  },
  inputs = list(X = X_pseudo, r = r, k = k),
  packages = "msPCA", reps = REPS, label = "msPCA (X)")

## --- elasticnet, Sigma input ---
b_enet_S <- bench_method(
  fun = function(Sigma, r, k) {
    elasticnet::spca(Sigma, K = r, sparse = "varnum",
                     para = rep(k, r), type = "Gram")$loadings
  },
  inputs = list(Sigma = pitprops, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (Sigma)")

## --- elasticnet, X input ---
b_enet_X <- bench_method(
  fun = function(X, r, k) {
    elasticnet::spca(X, K = r, sparse = "varnum",
                     para = rep(k, r), type = "predictor")$loadings
  },
  inputs = list(X = X_pseudo, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (X)")

## --- PMA ---
b_pma <- bench_method(
  fun = function(Sigma, r, sumabsv) {
    PMA::SPC(Sigma, sumabsv = sumabsv, K = r, orth = TRUE, trace = FALSE)$v
  },
  inputs = list(Sigma = pitprops, r = r, sumabsv = PMA_SUMABSV),
  packages = "PMA", reps = REPS, label = "PMA")

## --- sparsepca (normalised to unit norm) ---
b_spca <- bench_method(
  fun = function(Sigma, r, alpha) {
    unit_norm(sparsepca::spca(Sigma, k = r, alpha = alpha,
                              scale = TRUE, verbose = FALSE)$loadings)
  },
  inputs = list(Sigma = pitprops, r = r, alpha = SPARSEPCA_ALPHA),
  globals = list(unit_norm = unit_norm),
  packages = "sparsepca", reps = REPS, label = "sparsepca")

## --- amanpg ---
## type = 1: covariance matrix input; lambda1 is a vector of length r (one per PC)
b_amanpg <- bench_method(
  fun = function(Sigma, r, lambda1) {
    unit_norm(spca.amanpg(z = Sigma, lambda1 = lambda1, lambda2 = Inf,
                          k = r, type = 1, verbose = FALSE)$loadings)
  },
  inputs = list(Sigma = pitprops, r = r, lambda1 = AMANPG_LAMBDA1),
  globals = list(unit_norm = unit_norm),
  packages = "amanpg", reps = REPS, label = "amanpg")

## --- mixOmics::spca (keepX tuned to k per component; uses pseudo-data) ---
b_mixo <- bench_method(
  fun = function(X, r, k) {
    unit_norm(mixOmics::spca(X, ncomp = r, keepX = rep(k, r),
                             center = TRUE)$loadings$X)
  },
  inputs = list(X = X_pseudo, r = r, k = k),
  globals = list(unit_norm = unit_norm),
  packages = "mixOmics", reps = REPS, label = "mixOmics")

## --- nsprcomp (k nonzeros per PC; uses pseudo-data) ---
b_nspr <- bench_method(
  fun = function(X, r, k) {
    nsprcomp::nsprcomp(X, ncomp = r, k = rep(k, r),
                       nneg = FALSE, center = TRUE)$rotation
  },
  inputs = list(X = X_pseudo, r = r, k = k),
  packages = "nsprcomp", reps = REPS, label = "nsprcomp")

## --- nscumcomp (total nonzeros = r * k; uses pseudo-data) ---
## Retries on "Co-linear principal axes"; returns NULL if all attempts fail,
## which quality_table() propagates as NA rather than aborting the notebook.
b_nscum <- bench_method(
  fun = function(X, r, k, max_attempts) {
    for (attempt in seq_len(max_attempts)) {
      fit <- tryCatch(
        nsprcomp::nscumcomp(X, ncomp = r, k = r * k, nneg = FALSE, center = TRUE),
        error = function(e) e)
      if (!inherits(fit, "error")) return(fit$rotation)
    }
    NULL
  },
  inputs = list(X = X_pseudo, r = r, k = k, max_attempts = NSCUM_MAX_ATTEMPTS),
  packages = "nsprcomp", reps = REPS, label = "nscumcomp")

## --- Dense PCA ---
## pitprops is a correlation matrix, not a data matrix; use eigen() directly.
b_pca <- bench_method(
  fun = function(Sigma, r) {
    eigen(Sigma, symmetric = TRUE)$vectors[, seq_len(r)]
  },
  inputs = list(Sigma = pitprops, r = r),
  packages = character(), reps = REPS, label = "Dense PCA")


#### PYTHON METHODS ####

## --- sklearn SparsePCA ---
## Takes X_pseudo (raw data matrix). components_ is (r x p), transpose to (p x r).
## Initialising reticulate and importing sklearn costs on the order of 150 MB of
## NumPy/SciPy unrelated to the sparse-PCA solver, so it happens in setup():
## that cost lands in baseline_rss_mb and is excluded from mem_delta_mb.
b_sklearn <- bench_method(
  setup = function() {
    library(reticulate)
    use_python(Sys.which("python3"), required = TRUE)
    mod <- import("sklearn.decomposition")
    ## Throwaway fit: sklearn defers several imports (scipy.linalg, the
    ## dict-learning solver) until the first fit(), so warm them here rather
    ## than letting them land in the first timed repetition.
    mod$SparsePCA(n_components = 1L, alpha = 1)$fit(matrix(rnorm(20), 10, 2))
    mod
  },
  ## `seed` is supplied by the harness, one value per repetition; sklearn does
  ## not read R's RNG so it must be passed through as random_state.
  fun = function(X, r, alpha, seed, setup) {
    m <- setup$SparsePCA(n_components = r, alpha = alpha,
                         random_state = as.integer(seed))
    m$fit(X)
    unit_norm(t(m$components_))
  },
  inputs = list(X = X_pseudo, r = r, alpha = SKLEARN_ALPHA),
  globals = list(unit_norm = unit_norm),
  packages = character(), reps = REPS, label = "sklearn SparsePCA")


## ============================================================
## Console report and summary table
## ============================================================

benches <- list(b_mspca_S, b_mspca_X, b_enet_S, b_enet_X, b_pma, b_spca,
                b_amanpg, b_mixo, b_nspr, b_nscum, b_pca, b_sklearn)

cat("\n")
for (b in benches) report_bench(b, pitprops)

## FVE and orthogonality are evaluated on the published pitprops correlation
## matrix, not on the covariance of the pseudo-data.
results_pitprops <- cbind(bench_table(benches),
                          quality_table(benches, pitprops, r))

print(results_pitprops, row.names = FALSE)
write.csv(results_pitprops, "benchmarking/benchmarking_results_pitprops.csv",
					row.names = FALSE)
