############################################################
## Benchmarking notebook: msPCA vs competing packages on mtcars
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: datasets::mtcars (p = 11 variables, n = 32 observations)
## Task:    r = 3 sparse PCs with k = 4 nonzeros each
##
## Methods compared: msPCA (Sigma and X input), elasticnet (Sigma and X
##                   input), PMA, sparsepca, amanpg, mixOmics, nsprcomp,
##                   nscumcomp, dense PCA, scikit-learn SparsePCA
## For PMA (L1 penalty) and sparsepca (alpha penalty), the tuning
## parameter is set to achieve ~4 nonzeros per component.
## For mixOmics::spca(), keepX = rep(k, r) matches the cardinality
## budget exactly. For nsprcomp::nscumcomp(), k = r * k_per_pc
## gives the total nonzero budget across all components.
##
## Each method is timed AND memory-profiled by bench_utils.R, which runs
## it in a fresh R subprocess and records peak resident set size from
## getrusage(). See bench_utils.R for why R-level metrics (gc(),
## bench::mark) would bias the comparison in msPCA's favour.
##
## NOTE ON MEMORY: mtcars is p = 11, so the numerical working set of
## every method is a few kilobytes. msPCA, PMA, mixOmics, amanpg,
## nsprcomp, sklearn and dense PCA accordingly report close to zero.
## elasticnet, sparsepca and nscumcomp do not: they carry a fixed
## working-memory floor of roughly 10-25 MB that is present even at
## this size, on top of which cost scales with p. A single-repetition
## check confirmed that floor is genuine allocation rather than
## first-call package loading, which accounts for only about 4 MB of
## it. The datasets where memory scaling dominates the floor are
## notebook_breast.R (p = 500) and notebook_riboflavin.R (p = 4,088),
## where Sigma is 2 MB and 128 MB respectively.
############################################################

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("amanpg")
library("mixOmics")
library("nsprcomp")
library("datasets")
library("reticulate")

source("benchmarking/bench_utils.R")

py_bin <- Sys.which("python3")
use_python(py_bin, required = TRUE)
sk_decomp <- import("sklearn.decomposition")

## Helper: unit-normalise columns (matches existing pattern for sparsepca/mixOmics)
unit_norm <- function(L) {
  for (j in seq_len(ncol(L))) {
    nm <- sqrt(sum(L[, j]^2))
    if (nm > 0) L[, j] <- L[, j] / nm
  }
  L
}

S <- cor(datasets::mtcars)
X_mtcars <- scale(datasets::mtcars)
r <- 3L; k <- 4L

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
PMA_SUMABSV     <- 1.69
SPARSEPCA_ALPHA <- 0.004
AMANPG_LAMBDA1  <- c(10, 3, 0.01)
SKLEARN_ALPHA   <- 3.2


############################################################
## Tuning sweeps (untimed, unprofiled)
############################################################

if (RUN_TUNING_SWEEPS) {
  ## --- PMA (sumabsv tuned for ~k NNZ per PC) ---
  cat("PMA sumabsv sweep:\n")
  for (sv in seq(1.6, 1.8, by = 0.01)) {
    pma_try <- PMA::SPC(S, sumabsv = sv, K = r, orth = TRUE, trace = FALSE)
    cat("  sumabsv =", sv,
        "| NNZ:", paste(colSums(abs(pma_try$v) > 0), collapse = " "),
        "| FVE:", round(fraction_variance_explained(S, pma_try$v), 4), "\n")
  }

  ## --- sparsepca (alpha tuned for ~k NNZ per PC) ---
  cat("sparsepca alpha sweep:\n")
  for (a in seq(0.001, 0.01, by = 0.001)) {
    sp_try <- unit_norm(
      sparsepca::spca(S, k = r, alpha = a, verbose = FALSE, scale = FALSE)$loadings)
    cat("  alpha =", a,
        "| NNZ:", paste(colSums(abs(sp_try) > 0), collapse = " "),
        "| FVE:", round(fraction_variance_explained(S, sp_try), 4), "\n")
  }

  ## --- sklearn SparsePCA (alpha tuned for ~k NNZ per PC) ---
  cat("sklearn alpha sweep:\n")
  for (a in seq(3, 4, by = 0.1)) {
    m <- sk_decomp$SparsePCA(n_components = r, alpha = a, random_state = 43L)
    m$fit(X_mtcars)
    L <- unit_norm(t(m$components_))
    cat("  alpha =", a, "| NNZ:", paste(colSums(abs(L) > 0), collapse = " "), "\n")
  }
}


############################################################
## Measured runs -- one fresh subprocess per method.
## Each fun() returns a plain p x r loadings matrix.
############################################################

#### R METHODS ####

## --- msPCA, Sigma input (dense p x p covariance operator) ---
b_mspca_S <- bench_method(
  fun = function(S, r, k) {
    mspca(S, r = r, ks = rep(k, r), verbose = FALSE,
          feasibilityConstraintType = 0)$x_best
  },
  inputs = list(S = S, r = r, k = k),
  packages = "msPCA", reps = REPS, label = "msPCA (Sigma)")

## --- msPCA, X input (matrix-free Gram operator, never forms p x p) ---
## Same estimator, different internal covariance operator: Sigma is applied
## implicitly as X^T(X v)/(n-1). This is the row that demonstrates the O(np)
## rather than O(p^2) memory footprint; on mtcars the difference is nil, but
## the row must exist here for the tables to line up across datasets.
b_mspca_X <- bench_method(
  fun = function(X, r, k) {
    mspca(X, type = "X", r = r, ks = rep(k, r), verbose = FALSE,
          feasibilityConstraintType = 0)$x_best
  },
  inputs = list(X = X_mtcars, r = r, k = k),
  packages = "msPCA", reps = REPS, label = "msPCA (X)")

## --- elasticnet, Sigma input ---
b_enet_S <- bench_method(
  fun = function(S, r, k) {
    elasticnet::spca(S, K = r, sparse = "varnum",
                     para = rep(k, r), type = "Gram")$loadings
  },
  inputs = list(S = S, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (Sigma)")

## --- elasticnet, X input ---
b_enet_X <- bench_method(
  fun = function(X, r, k) {
    elasticnet::spca(X, K = r, sparse = "varnum",
                     para = rep(k, r), type = "predictor")$loadings
  },
  inputs = list(X = X_mtcars, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (X)")

## --- PMA ---
b_pma <- bench_method(
  fun = function(S, r, sumabsv) {
    PMA::SPC(S, sumabsv = sumabsv, K = r, orth = TRUE, trace = FALSE)$v
  },
  inputs = list(S = S, r = r, sumabsv = PMA_SUMABSV),
  packages = "PMA", reps = REPS, label = "PMA")

## --- sparsepca (normalised to unit norm) ---
b_spca <- bench_method(
  fun = function(S, r, alpha) {
    unit_norm(sparsepca::spca(S, k = r, alpha = alpha,
                              verbose = FALSE, scale = FALSE)$loadings)
  },
  inputs = list(S = S, r = r, alpha = SPARSEPCA_ALPHA),
  globals = list(unit_norm = unit_norm),
  packages = "sparsepca", reps = REPS, label = "sparsepca")

## --- amanpg ---
## type = 1: covariance matrix input; lambda1 is a vector of length r (one per PC)
b_amanpg <- bench_method(
  fun = function(S, r, lambda1) {
    unit_norm(spca.amanpg(z = S, lambda1 = lambda1, lambda2 = Inf,
                          k = r, type = 1, verbose = FALSE)$loadings)
  },
  inputs = list(S = S, r = r, lambda1 = AMANPG_LAMBDA1),
  globals = list(unit_norm = unit_norm),
  packages = "amanpg", reps = REPS, label = "amanpg")

## --- mixOmics::spca (keepX = rep(k, r) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
b_mixo <- bench_method(
  fun = function(X, r, k) {
    unit_norm(mixOmics::spca(X, ncomp = r, keepX = rep(k, r),
                             center = TRUE, scale = TRUE)$loadings$X)
  },
  inputs = list(X = X_mtcars, r = r, k = k),
  globals = list(unit_norm = unit_norm),
  packages = "mixOmics", reps = REPS, label = "mixOmics")

## --- nsprcomp (k nonzeros per PC; basic deflation) ---
b_nspr <- bench_method(
  fun = function(X, r, k) {
    nsprcomp::nsprcomp(X, ncomp = r, k = rep(k, r), nneg = FALSE,
                       center = TRUE, scale. = TRUE)$rotation
  },
  inputs = list(X = X_mtcars, r = r, k = k),
  packages = "nsprcomp", reps = REPS, label = "nsprcomp")

## --- nscumcomp (total budget r * k; joint optimisation) ---
b_nscum <- bench_method(
  fun = function(X, r, k) {
    nsprcomp::nscumcomp(X, ncomp = r, k = k * r, nneg = FALSE,
                        center = TRUE, scale. = TRUE)$rotation
  },
  inputs = list(X = X_mtcars, r = r, k = k),
  packages = "nsprcomp", reps = REPS, label = "nscumcomp")

## --- Dense PCA ---
b_pca <- bench_method(
  fun = function(X, r) prcomp(X, scale. = TRUE)$rotation[, seq_len(r)],
  inputs = list(X = datasets::mtcars, r = r),
  packages = character(), reps = REPS, label = "Dense PCA")


#### PYTHON METHODS ####

## --- sklearn SparsePCA ---
## Takes X (not Sigma). components_ is (r x p) so we transpose to (p x r).
## Initialising reticulate and importing sklearn costs on the order of 150 MB
## of NumPy/SciPy that has nothing to do with the sparse-PCA solver, so it is
## done in setup(): that cost lands in baseline_rss_mb and is excluded from
## mem_delta_mb, which then measures the fit itself. peak_rss_mb still shows
## the full footprint of the R + Python stack for readers who want it.
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
  inputs = list(X = X_mtcars, r = r, alpha = SKLEARN_ALPHA),
  globals = list(unit_norm = unit_norm),
  packages = character(), reps = REPS, label = "sklearn SparsePCA")


############################################################
## Console report and summary table
############################################################

benches <- list(b_mspca_S, b_mspca_X, b_enet_S, b_enet_X, b_pma, b_spca,
                b_amanpg, b_mixo, b_nspr, b_nscum, b_pca, b_sklearn)

cat("\n")
for (b in benches) report_bench(b, S)

## bench_table() supplies method, runtime_s (median of REPS), working_set_mb,
## input_mb, mem_delta_mb, peak_rss_mb, baseline_rss_mb, gc_max_mb.
## quality_table() supplies nnz_pc1..r, fve, orth from the returned loadings.
results <- cbind(bench_table(benches), quality_table(benches, S, r))

print(results, row.names = FALSE)
write.csv(results, "benchmarking/benchmarking_results_mtcars.csv", row.names = FALSE)
