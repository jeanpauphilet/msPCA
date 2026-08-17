############################################################
## Benchmarking notebook: msPCA vs competing packages on Breast data
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: PMA::breast_data (p = 500 genes, n = 89 samples)
## Task:    r = 3 sparse PCs with k = 20 nonzeros each
##
## Methods compared: msPCA (Sigma and X input), elasticnet (Sigma and X
##                   input), PMA, sparsepca, amanpg, mixOmics, nsprcomp,
##                   nscumcomp, dense PCA, scikit-learn SparsePCA
## For PMA (L1 penalty) and sparsepca (alpha penalty), the tuning
## parameter is set to achieve ~20 nonzeros per component.
## For mixOmics::spca(), keepX = rep(k, r) matches the cardinality
## budget exactly. For nsprcomp::nsprcomp(), k = rep(k, r)
## gives the per-component nonzero budget; nscumcomp() uses the
## total budget r * k.
##
## This is the canonical benchmark from Witten et al. (2009);
## including it allows direct comparison on equal footing with
## the original PMA paper.
##
## Each method is timed AND memory-profiled by bench_utils.R, which runs
## it in a fresh R subprocess and records peak resident set size from
## getrusage(). See bench_utils.R for why R-level metrics (gc(),
## bench::mark) would bias the comparison in msPCA's favour.
##
## NOTE ON MEMORY: p = 500, n = 89. The correlation matrix S is about
## 2 MB of doubles; the data matrix X500 is about 0.36 MB. This is the
## smallest dataset at which the Sigma-versus-X difference is visible
## above RSS granularity; riboflavin (p = 4,088) makes it decisive.
############################################################

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("amanpg")
library("mixOmics")
library("nsprcomp")
library("reticulate")

source("benchmarking/bench_utils.R")

## Python: the scikit-learn comparison runs through reticulate. The Python
## environment it uses is configured in one place -- python_setup.R -- which
## declares `scikit-learn` for reticulate's managed (uv) environment and hands
## back the interpreter that was actually resolved. Read that file for the
## RETICULATE_PYTHON / RETICULATE_USE_MANAGED_VENV escape hatches, and for why
## a system python3 with scikit-learn installed is no longer enough under
## reticulate >= 1.41.
##
## A missing scikit-learn is NOT fatal. The sklearn row is dropped and every
## other method still runs, in the same spirit as the case-study script when the
## raw returns file is absent.
source("python_setup.R")

py <- sklearn_setup()
HAVE_SKLEARN <- py$available
py_bin <- py$python   # resolved interpreter, passed to the benchmark subprocesses
if (HAVE_SKLEARN) {
  sk_decomp <- import("sklearn.decomposition")
  cat("Using Python at ", py_bin, " (scikit-learn ", py$version, ")\n", sep = "")
} else {
  sk_decomp <- NULL
  warning(sklearn_missing_message(py), call. = FALSE)
}

## Helper: unit-normalise columns
unit_norm <- function(L) {
  for (j in seq_len(ncol(L))) {
    nm <- sqrt(sum(L[, j]^2))
    if (nm > 0) L[, j] <- L[, j] / nm
  }
  L
}

## Load breast cancer data (89 samples x 19,672 genes).
## PMA::download_breast_data() downloads the RNA-seq expression data.
## Note: this is a network download and happens once, in the parent process;
## only the reduced X500 / S are shipped to the benchmark subprocesses.
X <- t(PMA::download_breast_data(url = "https://tibshirani.su.domains/PMA/breastdata.rda")$rna)
vars  <- apply(X, 2, var)
X500  <- X[, order(vars, decreasing = TRUE)[1:500]]
S     <- cor(X500)                  # 500 x 500 correlation matrix
r <- 3L; k <- 20L

## Repetitions per method. Each repetition runs in its OWN process, so every
## one yields an independent peak-RSS sample; runtime, memory and FVE are all
## reported as medians over repetitions. Five is the minimum that gives the
## memory median any resistance to outliers -- peak RSS varies substantially
## run to run for methods that allocate heavily on R's heap (see bench_utils.R).
REPS <- 5L

############################################################
## Automatic tuning of the penalty-parameterized methods
##
## msPCA, elasticnet, mixOmics and nsprcomp accept an exact cardinality per
## component, so they need no tuning at all. PMA, sparsepca, amanpg and
## scikit-learn expose only a penalty magnitude or an l1 bound: for each we
## sweep a grid and keep the value whose realised sparsity is closest to the
## target of k = 20 nonzeros per component, breaking ties on FVE.
##
## Nothing below is hand-picked. Selected values therefore cannot go stale when
## an input, a grid or a package version changes -- which is exactly what
## happened to PMA's `sumabsv` when the input to SPC() was corrected.
##
## Note on amanpg: spca.amanpg() takes a VECTOR of per-component l1 penalties,
## and it needs one. With lambda2 = Inf the method soft-thresholds the columns
## of Sigma %*% A, whose scale follows the eigenvalues, so a single shared value
## empties the trailing components as the spectrum decays -- it produced
## (5,0,0) here and (8,4,2,0,0,0) on Pitprops. tune_parameter_vector() therefore
## tunes the entries by coordinate descent. Every method whose interface offers
## per-component control gets it, so this is parity, not preferential treatment.
##
## Tuning runs in the parent process, before any measurement, so it cannot
## contaminate the peak-RSS readings.
############################################################

TARGET_NNZ <- rep(k, r)

PMA_SUMABSV <- tune_parameter(
  fit    = function(sv) PMA::SPC(X500, sumabsv = sv, K = r, orth = TRUE,
                          trace = FALSE)$v,
  grid   = seq(2.0, 8.0, by = 0.1),
  target = TARGET_NNZ, C = S, label = "PMA::SPC sumabsv")

SPARSEPCA_ALPHA <- tune_parameter(
  fit    = function(a) unit_norm(sparsepca::spca(S, k = r, alpha = a,
                          verbose = FALSE, scale = FALSE)$loadings),
  grid   = seq(0.0020, 0.0080, by = 0.0005),
  target = TARGET_NNZ, C = S, label = "sparsepca::spca alpha")

AMANPG_LAMBDA1 <- tune_parameter_vector(
  fit    = function(lam) unit_norm(spca.amanpg(z = X500, lambda1 = lam,
                          lambda2 = Inf, k = r, type = 0, verbose = FALSE,
                          normalize = TRUE)$loadings),
  grid   = seq(0.1, 12, by = 0.1),
  target = TARGET_NNZ, C = S, label = "amanpg::spca.amanpg lambda1")

SKLEARN_ALPHA <- if (!HAVE_SKLEARN) NA_real_ else tune_parameter(
  fit    = function(a) {
             m <- sk_decomp$SparsePCA(n_components = r, alpha = a, random_state = 43L)
             m$fit(X500)
             unit_norm(t(m$components_))
           },
  grid   = seq(5, 20, by = 0.5),
  target = TARGET_NNZ, C = S, label = "sklearn SparsePCA alpha")


## nscumcomp is the ONE competing function with a knob on non-redundancy:
## `gamma` penalises divergence from orthonormality of the loadings, and its
## default is 0 -- no penalty at all. Reporting its orthogonality violation
## while never asking it to be feasible would not be a fair comparison, so we
## tune gamma to the smallest violation msPCA's own tolerance allows, keeping
## the best FVE among the values that get there. Sparsity is unaffected:
## nscumcomp takes its cardinality budget separately.
##
## The grid runs to 1e8 deliberately. An earlier ceiling of 1e3 was selected at
## the boundary on both mtcars and Pitprops without reaching the 1e-4 target,
## which left it unknown whether nscumcomp plateaus above the tolerance or
## simply needed a larger penalty. Running well past the point of diminishing
## returns settles that; the tuning CSV flags `at_grid_edge` if it happens again.
NSCUM_GAMMA <- tune_for_feasibility(
  fit   = function(g) nsprcomp::nscumcomp(X500, ncomp = r, k = r * k, nneg = FALSE,
                        center = TRUE, gamma = g)$rotation,
  grid  = c(0, 10^seq(-3, 8, by = 0.5)),
  C     = S, ctype = 0, tol = 1e-4,
  label = "nsprcomp::nscumcomp gamma")


## Record the selected values next to the results, so the article's "Tuning:"
## lines can be read off a file rather than transcribed by hand. The table
## also carries the realised sparsity and FVE at each selected value, so the
## choice is auditable rather than just asserted.
write.csv(tuning_table(), "benchmarking/benchmarking_tuning_breast.csv",
          row.names = FALSE)

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
## Sigma is applied implicitly as X^T(X v)/(n-1), so the 2 MB correlation
## matrix is never allocated. Compare working_set_mb against the row above.
b_mspca_X <- bench_method(
  fun = function(X, r, k) {
    mspca(X, type = "X", r = r, ks = rep(k, r), verbose = FALSE,
          feasibilityConstraintType = 0)$x_best
  },
  inputs = list(X = X500, r = r, k = k),
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
  inputs = list(X = X500, r = r, k = k),
  packages = "elasticnet", reps = REPS, label = "elasticnet (X)")

## --- PMA ---
b_pma <- bench_method(
  fun = function(X, r, sumabsv) {
    PMA::SPC(X, sumabsv = sumabsv, K = r, orth = TRUE, trace = FALSE)$v
  },
  inputs = list(X = X500, r = r, sumabsv = PMA_SUMABSV),
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
## type = 0: data matrix input; lambda1 is a vector of length r (one per PC)
b_amanpg <- bench_method(
  fun = function(X, r, lambda1) {
    unit_norm(spca.amanpg(z = X, lambda1 = lambda1, lambda2 = Inf,
                          k = r, type = 0, verbose = FALSE,
                          normalize = TRUE)$loadings)
  },
  inputs = list(X = X500, r = r, lambda1 = AMANPG_LAMBDA1),
  globals = list(unit_norm = unit_norm),
  packages = "amanpg", reps = REPS, label = "amanpg")

## --- mixOmics::spca (keepX = rep(k, r) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
b_mixo <- bench_method(
  fun = function(X, r, k) {
    unit_norm(mixOmics::spca(X, ncomp = r, keepX = rep(k, r),
                             center = TRUE, scale = TRUE)$loadings$X)
  },
  inputs = list(X = X500, r = r, k = k),
  globals = list(unit_norm = unit_norm),
  packages = "mixOmics", reps = REPS, label = "mixOmics")

## --- nsprcomp (k nonzeros per PC; basic deflation) ---
b_nspr <- bench_method(
  fun = function(X, r, k) {
    nsprcomp::nsprcomp(X, ncomp = r, k = rep(k, r),
                       nneg = FALSE, center = TRUE)$rotation
  },
  inputs = list(X = X500, r = r, k = k),
  packages = "nsprcomp", reps = REPS, label = "nsprcomp")

## --- nscumcomp (total nonzeros = r * k; joint optimization) ---
## nscumcomp() may fail on high-correlation data ("Co-linear principal axes").
## Returning NULL is propagated as NA by quality_table() rather than aborting.
b_nscum <- bench_method(
  fun = function(X, r, k, gamma) {
    fit <- tryCatch(
      nsprcomp::nscumcomp(X, ncomp = r, k = r * k, nneg = FALSE, center = TRUE,
                          gamma = gamma),
      error = function(e) NULL)
    if (is.null(fit)) NULL else fit$rotation
  },
  inputs = list(X = X500, r = r, k = k, gamma = NSCUM_GAMMA),
  packages = "nsprcomp", reps = REPS, label = "nscumcomp")

## --- Dense PCA ---
b_pca <- bench_method(
  fun = function(X, r) prcomp(X, scale. = TRUE)$rotation[, seq_len(r)],
  inputs = list(X = X500, r = r),
  packages = character(), reps = REPS, label = "Dense PCA")


#### PYTHON METHODS ####

## --- sklearn SparsePCA ---
## Takes X500 (raw data matrix). components_ is (r x p), transpose to (p x r).
## Initialising reticulate and importing sklearn costs on the order of 150 MB of
## NumPy/SciPy unrelated to the sparse-PCA solver, so it happens in setup():
## that cost lands in baseline_rss_mb and is excluded from mem_delta_mb.
b_sklearn <- if (!HAVE_SKLEARN) NULL else bench_method(
  setup = function() {
    library(reticulate)
    if (nzchar(py_bin)) try(use_python(py_bin, required = TRUE), silent = TRUE)
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
  inputs = list(X = X500, r = r, alpha = SKLEARN_ALPHA),
  globals = list(unit_norm = unit_norm, py_bin = py_bin),
  packages = character(), reps = REPS, label = "sklearn SparsePCA")


############################################################
## Console report and summary table
############################################################

benches <- Filter(Negate(is.null), list(b_mspca_S, b_mspca_X, b_enet_S, b_enet_X, b_pma, b_spca,
                b_amanpg, b_mixo, b_nspr, b_nscum, b_pca, b_sklearn))

cat("\n")
for (b in benches) report_bench(b, S)

results <- cbind(bench_table(benches), quality_table(benches, S, r))

print(results, row.names = FALSE)
write.csv(results, "benchmarking/benchmarking_results_breast.csv", row.names = FALSE)
