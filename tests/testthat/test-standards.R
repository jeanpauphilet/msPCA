# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests written against the rOpenSci Statistical Software standards: correctness,
# parameter recovery, algorithm performance, edge conditions, noise
# susceptibility, and unsupervised-learning input/output properties.

#' @srrstats {G5.5} All tests that involve randomness call `set.seed()` first.

# A spiked covariance model: X = z v' + noise, with a 5-sparse unit-norm spike v.
spiked_data <- function(seed, n = 400L, p = 30L, k = 5L, strength = 3) {
  set.seed(seed)
  v <- numeric(p)
  v[seq_len(k)] <- 1 / sqrt(k)
  list(X = outer(rnorm(n) * strength, v) + matrix(rnorm(n * p), n, p), v = v)
}

# ---------- Correctness (G5.4) -----------------------------------------------

#' @srrstats {G5.0, G5.1} Tests use standard data sets with known properties:
#'   `mtcars` (base R) and the exported `snp500` data set.
#' @srrstats {G5.4, G5.4a} Correctness against a reference. msPCA is a new method
#'   with no prior R implementation of its multi-component formulation, so the
#'   single-component problem is checked against an exact reference: exhaustive
#'   enumeration of all supports, and the leading eigenvector when `k = p`.
test_that("single-component solutions match an exhaustive search over supports", {
  set.seed(5)
  p <- 7L
  k <- 3L
  C <- cor(matrix(rnorm(60 * p), 60, p))
  best <- max(apply(utils::combn(p, k), 2L, function(ix) max(eigen(C[ix, ix])$values)))
  set.seed(1)
  expect_equal(mspca(C, r = 1L, ks = k, verbose = FALSE)$variance_explained, best,
               tolerance = 1e-6)
  set.seed(1)
  expect_equal(tpm(C, k = k, verbose = FALSE)$objective_value, best, tolerance = 1e-6)
})

test_that("a dense component (k = p) recovers the leading eigenvalue", {
  set.seed(6)
  p <- 8L
  C <- cor(matrix(rnorm(80 * p), 80, p))
  set.seed(1)
  expect_equal(tpm(C, k = p, verbose = FALSE)$objective_value,
               max(eigen(C)$values), tolerance = 1e-6)
})

test_that("the Sigma and X inputs agree on the same problem", {
  X <- as.matrix(mtcars)
  set.seed(1)
  a <- mspca(cor(X), r = 2L, ks = c(4L, 4L), verbose = FALSE)
  set.seed(1)
  b <- mspca(X, r = 2L, ks = c(4L, 4L), type = "X", verbose = FALSE)
  expect_equal(a$variance_explained, b$variance_explained, tolerance = 1e-6)
})

# ---------- Parameter recovery (G5.6) ----------------------------------------

#' @srrstats {G5.6, G5.6a, G5.6b} Parameter recovery: the planted sparse spike is
#'   recovered (support exactly, direction up to a tolerance) under several
#'   random seeds.
test_that("a planted sparse spike is recovered across random seeds", {
  for (seed in 1:5) {
    d <- spiked_data(seed)
    set.seed(seed)
    fit <- mspca(d$X, r = 1L, ks = 5L, type = "X", verbose = FALSE)
    expect_setequal(which(fit$x_best[, 1] != 0), 1:5)
    expect_gt(abs(sum(fit$x_best[, 1] * d$v)), 0.99)
  }
})

# ---------- Algorithm performance (G5.7) -------------------------------------

#' @srrstats {G5.7} Algorithm performance: recovery of the planted direction
#'   improves as the number of observations grows.
test_that("recovery of the planted direction improves with more observations", {
  err <- function(n) {
    mean(vapply(1:4, function(seed) {
      d <- spiked_data(seed, n = n, strength = 1.5)
      set.seed(seed)
      fit <- mspca(d$X, r = 1L, ks = 5L, type = "X", verbose = FALSE)
      1 - abs(sum(fit$x_best[, 1] * d$v))
    }, numeric(1)))
  }
  expect_lt(err(2000L), err(60L))
})

# ---------- Edge conditions (G5.8) -------------------------------------------

#' @srrstats {G5.8, G5.8a} Zero-length data is rejected with an error.
test_that("zero-length input is rejected", {
  expect_error(mspca(matrix(numeric(0), 0, 0), r = 1L, ks = 1L, verbose = FALSE),
               "at least one row and column")
  expect_error(mspca(matrix(numeric(0), 0, 5), r = 1L, ks = 1L, type = "X",
                     verbose = FALSE), "at least 2 rows")
  expect_error(mspca(matrix(numeric(0), 5, 0), r = 1L, ks = 1L, type = "X",
                     verbose = FALSE), "at least one column")
  expect_error(tpm(matrix(numeric(0), 0, 0), k = 1L, verbose = FALSE),
               "at least one row and column")
})

#' @srrstats {G5.8b} Data of unsupported types (character, complex, list columns)
#'   is rejected.
test_that("unsupported data types are rejected", {
  expect_error(mspca(matrix("a", 3, 3), r = 1L, ks = 1L, verbose = FALSE),
               "numeric matrix")
  expect_error(mspca(matrix(1 + 0i, 3, 3), r = 1L, ks = 1L, verbose = FALSE),
               "numeric matrix")
  expect_error(mspca(data.frame(a = 1:5, b = letters[1:5]), r = 1L, ks = 1L,
                     type = "X", verbose = FALSE), "non-numeric")
  lst <- data.frame(a = 1:3)
  lst$b <- list(1, 2, 3)
  expect_error(mspca(lst, r = 1L, ks = 1L, type = "X", verbose = FALSE),
               "non-numeric")
  expect_error(mspca(list(1, 2), r = 1L, ks = 1L, verbose = FALSE), "numeric matrix")
})

#' @srrstats {G5.8c} Data with all-missing fields or identical columns.
test_that("all-NA and constant columns are handled with a clear condition", {
  expect_error(mspca(matrix(NA_real_, 6, 6), r = 1L, ks = 1L, verbose = FALSE),
               "non-finite")
  expect_error(mspca(matrix(NA_real_, 6, 6), r = 1L, ks = 1L, type = "X",
                     verbose = FALSE), "non-finite")
  expect_error(mspca(matrix(1, 10, 5), r = 1L, ks = 2L, type = "X", verbose = FALSE),
               "zero-variance")
  # constant up to rounding error is also caught
  Xc <- matrix(rnorm(50), 10, 5)
  Xc[, 2] <- 0.1 + 1e-17
  expect_error(mspca(Xc, r = 1L, ks = 2L, type = "X", verbose = FALSE),
               "zero-variance")
  # without scaling, a constant matrix is a legitimate (zero-variance) input:
  # the solver reports that no informative solution exists.
  expect_warning(
    mspca(matrix(1, 10, 5), r = 1L, ks = 2L, type = "X", scale = FALSE, verbose = FALSE),
    "feasible"
  )
})

#' @srrstats {G5.8d} Data outside the scope of the algorithm: a sparsity level
#'   larger than the dimension warns and is truncated; more variables than
#'   observations is supported (and is the case the raw-data path targets).
test_that("out-of-scope sparsity and p > n are handled", {
  expect_warning(mspca(cor(mtcars), r = 1L, ks = 20L, verbose = FALSE),
                 "exceeds the dimension")
  set.seed(2)
  X <- matrix(rnorm(5 * 12), 5, 12)
  fit <- mspca(X, r = 1L, ks = 3L, type = "X", verbose = FALSE)
  expect_equal(dim(fit$x_best), c(12L, 1L))
  expect_lte(sum(fit$x_best != 0), 3L)
})

# ---------- Noise susceptibility (G5.9) --------------------------------------

#' @srrstats {G5.9, G5.9a} Adding noise at the scale of `.Machine$double.eps` to
#'   the input does not change the solution.
test_that("trivial noise does not change the solution", {
  set.seed(2)
  S <- cor(matrix(rnorm(200 * 15), 200, 15))
  set.seed(3)
  E <- matrix(rnorm(15 * 15) * .Machine$double.eps, 15, 15)
  Sn <- S + (E + t(E)) / 2
  set.seed(1)
  a <- mspca(S, r = 2L, ks = c(4L, 4L), verbose = FALSE)
  set.seed(1)
  b <- mspca(Sn, r = 2L, ks = c(4L, 4L), verbose = FALSE)
  expect_equal(abs(a$x_best), abs(b$x_best), tolerance = 1e-6)
  expect_equal(a$variance_explained, b$variance_explained, tolerance = 1e-6)
})

#' @srrstats {G5.9b} Different random seeds give results of the same quality.
test_that("results are stable across random seeds", {
  set.seed(2)
  S <- cor(matrix(rnorm(200 * 15), 200, 15))
  obj <- vapply(1:4, function(s) {
    set.seed(s)
    mspca(S, r = 2L, ks = c(4L, 4L), verbose = FALSE)$objective_value
  }, numeric(1))
  expect_lt((max(obj) - min(obj)) / max(obj), 0.05)
})

# ---------- Output properties (G5.3, UL3.1, UL4.2) ---------------------------

#' @srrstats {G5.3} Returned numeric fields contain no missing or undefined values.
test_that("returned objects contain no NA, NaN or Inf in their numeric results", {
  set.seed(1)
  fit <- mspca(cor(mtcars), r = 2L, ks = c(4L, 4L), verbose = FALSE)
  for (nm in c("x_best", "variance_explained", "total_variance",
               "objective_value", "feasibility_violation")) {
    expect_true(all(is.finite(fit[[nm]])), info = nm)
  }
  fit_t <- tpm(cor(mtcars), k = 4L, verbose = FALSE)
  expect_true(all(is.finite(fit_t$x_best)))
  expect_true(is.finite(fit_t$objective_value))
})

#' @srrstats {UL3.1} Components are labelled in decreasing order of variance
#'   explained.
test_that("principal components are ordered by decreasing variance explained", {
  for (seed in 1:5) {
    set.seed(seed)
    X <- matrix(rnorm(100 * 12), 100, 12)
    fit <- mspca(X, r = 3L, ks = c(3L, 3L, 3L), type = "X", verbose = FALSE)
    expect_false(is.unsorted(rev(fit$variance_explained)))
  }
})

#' @srrstats {UL4.2} The returned object carries all parameters controlling the
#'   algorithm.
test_that("the returned objects record the algorithmic parameters", {
  set.seed(1)
  fit <- mspca(cor(mtcars), r = 2L, ks = c(4L, 3L), maxIter = 50, verbose = FALSE,
               timeLimitTPM = 3, maxRestartTPM = 4, minRestartTPM = 2)
  expect_equal(fit$control$ks, c(4L, 3L))
  expect_equal(fit$control$maxIter, 50)
  expect_equal(fit$control$maxRestartTPM, 4)
  expect_equal(fit$feasibilityConstraintType, 0L)
  fit_t <- tpm(cor(mtcars), k = 4L, maxIter = 30, timeLimit = 2, verbose = FALSE)
  expect_equal(fit_t$control, list(k = 4L, maxIter = 30, timeLimit = 2))
})

# ---------- Input handling (G2.7, G2.11, G2.12, G3.1, UL1.3, UL7.x) ----------

#' @srrstats {G2.7, G2.8} Numeric data.frames are accepted alongside matrices and
#'   converted to a matrix once, at the entry point.
test_that("an all-numeric data.frame gives the same fit as the equivalent matrix", {
  set.seed(1)
  a <- mspca(mtcars, r = 2L, ks = c(4L, 4L), type = "X", verbose = FALSE)
  set.seed(1)
  b <- mspca(as.matrix(mtcars), r = 2L, ks = c(4L, 4L), type = "X", verbose = FALSE)
  expect_equal(a$x_best, b$x_best)
  expect_equal(rownames(a$x_best), colnames(mtcars))
  set.seed(1)
  s <- mspca(as.data.frame(cor(mtcars)), r = 1L, ks = 3L, verbose = FALSE)
  expect_equal(dim(s$x_best), c(ncol(mtcars), 1L))
})

#' @srrstats {G3.1, G3.1a} Any covariance estimator can be used, by supplying its
#'   output as `type = "Sigma"`; this is documented in the examples.
test_that("alternative covariance estimators can be supplied through type = 'Sigma'", {
  set.seed(1)
  fit <- mspca(cor(mtcars, method = "spearman"), r = 2L, ks = c(4L, 4L),
               verbose = FALSE)
  expect_s3_class(fit, "mspca")
  X <- as.matrix(mtcars)
  X[3, 2] <- NA
  S <- cov(X, use = "pairwise.complete.obs")
  expect_error(mspca(X, r = 1L, ks = 2L, type = "X", verbose = FALSE), "non-finite")
  expect_s3_class(mspca(S, r = 1L, ks = 2L, verbose = FALSE), "mspca")
})

#' @srrstats {UL7.0} Inappropriate types of input are rejected with informative
#'   errors (see also `test-validators.R`).
#' @srrstats {UL7.1} Violating an assumption on the input gives invalid output,
#'   which the tests show explicitly: with the PSD check disabled, a matrix that
#'   is not a covariance matrix yields a negative "variance explained".
test_that("a non-PSD input with checkPSD = FALSE yields an invalid, negative variance", {
  set.seed(3)
  S <- crossprod(matrix(rnorm(8 * 40), 40, 8))
  fit <- suppressWarnings(mspca(-S, r = 1L, ks = 3L, checkPSD = FALSE, verbose = FALSE))
  expect_lt(fit$variance_explained, 0)
})

#' @srrstats {UL1.3, UL7.3} Variable names on the input are propagated to the
#'   rows of the loadings, in both input modes.
test_that("variable names are propagated to the loadings", {
  set.seed(1)
  a <- mspca(cor(mtcars), r = 1L, ks = 3L, verbose = FALSE)
  b <- mspca(as.matrix(mtcars), r = 1L, ks = 3L, type = "X", verbose = FALSE)
  expect_equal(rownames(a$x_best), colnames(mtcars))
  expect_equal(rownames(b$x_best), colnames(mtcars))
})

#' @srrstats {G2.3, G2.3a, G2.3b} Character parameters are matched with
#'   `match.arg()` after `tolower()`, so they are case-insensitive.
test_that("character options are case-insensitive", {
  set.seed(1)
  X <- as.matrix(mtcars)
  a <- mspca(X, r = 1L, ks = 3L, type = "X", divisor = "N", verbose = FALSE)
  b <- mspca(X, r = 1L, ks = 3L, type = "x", divisor = "n", verbose = FALSE)
  expect_equal(a$divisor, "n")
  expect_equal(b$divisor, "n")
  expect_error(mspca(X, r = 1L, ks = 3L, type = "X", divisor = "m", verbose = FALSE),
               "should be one of")
})
