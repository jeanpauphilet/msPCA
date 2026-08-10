# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for mspca() -- both the Sigma and X input paths.

set.seed(42)
p <- 8L
Sigma <- crossprod(matrix(rnorm(p * 20), 20, p))   # PSD p x p
Sigma <- Sigma / Sigma[1, 1]                        # scale so trace is informative
X     <- matrix(rnorm(50 * p), 50, p,
                dimnames = list(NULL, paste0("v", seq_len(p))))

# ---------- Sigma path -------------------------------------------------------

test_that("mspca (Sigma) returns an object of class 'mspca'", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_s3_class(res, "mspca")
})

test_that("mspca (Sigma) x_best has correct dimensions", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_equal(dim(res$x_best), c(p, 2L))
})

test_that("mspca (Sigma) each PC has at most k non-zero loadings", {
  ks  <- c(3L, 3L)
  res <- mspca(Sigma, r = 2L, ks = ks, verbose = FALSE)
  for (j in seq_len(2L)) {
    expect_lte(sum(res$x_best[, j] != 0), ks[j])
  }
})

test_that("mspca (Sigma) stores required output fields", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_true(all(c("x_best", "objective_value", "feasibility_violation",
                    "runtime", "variance_explained", "total_variance",
                    "inputType") %in% names(res)))
  expect_equal(res$inputType, "Sigma")
})

test_that("mspca (Sigma) loadings are unit-norm", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  norms <- colSums(res$x_best^2)
  expect_equal(norms, rep(1, 2L), tolerance = 1e-6)
})

test_that("mspca (Sigma) variance_explained is consistent with total_variance", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_true(all(res$variance_explained >= 0))
  expect_lte(sum(res$variance_explained), res$total_variance + 1e-9)
})

test_that("mspca (Sigma) objective_value is non-negative", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_gte(res$objective_value, 0)
})

test_that("mspca (Sigma) runtime is a non-negative number", {
  res <- mspca(Sigma, r = 2L, ks = c(3L, 3L), verbose = FALSE)
  expect_true(is.numeric(res$runtime) && res$runtime >= 0)
})

# ---------- X path -----------------------------------------------------------

test_that("mspca (X) returns an object of class 'mspca'", {
  res <- mspca(X, r = 2L, ks = c(3L, 3L), type = "X", verbose = FALSE)
  expect_s3_class(res, "mspca")
})

test_that("mspca (X) x_best has correct dimensions and rownames", {
  res <- mspca(X, r = 2L, ks = c(3L, 3L), type = "X", verbose = FALSE)
  expect_equal(dim(res$x_best), c(p, 2L))
  expect_equal(rownames(res$x_best), colnames(X))
})

test_that("mspca (X) stores X-specific metadata fields", {
  res <- mspca(X, r = 2L, ks = c(3L, 3L), type = "X", verbose = FALSE)
  expect_equal(res$inputType, "X")
  expect_equal(res$nObs, nrow(X))
  expect_equal(res$p, ncol(X))
  expect_true("center" %in% names(res))
  expect_true("scale"  %in% names(res))
  expect_true("divisor" %in% names(res))
})

test_that("mspca (X) each PC has at most k non-zero loadings", {
  ks  <- c(4L, 4L)
  res <- mspca(X, r = 2L, ks = ks, type = "X", verbose = FALSE)
  for (j in seq_len(2L)) {
    expect_lte(sum(res$x_best[, j] != 0), ks[j])
  }
})

test_that("mspca (X) with scale=FALSE and divisor='n' runs without error", {
  res <- mspca(X, r = 1L, ks = 3L, type = "X", scale = FALSE,
               divisor = "n", verbose = FALSE)
  expect_s3_class(res, "mspca")
})

# ---------- Real-data smoke test (mtcars) ------------------------------------

test_that("mspca (Sigma) produces stable results on mtcars correlation matrix", {
  Cm  <- cor(datasets::mtcars)
  res <- mspca(Cm, r = 2L, ks = c(4L, 4L), verbose = FALSE)
  expect_equal(dim(res$x_best), c(11L, 2L))
  # Fraction of variance explained should be positive and <=1
  fve <- sum(res$variance_explained) / res$total_variance
  expect_gt(fve, 0)
  expect_lte(fve, 1 + 1e-9)
})

test_that("mspca (X) produces stable results on mtcars raw matrix", {
  Xm  <- as.matrix(datasets::mtcars)
  res <- mspca(Xm, r = 2L, ks = c(4L, 4L), type = "X", verbose = FALSE)
  expect_equal(dim(res$x_best), c(11L, 2L))
})
