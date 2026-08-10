# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for tpm() -- both the Sigma and X input paths.

set.seed(7)
p <- 8L
Sigma <- crossprod(matrix(rnorm(p * 20), 20, p))
Sigma <- Sigma / Sigma[1, 1]
X     <- matrix(rnorm(40 * p), 40, p,
                dimnames = list(NULL, paste0("w", seq_len(p))))

# ---------- Sigma path -------------------------------------------------------

test_that("tpm (Sigma) returns an object of class 'tpm'", {
  res <- tpm(Sigma, k = 3L, verbose = FALSE)
  expect_s3_class(res, "tpm")
})

test_that("tpm (Sigma) x_best is a p x 1 matrix", {
  res <- tpm(Sigma, k = 3L, verbose = FALSE)
  expect_equal(dim(res$x_best), c(p, 1L))
})

test_that("tpm (Sigma) loading has at most k non-zero entries", {
  k   <- 3L
  res <- tpm(Sigma, k = k, verbose = FALSE)
  expect_lte(sum(res$x_best != 0), k)
})

test_that("tpm (Sigma) loading is unit-norm", {
  res <- tpm(Sigma, k = 4L, verbose = FALSE)
  expect_equal(sum(res$x_best^2), 1, tolerance = 1e-6)
})

test_that("tpm (Sigma) stores required fields", {
  res <- tpm(Sigma, k = 3L, verbose = FALSE)
  expect_true(all(c("x_best", "objective_value", "runtime", "inputType") %in% names(res)))
  expect_equal(res$inputType, "Sigma")
})

test_that("tpm (Sigma) objective_value matches x^T Sigma x", {
  res <- tpm(Sigma, k = 3L, verbose = FALSE)
  v   <- res$x_best
  expected <- drop(t(v) %*% Sigma %*% v)
  expect_equal(res$objective_value, expected, tolerance = 1e-5)
})

# ---------- X path -----------------------------------------------------------

test_that("tpm (X) returns an object of class 'tpm'", {
  res <- tpm(X, k = 3L, type = "X", verbose = FALSE)
  expect_s3_class(res, "tpm")
})

test_that("tpm (X) x_best has correct dimensions and rownames", {
  res <- tpm(X, k = 3L, type = "X", verbose = FALSE)
  expect_equal(dim(res$x_best), c(p, 1L))
  expect_equal(rownames(res$x_best), colnames(X))
})

test_that("tpm (X) stores X-specific metadata fields", {
  res <- tpm(X, k = 3L, type = "X", verbose = FALSE)
  expect_equal(res$inputType, "X")
  expect_equal(res$nObs, nrow(X))
  expect_equal(res$p, ncol(X))
})

test_that("tpm (X) k = p (dense) runs without error", {
  res <- tpm(X, k = p, type = "X", verbose = FALSE)
  expect_s3_class(res, "tpm")
})
