# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for input validation: .validate_sigma_matrix() and .validate_x_matrix()
# exercised indirectly through mspca() and tpm().

set.seed(99)
p <- 6L

# A valid PSD Sigma for reference
Sigma_ok <- crossprod(matrix(rnorm(p * 20), 20, p))

# ---------- Sigma-path validators --------------------------------------------

test_that("non-numeric M with type='Sigma' raises an error", {
  bad <- matrix(as.character(1:(p * p)), p, p)
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "numeric matrix")
})

test_that("non-square M with type='Sigma' raises an error", {
  bad <- matrix(rnorm((p - 1) * p), p - 1, p)
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "square")
})

test_that("M with NA values with type='Sigma' raises an error", {
  bad <- Sigma_ok
  bad[2, 2] <- NA_real_
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "non-finite")
})

test_that("M with Inf values with type='Sigma' raises an error", {
  bad <- Sigma_ok
  bad[1, 3] <- Inf
  bad[3, 1] <- Inf
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "non-finite")
})

test_that("non-symmetric M with type='Sigma' raises an error", {
  bad      <- Sigma_ok
  bad[1, 2] <- bad[1, 2] + 1
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "symmetric")
})

test_that("non-PSD M with type='Sigma' raises an error by default", {
  bad      <- Sigma_ok
  diag(bad) <- diag(bad) - max(diag(bad)) - 1   # force negative eigenvalues
  bad      <- (bad + t(bad)) / 2                 # keep symmetric
  expect_error(mspca(bad, r = 1L, ks = 2L, verbose = FALSE),
               regexp = "positive semidefinite")
})

test_that("non-PSD M is accepted when checkPSD = FALSE", {
  bad      <- -Sigma_ok   # all eigenvalues negative
  bad      <- (bad + t(bad)) / 2
  expect_error(
    mspca(bad, r = 1L, ks = 2L, verbose = FALSE, checkPSD = FALSE),
    NA   # no error expected
  )
})

# ---------- X-path validators ------------------------------------------------

test_that("non-numeric M with type='X' raises an error", {
  bad <- matrix(as.character(1:(10 * p)), 10, p)
  expect_error(mspca(bad, r = 1L, ks = 2L, type = "X", verbose = FALSE),
               regexp = "numeric matrix")
})

test_that("M with NA values with type='X' raises an error", {
  bad <- matrix(rnorm(30 * p), 30, p)
  bad[5, 3] <- NA_real_
  expect_error(mspca(bad, r = 1L, ks = 2L, type = "X", verbose = FALSE),
               regexp = "non-finite")
})

test_that("M with too few rows (n < 2) with type='X' raises an error", {
  bad <- matrix(rnorm(p), 1L, p)
  expect_error(mspca(bad, r = 1L, ks = 2L, type = "X", verbose = FALSE),
               regexp = "at least 2")
})

test_that("zero-variance column with scale=TRUE raises an error", {
  bad        <- matrix(rnorm(30 * p), 30, p)
  bad[, 1]   <- 0
  expect_error(mspca(bad, r = 1L, ks = 2L, type = "X",
                     scale = TRUE, verbose = FALSE),
               regexp = "zero-variance")
})

# ---------- tpm() uses the same validators -----------------------------------

test_that("tpm() (Sigma) rejects non-square input", {
  bad <- matrix(rnorm((p - 1) * p), p - 1, p)
  expect_error(tpm(bad, k = 2L, verbose = FALSE), regexp = "square")
})

test_that("tpm() (X) rejects non-finite input", {
  bad <- matrix(rnorm(20 * p), 20, p)
  bad[1, 1] <- NaN
  expect_error(tpm(bad, k = 2L, type = "X", verbose = FALSE),
               regexp = "non-finite")
})
