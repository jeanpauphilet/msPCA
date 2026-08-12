# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for X-input preprocessing shared by mspca() and tpm().

set.seed(123)
n_obs <- 30L
p_pre <- 5L
X_pre <- matrix(rnorm(n_obs * p_pre), n_obs, p_pre,
                dimnames = list(NULL, paste0("x", seq_len(p_pre))))

# ---------- mspca() ----------------------------------------------------------

test_that("mspca (X) records center=FALSE and scale=TRUE", {
  res <- mspca(X_pre, r = 1L, ks = 2L, type = "X",
               center = FALSE, scale = TRUE, verbose = FALSE)

  expect_s3_class(res, "mspca")
  expect_false(res$center)
  expect_true(res$scale)
  expect_equal(res$divisor, "n-1")
  expect_equal(dim(res$x_best), c(p_pre, 1L))
  expect_true(all(is.finite(res$x_best)))
  expect_equal(sum(res$x_best^2), 1, tolerance = 1e-6)
})

test_that("mspca (X) divisor changes variance scale but not loadings", {
  set.seed(321)
  res_n1 <- mspca(X_pre, r = 2L, ks = c(2L, 2L), type = "X",
                  divisor = "n-1", verbose = FALSE)
  set.seed(321)
  res_n <- mspca(X_pre, r = 2L, ks = c(2L, 2L), type = "X",
                 divisor = "n", verbose = FALSE)

  loading_alignment <- abs(diag(crossprod(res_n1$x_best, res_n$x_best)))
  expect_equal(loading_alignment, c(1, 1), tolerance = 1e-6)
  expect_equal(res_n$total_variance / res_n1$total_variance,
               (n_obs - 1) / n_obs, tolerance = 1e-8)
  expect_equal(sum(res_n$variance_explained) / res_n$total_variance,
               sum(res_n1$variance_explained) / res_n1$total_variance,
               tolerance = 1e-8)
})

test_that("mspca (X) centered and scaled data matches the correlation scale", {
  set.seed(456)
  res_x <- mspca(X_pre, r = 1L, ks = 2L, type = "X",
                 center = TRUE, scale = TRUE, divisor = "n-1",
                 verbose = FALSE)
  set.seed(456)
  res_sigma <- mspca(stats::cor(X_pre), r = 1L, ks = 2L,
                     verbose = FALSE)

  expect_equal(dim(res_x$x_best), dim(res_sigma$x_best))
  expect_equal(colSums(res_x$x_best^2), c(1), tolerance = 1e-6)
  expect_equal(res_x$variance_explained, res_sigma$variance_explained,
               tolerance = 1e-6)
  expect_equal(res_x$total_variance, res_sigma$total_variance,
               tolerance = 1e-6)
})

# ---------- tpm() ------------------------------------------------------------

test_that("tpm (X) records center=FALSE and scale=TRUE", {
  res <- tpm(X_pre, k = 2L, type = "X", center = FALSE, scale = TRUE,
             verbose = FALSE)

  expect_s3_class(res, "tpm")
  expect_false(res$center)
  expect_true(res$scale)
  expect_equal(res$divisor, "n-1")
  expect_equal(dim(res$x_best), c(p_pre, 1L))
  expect_true(all(is.finite(res$x_best)))
  expect_equal(sum(res$x_best^2), 1, tolerance = 1e-6)
})

test_that("tpm (X) divisor changes objective scale but not loading", {
  set.seed(654)
  res_n1 <- tpm(X_pre, k = 2L, type = "X", divisor = "n-1",
                verbose = FALSE)
  set.seed(654)
  res_n <- tpm(X_pre, k = 2L, type = "X", divisor = "n",
               verbose = FALSE)

  loading_alignment <- abs(drop(crossprod(res_n1$x_best, res_n$x_best)))
  expect_equal(loading_alignment, 1, tolerance = 1e-6)
  expect_equal(res_n$objective_value / res_n1$objective_value,
               (n_obs - 1) / n_obs, tolerance = 1e-8)
})

test_that("tpm supports one-sparse and dense components for both input types", {
  tpm_sigma_one <- tpm(crossprod(matrix(rnorm(p_pre * 10), 10, p_pre)),
                       k = 1L, verbose = FALSE)
  tpm_sigma_dense <- tpm(crossprod(matrix(rnorm(p_pre * 10), 10, p_pre)),
                         k = p_pre, verbose = FALSE)
  tpm_x_one <- tpm(X_pre, k = 1L, type = "X", verbose = FALSE)
  tpm_x_dense <- tpm(X_pre, k = p_pre, type = "X", verbose = FALSE)

  expect_equal(sum(tpm_sigma_one$x_best != 0), 1L)
  expect_equal(sum(tpm_x_one$x_best != 0), 1L)
  expect_equal(sum(tpm_sigma_one$x_best^2), 1, tolerance = 1e-6)
  expect_equal(sum(tpm_x_one$x_best^2), 1, tolerance = 1e-6)
  expect_lte(sum(tpm_sigma_dense$x_best != 0), p_pre)
  expect_lte(sum(tpm_x_dense$x_best != 0), p_pre)
  expect_true(all(is.finite(c(tpm_sigma_dense$objective_value,
                              tpm_x_dense$objective_value))))
})
