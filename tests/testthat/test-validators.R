# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for input validation: the matrix validators
# (.validate_sigma_matrix() / .validate_x_matrix()) and the scalar/vector
# argument validators (.validate_solver_args() and friends), all exercised
# indirectly through mspca() and tpm().

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

# ---------- Cardinality (`ks` / `k`) -----------------------------------------
#
# These values reach C++ as an `int` and end up in `truncateVector()`, which forms
# the iterator `idx.begin() + k`. Zero gives an all-zero "component"; negative and
# NA (which becomes NA_INTEGER, i.e. INT_MIN) are undefined behaviour. R must
# reject them before the call, so each case below asserts an error rather than
# whatever the solver happens to do with it.

test_that("k = 0 is rejected", {
  expect_error(tpm(Sigma_ok, k = 0, verbose = FALSE), regexp = "must be >= 1")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 0, verbose = FALSE), regexp = "must be >= 1")
})

test_that("negative k is rejected", {
  expect_error(tpm(Sigma_ok, k = -2, verbose = FALSE), regexp = "must be >= 1")
  expect_error(mspca(Sigma_ok, r = 2L, ks = c(2, -1), verbose = FALSE),
               regexp = "must be >= 1")
})

test_that("fractional k is rejected rather than silently truncated", {
  expect_error(tpm(Sigma_ok, k = 2.5, verbose = FALSE), regexp = "whole number")
  expect_error(mspca(Sigma_ok, r = 2L, ks = c(2, 3.7), verbose = FALSE),
               regexp = "whole number")
})

test_that("NA, NaN and Inf cardinalities are rejected", {
  expect_error(tpm(Sigma_ok, k = NA, verbose = FALSE), regexp = "numeric vector|NA")
  expect_error(tpm(Sigma_ok, k = NA_real_, verbose = FALSE), regexp = "NA, NaN or Inf")
  expect_error(tpm(Sigma_ok, k = NaN, verbose = FALSE), regexp = "NA, NaN or Inf")
  expect_error(tpm(Sigma_ok, k = Inf, verbose = FALSE), regexp = "NA, NaN or Inf")
  expect_error(mspca(Sigma_ok, r = 2L, ks = c(2, NA_real_), verbose = FALSE),
               regexp = "NA, NaN or Inf")
})

test_that("non-numeric and empty cardinalities are rejected", {
  expect_error(tpm(Sigma_ok, k = "2", verbose = FALSE), regexp = "numeric vector")
  expect_error(mspca(Sigma_ok, r = 1L, ks = integer(0), verbose = FALSE),
               regexp = "at least one element")
})

test_that("tpm() requires a single k", {
  expect_error(tpm(Sigma_ok, k = c(2L, 3L), verbose = FALSE),
               regexp = "single number")
})

test_that("the error names the offending element of ks", {
  expect_error(mspca(Sigma_ok, r = 3L, ks = c(2, 2, 0), verbose = FALSE),
               regexp = "element\\(s\\) 3")
})

# ---------- Number of components (`r`) ---------------------------------------

test_that("invalid r is rejected", {
  expect_error(mspca(Sigma_ok, r = 0L, ks = 2L, verbose = FALSE), regexp = "must be >= 1")
  expect_error(mspca(Sigma_ok, r = -1L, ks = 2L, verbose = FALSE), regexp = "must be >= 1")
  expect_error(mspca(Sigma_ok, r = 2.5, ks = c(2L, 2L), verbose = FALSE),
               regexp = "whole number")
  expect_error(mspca(Sigma_ok, r = NA_integer_, ks = 2L, verbose = FALSE),
               regexp = "finite number")
  expect_error(mspca(Sigma_ok, r = c(1L, 2L), ks = 2L, verbose = FALSE),
               regexp = "single number")
})

test_that("whole-number doubles are accepted for r and ks", {
  # The documentation and examples pass r = 2 and ks = c(4, 4), which are doubles;
  # requiring an integer type would reject the package's own usage.
  expect_error(mspca(Sigma_ok, r = 2, ks = c(2, 2), verbose = FALSE), NA)
  expect_error(tpm(Sigma_ok, k = 2, verbose = FALSE), NA)
})

test_that("more sparsity levels than components is an error", {
  expect_error(mspca(Sigma_ok, r = 2L, ks = c(2L, 2L, 2L), verbose = FALSE),
               regexp = "only 2 principal component")
})

test_that("fewer sparsity levels than components still warns and shrinks r", {
  # Documented behaviour, deliberately left to the solver rather than turned into
  # an error by the validator.
  expect_warning(
    res <- mspca(Sigma_ok, r = 3L, ks = c(2L, 2L), verbose = FALSE),
    regexp = "only provided"
  )
  expect_equal(ncol(res$x_best), 2L)
})

# ---------- Budgets and tolerances -------------------------------------------

test_that("nonpositive maxIter is rejected", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, maxIter = 0L, verbose = FALSE),
               regexp = "`maxIter` must be >= 1")
  expect_error(tpm(Sigma_ok, k = 2L, maxIter = -5L, verbose = FALSE),
               regexp = "`maxIter` must be >= 1")
})

test_that("negative restart budgets are rejected but zero is allowed", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, maxRestartTPM = -1L, verbose = FALSE),
               regexp = "`maxRestartTPM` must be >= 0")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, minRestartTPM = -1L, verbose = FALSE),
               regexp = "`minRestartTPM` must be >= 0")
  # Zero restarts is a legitimate configuration: seeded run only.
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, maxRestartTPM = 0L,
                     minRestartTPM = 0L, verbose = FALSE), NA)
})

test_that("negative and missing tolerances are rejected", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityTolerance = -1e-4,
                     verbose = FALSE), regexp = "`feasibilityTolerance` must be >= 0")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, stallingTolerance = -1,
                     verbose = FALSE), regexp = "`stallingTolerance` must be >= 0")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityTolerance = NA_real_,
                     verbose = FALSE), regexp = "must not be NA")
})

test_that("Inf tolerances remain supported", {
  # feasibilityTolerance = Inf accepts any solution as feasible; the C++ layer
  # takes the tolerances as doubles, so Inf passes through meaningfully.
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityTolerance = Inf,
                     verbose = FALSE), NA)
})

test_that("invalid time limits are rejected, including Inf", {
  # The C++ parameter is an int, so Inf would arrive as NA_INTEGER (INT_MIN) and
  # silently disable the restarts instead of removing the limit.
  expect_error(tpm(Sigma_ok, k = 2L, timeLimit = Inf, verbose = FALSE),
               regexp = "finite number of seconds")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, timeLimitTPM = Inf, verbose = FALSE),
               regexp = "finite number of seconds")
  expect_error(tpm(Sigma_ok, k = 2L, timeLimit = -1, verbose = FALSE),
               regexp = "`timeLimit` must be >= 0")
  expect_error(tpm(Sigma_ok, k = 2L, timeLimit = NA_real_, verbose = FALSE),
               regexp = "finite number of seconds")
})

test_that("invalid feasibilityConstraintType is rejected", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityConstraintType = 2L,
                     verbose = FALSE), regexp = "must be 0")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityConstraintType = NA,
                     verbose = FALSE), regexp = "must be 0")
})

test_that("fractional feasibilityConstraintType is rejected, not truncated", {
  # The old order of operations -- as.integer() first, membership test second --
  # accepted anything in (-1, 2): 0.5, 0.9 and -0.5 all truncate to 0 and 1.7
  # truncates to 1, so the solver silently enforced a constraint that was never
  # requested. Each value below must error rather than be rounded into range.
  for (bad_val in list(0.5, 0.9, 1.5, 1.7, -0.5)) {
    expect_error(mspca(Sigma_ok, r = 1L, ks = 2L,
                       feasibilityConstraintType = bad_val, verbose = FALSE),
                 regexp = "must be 0")
  }
})

test_that("non-numeric, missing and non-finite constraint types are rejected", {
  for (bad_val in list("0", "orthogonality", NA, NA_real_, NaN, Inf, -Inf,
                       NULL, c(0L, 1L), TRUE)) {
    expect_error(mspca(Sigma_ok, r = 1L, ks = 2L,
                       feasibilityConstraintType = bad_val, verbose = FALSE),
                 regexp = "must be 0")
  }
})

test_that("valid constraint types are accepted in either numeric type", {
  for (good_val in list(0L, 1L, 0, 1)) {
    expect_error(mspca(Sigma_ok, r = 2L, ks = c(2L, 2L),
                       feasibilityConstraintType = good_val, verbose = FALSE), NA)
  }
})

test_that("summary.mspca applies the same rule, and still accepts NULL", {
  res <- mspca(Sigma_ok, r = 2L, ks = c(2L, 2L), verbose = FALSE)
  for (bad_val in list(0.5, 1.5, -0.5, NA, Inf, "1", c(0L, 1L))) {
    expect_error(
      capture.output(summary(res, feasibilityConstraintType = bad_val)),
      regexp = "must be 0"
    )
  }
  # NULL means "report under the type used to fit", and must keep working.
  expect_error(capture.output(summary(res, feasibilityConstraintType = NULL)), NA)
  expect_error(capture.output(summary(res)), NA)
})

test_that("the constraint-type error mentions NULL only where NULL is allowed", {
  res <- mspca(Sigma_ok, r = 2L, ks = c(2L, 2L), verbose = FALSE)
  e_fit <- tryCatch(mspca(Sigma_ok, r = 1L, ks = 2L, feasibilityConstraintType = 0.5,
                          verbose = FALSE),
                    error = conditionMessage)
  e_sum <- tryCatch(capture.output(summary(res, feasibilityConstraintType = 0.5)),
                    error = conditionMessage)
  expect_false(grepl("NULL", e_fit))
  expect_true(grepl("or NULL", e_sum))
})

# ---------- On/off switches ---------------------------------------------------
#
# `isTRUE()` is deliberately strict: isTRUE(1) is FALSE. Using it to read a
# user-supplied switch meant that `checkPSD = 1` silently skipped the PSD check,
# i.e. failed in the unsafe direction with no message. .validate_flag() coerces
# an unambiguous 0/1 and rejects everything else.

test_that("checkPSD = 1 runs the PSD check rather than silently skipping it", {
  bad      <- Sigma_ok
  diag(bad) <- diag(bad) - max(diag(bad)) - 1
  bad      <- (bad + t(bad)) / 2            # force negative eigenvalues
  expect_error(mspca(bad, r = 1L, ks = 2L, checkPSD = 1, verbose = FALSE),
               regexp = "positive semidefinite")
  expect_error(tpm(bad, k = 2L, checkPSD = 1, verbose = FALSE),
               regexp = "positive semidefinite")
  # and 0 still means "skip", as FALSE does
  expect_error(mspca(bad, r = 1L, ks = 2L, checkPSD = 0, verbose = FALSE), NA)
})

test_that("scale = 1 standardizes rather than being read as FALSE", {
  # With scale read as FALSE the zero-variance guard never fires; with it read as
  # TRUE the column is caught. The error is the observable difference.
  bad      <- matrix(rnorm(30 * p), 30, p)
  bad[, 1] <- 0
  expect_error(mspca(bad, r = 1L, ks = 2L, type = "X", scale = 1, verbose = FALSE),
               regexp = "zero-variance")
})

test_that("0/1 are accepted wherever a Boolean is expected", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, verbose = 0), NA)
  expect_error(tpm(Sigma_ok, k = 2L, verbose = 0), NA)
  expect_error(mspca(matrix(rnorm(30 * p), 30, p), r = 1L, ks = 2L, type = "X",
                     center = 1, scale = 0, verbose = 0), NA)
})

test_that("ambiguous or missing flag values are rejected", {
  for (bad_val in list(2, -1, NA, "TRUE", c(TRUE, TRUE), NULL)) {
    expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, checkPSD = bad_val, verbose = FALSE),
                 regexp = "`checkPSD` must be a single TRUE or FALSE")
  }
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, verbose = NA),
               regexp = "`verbose` must be a single TRUE or FALSE")
})

test_that("the flag error names the argument and reports what was passed", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, verbose = "yes"),
               regexp = "`verbose`.*Got: character yes")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, verbose = c(TRUE, FALSE)),
               regexp = "Got: logical of length 2")
})

# ---------- Matrix-check tolerances -------------------------------------------

test_that("invalid symmetry/PSD tolerances are rejected", {
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, symTolerance = -1, verbose = FALSE),
               regexp = "`symTolerance` must be >= 0")
  expect_error(mspca(Sigma_ok, r = 1L, ks = 2L, psdTolerance = NA_real_, verbose = FALSE),
               regexp = "`psdTolerance` must not be NA")
  expect_error(tpm(Sigma_ok, k = 2L, symTolerance = NA_real_, verbose = FALSE),
               regexp = "`symTolerance` must not be NA")
})

# ---------- Arguments are validated before the matrix is touched -------------

test_that("argument validation does not depend on a valid matrix", {
  # Cheap checks should fire regardless; in particular the O(p^3) PSD eigen
  # decomposition should not run before k has been rejected.
  bad <- matrix(as.character(1:(p * p)), p, p)
  expect_error(tpm(bad, k = -1, verbose = FALSE), regexp = "must be >= 1")
})

# ---------- Clamping above the dimension -------------------------------------

test_that("mspca() and tpm() clamp k above the dimension", {
  expect_warning(
    mspca_res <- mspca(Sigma_ok, r = 1L, ks = p + 1L, verbose = FALSE),
    regexp = "exceeds the dimension"
  )
  expect_warning(
    tpm_res <- tpm(Sigma_ok, k = p + 1L, verbose = FALSE),
    regexp = "exceeds the dimension"
  )

  expect_lte(sum(mspca_res$x_best != 0), p)
  expect_lte(sum(tpm_res$x_best != 0), p)
})
