# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet

# Tests for S3 methods: print.mspca, summary.mspca, print.summary.mspca

set.seed(11)
p    <- 8L
Cm   <- cor(datasets::mtcars)
res  <- mspca(Cm, r = 2L, ks = c(4L, 4L), verbose = FALSE)
resX <- mspca(as.matrix(datasets::mtcars), r = 2L, ks = c(4L, 4L),
              type = "X", verbose = FALSE)

# ---------- print.mspca ------------------------------------------------------

test_that("print.mspca returns x invisibly", {
  out <- withVisible(print(res, Cm))
  expect_false(out$visible)
  expect_identical(out$value, res)
})

test_that("print.mspca produces output containing 'msPCA solution'", {
  txt <- capture.output(print(res, Cm))
  expect_true(any(grepl("msPCA solution", txt)))
})

test_that("print.mspca works with type='X' result without C argument", {
  expect_no_error(capture.output(print(resX)))
})

test_that("print.mspca works with digits argument", {
  expect_no_error(capture.output(print(res, Cm, digits = 3L)))
})

test_that("print.mspca errors when C is NULL for a Sigma-fit result with no stored variance", {
  # Construct a minimal mspca object missing variance_explained / total_variance
  bad <- structure(list(x_best = res$x_best, inputType = "Sigma"),
                   class = "mspca")
  expect_error(print(bad), regexp = "variance_explained")
})

# ---------- summary.mspca ----------------------------------------------------

test_that("summary.mspca returns invisibly an object of class 'summary.mspca'", {
  out <- withVisible(summary(res, Cm))
  expect_false(out$visible)
  expect_s3_class(out$value, "summary.mspca")
})

test_that("summary.mspca contains all required list fields", {
  s <- summary(res, Cm)
  expect_true(all(c("table", "feasibility_mat", "feasibility",
                    "runtime", "r", "inputType") %in% names(s)))
})

test_that("summary.mspca table has correct structure", {
  s   <- summary(res, Cm)
  tbl <- s$table
  expect_s3_class(tbl, "data.frame")
  expect_equal(nrow(tbl), 2L)
  expect_true(all(c("PC", "nonzero", "variance", "fve", "cumulative_fve") %in%
                    names(tbl)))
})

test_that("summary.mspca FVE values are in [0, 1]", {
  s <- summary(res, Cm)
  expect_true(all(s$table$fve >= 0 & s$table$fve <= 1 + 1e-9))
  expect_true(all(s$table$cumulative_fve <= 1 + 1e-9))
})

test_that("summary.mspca feasibility_mat is an r x r matrix", {
  s <- summary(res, Cm)
  expect_equal(dim(s$feasibility_mat), c(2L, 2L))
})

test_that("summary.mspca scalar feasibility is non-negative", {
  s <- summary(res, Cm)
  expect_gte(s$feasibility, 0)
})

test_that("summary.mspca works with type='X' result without C argument", {
  expect_no_error(capture.output(summary(resX)))
})

test_that("summary.mspca feasibilityConstraintType=1 runs without error", {
  expect_no_error(capture.output(
    suppressWarnings(summary(res, Cm, feasibilityConstraintType = 1L))))
})

# ---------- constraint type is stored and reused (regression) ----------------
# Regression tests for the referee's report: summary() defaulted to
# orthogonality regardless of how the model was fitted.

resU <- mspca(Cm, r = 2L, ks = c(4L, 4L),
              feasibilityConstraintType = 1L, verbose = FALSE)

test_that("Sigma-mode fits carry the variable names on x_best", {
  expect_identical(rownames(res$x_best), colnames(Cm))
  # print() must therefore label rows correctly with no C argument.
  txt <- capture.output(print(res))
  expect_true(any(grepl(colnames(Cm)[which(res$x_best[, 1L] != 0)[1L]], txt,
                        fixed = TRUE)))
})

test_that("print() output is the same with and without C", {
  expect_identical(capture.output(print(res)), capture.output(print(res, Cm)))
})

test_that("mspca stores the constraint type it was fitted with", {
  expect_identical(res$feasibilityConstraintType,  0L)
  expect_identical(resU$feasibilityConstraintType, 1L)
})

test_that("mspca stores both non-redundancy matrices, correctly oriented", {
  for (o in list(res, resU)) {
    nr <- o$nonredundancy
    expect_true(all(c("orthogonality", "uncorrelatedness") %in% names(nr)))
    for (M in nr) {
      expect_equal(dim(M), c(2L, 2L))
      expect_true(all(is.na(diag(M))))
      expect_true(all(is.na(M[lower.tri(M)])))
      expect_true(all(M[upper.tri(M)] >= 0))
      expect_identical(colnames(M), c("PC1", "PC2"))
    }
  }
})

# Reference implementations of the two pairwise statistics, written out in
# full so the tests do not simply mirror the internal helpers.
ref_ortho <- function(v, i, j) abs(sum(v[, i] * v[, j]))
ref_uncor <- function(v, C, i, j) {
  abs(sum(v[, i] * (C %*% v[, j]))) / (sum(diag(C)) / nrow(C))
}

test_that("stored matrices match a direct recomputation from C", {
  nr <- resU$nonredundancy
  v  <- resU$x_best
  expect_equal(nr$orthogonality[1L, 2L],    ref_ortho(v, 1L, 2L))
  expect_equal(nr$uncorrelatedness[1L, 2L], ref_uncor(v, Cm, 1L, 2L))
})

test_that("the X path reproduces the Sigma-path non-redundancy values", {
  Xm  <- as.matrix(datasets::mtcars)
  nrX <- resX$nonredundancy
  v   <- resX$x_best
  Cx  <- stats::cor(Xm)   # resX was fitted with center = TRUE, scale = TRUE
  expect_equal(nrX$uncorrelatedness[1L, 2L], ref_uncor(v, Cx, 1L, 2L),
               tolerance = 1e-8)
  expect_equal(nrX$orthogonality[1L, 2L], ref_ortho(v, 1L, 2L),
               tolerance = 1e-8)
})

# ---------- normalization of the uncorrelatedness measure --------------------

test_that("uncorrelatedness violations are normalized by tr(Sigma)/p", {
  Sc <- stats::cov(as.matrix(datasets::mtcars))   # tr(Sc)/p is far from 1
  expect_false(isTRUE(all.equal(sum(diag(Sc)) / nrow(Sc), 1)))
  resS <- mspca(Sc, r = 2L, ks = c(4L, 4L), verbose = FALSE)
  v    <- resS$x_best
  expect_equal(resS$nonredundancy$uncorrelatedness[1L, 2L],
               ref_uncor(v, Sc, 1L, 2L))
  # The unnormalized quantity is a different number, so the normalization is
  # demonstrably being applied rather than silently equal to 1.
  expect_false(isTRUE(all.equal(resS$nonredundancy$uncorrelatedness[1L, 2L],
                               abs(crossprod(v, Sc %*% v)[1L, 2L]))))
  Mu <- resS$nonredundancy$uncorrelatedness
  expect_equal(feasibility_violation_off(Sc, v, 1L), sum(Mu[upper.tri(Mu)]))
})

test_that("the uncorrelatedness measure is invariant to a rescaling of Sigma", {
  Sc <- stats::cov(as.matrix(datasets::mtcars))
  v  <- mspca(Sc, r = 2L, ks = c(4L, 4L), verbose = FALSE)$x_best
  expect_equal(feasibility_violation_off(Sc,       v, 1L),
               feasibility_violation_off(1000 * Sc, v, 1L))
  # On a correlation matrix tr(C)/p = 1, so the measure is unchanged from
  # earlier versions of the package.
  expect_equal(sum(diag(Cm)) / nrow(Cm), 1)
})

test_that("a degenerate zero-trace input does not produce NaN", {
  Z <- matrix(0, 4L, 4L)
  U <- cbind(c(1, 0, 0, 0), c(0, 1, 0, 0))
  expect_equal(feasibility_violation_off(Z, U, 1L), 0)
  expect_true(is.finite(feasibility_violation_off(Z, U, 1L)))
})

test_that("summary() defaults to the fitted constraint type, not to orthogonality", {
  sU <- summary(resU)
  expect_identical(sU$feasibilityConstraintType, 1L)
  expect_identical(sU$fittedConstraintType,      1L)
  expect_equal(sU$feasibility,
               feasibility_violation_off(Cm, resU$x_best, 1L), tolerance = 1e-10)

  s0 <- summary(res)
  expect_identical(s0$feasibilityConstraintType, 0L)
  expect_equal(s0$feasibility,
               feasibility_violation_off(Cm, res$x_best, 0L), tolerance = 1e-10)
})

test_that("summary() reports the uncorrelatedness figure, not the orthogonality one", {
  # Deterministic fixture with hand-computable values, so the assertion cannot
  # pass by accident. With C = diag(4, 1, 1), u1 = (1, 0, 0), u2 = (0.6, 0.8, 0):
  #   orthogonality    |u1'u2|                  = 0.6
  #   uncorrelatedness |u1'C u2| / (tr(C)/p)    = 2.4 / 2 = 1.2
  # Under the old default, summary() of an uncorrelatedness fit reported 0.6.
  C  <- diag(c(4, 1, 1))
  v  <- cbind(c(1, 0, 0), c(0.6, 0.8, 0))
  expect_equal(feasibility_violation_off(C, v, 0L), 0.6)
  expect_equal(feasibility_violation_off(C, v, 1L), 1.2)

  obj <- structure(
    list(x_best = v,
         variance_explained = as.numeric(diag(crossprod(v, C %*% v))),
         total_variance = sum(diag(C)),
         inputType = "Sigma",
         runtime = 0,
         feasibilityConstraintType = 1L,
         nonredundancy = .nonredundancy(v, C %*% v, .avg_variance(C))),
    class = "mspca")

  s <- summary(obj)
  expect_equal(s$feasibility, feasibility_violation_off(C, v, 1L))
  expect_false(isTRUE(all.equal(s$feasibility,
                               feasibility_violation_off(C, v, 0L))))

  # Explicitly requesting the other type returns the other number.
  s0 <- suppressWarnings(summary(obj, feasibilityConstraintType = 0L))
  expect_equal(s0$feasibility, feasibility_violation_off(C, v, 0L))
})

test_that("summary() prints the constraint definition it is reporting", {
  txt <- capture.output(summary(resU))
  expect_true(any(grepl("uncorrelatedness", txt)))
  expect_true(any(grepl("as fitted", txt)))
})

test_that("summary() warns when asked for a type other than the fitted one", {
  expect_warning(capture.output(summary(resU, feasibilityConstraintType = 0L)),
                 regexp = "fitted under uncorrelatedness")
  expect_no_warning(capture.output(summary(resU)))
  expect_no_warning(capture.output(summary(resU, feasibilityConstraintType = 1L)))
})

test_that("summary() with an explicit type returns that type's figures", {
  s <- suppressWarnings(summary(resU, feasibilityConstraintType = 0L))
  expect_identical(s$feasibilityConstraintType, 0L)
  expect_identical(s$fittedConstraintType,      1L)
  expect_equal(s$feasibility,
               feasibility_violation_off(Cm, resU$x_best, 0L), tolerance = 1e-10)
})

test_that("summary() needs no C argument, for either constraint type", {
  # Previously type 1 without C silently substituted the identity for Sigma
  # in the scalar, and errored on NULL in the pairwise matrix.
  expect_no_error(capture.output(summary(resU)))
  s_noC <- suppressWarnings(summary(resU))
  s_C   <- suppressWarnings(summary(resU, Cm))
  expect_equal(s_noC$feasibility, s_C$feasibility)
  expect_equal(s_noC$feasibility_mat, s_C$feasibility_mat)
})

test_that("summary() of an X-mode fit under uncorrelatedness needs no C", {
  resXU <- mspca(as.matrix(datasets::mtcars), r = 2L, ks = c(4L, 4L),
                 type = "X", feasibilityConstraintType = 1L, verbose = FALSE)
  s <- summary(resXU)
  expect_identical(s$feasibilityConstraintType, 1L)
  expect_equal(s$feasibility,
               feasibility_violation_off(stats::cor(as.matrix(datasets::mtcars)),
                                         resXU$x_best, 1L),
               tolerance = 1e-8)
})

# ---------- per-PC violations ------------------------------------------------

test_that("feasibility_perPC is the row-wise max of the symmetrised matrix", {
  s <- summary(res)
  M <- s$feasibility_mat
  expect_length(s$feasibility_perPC, 2L)
  expect_equal(unname(s$feasibility_perPC), rep(M[1L, 2L], 2L))
  expect_equal(s$table$max_violation, unname(s$feasibility_perPC))
})

test_that("feasibility_perPC picks the largest pair, not a sum", {
  res3 <- mspca(Cm, r = 3L, ks = c(4L, 4L, 4L), verbose = FALSE)
  s3   <- summary(res3)
  M    <- s3$feasibility_mat
  expect_equal(unname(s3$feasibility_perPC[1L]), max(M[1L, 2L], M[1L, 3L]))
  expect_equal(unname(s3$feasibility_perPC[3L]), max(M[1L, 3L], M[2L, 3L]))
  expect_equal(s3$feasibility, sum(M[upper.tri(M)]))
})

test_that("r = 1 gives zero total and zero per-PC violation", {
  res1 <- mspca(Cm, r = 1L, ks = 4L, verbose = FALSE)
  s1   <- summary(res1)
  expect_equal(s1$feasibility, 0)
  expect_equal(unname(s1$feasibility_perPC), 0)
})

# ---------- backward compatibility -------------------------------------------

test_that("legacy objects without stored diagnostics still summarize with C", {
  legacy <- res
  legacy$nonredundancy <- NULL
  legacy$feasibilityConstraintType <- NULL
  s <- summary(legacy, Cm)
  expect_true(is.na(s$fittedConstraintType))
  expect_identical(s$feasibilityConstraintType, 0L)   # historical default
  expect_equal(s$feasibility, feasibility_violation_off(Cm, res$x_best, 0L))
  expect_equal(s$feasibility_mat, res$nonredundancy$orthogonality)
})

test_that("legacy objects error, rather than mislead, for type 1 without C", {
  legacy <- res
  legacy$nonredundancy <- NULL
  legacy$feasibilityConstraintType <- NULL
  expect_error(summary(legacy, feasibilityConstraintType = 1L),
               regexp = "predates msPCA 0.5.1")
})

test_that("mspca rejects an invalid feasibilityConstraintType", {
  expect_error(mspca(Cm, r = 2L, ks = c(4L, 4L), feasibilityConstraintType = 2L,
                     verbose = FALSE),
               regexp = "must be 0")
  expect_error(summary(res, feasibilityConstraintType = 3L), regexp = "must be 0")
})

test_that("the solver's feasibility_violation is not the summary statistic", {
  # They differ by the diagonal normalization terms; documented in ?mspca.
  s <- summary(res)
  expect_gte(res$feasibility_violation, 0)
  expect_equal(s$feasibility, feasibility_violation_off(Cm, res$x_best, 0L))
})

test_that("summary.mspca r=1 (single PC) runs without error", {
  res1 <- mspca(Cm, r = 1L, ks = 4L, verbose = FALSE)
  expect_no_error(capture.output(summary(res1, Cm)))
})

# ---------- helper functions -------------------------------------------------

test_that("fraction_variance_explained returns a scalar in (0, 1]", {
  fve <- fraction_variance_explained(Cm, res$x_best)
  expect_length(fve, 1L)
  expect_gt(fve, 0)
  expect_lte(fve, 1 + 1e-9)
})

test_that("fraction_variance_explained_perPC returns a vector of length r", {
  fve_pc <- fraction_variance_explained_perPC(Cm, res$x_best)
  expect_length(fve_pc, 2L)
  expect_true(all(fve_pc >= 0))
})

test_that("feasibility_violation_off is 0 for a single PC", {
  res1 <- mspca(Cm, r = 1L, ks = 4L, verbose = FALSE)
  fv   <- feasibility_violation_off(Cm, res1$x_best, 0L)
  expect_equal(fv, 0, tolerance = 1e-9)
})
