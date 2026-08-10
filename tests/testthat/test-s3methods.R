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
  expect_no_error(capture.output(summary(res, Cm, feasibilityConstraintType = 1L)))
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
