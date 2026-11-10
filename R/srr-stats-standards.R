#' srr_stats
#'
#' Statistical software review standards (general and unsupervised learning
#' categories). Each standard is addressed here or at the location named in the
#' tag text; standards that are addressed by tests are tagged in `tests/testthat/`.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstats {G1.0} Primary reference: Cory-Wright and Pauphilet, 'Sparse PCA with Multiple Components', Operations Research (doi:10.1287/opre.2023.0598), cited in DESCRIPTION, the README and `citation("msPCA")`; the truncated power method of Yuan and Zhang (2013) is cited in `?tpm`.
#' @srrstats {G1.1} The multi-component deflation/penalty algorithm is a novel algorithm (the first implementation, introduced in the reference paper); the single-component `tpm()` is a new implementation of the existing truncated power method of Yuan and Zhang (2013). Both statements are made in the 'Algorithm and implementation notes' vignette and the `?mspca` and `?tpm` documentation.
#' @srrstats {G1.2} A Life Cycle Statement (actively maintained; stable API; anticipated developments) is given in the README, section 'Life cycle'.
#' @srrstats {G1.3} Terminology (sparsity level, orthogonality and uncorrelatedness constraints, feasibility violation, fraction of variance explained) is defined in the function documentation and the 'Algorithm and implementation notes' vignette.
#' @srrstats {G1.4} All exported functions and the exported data set are documented with roxygen2.
#' @srrstats {G1.4a} All internal (non-exported) functions in `R/main.R` carry roxygen2 comment blocks (tagged `@noRd`) describing their purpose, inputs and return value.
#' @srrstats {G1.5} The `replication/` directory holds all scripts needed to reproduce the benchmarking and case-study results reported in the paper and the vignettes.
#' @srrstats {G1.6} The website article 'Benchmarking against other sparse PCA packages' (`vignettes/articles/`) and the scripts in `replication/benchmarking/` compare msPCA with seven alternative R implementations.
#' @srrstats {G2.0a} Expectations on input lengths are stated in the `@param` documentation (for example 'a single whole number', '`length(ks)` must not exceed `r`').
#' @srrstats {G2.1a} Expected data types of all inputs are documented in the `@param` entries of `mspca()` and `tpm()`.
#' @srrstats {G2.4} Conversions between types are explicit and happen once at the entry points.
#' @srrstats {G2.14} The documented policy for missing data is to error; users who need another policy can compute a covariance matrix with `cov(use = ...)` and pass it with `type = "Sigma"` (documented in `?mspca`, tested in `test-standards.R`).
#' @srrstats {G2.14a} Missing data produces an informative error.
#' @srrstats {G5.2a} Condition messages are unique; messages produced by the shared validators are parametrized by the name of the offending argument, and the two matrix validators use distinct messages.
#' @srrstats {UL1.0} The expected format of the input (numeric matrix or all-numeric data.frame; covariance/correlation matrix or raw data) and the rejected formats are documented in the `@param M` entry of `mspca()` and `tpm()`.
#' @srrstats {UL1.3a} Attributes of the input other than the column names (which label the rows of `x_best`) are not transferred to the result; this is stated in the `@return` documentation of `mspca()`.
#' @srrstats {UL1.4} Assumptions on the input (symmetric positive semidefinite covariance matrix, finite entries, at least two observations) are documented in `@param M`, `@param checkPSD` and in the vignette.
#' @srrstats {UL1.4a} The sensitivity of the solution to the scale of the variables, and the role of `scale`, are documented in `@param scale` and illustrated in the examples of `?mspca` (covariance versus correlation).
#' @srrstats {UL1.4b} The `?mspca` examples explain why `scale = TRUE` is used and contrast it with `scale = FALSE`.
#' @srrstats {UL4.0} `mspca()` and `tpm()` return model objects of S3 classes `"mspca"` and `"tpm"`.
#' @srrstats {UL4.3} `print.mspca()` provides an on-screen summary of the model (number of components, variance explained, sparsity, active loadings).
#' @srrstats {UL4.3a} `print.mspca()` prints only the rows (variables) with a non-zero loading, at most the number of variables selected by the user through `ks`.
#' @srrstats {UL4.4} `summary.mspca()` summarizes the per-component statistics and the pairwise non-redundancy violations.
#' @noRd
NULL

#' NA_standards
#'
#' Standards that do not apply to msPCA, with the reason.
#'
#' @srrstatsNA {G2.4c} The package has no character input to convert.
#' @srrstatsNA {G2.4d} The package has no factor input to convert.
#' @srrstatsNA {G2.4e} The package has no factor input to convert.
#' @srrstatsNA {G2.5} The package accepts no factor inputs.
#' @srrstatsNA {G2.6} The package accepts no one-dimensional input; the data is a matrix.
#' @srrstatsNA {G2.10} The package does not extract single columns from tabular inputs; data.frames are converted to a matrix as a whole.
#' @srrstatsNA {G2.14b} Ignoring missing values has no well-defined meaning for the algorithm, which operates on a covariance matrix; missing values are an error (see G2.14).
#' @srrstatsNA {G2.14c} Imputation is outside the scope of the package; see G2.14 for how to supply a covariance matrix estimated with missing data.
#' @srrstatsNA {G4.0} The package does not write outputs to files.
#' @srrstatsNA {G5.4b} msPCA is a new method; no prior implementation exists for its multi-component formulation. Correctness of the single-component problem is checked against exact references (G5.4), and the quality of the solutions against other packages is shown in the benchmarking article (G1.6).
#' @srrstatsNA {G5.4c} No published stored values exist for this method on shared data sets; see G5.4b.
#' @srrstatsNA {G5.10} All tests run in under two minutes; no extended test suite is required.
#' @srrstatsNA {G5.11} No test requires a large external data set.
#' @srrstatsNA {G5.11a} No test downloads data.
#' @srrstatsNA {G5.12} There are no extended tests.
#' @srrstatsNA {UL1.2} Row and column names are propagated when present (see UL1.3) but are not required: unnamed inputs yield unnamed loadings, so no assertion or message is issued.
#' @srrstatsNA {UL2.2} The package does not accept missing values; see G2.14.
#' @srrstatsNA {UL2.3} Collinearity is not a problem for the algorithm: the covariance matrix may be rank deficient (for example when n < p), and this case is explicitly supported.
#' @srrstatsNA {UL3.0} The package does not apply sequential group labels (it is not a clustering or partitioning algorithm).
#' @srrstatsNA {UL3.2} The input data has no labels beyond the column names, which are propagated as described in UL1.3.
#' @srrstatsNA {UL3.3} The fitted object contains the loadings; scores for new observations are obtained by projecting them on `x_best` (illustrated in the examples of `?mspca`). A dedicated predict method is not provided.
#' @srrstatsNA {UL3.4} The package does not assign data to discrete groups.
#' @srrstatsNA {UL4.1} Creating an unfitted model object is not needed: fitting is fast and non-iterative batch control is not provided.
#' @srrstatsNA {UL6.0} The sparse loadings are communicated through `print()` and `summary()`; a plot method is not provided in this version.
#' @srrstatsNA {UL6.1} No plot method is provided; see UL6.0.
#' @srrstatsNA {UL6.2} No plot method is provided; see UL6.0.
#' @srrstatsNA {UL7.2} The package does not apply group labels; see UL3.0.
#' @srrstatsNA {UL7.4} No predict routine is provided; see UL3.3.
#' @srrstatsNA {UL7.5} The package has no batch-processing routines.
#' @srrstatsNA {UL7.5a} The package has no batch-processing routines.
#' @noRd
NULL
