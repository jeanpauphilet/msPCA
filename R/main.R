## This files contains all the functions of the msPCA package.
## First, it interfaces all the functions defined through RcppExports.R. In particular, it provides user-friendly names and documentation.
## Second, it defines a series of useful function for model evaluation/inspection.


## 0 - Internal helpers (not exported) ------------------------------------------

# Validate a covariance/correlation matrix (type = "Sigma"): square, finite,
# symmetric, (optionally) PSD.
.validate_sigma_matrix <- function(M, symTolerance = 1e-8, psdTolerance = 1e-8, checkPSD = TRUE) {
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("With `type = \"Sigma\"`, `M` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(M) != ncol(M)) {
    stop("With `type = \"Sigma\"`, `M` must be a square matrix.", call. = FALSE)
  }
  if (any(!is.finite(M))) {
    stop("`M` contains non-finite values (NA, NaN or Inf).", call. = FALSE)
  }
  if (max(abs(M - t(M))) > symTolerance) {
    stop("With `type = \"Sigma\"`, `M` is not symmetric within `symTolerance`.", call. = FALSE)
  }
  if (isTRUE(checkPSD)) {
    smallest <- min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)
    if (smallest < -psdTolerance) {
      stop(sprintf(paste0("With `type = \"Sigma\"`, `M` is not positive semidefinite (smallest eigenvalue %.3e < -psdTolerance). ",
                          "Pass a valid covariance/correlation matrix or set `checkPSD = FALSE` to override."),
                  smallest), call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Validate a raw data matrix (type = "X"): n observations x p variables.
.validate_x_matrix <- function(M, scale = FALSE) {
  if (!is.matrix(M) || !is.numeric(M)) {
    stop("With `type = \"X\"`, `M` must be a numeric matrix (rows = observations, columns = variables).", call. = FALSE)
  }
  if (any(!is.finite(M))) {
    stop("`M` contains non-finite values (NA, NaN or Inf).", call. = FALSE)
  }
  if (nrow(M) < 2L) {
    stop("With `type = \"X\"`, `M` must have at least 2 rows (observations).", call. = FALSE)
  }
  if (isTRUE(scale)) {
    col_var <- apply(M, 2, stats::var)
    if (any(col_var == 0)) {
      stop("`M` has zero-variance column(s); cannot scale to a correlation matrix. Set `scale = FALSE` or drop these columns.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# Center/scale the data matrix (cheap O(np) operations on the data, never the
# p x p matrix) and return the processed matrix together with the covariance
# divisor 1/(n-1) or 1/n. With center = TRUE, scale = FALSE, divisor = "n-1" the
# implicit covariance invDivisor * Xproc^T Xproc equals cov(M); with scale = TRUE
# it equals cor(M).
.preprocess_x <- function(M, center = TRUE, scale = FALSE, divisor = c("n-1", "n")) {
  divisor <- match.arg(divisor)
  nObs <- nrow(M)
  Xproc <- base::scale(M, center = center, scale = scale)
  # Drop the attributes added by scale() and ensure a plain numeric matrix.
  Xproc <- matrix(as.numeric(Xproc), nrow = nObs, ncol = ncol(M),
                  dimnames = list(rownames(M), colnames(M)))
  invDivisor <- if (divisor == "n-1") 1 / (nObs - 1) else 1 / nObs
  list(Xproc = Xproc, invDivisor = invDivisor, divisor = divisor)
}


## 1 - Interface with RcppExports ----------------------------------------------
#' Multiple Sparse PCA
#'
#' Returns multiple sparse principal components of a dataset using an iterative
#' deflation heuristic. As in the `elasticnet` package, the data is passed as a
#' single argument `M` whose interpretation is set by `type`: `"Sigma"` (the
#' default) treats `M` as a covariance/correlation matrix (p x p) and `"X"` treats
#' `M` as a raw data matrix (n observations x p variables). With `type = "X"` the
#' algorithm operates on the data directly via the products \eqn{X^\top(X\beta)}
#' and never forms the p x p matrix, which is substantially more scalable when
#' \eqn{n \ll p}.
#' @param M A matrix. The data, interpreted according to `type`: a covariance/
#'   correlation matrix (p x p) when `type = "Sigma"`, or a raw data matrix
#'   (n x p) when `type = "X"`.
#' @param r An integer. Number of principal components (PCs) to be computed.
#' @param ks An integer vector. Target sparsity of each PC.
#' @param type (optional) Either "Sigma" (default; `M` is a covariance/correlation matrix) or "X" (`M` is a raw data matrix).
#' @param feasibilityConstraintType (optional) An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints. Default 0.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param feasibilityTolerance (optional) A float. Tolerance for constraint violation (orthogonality/uncorrelatedness, according to `feasibilityConstraintType`). Default 1e-4.
#' @param stallingTolerance (optional) A float. Controls the objective improvement below which the algorithm is considered to have stalled. Default 1e-8.
#' @param timeLimitTPM (optional) An integer. Maximum time in seconds for the truncated power method (inner iteration). Default 20.
#' @param maxRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for the first outer iteration. Default 30.
#' @param minRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for outer iterations >= 2. Default 20.
#' @param center (optional, type = "X") A Boolean. Center the columns of `M` before computing the covariance. Default TRUE.
#' @param scale (optional, type = "X") A Boolean. Scale the columns of `M` to unit variance, i.e. operate on the correlation matrix. Default TRUE.
#' @param divisor (optional, type = "X") Either "n-1" (default, sample covariance, matches `cov`/`cor`) or "n" (population covariance). Default "n-1".
#' @param checkPSD (optional, type = "Sigma") A Boolean. Verify that `M` is positive semidefinite. Default TRUE.
#' @param symTolerance (optional, type = "Sigma") A float. Tolerance for the symmetry  check on `M`. Default 1e-8.
#' @param psdTolerance (optional, type = "Sigma") A float. Tolerance (on the smallest eigenvalue) for the PSD check on `M`. Default 1e-8.
#' @return An object of class `"mspca"` (a list) with fields: `x_best` (p x r
#'   matrix of sparse PC loadings), `objective_value`, `feasibility_violation`,
#'   `runtime`, `variance_explained` (per-PC explained variance), and
#'   `total_variance` (trace of the covariance matrix). With `type = "X"` it
#'   additionally records `inputType`, `center`, `scale`, `divisor`, `nObs`,
#'   and `p`. Use [print()] to display the sparse loadings and [summary()] for
#'   a full per-PC breakdown.
#' @examples
#' # From a covariance/correlation matrix (the default type):
#' TestMat <- cor(mtcars)
#' res <- mspca(TestMat, r = 2, ks = c(4, 4), verbose = FALSE)
#' print(res, TestMat)
#' # Equivalent call from the raw data matrix (C need not be passed to print):
#' res_X <- mspca(as.matrix(mtcars), r = 2, ks = c(4, 4), type = "X", verbose = FALSE)
#' print(res_X)
mspca <- function(M, r, ks, type = c("Sigma", "X"),
                  feasibilityConstraintType = 0, verbose = TRUE, maxIter = 200,
                  feasibilityTolerance = 1e-4, stallingTolerance = 1e-8,
                  timeLimitTPM = 20, maxRestartTPM = 30, minRestartTPM = 20,
                  center = TRUE, scale = TRUE, divisor = c("n-1", "n"),
                  checkPSD = TRUE, symTolerance = 1e-8, psdTolerance = 1e-8) {
  type <- match.arg(type)
  if (type == "Sigma") {
    .validate_sigma_matrix(M, symTolerance, psdTolerance, checkPSD)
    res <- iterativeDeflationHeuristic(M, r, ks, feasibilityConstraintType,
                                       verbose, maxIter, feasibilityTolerance,
                                       stallingTolerance, timeLimitTPM,
                                       maxRestartTPM, minRestartTPM)
    res$inputType <- "Sigma"
  } else {
    .validate_x_matrix(M, scale)
    pp <- .preprocess_x(M, center, scale, divisor)
    res <- iterativeDeflationHeuristicX(pp$Xproc, pp$invDivisor, r, ks,
                                        feasibilityConstraintType, verbose, maxIter,
                                        feasibilityTolerance, stallingTolerance,
                                        timeLimitTPM, maxRestartTPM, minRestartTPM)
    if (!is.null(colnames(M))) rownames(res$x_best) <- colnames(M)
    res$inputType <- "X"
    res$center <- center
    res$scale <- scale
    res$divisor <- pp$divisor
    res$nObs <- nrow(M)
    res$p <- ncol(M)
  }
  class(res) <- "mspca"
  res
}

#' Truncated Power Method
#'
#' Returns the leading sparse principal component of a dataset using the truncated
#' power method. As in [mspca()], the data is passed as a single argument `M` whose
#' interpretation is set by `type`: `"Sigma"` (default) for a covariance/correlation
#' matrix (p x p) or `"X"` for a raw data matrix (n x p). See [mspca()] for the
#' raw-data preprocessing controls.
#' @param M A matrix. The data, interpreted according to `type`: a covariance/
#'   correlation matrix (p x p) when `type = "Sigma"`, or a raw data matrix
#'   (n x p) when `type = "X"`.
#' @param k An integer. Target sparsity of the PC.
#' @param type (optional) Either "Sigma" (default; `M` is a covariance/correlation
#'   matrix) or "X" (`M` is a raw data matrix).
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param timeLimit (optional) An integer. Maximum time in seconds. Default 10.
#' @param center (optional, type = "X") A Boolean. Center the columns of `M`. Default TRUE.
#' @param scale (optional, type = "X") A Boolean. Scale the columns of `M` to unit variance. Default FALSE.
#' @param divisor (optional, type = "X") Either "n-1" (default) or "n".
#' @param checkPSD (optional, type = "Sigma") A Boolean. Verify `M` is PSD. Default TRUE.
#' @param symTolerance (optional, type = "Sigma") A float. Symmetry-check tolerance. Default 1e-8.
#' @param psdTolerance (optional, type = "Sigma") A float. PSD-check tolerance. Default 1e-8.
#' @return An object of class `"tpm"` (a list) with fields: `x_best` (p x 1
#'   matrix containing the sparse PC loading), `objective_value`, and `runtime`.
#'   With `type = "X"` it additionally records `inputType`, `center`, `scale`,
#'   `divisor`, `nObs`, and `p`.
#' @references Yuan, X. T., & Zhang, T. (2013). Truncated power method for
#'   sparse eigenvalue problems. \emph{The Journal of Machine Learning
#'   Research}, \bold{14}(1), 899--925.
#' @examples
#' TestMat <- cor(mtcars)
#' tpm(TestMat, 4)
tpm <- function(M, k, type = c("Sigma", "X"),
                maxIter = 200, verbose = TRUE, timeLimit = 10,
                center = TRUE, scale = FALSE, divisor = c("n-1", "n"),
                checkPSD = TRUE, symTolerance = 1e-8, psdTolerance = 1e-8) {
  type <- match.arg(type)
  if (type == "Sigma") {
    .validate_sigma_matrix(M, symTolerance, psdTolerance, checkPSD)
    res <- truncatedPowerMethod(M, k, maxIter, verbose, timeLimit)
    res$inputType <- "Sigma"
  } else {
    .validate_x_matrix(M, scale)
    pp <- .preprocess_x(M, center, scale, divisor)
    res <- truncatedPowerMethodX(pp$Xproc, pp$invDivisor, k, maxIter, verbose, timeLimit)
    if (!is.null(colnames(M))) rownames(res$x_best) <- colnames(M)
    res$inputType <- "X"
    res$center <- center
    res$scale <- scale
    res$divisor <- pp$divisor
    res$nObs <- nrow(M)
    res$p <- ncol(M)
  }
  class(res) <- "tpm"
  res
}


## 2 - Useful functions
#' Variance Explained Per PC
#'
#' Computes the variance explained by each PC.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return A numeric vector of length `ncol(U)`.
variance_explained_perPC <- function(C, U){
  colSums(U * (C %*% U))
}

#' Fraction of Variance Explained Per PC
#'
#' Computes the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix) by each PC.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return A numeric vector of length `ncol(U)`.
fraction_variance_explained_perPC <- function(C, U){
  fve <- variance_explained_perPC(C, U)
  fve <- fve / sum(diag(C))
  fve
}

#' Fraction of Variance Explained
#'
#' Computes the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix) by a set of PCs.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return A float.
#' @examples
#' TestMat <- cor(mtcars)
#' mspcares <- mspca(TestMat, r = 2, ks = c(4, 4), verbose = FALSE)
#' fraction_variance_explained(TestMat, mspcares$x_best)
fraction_variance_explained <- function(C, U){
  sum(U * (C %*% U)) / sum(diag(C))
}

#' Feasibility Violation
#'
#' Computes the feasibility violation defined as
#' \eqn{\sum_{t > s} |u_{t}^\top u_{s}|} if orthogonality constraints are
#' enforced (`feasibilityConstraintType = 0`) and
#' \eqn{\sum_{t > s} |u_{t}^\top C u_{s}|} if zero-correlation constraints are
#' enforced (`feasibilityConstraintType = 1`).
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. Each column corresponds to a p-dimensional PC.
#' @param feasibilityConstraintType An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints.
#' @return A float.
#' @examples
#' TestMat <- cor(mtcars)
#' mspcares <- mspca(TestMat, r = 2, ks = c(4, 4), verbose = FALSE)
#' feasibility_violation_off(TestMat, mspcares$x_best, 0)
feasibility_violation_off <- function(C, U, feasibilityConstraintType){
  M = if (feasibilityConstraintType == 0) {
    crossprod(U)
  } else {
   crossprod(U, C %*% U)
  }
  sum(abs(M[upper.tri(M, diag=FALSE)]))
}

#' Print an mspca Object
#'
#' S3 print method for objects of class `"mspca"` returned by [mspca()].
#' Displays the sparse loading matrix (restricted to the union of non-zero
#' rows) together with the percentage of variance explained and the number of
#' non-zero loadings per component.
#'
#' When the model was fit from a covariance/correlation matrix
#' (`type = "Sigma"`), pass that matrix as `C` so that per-PC variance figures
#' can be computed; when it was fit from a raw data matrix (`type = "X"`), `C`
#' may be omitted because the figures are stored inside the object.
#'
#' @param x An object of class `"mspca"`, as returned by [mspca()].
#' @param C (optional) A numeric matrix (p x p). The covariance or correlation
#'   matrix used when fitting. May be omitted for `type = "X"` results.
#' @param digits An integer or `NULL`. Number of significant digits for
#'   display. When `NULL` (the default), \code{getOption("digits")} is used,
#'   so the output respects \code{options(digits = ...)}.
#' @param ... Further arguments required by the `print()` generic; not used by
#'   this method.
#' @return Invisibly returns `x`.
#' @method print mspca
#' @export
#' @examples
#' TestMat <- cor(mtcars)
#' res <- mspca(TestMat, r = 2, ks = c(4, 4), verbose = FALSE)
#' print(res, TestMat)
print.mspca <- function(x, C = NULL, digits = NULL, ...) {
  dg <- if (is.null(digits)) getOption("digits") else digits
  v  <- x$x_best
  r  <- ncol(v)

  if (!is.null(C)) {
    fve <- fraction_variance_explained_perPC(C, v)
    rn  <- rownames(C)
  } else if (!is.null(x$variance_explained) && !is.null(x$total_variance)) {
    fve <- x$variance_explained / x$total_variance
    rn  <- rownames(v)
  } else {
    stop(paste("Provide the covariance/correlation matrix `C`,",
               "or a result object that contains `variance_explained`",
               "and `total_variance`."), call. = FALSE)
  }
  rownames(v) <- if (!is.null(rn)) rn else seq_len(nrow(v))

  cat("\nmsPCA solution:", r, "sparse PCs\n")
  cat("Pct. variance explained:", format(fve * 100, digits = dg), "\n")
  cat("Non-zero loadings per PC:", colSums(v != 0), "\n")
  cat("\nSparse PCs\n")
  print(v[rowSums(v != 0) > 0, , drop = FALSE], digits = dg)
  invisible(x)
}

#' Summarize an mspca Object
#'
#' S3 summary method for objects of class `"mspca"` returned by [mspca()].
#' Returns (and prints) a per-PC summary table (number of non-zero loadings,
#' variance explained, FVE, and cumulative FVE) together with the pairwise
#' feasibility violation matrix and the total solver runtime.
#'
#' @param object An object of class `"mspca"`, as returned by [mspca()].
#' @param C (optional) A numeric matrix (p x p). The covariance or correlation
#'   matrix used when fitting. May be omitted for `type = "X"` results, where
#'   the figures are stored inside the object.
#' @param feasibilityConstraintType An integer. Type of constraint used to
#'   compute the feasibility violation reported in the summary. `0` (default)
#'   for orthogonality; `1` for zero pairwise correlation.
#' @param digits An integer or `NULL`. Number of significant digits for
#'   display. When `NULL` (the default), \code{getOption("digits")} is used.
#' @param ... Further arguments required by the `summary()` generic; not used
#'   by this method.
#' @return Invisibly returns a list of class `"summary.mspca"` with fields:
#'   \describe{
#'     \item{`table`}{Data frame with columns `PC`, `nonzero`, `variance`,
#'       `fve`, and `cumulative_fve`.}
#'     \item{`feasibility_mat`}{r x r matrix of pairwise feasibility violations
#'       (\eqn{|u_i^\top u_j|} or \eqn{|u_i^\top C u_j|}). Diagonal and lower
#'       triangle are `NA`.}
#'     \item{`feasibility`}{Scalar total feasibility violation (sum of the
#'       upper triangle of `feasibility_mat`).}
#'     \item{`runtime`}{Solver runtime in seconds (if stored in the object).}
#'     \item{`r`}{Number of sparse PCs.}
#'     \item{`inputType`}{`"Sigma"` or `"X"`.}
#'   }
#' @method summary mspca
#' @export
#' @examples
#' TestMat <- cor(mtcars)
#' res <- mspca(TestMat, r = 2, ks = c(4, 4), verbose = FALSE)
#' summary(res, TestMat)
summary.mspca <- function(object, C = NULL,
                          feasibilityConstraintType = 0L,
                          digits = NULL, ...) {
  dg <- if (is.null(digits)) getOption("digits") else digits
  v  <- object$x_best
  r  <- ncol(v)

  ## -- Variance figures --
  if (!is.null(C)) {
    var_pc    <- variance_explained_perPC(C, v)
    total_var <- sum(diag(C))
  } else if (!is.null(object$variance_explained) && !is.null(object$total_variance)) {
    var_pc    <- object$variance_explained
    total_var <- object$total_variance
  } else {
    stop(paste("Provide the covariance/correlation matrix `C`,",
               "or a result object that contains `variance_explained`",
               "and `total_variance`."), call. = FALSE)
  }
  fve     <- var_pc / total_var
  cum_fve <- cumsum(fve)

  ## -- Pairwise feasibility matrix (r x r, upper triangle meaningful) --
  pc_labels <- paste0("PC", seq_len(r))
  if (r >= 2L) {
    if (feasibilityConstraintType == 0L) {
      G <- crossprod(v)
    } else {
      G <- crossprod(v, C %*% v)
    }
    feas_mat <- abs(G)
    diag(feas_mat) <- NA_real_         # diagonal is trivially 1; suppress it
    feas_mat[lower.tri(feas_mat)] <- NA_real_
  } else {
    feas_mat <- matrix(NA_real_, 1L, 1L)
  }
  dimnames(feas_mat) <- list(pc_labels, pc_labels)

  tbl <- data.frame(
    PC             = pc_labels,
    nonzero        = colSums(v != 0),
    variance       = var_pc,
    fve            = fve,
    cumulative_fve = cum_fve,
    row.names      = NULL
  )

  total_feas <- feasibility_violation_off(
    if (!is.null(C)) C else diag(nrow(v)),
    v, feasibilityConstraintType)

  out <- structure(
    list(
      table            = tbl,
      feasibility_mat  = feas_mat,
      feasibility      = total_feas,
      runtime          = object$runtime,
      r                = r,
      inputType        = object$inputType
    ),
    class = "summary.mspca"
  )

  ## -- Print --
  cat("\nmsPCA summary:", r, "sparse PC(s)\n")
  if (!is.null(object$inputType))
    cat("Input type   :", object$inputType, "\n")
  if (!is.null(object$runtime))
    cat("Runtime (s)  :", format(object$runtime, digits = dg), "\n")
  cat("\nPer-component statistics:\n")
  fmt_tbl                <- tbl
  fmt_tbl$variance       <- format(tbl$variance,       digits = dg)
  fmt_tbl$fve            <- format(tbl$fve,            digits = dg)
  fmt_tbl$cumulative_fve <- format(tbl$cumulative_fve, digits = dg)
  print(fmt_tbl, row.names = FALSE)
  if (r >= 2L) {
    cat("\nPairwise feasibility violations (upper triangle):\n")
    print(feas_mat, digits = dg, na.print = ".")
    cat("Total:", format(total_feas, digits = dg, scientific = TRUE), "\n")
  }
  invisible(out)
}
