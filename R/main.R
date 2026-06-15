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
#' @param ks A list of integers. Target sparsity of each PC.
#' @param type (optional) Either "Sigma" (default; `M` is a covariance/correlation matrix) or "X" (`M` is a raw data matrix).
#' @param feasibilityConstraintType (optional) An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints. Default 0.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param feasibilityTolerance (optional) A float. Tolerance for the violation of the orthogonality constraints. Default 1e-4
#' @param stallingTolerance (optional) A float. Controls the objective improvement below which the algorithm is considered to have stalled. Default 1e-8
#' @param timeLimitTPM (optional) An integer. Maximum time in seconds for the truncated power method (inner iteration). Default 20.
#' @param maxRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for the first outer iteration. Default 20.
#' @param minRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for outer iterations >= 2. Default 20.
#' @param center (optional, type = "X") A Boolean. Center the columns of `M` before computing the covariance. Default TRUE.
#' @param scale (optional, type = "X") A Boolean. Scale the columns of `M` to unit variance, i.e. operate on the correlation matrix. Default TRUE.
#' @param divisor (optional, type = "X") Either "n-1" (default, sample covariance, matches `cov`/`cor`) or "n" (population covariance). Default "n-1".
#' @param checkPSD (optional, type = "Sigma") A Boolean. Verify that `M` is positive semidefinite. Default TRUE.
#' @param symTolerance (optional, type = "Sigma") A float. Tolerance for the symmetry  check on `M`. Default 1e-8.
#' @param psdTolerance (optional, type = "Sigma") A float. Tolerance (on the smallest eigenvalue) for the PSD check on `M`. Default 1e-8.
#' @return An object with fields: `x_best` (p x r array containing the sparse PCs),
#'   `objective_value`, `feasibility_violation`, `runtime`, `variance_explained`
#'   (per-PC explained variance) and `total_variance` (trace of the covariance).
#'   With `type = "X"` it additionally records `inputType`, `center`, `scale`,
#'   `divisor`, `nObs` and `p`.
#' @examples
#' # From a covariance/correlation matrix (the default type):
#' TestMat <- cor(datasets::mtcars)
#' mspca(TestMat, 2, c(4,4))
#' # Equivalent call from the raw data matrix:
#' mspca(scale(datasets::mtcars), r = 2, ks = c(4,4), type = "X", verbose = FALSE)
mspca <- function(M, r, ks, type = c("Sigma", "X"),
                  feasibilityConstraintType = 0, verbose = TRUE, maxIter = 200,
                  feasibilityTolerance = 1e-4, stallingTolerance = 1e-8,
                  timeLimitTPM = 20, maxRestartTPM = 20, minRestartTPM = 20,
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
#' @return An object with fields: `x_best` (p x 1 array containing the sparse PC),
#'   `objective_value`, `runtime`. With `type = "X"` it additionally records
#'   `inputType`, `center`, `scale`, `divisor`, `nObs` and `p`.
#' @references Yuan, X. T., & Zhang, T. (2013). Truncated power method for sparse eigenvalue problems. The Journal of Machine Learning Research, 14(1), 899-925.
#' @examples
#' TestMat <- cor(datasets::mtcars)
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
  res
}


## 2 - Useful functions
#' Variance Explained per PC
#'
#' Computes the variance explained by each PC.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return An array.
variance_explained_perPC <- function(C, U){
  colSums(U * (C %*% U))
}

#' Fraction of Variance Explained per PC
#'
#' Computes the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix) by each PC.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return An array.
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
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4), type = "Sigma")
#' fraction_variance_explained(TestMat, mspcares$x_best)
fraction_variance_explained <- function(C, U){
  sum(U * (C %*% U)) / sum(diag(C))
}

#' Feasibility Violation
#'
#' Computes the feasibility violation defined as \eqn{\sum_{t > s} u_{t}^\top u_{s}} if orthogonality constraints are enforced (feasibilityConstraintType = 0) and \eqn{\sum_{t > s} u_{t}^\top C u_{s}} if zero-correlation constraints are enforced (feasibilityConstraintType = 1).
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. Each column correspond to an p-dimensional PC.
#' @param feasibilityConstraintType An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints.
#' @return A float.
#' @examples
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4), type = "Sigma")
#' feasibility_violation_off(TestMat, mspcares$x_best, 0)
feasibility_violation_off <- function(C, U, feasibilityConstraintType){
  M = if (feasibilityConstraintType == 0) {
    crossprod(U)
  } else {
   crossprod(U, C %*% U)
  }
  sum(abs(M[upper.tri(M, diag=FALSE)]))
}

#' Print msPCA Output
#'
#' Displays the output of the msPCA algorithm. When the model was fit from a
#' covariance/correlation matrix (`type = "Sigma"`), pass that matrix as `C`; when
#' it was fit from a raw data matrix (`type = "X"`), `C` can be omitted, in which
#' case the explained-variance figures stored in the result object are used.
#' @param sol_object A list. The output of the mspca or tpm function.
#' @param C (optional) A matrix. The correlation or covariance matrix (p x p). May
#'   be omitted when `sol_object` already carries `variance_explained` /
#'   `total_variance` (i.e. a `type = "X"` result).
#' @param digits An integer. Number of digits used for rounded display. Default 3.
#' @return None. Prints output to console.
#' @examples
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4), type = "Sigma")
#' print_mspca(mspcares, TestMat, digits = 3)
print_mspca <- function(sol_object, C = NULL, digits = 3){
  cat("\nmsPCA solution:\n")
  v <- sol_object$x_best
  r <- ncol(v)
  cat(paste(r,"sparse PCs",sep=" "), "\n")

  if (!is.null(C)) {
    fve <- fraction_variance_explained_perPC(C, v)
    rn <- row.names(C)
  } else if (!is.null(sol_object$variance_explained) && !is.null(sol_object$total_variance)) {
    fve <- sol_object$variance_explained / sol_object$total_variance
    rn <- row.names(v)
  } else {
    stop("Provide the covariance/correlation matrix `C`, or a result object that contains `variance_explained` and `total_variance`.", call. = FALSE)
  }

  cat("Pct. of variance explained:", format(round(fve, digits) * 100), "\n")

  k <- colSums(v != 0)
  if (!is.null(rn)) row.names(v) <- rn

  union_of_supports <- rowSums(v != 0) > 0

  cat("Num. of non-zero loadings : ", k, "\n")
  cat("Sparse PCs \n")
  print(round(v[union_of_supports, , drop = FALSE], digits))

}
