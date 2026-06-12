## This files contains all the functions of the msPCA package. 
## First, it interfaces all the functions defined through RcppExports.R. In particular, it provides user-friendly names and documentation.
## Second, it defines a series of useful function for model evaluation/inspection.


## 1 - Interface with RcppExports
#' Multiple Sparse PCA
#'
#' Returns multiple sparse principal component of a matrix using an iterative deflation heuristic.
#' @param Sigma A matrix. The correlation or covariance matrix, whose sparse PCs will be computed.
#' @param r An integer. Number of principal components (PCs) to be computed.
#' @param ks A list of integers. Target sparsity of each PC.
#' @param feasibilityConstraintType (optional) An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints. Default 0.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param feasibilityTolerance (optional) A float. Tolerance for the violation of the orthogonality constraints. Default 1e-4
#' @param stallingTolerance (optional) A float. Controls the objective improvement below which the algorithm is considered to have stalled. Default 1e-8
#' @param timeLimitTPM (optional) An integer. Maximum time in seconds for the truncated power method (inner iteration). Default 20.
#' @param maxRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for the first outer iteration. Default 20.
#' @param minRestartTPM (optional) An integer. Number of random restarts of the truncated power method (inner iteration) for outer iterations >= 2. Default 10.
#' @return An object with 4 fields: `x_best` (p x r array containing the sparse PCs), `objective_value`, `feasibility_violation`, `runtime`.
#' @examples
#' TestMat <- cor(datasets::mtcars)
#' mspca(TestMat, 2, c(4,4))
mspca <- iterativeDeflationHeuristic

#' Truncated Power Method
#'
#' Returns the leading sparse principal component of a matrix using the truncated power method.
#' @param Sigma A matrix. The correlation or covariance matrix, whose sparse PCs will be computed.
#' @param k An integer. Target sparsity of the PC.
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param timeLimit (optional) An integer. Maximum time in seconds. Default 10.
#' @return An object with 3 fields: `x_best` (p x 1 array containing the sparse PC), `objective_value`, `runtime`.
#' @references Yuan, X. T., & Zhang, T. (2013). Truncated power method for sparse eigenvalue problems. The Journal of Machine Learning Research, 14(1), 899-925.
#' @examples
#' TestMat <- cor(datasets::mtcars)
#' tpm(TestMat, 4)
tpm <- truncatedPowerMethod


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
#' mspcares <- mspca(TestMat, 2, c(4,4))
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
#' mspcares <- mspca(TestMat, 2, c(4,4))
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
#' Displays the output of the msPCA algorithm.
#' @param sol_object A list. The output of the mspca or tpm function.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param digits An integer. Number of digits used for rounded display. Default 3.
#' @return None. Prints output to console.
#' @examples
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4))
#' print_mspca(mspcares, TestMat, digits = 3)
print_mspca <- function(sol_object, C, digits = 3){
  cat("\nmsPCA solution:\n")
  v <- sol_object$x_best
  r <- ncol(v)
  cat(paste(r,"sparse PCs",sep=" "), "\n")
  
  fve <- fraction_variance_explained_perPC(C, v)
  cat("Pct. of variance explained:", format(round(fve, digits) * 100), "\n")
  
  k <- colSums(v != 0)
  row.names(v) <- row.names(C)
  
  union_of_supports <- rowSums(v != 0) > 0

  cat("Num. of non-zero loadings : ", k, "\n")
  cat("Sparse PCs \n")
  print(round(v[union_of_supports, , drop = FALSE], digits))
  
}
