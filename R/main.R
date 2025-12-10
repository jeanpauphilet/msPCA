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
#' @param maxIter (optional) An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param verbose (optional) A Boolean. Controls console output. Default TRUE.
#' @param feasibilityConstraintType (optional) An integer. Type of feasibility constraints to be enforced. 0: orthogonality constraints; 1: uncorrelatedness constraints. Default 0.
#' @param violationTolerance (optional) A float. Tolerance for the violation of the orthogonality constraints. Default 1e-4
#' @param stallingTolerance (optional) A float. Controls the objective improvement below which the algorithm is considered to have stalled. Default 1e-8
#' @param maxIterTPW (optional) An integer. Maximum number of iterations of the truncated power method (inner iteration). Default 200.
#' @param  timeLimitTPW (optional) An integer. Maximum time in seconds for the truncated power method (inner iteration). Default 20.
#' @return An object with 4 fields: `x_best` (p x r array containing the sparse PCs), `objective_value`, `feasibility_violation`, `runtime`.
#' @examples
#' library(datasets)
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
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' tpw(TestMat, 4)
tpw <- truncatedPowerMethod


## 2 - Useful functions
#' Variance explained per PC
#'
#' Computes the variance explained by each PC.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return An array.
variance_explained_perPC <- function(C, U){
  r <- dim(U)[2]
  ve <- numeric(r)
  for (i in 1:r){
    ve[i] <- t(U[, i]) %*% C %*% U[, i]
  }
  ve
}

#' Fraction of variance explained per PC
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

#' Fraction of variance explained
#'
#' Computes the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix) by a set of PCs.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return A float.
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4))
#' fraction_variance_explained(TestMat, mspcares$x_best)
fraction_variance_explained <- function(C, U){
  sum(fraction_variance_explained_perPC(C, U))
}

#' Orthogonality constraint violation
#'
#' Computes the orthogonality constraint violation defined as the distance (infinity norm) between \eqn{U^\top U} and the identity matrix.
#' @param U A matrix. Each column correspond to an p-dimensional PC.
#' @return A float.
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4))
#' orthogonality_violation(mspcares$x_best)
orthogonality_violation <- function(U){
  sum(abs(t(U) %*% U - diag(dim(U)[2])))
}


#' Pairwsise correlation
#'
#' Computes the pairwise correlations between PCs defined as \eqn{u_{t}^\top C u_{s}}.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. Each column correspond to an p-dimensional PC.
#' @return A float matrix (r x r).
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4))
#' pairwise_correlation(TestMat, mspcares$x_best)
pairwise_correlation <- function(C, U){
  t(U) %*% C %*% U
}
#' Print mspca output
#'
#' Displays the output of the msPCA algorithm.
#' @param sol_object A list. The output of the mspca or twp function.
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @return None. Prints output to console.
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' mspcares <- mspca(TestMat, 2, c(4,4))
#' print_mspca(mspcares, TestMat)
print_mspca <- function(sol_object, C){
  cat("\nmsPCA solution:\n")
  r <- dim(sol_object$x_best)[2] 
  cat(paste(r,"sparse PCs",sep=" "), "\n")
  
  fve <- fraction_variance_explained_perPC(C, sol_object$x_best)
  cat("Pct. of variance explained:", format(round(fve, 3)*100), "\n")
  
  v <- sol_object$x_best
  k <- 1:r
  for (j in 1:r){
    k[j] <- sum(v[,j] != 0)
  }
  row.names(v) <- row.names(C)
  
  union_of_supports <- rowSums(abs(v) > 0) > 0

  cat("Num. of non-zero loadings : ", k, "\n")
  cat("Sparse PCs \n")
  print(round(v[union_of_supports,],3))
  
}
