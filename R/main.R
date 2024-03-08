## This files contains all the functions of the msPCA package. 
## First, it interfaces all the functions defined through RcppExports.R. In particular, it provides user-friendly names and documentation
## Second, it defines a series of useful function for model evaluation/inspection


## 1 - Interface with RcppExports
#' Returns multiple sparse principal component of a matrix using an iterative deflation heuristic
#'
#' @param Sigma A matrix. The correlation or covariance matrix, whose sparse PCs will be computed.
#' @param r An integer. Number of principal components (PCs) to be computed.
#' @param ks A l ist of integers. Target sparsity of each PC.
#' @param maxIter An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param verbose A Boolean. Controls console output. Default TRUE.
#' @param violationTolerance A float. Tolerance for the violation of the orthogonality constraints. Default 1e-4
#' @param stallingTolerance A float. Controls the objective improvement below which the algorithm is considered to have stalled. Default 1e-8
#' @param maxIterTPW An integer. Maximum number of iterations of the truncated power method (inner iteration). Default 200.
#' @param  timeLimitTPW An integer. Maximum time in seconds for the truncated power method (inner iteration). Default 20.
#' @return An object with 4 fields: x_best (array containing the sparse PCs), objective_value, orthogonality_violation, runtime.
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' mspca(TestMat, 2, c(4,4))
mspca <- iterativeDeflationHeuristic 

#' Returns the leading sparse principal component of a matrix using the truncated power method
#'
#' @param Sigma A matrix. The correlation or covariance matrix, whose sparse PCs will be computed.
#' @param k An integer. Target sparsity of the PC.
#' @param maxIter An integer. Maximum number of iterations of the algorithm. Default 200.
#' @param verbose A Boolean. Controls console output. Default TRUE.
#' @param timeLimit An integer. Maximum time in seconds. Default 10.
#' @return An object with 3 fields: x_best (array containing the sparse PCs), objective_value, runtime.
#' @examples
#' library(datasets)
#' TestMat <- cor(datasets::mtcars)
#' tpw(TestMat, 4)
tpw <- truncatedPowerMethod 


## 2 - Useful functions
variance_explained <- function(C, U){
  fve = 0
  k = dim(U)[2]
  for (i in 1:k){
    fve = fve + U[,i] %*% C %*% U[,i]
  }
  fve
}

#' Compute the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix)
#' @param Sigma A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the PCs (p x k). 
#' @return A float. 
fraction_variance_explained <- function(C, U){
  fve = variance_explained(C,U)
  fve = fve / sum(diag(C))
  fve
}

#' Computes the orthogonality constraint violation defined as the distance (L infinity norm) between $U^T U$ and the identity matrix
#' @param U A matrix. Each column correspond to an n-dimensional PC. 
#' @return A float.
orthogonality_violation <- function(U){
  sum(abs(t(U) %*% U - diag(dim(U)[2])))
}

#' Compresses a dataset according to new coordinates
#'
#' @param df A matrix. The original dataframe. Each row corresponds to an observation (n rows, p columns) 
#' @param U A matrix. The new coordinates (PC). Each column corresponds to a PC (p rows, k columns)
#' @return A matrix with n rows and k columns. 
compress_data <- function(df, U){
  data.frame(as.matrix(df) %*% U)
}

