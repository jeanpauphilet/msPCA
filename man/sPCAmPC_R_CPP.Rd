\name{msPCA-Functions}
\alias{spca}
\title{Set of functions used in msPCA package}
\description{
  TO DO: Add description
}
\usage{
spca(Sigma, r, ks, numIters = 200L, verbose = TRUE, violation_tolerance = 1e-4)
}

\arguments{
  \item{Sigma}{Matrix. The correlation or covariance matrix, whose sparse PCs will be computed.}
  \item{r}{Integer. Number of principal components (PCs) to be computed.}
  \item{ks}{List of integers. Target sparsity of each PC.}
  \item{numIters}{Integer. Maximum number of iterations of the algorithm. Default 200.}
  \item{verbose}{Boolean. Output to be displayed on the console. Default TRUE.}
  \item{violation_tolerance}{Float. Tolerance for the violation of the orthogonality constraints.}
}

\value{
  \code{spca(Sigma, r, ks, numIters, verbose, violation_tolerance)} returns an object with the following fields:
  \item{x_best}{Matrix. An n times r matrix containing the r sparse components.}
  \item{objective_valute}{Float. Value of x_best.T*Sigma*x_best.}
  \item{orthogonality_violation}{Float. Distance (L1) between x_best.T*x_best and the identity matrix.}
  \item{runtime}{Float. Time taken by the algorithm (in seconds).}
  
}

\references{
This function implements the iterative deflation heuristic (Algorithm 2) of:
Cory-Wright, R., & Pauphilet, J. (2022). Sparse PCA With Multiple Components. arXiv preprint arXiv:2209.14790.
}

\examples{
  ## these are random numbers
library(datasets)
TestMat <- cor(datasets::mtcars)
TestRound <- 2
TestKS <- c(4, 4)

spca(TestMat, TestRound, TestKS)
}

