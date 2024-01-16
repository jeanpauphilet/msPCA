\name{sPCAmPC-Functions}
\alias{spca}
\title{Set of functions used in sPCAmPC package}
\description{
  I only wrote the description for cpp_findmultPCs_deflation.
  The description of other functions can also be added if needed.
  Notice only cpp_findmultPCs_deflation is the interface for usage in R.
}
\usage{
spca(Sigma, r, ks, numIters = 200L, verbose = TRUE, violation_tolerance = 1e-4)
}
\arguments{
  \item{Sigma}{Matrix. The correlation or covariance matrix, whose sparse PCs will be computed.}
  \item{r}{Integer. Number of principal components (PCs) to be computed.}
  \item{ks}{List of integers. Target sparsity of each PC.}
  \item{numIters}{Integer. Maximum number of iterations of the algorithm. Default 200.}
  \item{verbose}{Boolean. Default TRUE.}
  \item{violation_tolerance}{Float. Tolerance for the violation of the orthogonality constraints.}
}
\value{
  \code{spca(Sigma, r, ks, numIters, verbose, violation_tolerance)} will return a value to indicate something, and may output some message on the console.

}
\details{
  Introduce these functions more.
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
\author{Dr. Jean Pauphilet}

