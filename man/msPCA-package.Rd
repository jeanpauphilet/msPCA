\name{msPCA-package}
\alias{msPCA}
\docType{package}
\title{
\packageTitle{msPCA}
}
\description{
\packageDescription{msPCA}
}
\details{
\packageIndices{msPCA}

The DESCRIPTION file:
\packageDESCRIPTION{msPCA}
}
\author{
\packageAuthor{msPCA}

Maintainer: \packageMaintainer{msPCA}
}
\references{
This package implements one of the methods proposed by:
Cory-Wright, R., & Pauphilet, J. (2022). Sparse PCA With Multiple Components. arXiv preprint arXiv:2209.14790.

If you find this package useful for your work, please cite
Cory-Wright, R., & Pauphilet, J. (2024). msPCA: An R Package for Sparse PCA with Multiple Principal Components.

}
\keyword{principal component analysis; sparsity; interpretability; penalized optimization}
\seealso{
There are no other man pages.
}
\examples{
library(datasets)
TestMat <- cor(datasets::mtcars)
TestRound <- 2
TestKS <- c(4, 4)

iterativeDeflationHeuristic(TestMat, TestRound, TestKS)
}
