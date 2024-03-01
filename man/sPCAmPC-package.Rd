\name{msPCA-package}
\alias{msPCA}
\alias{msPCA-package}
\docType{package}
\title{
\packageTitle{msPCA}
}
\description{
\packageDescription{msPCA}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{msPCA}
\packageIndices{msPCA}
 An overview of how to use the package, including the most important functions
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

spca(TestMat, TestRound, TestKS)
}
