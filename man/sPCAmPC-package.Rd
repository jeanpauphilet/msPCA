\name{sPCAmPC-package}
\alias{sPCAmPC}
\docType{package}
\title{
\packageTitle{sPCAmPC}
}
\description{
\packageDescription{sPCAmPC}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{sPCAmPC}
\packageIndices{sPCAmPC}
 An overview of how to use the package, including the most important functions
}
\author{
\packageAuthor{sPCAmPC}

Maintainer: \packageMaintainer{sPCAmPC}
}
\references{
This package implements one of the methods proposed by:
Cory-Wright, R., & Pauphilet, J. (2022). Sparse PCA With Multiple Components. arXiv preprint arXiv:2209.14790.

If you find this package useful for your work, please cite
Cory-Wright, R., & Pauphilet, J. (2024). sPCAmPC: An R Package for Sparse PCA with Multiple Principal Components.

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
