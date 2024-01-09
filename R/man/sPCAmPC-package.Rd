\name{sPCAmPC-package}
\alias{sPCAmPC-package}
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
Literature or other references for background information
}
\keyword{ package }
\seealso{
There are no other man pages.
}
\examples{
## these are random numbers
TestMat <- matrix(
  c(0.907247, -0.00283434, 0.0259237, -0.0184396, -0.0179678,
    -0.0593399, 0.0556192, -0.0373431, -0.0370315, -0.00805812,
    -0.00283434, 1.24395, -0.310638, -0.000514658, -0.00193052,
    -0.0119863, 0.0160602, 0.215909, -0.0671475, 0.272086,
    0.0259237, -0.310638, 1.24793, -0.00777525, 0.0349503,
    0.0351923, -0.0351558, -0.191682, 0.0315086, -0.290368,
    -0.0184396, -0.000514658, -0.00777525, 1.26401, -0.0102385,
    0.264942, 0.257991, -0.0007677, -0.226642, 0.0608659,
    -0.0179678, -0.00193052, 0.0349503, -0.0102385, 0.961997,
    0.00842352, 0.0183949, 0.00120234, 0.0351232, 0.036898,
    -0.0593399, -0.0119863, 0.0351923, 0.264942, 0.00842352,
    1.21353, 0.287146, -0.00727776, -0.205039, -0.0587133,
    0.0556192, 0.0160602, -0.0351558, 0.257991, 0.0183949,
    0.287146, 1.28508, 0.0155203, -0.260995, -0.0150086,
    -0.0373431, 0.215909, -0.191682, -0.0007677, 0.00120234,
    -0.00727776, 0.0155203, 1.28349, 0.0283306, 0.273234,
    -0.0370315, -0.0671475, 0.0315086, -0.226642, 0.0351232,
    -0.205039, -0.260995, 0.0283306, 1.21418, -0.0375033,
    -0.00805812, 0.272086, -0.290368, 0.0608659, 0.036898,
    -0.0587133, -0.0150086, 0.273234, -0.0375033, 1.18052),
  nrow = 10,
  ncol = 10,
  byrow = TRUE
)
TestRound <- 2
TestKS <- matrix( c(4, 4), nrow = 1, ncol = 2, byrow = TRUE )
TestNumIters <- 20

cpp_findmultPCs_deflation(TestMat, TestRound, TestKS, TestNumIters)
}
