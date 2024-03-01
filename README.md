# msPCA
An R Package for Sparse PCA with Multiple Principal Components

## Installation 
This package can be installed in R directly from the Github repository. To do so, you will first need to install `devtools`:
<br />`install.packages("devtools")`<br />
And then run the following commands: 
<br /> `library(devtools)`<br />
`install_github('jeanpauphilet/msPCA')`<br />

## Getting started
The package consists of one main function, `msPCA`, which takes as input:
- a data matrix (either the correlation or covariance matrix of the dataset),
- the number of principal components (PCs) to be computed, r,
- a list of r integers corresponding to the sparsity of each PC.
  
It returns...

Here is a short example demonstrating how to use the package. First, you need to load the library. 
<br />
`library(msPCA)`<br />
Then, define the input variables.
<br />
`TestMat <- matrix(`<br />
`  c(0.907247, -0.00283434, 0.0259237, -0.0184396, -0.0179678,`<br />
`    -0.0593399, 0.0556192, -0.0373431, -0.0370315, -0.00805812,`<br />
`    -0.00283434, 1.24395, -0.310638, -0.000514658, -0.00193052,`<br />
`    -0.0119863, 0.0160602, 0.215909, -0.0671475, 0.272086,`<br />
`    0.0259237, -0.310638, 1.24793, -0.00777525, 0.0349503,`<br />
`    0.0351923, -0.0351558, -0.191682, 0.0315086, -0.290368,`<br />
`    -0.0184396, -0.000514658, -0.00777525, 1.26401, -0.0102385,`<br />
`    0.264942, 0.257991, -0.0007677, -0.226642, 0.0608659,`<br />
`    -0.0179678, -0.00193052, 0.0349503, -0.0102385, 0.961997,`<br />
`    0.00842352, 0.0183949, 0.00120234, 0.0351232, 0.036898,`<br />
`    -0.0593399, -0.0119863, 0.0351923, 0.264942, 0.00842352,`<br />
`    1.21353, 0.287146, -0.00727776, -0.205039, -0.0587133,`<br />
`    0.0556192, 0.0160602, -0.0351558, 0.257991, 0.0183949,`<br />
`    0.287146, 1.28508, 0.0155203, -0.260995, -0.0150086,`<br />
`    -0.0373431, 0.215909, -0.191682, -0.0007677, 0.00120234,`<br />
`    -0.00727776, 0.0155203, 1.28349, 0.0283306, 0.273234,`<br />
`    -0.0370315, -0.0671475, 0.0315086, -0.226642, 0.0351232,`<br />
`    -0.205039, -0.260995, 0.0283306, 1.21418, -0.0375033,`<br />
`    -0.00805812, 0.272086, -0.290368, 0.0608659, 0.036898,`<br />
`    -0.0587133, -0.0150086, 0.273234, -0.0375033, 1.18052),`<br />
`  nrow = 10,`<br />
`  ncol = 10,`<br />
`  byrow = TRUE`<br />
`)`<br />
`TestRound <- 2`<br />
`TestKS <- c(4, 4)`<br />
And then simmply call the function
<br />
`cpp_findmultPCs_deflation(TestMat, TestRound, TestKS)`
<a id="Files"></a>

## Development 
Here, we provide more information about the code structure and organization to help developers that would like to improve the method or build up on it. 

## Notes
- You can change the variable names.
- You can decide if you want to assign names to plain numbers or add comments for them.
- The help menu of the R package needs to change.

## Files
- R
  - RcppExports.R<br />
      It offers the R interface, which will call the corresponding C++ interface. Regenerate or change it manually if needed (e.g., if the interface changes). It can be generated using `Rcpp::compileAttributes()`.
- man/ contains the pages of the manual: one page for the package and one per function. 
- src/ contains the source files of the algorithm, in C++. 
  - ConstantArguments.h<br />
      It contains some parameters of the algorithm that are not directly tuneable by the end user.
  - msPCA_R_CPP.cpp<br />
      It contains the implementation of the algorithm.
  - RcppExports.cpp<br />
      It contains the converted function that can be used by R. Regenerate or change it manually if needed (e.g., if the interface changes). It can be generated using `Rcpp::compileAttributes()`.
  - Makevars<br />
      This is not currently used. Use it to set attributes, such as the version of C++ for compilation.
  - Makevars.win<br />
      This is not currently used. Use it to set attributes, such as the version of C++ for compilation.
- NAMESPACE<br />
    It is used to build this package. Change it if needed (e.g., if the interface changes).
- DESCRIPTION<br />
    It contains the description of this package.
- LICENSE<br />
    It contains the license information.
- msPCA.Rproj<br />
    It contains the settings of this R project. It is used by RStudio and often does not need to be changed.
## Guidance to future developers
- The essence of this algorithm is in the file "msPCA_R_CPP.cpp" and the file "ConstantArguments.h", where "msPCA_R_CPP.cpp" handles the computation and "ConstantArguments.h" lists all arguments. (This needs to change depending on your decisions on arguments.)
