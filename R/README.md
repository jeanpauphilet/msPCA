# R Package
## Tools
IDE: RStudio
## Implementation
- This code is based on the C++ code in sPCAmPC.Julia_and_CPP, whose details can be found in its readme file.
- As said in the readme file in sPCAmPC.Julia_and_CPP, I do not regard this C++ code as perfect, or even passable. However, before further discussion, perhaps it is good enough from some point of view.
  Also, the focus is to make the R package based on the given C++ code, not to judge it or make it better beyond the correctness. Therefore, I used the C++ code in sPCAmPC.Julia_and_CPP, instead of the one in sPCAmPC.CPP.Better.
- Some code changes need to be made to make it suitable for the R package, though fundamentally they are the same.
## R Interface
The interface in R is the same as the one used in C++ code in sPCAmPC.Julia_and_CPP.<br />
Though better or more interfaces can be made, the number of arguments is so big. Before further discussion, the correct decision seemed to be just using the old interface.
## Test Template
`TestMat <- matrix(`<br />
`  c(1.0059, -0.0089081, 0.0418695, -0.000217439, 0.0453617, 0.00928114, 0.0536835, -0.0121021, -0.0217993, -0.0168165,`<br />
`    -0.0089081, 1.21568, 0.0421938, -0.190951, 0.0292236, 0.168849, 0.278629, -0.0212905, 0.00591782, -0.0409329,`<br />
`    0.0418695, 0.0421938, 1.01591, 0.0305831, 0.0125232, 0.0125612, -0.00797738, 0.0221298, -0.0549147, -0.00248777,`<br />
`    -0.000217439, -0.190951, 0.0305831, 1.23649, 0.0541234, -0.228092, -0.184682, -0.0682807, -0.0321624, 0.00907545,`<br />
`    0.0453617, 0.0292236, 0.0125232, 0.0541234, 1.17742, 0.0524149, -0.016881, -0.283933, 0.240862, 0.247478,`<br />
`    0.00928114, 0.168849, 0.0125612, -0.228092, 0.0524149, 1.21963, 0.244194, -0.0221961, -0.0391619, 0.0413205,`<br />
`    0.0536835, 0.278629, -0.00797738, -0.184682, -0.016881, 0.244194, 1.31736, 0.0464585, 0.000602352, 0.013095,`<br />
`    -0.0121021, -0.0212905, 0.0221298, -0.0682807, -0.283933, -0.0221961, 0.0464585, 1.28139, -0.284837, -0.346121,`<br />
`    -0.0217993, 0.00591782, -0.0549147, -0.0321624, 0.240862, -0.0391619, 0.000602352, -0.284837, 1.2565, 0.219123,`<br />
`    -0.0168165, -0.0409329, -0.00248777, 0.00907545, 0.247478, 0.0413205, 0.013095, -0.346121, 0.219123, 1.23185),`<br />
`  nrow = 10,`<br />
`  ncol = 10,`<br />
`  byrow = TRUE`<br />
`)`<br />
`TestRound <- 2`<br />
`TestKS <- matrix(`<br />
`  c(4, 4),`<br />
`  nrow = 1,`<br />
`  ncol = 2,`<br />
`  byrow = TRUE`<br />
`)`<br />
`TestNumIters <- 20`<br />
<br />
`cpp_findmultPCs_deflation(TestMat, TestRound, TestKS, TestNumIters)`
