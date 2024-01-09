# sPCAmPC
An R Package for Sparse PCA with Multiple Principal Components
## Notes
- You can change the variable names.
- You can decide if you want to assign names to plain numbers or add comments for them.
- The help menu of the R package needs to change.
## Tools
IDE: RStudio<br />
(It can be used without RStudio, though I have not yet tried this.) 
## Usage
If you use RStudio, you can use this library by typing the following into the console,<br />
`library(devtools)`<br />
I guess devtools needs to be installed if it has never been installed. To install, use `install.packages("devtools")`.<br />
`install_github('jeanpauphilet/sPCAmPC/R')`<br />
wait until the installation is finished<br />
`library(sPCAmPC)`<br />
Then, input the variables and call the function.<br />
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
