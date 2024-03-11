source("utils.R")

############################################################
## In this notebook, we illustrate the use of the msPCA package on the mtcars dataset
## Our goal is to find 2 4-sparse PCs that explain most of the variance.
##
## We compare the methods from the msPCA package with other packages.
##
## For packages that cannot constraint the sparsity of each PC directly (sparsepca, PMA),
## we tuned the penalty on the sparsity/L1 norm in order to achieve a sparsity pattern
## as close as possible (sparsity of 4 for each PC).
############################################################
library(datasets)
df <- datasets::mtcars
S <- cor(df)

############################################################
# msPCA package
############################################################
library(msPCA)

## First method: Truncated Power Method for a single sparse PC
tpw_results <- msPCA::tpw(S, 4, maxIter=100)
U <- as.matrix(tpw_results$x_best)
#Sparsity
colSums(abs(U) > 0)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,U)


## Second method: Iterative Deflation Heuristic for multiple sparse PCs
mspca_results <- msPCA::mspca(S, 2, c(4,4), verbose=TRUE)
#Sparsity
colSums(abs(mspca_results$x_best) > 0)
#Orthogonality
msPCA:::orthogonality_violation(mspca_results$x_best)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,mspca_results$x_best)

#' Compute the fraction of variance explained by each PC
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return An array.
variance_explained_perPC <- function(C, U){
  r = dim(U)[2]
  ve = numeric(r)
  for (i in 1:r){
    ve[i] <- sum(U[,i] %*% C %*% U[,i])
  }
  ve
}

#' Compute the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix) by each PC
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return An array.
fraction_variance_explained_perPC <- function(C, U){
  fve <- variance_explained_perPC(C,U)
  fve <- fve / sum(diag(C))
  fve
}

#' Compute the fraction of variance explained (variance explained normalized by the trace of the covariance/correlation matrix)
#' @param C A matrix. The correlation or covariance matrix (p x p).
#' @param U A matrix. The matrix containing the r PCs (p x r).
#' @return A float.
fraction_variance_explained <- function(C, U){
  sum(fraction_variance_explained_perPC(C,U))
}

print.mspca<-function(sol_object, C){
  cat("\n msPCA solution:\n")
  r <- dim(sol_object$x_best)[2]
  cat("\n")
  cat(paste(r,"sparse PCs",sep=" "), "\n")

  fve <- fraction_variance_explained_perPC(C, sol_object$x_best)
  cat("Pct. of exp. var. :", format(round(fve, 3)*100), "\n")

  v <- sol_object$x_best
  k <- 1:r
  for (j in 1:r){
    k[j] <- sum(v[,j] != 0)
  }
  cat("Num. of non-zero loadings : ", k, "\n")
  cat("Sparse loadings \n")
  print(round(v,3))

}


variance_explained_perPC(S,mspca_results$x_best)
fraction_variance_explained_perPC(S,mspca_results$x_best)
fraction_variance_explained(S,mspca_results$x_best)

############################################################
# Zou et al. (2008) (elasticnet)
############################################################
library(elasticnet)
enet_results <- elasticnet::spca(S, 2, sparse="varnum", para = c(4,4), type="Gram")
#Sparsity
colSums(abs(enet_results$loadings) > 0)
#Orthogonality
msPCA:::orthogonality_violation(enet_results$loadings)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,enet_results$loadings)



############################################################
# Witten et al. (2009) (PMA)
############################################################
library(PMA)
pma_results <- PMA::SPC(as.matrix(df), sumabsv = sqrt(1.8291), niter = 20, K = 2, orth = FALSE, trace = TRUE)
#Sparsity
colSums(abs(pma_results$v) > 0)
#Orthogonality
msPCA:::orthogonality_violation(pma_results$v)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,pma_results$v)

pma_results <- PMA::SPC(as.matrix(df), sumabsv = sqrt(1.83), niter = 20, K = 2, orth = TRUE, trace = TRUE)
#Sparsity
colSums(abs(pma_results$v) > 0)
#Orthogonality
msPCA:::orthogonality_violation(pma_results$v)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,pma_results$v)



############################################################
###Erichson et al. (2018) (sparsepca)
############################################################
library(sparsepca)

sparsepca <- sparsepca::spca(S, k=2, alpha=0.01, verbose=F)
sparsepca$loadings[,1] <- sparsepca$loadings[,1] / sum(sparsepca$loadings[,1]**2)
sparsepca$loadings[,2] <- sparsepca$loadings[,2] / sum(sparsepca$loadings[,2]**2)
#Sparsity
colSums(abs(sparsepca$loadings) > 0)
#Orthogonality
msPCA:::orthogonality_violation(sparsepca$loadings)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,sparsepca$loadings)


#randomized PCA -> rpca
sparsepca <- sparsepca::rspca(S, k=2, alpha=0.01, verbose=F)
sparsepca$loadings[,1] <- sparsepca$loadings[,1] / sum(sparsepca$loadings[,1]**2)
sparsepca$loadings[,2] <- sparsepca$loadings[,2] / sum(sparsepca$loadings[,2]**2)
#Sparsity
colSums(abs(sparsepca$loadings) > 0)
#Orthogonality
msPCA:::orthogonality_violation(sparsepca$loadings)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,sparsepca$loadings)


#robust PCA -> robpca
sparsepca <- sparsepca::robspca(S, k=2, alpha=0.01, verbose=F)
sparsepca$loadings[,1] <- sparsepca$loadings[,1] / sum(sparsepca$loadings[,1]**2)
sparsepca$loadings[,2] <- sparsepca$loadings[,2] / sum(sparsepca$loadings[,2]**2)
#Sparsity
colSums(abs(sparsepca$loadings) > 0)
#Orthogonality
msPCA:::orthogonality_violation(sparsepca$loadings)
#Fraction of variance explained
msPCA:::fraction_variance_explained(S,sparsepca$loadings)


