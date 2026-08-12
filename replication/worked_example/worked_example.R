## ============================================================
## Section 4 Worked Example: msPCA on the mtcars dataset
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: datasets::mtcars (p = 11 variables, n = 32 observations)
## Goal: extract r = 2 sparse PCs with k = 4 nonzeros each
## ============================================================

# install.packages("msPCA")
library("msPCA")

## Compute the correlation matrix
Sigma <- cor(datasets::mtcars)

## Fit two 4-sparse PCs with the orthogonality constraint
set.seed(42)
res <- mspca(Sigma, r = 2, ks = c(4, 4), verbose = FALSE)
print_mspca(res, Sigma)

## Inspect feasibility and fraction of variance explained
feasibility_violation_off(Sigma, res$x_best, 0)
fraction_variance_explained(Sigma, res$x_best)

## Compare with dense PCA (prcomp)
pca_res <- prcomp(datasets::mtcars, scale. = TRUE)
fraction_variance_explained(Sigma, pca_res$rotation[, 1:2])

