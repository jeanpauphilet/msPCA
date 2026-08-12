############################################################
## Benchmarking notebook: msPCA vs competing packages on Riboflavin
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: hdi::riboflavin (p = 4,088 genes, n = 71 bacterial cultures)
## Task:    r = 2 sparse PCs with k = 20 nonzeros each
##
## Methods compared: msPCA, elasticnet, PMA, sparsepca,
##                   mixOmics, nsprcomp, nscumcomp, dense PCA
## For PMA (L1 penalty) and sparsepca (alpha penalty), the tuning
## parameter is set to achieve ~20 nonzeros per component.
## For mixOmics::spca(), keepX = c(k, k) matches the cardinality
## budget exactly. For nsprcomp::nsprcomp(), k = c(k, k)
## gives the per-component nonzero budget; nscumcomp() uses the
## total budget r * k.
##
## NOTE: p >> n (4,088 >> 71). The correlation matrix S is rank-
## deficient (rank <= 70). Methods that operate on S directly
## (msPCA, elasticnet, sparsepca) may produce fewer meaningful
## components; results should be interpreted accordingly.
############################################################

library("msPCA")
library("elasticnet")
library("mixOmics")
library("nsprcomp")

## Load Riboflavin production data.
## hdi::riboflavin is a list with $x (71 x 4088 log-expression matrix)
## and $y (71-element log-riboflavin-production response vector).
data("riboflavin", package = "hdi")
X_ribo <- riboflavin$x               # n = 71, p = 4088
X_ribo <- scale(X_ribo)              # centre and scale

## Build correlation matrix (rank-deficient: rank <= n-1 = 70).
cat("Computing correlation matrix (p =", ncol(X_ribo), ")... ")
S <- cor(riboflavin$x)
cat("done.\n")

r <- 2L; k <- 20L

## --- msPCA ---
cat("Running msPCA...\n")
set.seed(43)
run_mspca <- system.time( 
  mspca_res <- mspca(S, r = r, ks = rep(k, r), verbose = TRUE,
          feasibilityConstraintType = 0,
          maxRestartTPM = 5, minRestartTPM=1, maxIter=60)
)
cat("msPCA | NNZ:", colSums(abs(mspca_res$x_best) > 0),
           "| FVE:", round(fraction_variance_explained(S, mspca_res$x_best), 4),
           "| Orth:", format(feasibility_violation_off(S, mspca_res$x_best, 0),
                             scientific = TRUE, digits = 3),
           "| Time:", round(run_mspca["elapsed"], 3), "s\n") 

## --- elasticnet ---
cat("Running elasticnet...\n")
run_enet <- system.time(
    enet_res <- elasticnet::spca(X_ribo, K = r, sparse = "varnum", trace=T,
                                 para = rep(k, r), type = "predictor")
  )
cat("elasticnet | NNZ:", colSums(abs(enet_res$loadings) > 0),
           "| FVE:", round(fraction_variance_explained(S, enet_res$loadings), 4),
           "| Orth:", format(feasibility_violation_off(S, enet_res$loadings, 0),
                             scientific = TRUE, digits = 3),
           "| Time:", round(run_enet["elapsed"], 3), "s\n")

## --- mixOmics::spca (keepX = c(k, k) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
cat("Running mixOmics...\n")
run_mixo <- system.time(
   mixo_res  <- mixOmics::spca(X_ribo, ncomp = r,
                   keepX = rep(k, r), center = TRUE)
    )
# mixo_res  <- run_mixo$result
mixo_load <- mixo_res$loadings$X
for (j in seq_len(r)) {
    nm <- sqrt(sum(mixo_load[, j]^2))
    if (nm > 0) mixo_load[, j] <- mixo_load[, j] / nm
  }
cat("mixOmics | NNZ:", colSums(abs(mixo_load) > 0),
        "| FVE:", round(fraction_variance_explained(S, mixo_load), 4),
        "| Orth:", format(feasibility_violation_off(S, mixo_load, 0),
                          scientific = TRUE, digits = 3),
        "| Time:", round(run_mixo["elapsed"], 3), "s\n")


## --- nsprcomp (k nonzeros per PC; basic deflation) ---
t_nspr <- system.time(
    nspr_res <- nsprcomp::nsprcomp(X_ribo, ncomp = r, k = rep(k, r), nneg = FALSE, center = TRUE)
  )
nspr_load <- nspr_res$rotation
cat("nsprcomp | NNZ:", colSums(abs(nspr_load) > 0),
       "| FVE:", round(fraction_variance_explained(S, nspr_load), 4),
  "| Orth:", format(feasibility_violation_off(S, nspr_load, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_nspr["elapsed"], 3), "s\n")

## --- Dense PCA ---
t_pca <- system.time(
    pca_res <- prcomp(riboflavin$x, scale. = TRUE)
  )
cat("Dense PCA | FVE:",
  round(fraction_variance_explained(S, pca_res$rotation[, 1:r]), 4),
  "| Time:", round(t_pca["elapsed"], 3), "s\n")


## --- Summary data frame ---
## Helper to extract NNZ safely (returns NA when method timed out / failed).
nnz <- function(mat, j) if (!is.null(mat)) colSums(abs(mat) > 0)[j] else NA_integer_
fve_safe  <- function(mat) if (!is.null(mat)) fraction_variance_explained(S, mat) else NA_real_
orth_safe <- function(mat) if (!is.null(mat)) feasibility_violation_off(S, mat, 0) else NA_real_

results <- data.frame(
    method  = c("msPCA", "elasticnet", "mixOmics", "nsprcomp", "Dense PCA"),
    nnz_pc1 = c(nnz(mspca_res$x_best, 1L),
                nnz(enet_res$loadings,  1L),
                nnz(mixo_load, 1L),
                nnz(nspr_load, 1L),
                ncol(S)),
    nnz_pc2 = c(nnz(mspca_res$x_best, 2L),
                nnz(enet_res$loadings,  2L),
                nnz(mixo_load, 2L),
                nnz(nspr_load, 2L),
                ncol(S)),
    fve     = round(c(fve_safe(mspca_res$x_best),
                fve_safe(enet_res$loadings),
                fve_safe(mixo_load),
                fve_safe(nspr_load),
                fraction_variance_explained(S, pca_res$rotation[, 1:r])), 3),
    orth    = c(orth_safe(mspca_res$x_best),
                orth_safe(enet_res$loadings),
                orth_safe(mixo_load),
                orth_safe(nspr_load),
                orth_safe(pca_res$rotation[, 1:r])),
    runtime = round(c(run_mspca["elapsed"], run_enet["elapsed"],
          run_mixo["elapsed"], t_nspr["elapsed"], t_pca["elapsed"]), 3)
  )
write.csv(results, "benchmarking/benchmarking_results_riboflavin.csv", row.names = FALSE)
