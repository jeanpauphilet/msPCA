############################################################
## Benchmarking notebook: msPCA vs competing packages on mtcars
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: datasets::mtcars (p = 11 variables, n = 32 observations)
## Task:    r = 2 sparse PCs with k = 4 nonzeros each
##
## Methods compared: msPCA, elasticnet, PMA, sparsepca,
##                   mixOmics, nsprcomp, dense PCA
## For PMA (L1 penalty) and sparsepca (alpha penalty), the tuning
## parameter is set to achieve ~4 nonzeros per component.
## For mixOmics::spca(), keepX = c(k, k) matches the cardinality
## budget exactly. For nsprcomp::nscumcomp(), k = r * k_per_pc
## gives the total nonzero budget across all components.
############################################################

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("mixOmics")
library("nsprcomp")
library("datasets")

S <- cor(datasets::mtcars)
r <- 3L; k <- 4L

## --- msPCA ---
set.seed(43)
t_mspca <- system.time(
    mspca_res <- mspca(S, r = r, ks = rep(k, r), verbose = FALSE,
                       feasibilityConstraintType = 0)
  )
cat("msPCA | NNZ:", colSums(abs(mspca_res$x_best) > 0),
       "| FVE:", round(fraction_variance_explained(S, mspca_res$x_best), 4),
       "| Orth:", format(feasibility_violation_off(S, mspca_res$x_best, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_mspca["elapsed"], 3), "s\n")

## --- elasticnet ---
t_enet <- system.time(
    enet_res <- elasticnet::spca(S, K = r, sparse = "varnum",
                                 para = rep(k, r), type = "Gram")
  )
cat("elasticnet | NNZ:", colSums(abs(enet_res$loadings) > 0),
       "| FVE:", round(fraction_variance_explained(S, enet_res$loadings), 4),
       "| Orth:", format(feasibility_violation_off(S, enet_res$loadings, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_enet["elapsed"], 3), "s\n")

## --- PMA (sumabsv tuned for ~3 NNZ per PC) ---
for (sv in seq(1.6, 1.8, by = 0.01)) {
	pma_try <- PMA::SPC(S, sumabsv = sv, K = r, orth = TRUE, trace = FALSE)
	cat("  sumabsv =", sv,
              "| NNZ:", paste(colSums(abs(pma_try$v) > 0), collapse = " "), 
              "| FVE:", round(fraction_variance_explained(S, pma_try$v), 4),
"\n")
}
t_pma <- system.time(
    pma_res <- PMA::SPC(S, sumabsv = 1.69, K = r, orth = TRUE, trace = FALSE)
  )
cat("PMA | NNZ:", colSums(abs(pma_res$v) > 0),
       "| FVE:", round(fraction_variance_explained(S, pma_res$v), 4),
       "| Orth:", format(feasibility_violation_off(S, pma_res$v, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_pma["elapsed"], 3), "s\n")

## --- sparsepca (alpha tuned for ~3 NNZ per PC; normalize to unit norm) ---
for (a in seq(0.001, 0.01, by = 0.001)) {
	sparsepca_try <- sparsepca::spca(S, k = r, alpha = a, verbose = FALSE, scale = FALSE)
       for (j in seq_len(r)) {
              nm <- sqrt(sum(sparsepca_try$loadings[, j]^2))
              if (nm > 0) sparsepca_try$loadings[, j] <- sparsepca_try$loadings[, j] / nm
       }       
	cat("  alpha =", a,
              "| NNZ:", paste(colSums(abs(sparsepca_try$loadings) > 0), collapse = " "), 
              "| FVE:", round(fraction_variance_explained(S, sparsepca_try$loadings), 4),
"\n")
}
t_spca <- system.time(
    spca_res <- sparsepca::spca(S, k = r, alpha = 0.004, verbose = FALSE, scale = FALSE)
  )
for (j in seq_len(r)) {
    nm <- sqrt(sum(spca_res$loadings[, j]^2))
    if (nm > 0) spca_res$loadings[, j] <- spca_res$loadings[, j] / nm
  }
cat("sparsepca | NNZ:", colSums(abs(spca_res$loadings) > 0),
       "| FVE:", round(fraction_variance_explained(S, spca_res$loadings), 4),
       "| Orth:", format(feasibility_violation_off(S, spca_res$loadings, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_spca["elapsed"], 3), "s\n")


## --- mixOmics::spca (keepX = c(k, k) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
X_mtcars <- scale(datasets::mtcars)
t_mixo <- system.time(
    mixo_res <- mixOmics::spca(X_mtcars, ncomp = r,
                               keepX = rep(k, r), center = TRUE, scale = TRUE)
  )
mixo_load <- mixo_res$loadings$X
for (j in seq_len(r)) {
    nm <- sqrt(sum(mixo_load[, j]^2))
    if (nm > 0) mixo_load[, j] <- mixo_load[, j] / nm
  }
cat("mixOmics | NNZ:", colSums(abs(mixo_load) > 0),
       "| FVE:", round(fraction_variance_explained(S, mixo_load), 4),
       "| Orth:", format(feasibility_violation_off(S, mixo_load, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_mixo["elapsed"], 3), "s\n")

## --- nsprcomp (k nonzeros per PC) ---
# Basic deflation
t_nscomp <- system.time(
              nscomp_load <- nsprcomp::nsprcomp(X_mtcars, ncomp = r, k = rep(k, r), nneg = FALSE, center = TRUE, scale. = TRUE)
       )
nscomp_load <- nscomp_load$rotation
cat("nsprcomp | NNZ:", colSums(abs(nscomp_load) > 0),
       "| FVE:", round(fraction_variance_explained(S, nscomp_load), 4),
       "| Orth:", format(feasibility_violation_off(S, nscomp_load, 1),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_nscomp["elapsed"], 3), "s\n")

# Joint optimization
t_nscum <- system.time(
              nscum_res <- nsprcomp::nscumcomp(X_mtcars, ncomp = r, k = k*r, nneg = FALSE,
                                              center = TRUE, scale. = TRUE)
       )
nscum_load <- nscum_res$rotation
cat("nscumcomp | NNZ:", colSums(abs(nscum_load) > 0),
       "| FVE:", round(fraction_variance_explained(S, nscum_load), 4),
       "| Orth:", format(feasibility_violation_off(S, nscum_load, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_nscum["elapsed"], 3), "s\n")

## --- Dense PCA ---
t_pca <- system.time(
       pca_res <- prcomp(datasets::mtcars, scale. = TRUE)
)
cat("Dense PCA | FVE:",
       round(fraction_variance_explained(S, pca_res$rotation[, 1:2]), 4), "\n")

## --- Summary data frame ---
results <- data.frame(
    method  = c("msPCA", "elasticnet", "PMA", "sparsepca",
                "mixOmics", "nsprcomp", "nscumcomp", "Dense PCA"),
    nnz_pc1 = c(colSums(abs(mspca_res$x_best) > 0)[1],
                colSums(abs(enet_res$loadings)  > 0)[1],
                colSums(abs(pma_res$v)           > 0)[1],
                colSums(abs(spca_res$loadings)   > 0)[1],
                colSums(abs(mixo_load)           > 0)[1],
                colSums(abs(nscomp_load)          > 0)[1],
                colSums(abs(nscum_load)          > 0)[1],
                11L),
    nnz_pc2 = c(colSums(abs(mspca_res$x_best) > 0)[2],
                colSums(abs(enet_res$loadings)  > 0)[2],
                colSums(abs(pma_res$v)           > 0)[2],
                colSums(abs(spca_res$loadings)   > 0)[2],
                colSums(abs(mixo_load)           > 0)[2],
                colSums(abs(nscomp_load)          > 0)[2],
                colSums(abs(nscum_load)          > 0)[2],
                11L),
    nnz_pc3 = c(colSums(abs(mspca_res$x_best) > 0)[3],
                colSums(abs(enet_res$loadings)  > 0)[3],
                colSums(abs(pma_res$v)           > 0)[3],
                colSums(abs(spca_res$loadings)   > 0)[3],
                colSums(abs(mixo_load)           > 0)[3],
                colSums(abs(nscomp_load)          > 0)[3],
                colSums(abs(nscum_load)          > 0)[3],
                11L),
    fve     = round(c(fraction_variance_explained(S, mspca_res$x_best),
                fraction_variance_explained(S, enet_res$loadings),
                fraction_variance_explained(S, pma_res$v),
                fraction_variance_explained(S, spca_res$loadings),
                fraction_variance_explained(S, mixo_load),
                fraction_variance_explained(S, nscomp_load),
                fraction_variance_explained(S, nscum_load),
                fraction_variance_explained(S, pca_res$rotation[, 1:r])
                ),3),
    orth    = c(feasibility_violation_off(S, mspca_res$x_best, 0),
                feasibility_violation_off(S, enet_res$loadings, 0),
                feasibility_violation_off(S, pma_res$v, 0),
                feasibility_violation_off(S, spca_res$loadings, 0),
                feasibility_violation_off(S, mixo_load, 0),
                feasibility_violation_off(S, nscomp_load, 0),
                feasibility_violation_off(S, nscum_load, 0),
                feasibility_violation_off(S, pca_res$rotation[, 1:r], 0)
              ),
    runtime = round(c(t_mspca["elapsed"], t_enet["elapsed"],
                t_pma["elapsed"],   t_spca["elapsed"],
                t_mixo["elapsed"],  t_nscomp["elapsed"], t_nscum["elapsed"], t_pca["elapsed"]), 3)
)

write.csv(results, "benchmarking/benchmarking_results_mtcars.csv", row.names = FALSE)


