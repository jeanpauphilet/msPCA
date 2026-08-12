############################################################
## Benchmarking notebook: msPCA vs competing packages on Breast data
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: PMA::breast_data (p = 500 genes, n = 89 samples)
## Task:    r = 3 sparse PCs with k = 10 nonzeros each
##
## Methods compared: msPCA, elasticnet, PMA, sparsepca,
##                   mixOmics, nsprcomp, nscumcomp, dense PCA
## For PMA (L1 penalty) and sparsepca (alpha penalty), the tuning
## parameter is set to achieve ~10 nonzeros per component.
## For mixOmics::spca(), keepX = rep(k, r) matches the cardinality
## budget exactly. For nsprcomp::nsprcomp(), k = rep(k, r)
## gives the per-component nonzero budget; nscumcomp() uses the
## total budget r * k.
##
## This is the canonical benchmark from Witten et al. (2009);
## including it allows direct comparison on equal footing with
## the original PMA paper.
############################################################

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("mixOmics")
library("nsprcomp")

## Load breast cancer data (89 samples x 19,672 genes).
## PMA::download_breast_data() downloads the RNA-seq expression data.
X <- t(PMA::download_breast_data(url = "https://tibshirani.su.domains/PMA/breastdata.rda")$rna)
vars  <- apply(X, 2, var)
X500  <- X[, order(vars, decreasing = TRUE)[1:500]]
S     <- cor(X500)                  # 500 x 500 correlation matrix
r <- 3L; k <- 20L

## --- msPCA ---
set.seed(43)
t_mspca <- system.time(
    mspca_res <- mspca(S, r = r, ks = rep(k, r), verbose = TRUE,
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

## --- PMA (sumabsv tuned for ~20 NNZ per PC on Breast data) ---
## Run the sweep first to find the value closest to k = 20 nonzeros per PC;
cat("PMA sumabsv sweep:\n")
for (sv in seq(2.9, 3, by = 0.01)) {
	pma_try <- PMA::SPC(X500, sumabsv = sv, K = r, orth = TRUE, trace = FALSE)
	cat("  sumabsv =", sv,
              "| NNZ:", paste(colSums(abs(pma_try$v) > 0), collapse = " "), 
              "| FVE:", round(fraction_variance_explained(S, pma_try$v), 4),
"\n")
}
## Set sumabsv to the value that achieves ~20 NNZ per PC above, then run:
t_pma <- system.time(
    pma_res <- PMA::SPC(X500, sumabsv = 2.92, K = r, orth = TRUE, trace = FALSE)
  )
cat("PMA | NNZ:", colSums(abs(pma_res$v) > 0),
       "| FVE:", round(fraction_variance_explained(S, pma_res$v), 4),
       "| Orth:", format(feasibility_violation_off(S, pma_res$v, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_pma["elapsed"], 3), "s\n")

## --- sparsepca (alpha tuned for ~20 NNZ per PC; normalize to unit norm) ---
for (a in seq(0.003, 0.005, by = 0.0005)) {
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

## --- mixOmics::spca (keepX = rep(k, r) matches cardinality budget exactly) ---
## mixOmics::spca() requires a raw data matrix; we pass the centred/scaled data.
t_mixo <- system.time(
    mixo_res <- mixOmics::spca(X500, ncomp = r,
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

## --- nsprcomp (k nonzeros per PC; basic deflation) ---
t_nspr <- system.time(
    nspr_res <- nsprcomp::nsprcomp(X500, ncomp = r, k = rep(k, r), nneg = FALSE, center = TRUE)
  )
nspr_load <- nspr_res$rotation
cat("nsprcomp | NNZ:", colSums(abs(nspr_load) > 0),
       "| FVE:", round(fraction_variance_explained(S, nspr_load), 4),
       "| Orth:", format(feasibility_violation_off(S, nspr_load, 0),
                         scientific = TRUE, digits = 3),
       "| Time:", round(t_nspr["elapsed"], 3), "s\n")

## --- nscumcomp (total nonzeros = r * k; joint optimization) ---
## nscumcomp() may fail on high-correlation data ("Co-linear principal axes").
## We use tryCatch() for graceful failure reporting.
t_nscum <- system.time(
    nscum_res <- nsprcomp::nscumcomp(X500, ncomp = r, k = r * k, nneg = FALSE, center = TRUE)
  )
nscum_load <- nscum_res$rotation
cat("nscumcomp | NNZ:", colSums(abs(nscum_load) > 0),
           "| FVE:", round(fraction_variance_explained(S, nscum_load), 4),
           "| Orth:", format(feasibility_violation_off(S, nscum_load, 0),
                             scientific = TRUE, digits = 3),
           "| Time:", round(t_nscum["elapsed"], 3), "s\n") 

## --- Dense PCA ---
t_pca <- system.time(
    pca_res <- prcomp(X500, scale. = TRUE)
  )
cat("Dense PCA | FVE:",
       round(fraction_variance_explained(S, pca_res$rotation[, 1:r]), 4),
       "| Time:", round(t_pca["elapsed"], 3), "s\n")

## --- Summary data frame ---
results <- data.frame(
    method  = c("msPCA", "elasticnet", "PMA", "sparsepca",
                "mixOmics", "nsprcomp", "nscumcomp", "Dense PCA"),
    nnz_pc1 = c(colSums(abs(mspca_res$x_best) > 0)[1],
                colSums(abs(enet_res$loadings)  > 0)[1],
                colSums(abs(pma_res$v)           > 0)[1],
                colSums(abs(spca_res$loadings)   > 0)[1],
                colSums(abs(mixo_load)           > 0)[1],
                colSums(abs(nspr_load)           > 0)[1],
                if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[1] else NA_integer_,
                500L),
    nnz_pc2 = c(colSums(abs(mspca_res$x_best) > 0)[2],
                colSums(abs(enet_res$loadings)  > 0)[2],
                colSums(abs(pma_res$v)           > 0)[2],
                colSums(abs(spca_res$loadings)   > 0)[2],
                colSums(abs(mixo_load)           > 0)[2],
                colSums(abs(nspr_load)           > 0)[2],
                if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[2] else NA_integer_,
                500L),
    nnz_pc3 = c(colSums(abs(mspca_res$x_best) > 0)[3],
                colSums(abs(enet_res$loadings)  > 0)[3],
                colSums(abs(pma_res$v)           > 0)[3],
                colSums(abs(spca_res$loadings)   > 0)[3],
                colSums(abs(mixo_load)           > 0)[3],
                colSums(abs(nspr_load)           > 0)[3],
                if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[3] else NA_integer_,
                500L),
    fve     = round(c(fraction_variance_explained(S, mspca_res$x_best),
                fraction_variance_explained(S, enet_res$loadings),
                fraction_variance_explained(S, pma_res$v),
                fraction_variance_explained(S, spca_res$loadings),
                fraction_variance_explained(S, mixo_load),
                fraction_variance_explained(S, nspr_load),
                if (!is.null(nscum_load)) fraction_variance_explained(S, nscum_load) else NA_real_,
                fraction_variance_explained(S, pca_res$rotation[, 1:r])), 3),
    orth    = c(feasibility_violation_off(S, mspca_res$x_best, 0),
                feasibility_violation_off(S, enet_res$loadings, 0),
                feasibility_violation_off(S, pma_res$v, 0),
                feasibility_violation_off(S, spca_res$loadings, 0),
                feasibility_violation_off(S, mixo_load, 0),
                feasibility_violation_off(S, nspr_load, 0),
                if (!is.null(nscum_load)) feasibility_violation_off(S, nscum_load, 0) else NA_real_,
                feasibility_violation_off(S, pca_res$rotation[, 1:r], 0)),
    runtime = round(c(t_mspca["elapsed"], t_enet["elapsed"],
                t_pma["elapsed"],   t_spca["elapsed"],
                t_mixo["elapsed"],  t_nspr["elapsed"],
                t_nscum["elapsed"], t_pca["elapsed"]), 3)
  )
write.csv(results, "benchmarking/benchmarking_results_breast.csv", row.names = FALSE)
