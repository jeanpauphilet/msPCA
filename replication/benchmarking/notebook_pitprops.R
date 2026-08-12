## ============================================================
## Pitprops benchmark: msPCA vs competing packages
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: Pitprops correlation matrix (Jeffers 1967, p = 13)
##   13 physical measurements on 180 timber pit-prop specimens.
##   Classic sparse PCA benchmark; reference results in
##   Zou, Hastie & Tibshirani (2006, Table 3).
##
## Task: r = 6 sparse PCs with k = 4 nonzeros each
##
## Note: mixOmics::spca() and nsprcomp::nscumcomp() require a raw
## data matrix. Since only the correlation matrix is available for
## this dataset, we generate pseudo-data via MASS::mvrnorm() with
## the Pitprops matrix as the population covariance. FVE and
## orthogonality are then evaluated on the original pitprops matrix.
## ============================================================

library("msPCA")
library("elasticnet")
library("PMA")
library("sparsepca")
library("mixOmics")
library("nsprcomp")
library("MASS")

## Pitprops correlation matrix from Jeffers (1967), Table 1.
## Variables: topdiam, length, moist, testsg, ovensg, ringtop,
##            ringbut, bowmax, bowdist, whorls, clear, knots, diaknot
pitprops <- matrix(c(
	1.000,  0.954,  0.364,  0.342, -0.129,  0.313,  0.496,  0.424,
	0.592,  0.545,  0.084, -0.019,  0.134,
	0.954,  1.000,  0.297,  0.284, -0.118,  0.291,  0.503,  0.419,
	0.648,  0.569,  0.076, -0.036,  0.144,
	0.364,  0.297,  1.000,  0.882, -0.148,  0.153, -0.029, -0.054,
	0.125, -0.081,  0.162,  0.220,  0.126,
	0.342,  0.284,  0.882,  1.000,  0.220,  0.381,  0.174, -0.059,
	0.137, -0.014,  0.097,  0.169,  0.015,
 -0.129, -0.118, -0.148,  0.220,  1.000,  0.364,  0.296,  0.004,
 -0.039,  0.037, -0.091, -0.145, -0.208,
	0.313,  0.291,  0.153,  0.381,  0.364,  1.000,  0.813,  0.090,
	0.211,  0.274, -0.036,  0.024, -0.329,
	0.496,  0.503, -0.029,  0.174,  0.296,  0.813,  1.000,  0.372,
	0.465,  0.679, -0.113, -0.232, -0.424,
	0.424,  0.419, -0.054, -0.059,  0.004,  0.090,  0.372,  1.000,
	0.482,  0.557,  0.061, -0.357, -0.172,
	0.592,  0.648,  0.125,  0.137, -0.039,  0.211,  0.465,  0.482,
	1.000,  0.526,  0.085, -0.127, -0.175,
	0.545,  0.569, -0.081, -0.014,  0.037,  0.274,  0.679,  0.557,
	0.526,  1.000, -0.319, -0.368, -0.545,
	0.084,  0.076,  0.162,  0.097, -0.091, -0.036, -0.113,  0.061,
	0.085, -0.319,  1.000,  0.029,  0.054,
 -0.019, -0.036,  0.220,  0.169, -0.145,  0.024, -0.232, -0.357,
 -0.127, -0.368,  0.029,  1.000,  0.220,
	0.134,  0.144,  0.126,  0.015, -0.208, -0.329, -0.424, -0.172,
 -0.175, -0.545,  0.054,  0.220,  1.000
), nrow = 13, byrow = TRUE)
varnames <- c("topdiam", "length", "moist", "testsg", "ovensg",
							"ringtop", "ringbut", "bowmax", "bowdist", "whorls",
							"clear", "knots", "diaknot")
rownames(pitprops) <- colnames(pitprops) <- varnames
p <- 13L; r <- 6L; k <- 4L

## Pseudo-data for methods requiring a raw data matrix (mixOmics, nsprcomp).
## Generated from the pitprops correlation matrix; n = 500 gives stable loadings.
set.seed(42)
X_pseudo <- MASS::mvrnorm(n = 500L, mu = rep(0, p), Sigma = pitprops)
colnames(X_pseudo) <- varnames

## --- msPCA ---
set.seed(43)
t_mspca <- system.time(
	mspca_res <- mspca(pitprops, r = r, ks = rep(k, r), 
										 verbose = FALSE, feasibilityConstraintType = 0)
)
cat("msPCA | NNZ:", colSums(abs(mspca_res$x_best) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, mspca_res$x_best), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, mspca_res$x_best, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_mspca["elapsed"], 4), "s\n")

## --- elasticnet ---
t_enet <- system.time(
	enet_res <- elasticnet::spca(pitprops, K = r, sparse = "varnum",
															 para = rep(k, r), type = "Gram")
)
cat("elasticnet | NNZ:", colSums(abs(enet_res$loadings) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, enet_res$loadings), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, enet_res$loadings, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_enet["elapsed"], 4), "s\n")

## --- PMA (sumabsv tuned for ~4 NNZ per	 PC on Pitprops) ---
## Sweep to find appropriate sumabsv
for (sv in seq(1.7, 1.8, by = 0.001)) {
	pma_try <- PMA::SPC(pitprops, sumabsv = sv, K = r, orth = TRUE, trace = FALSE)
	cat("  sumabsv =", sv,
              "| NNZ:", paste(colSums(abs(pma_try$v) > 0), collapse = " "), 
              "| FVE:", round(fraction_variance_explained(pitprops, pma_try$v), 4),
"\n")
}
## Set sumabsv to the value that achieves ~4 NNZ per PC above, then run:
t_pma <- system.time(
	pma_res <- PMA::SPC(pitprops, sumabsv = 1.727, K = r,
											orth = TRUE, trace = FALSE)
)
cat("PMA | NNZ:", colSums(abs(pma_res$v) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, pma_res$v), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, pma_res$v, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_pma["elapsed"], 4), "s\n")

## --- sparsepca ---
## --- sparsepca (alpha = 0.004 tuned for ~4 NNZ per PC; normalize to unit norm) ---
for (a in seq(0.001, 0.005, by = 0.0005)) {
	sparsepca_try <- sparsepca::spca(pitprops, k = r, alpha = a, verbose = FALSE, scale = TRUE)
       for (j in seq_len(r)) {
              nm <- sqrt(sum(sparsepca_try$loadings[, j]^2))
              if (nm > 0) sparsepca_try$loadings[, j] <- sparsepca_try$loadings[, j] / nm
       }       
	cat("  alpha =", a,
              "| NNZ:", paste(colSums(abs(sparsepca_try$loadings) > 0), collapse = " "), 
              "| FVE:", round(fraction_variance_explained(pitprops, sparsepca_try$loadings), 4),
"\n")
}
t_spca <- system.time(
	spca_res <- sparsepca::spca(pitprops, k = r, alpha = 0.004, scale = TRUE, verbose = FALSE)
)
for (j in seq_len(r)) {
	nm <- sqrt(sum(spca_res$loadings[, j]^2))
	if (nm > 0) spca_res$loadings[, j] <- spca_res$loadings[, j] / nm
}
cat("sparsepca | NNZ:", colSums(abs(spca_res$loadings) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, spca_res$loadings), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, spca_res$loadings, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_spca["elapsed"], 4), "s\n")

## --- mixOmics::spca (keepX tuned to k per component; uses pseudo-data) ---
t_mixo <- system.time(
	mixo_res <- mixOmics::spca(X_pseudo, ncomp = r, keepX = rep(k, r), center = TRUE)
)
mixo_load <- mixo_res$loadings$X
for (j in seq_len(r)) {
	nm <- sqrt(sum(mixo_load[, j]^2))
	if (nm > 0) mixo_load[, j] <- mixo_load[, j] / nm
}
cat("mixOmics | NNZ:", colSums(abs(mixo_load) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, mixo_load), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, mixo_load, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_mixo["elapsed"], 4), "s\n")

## --- nsprcomp (k nonzeros per PC; uses pseudo-data) ---
## Basic deflation
t_nspr <- system.time(
	nspr_res <- nsprcomp::nsprcomp(X_pseudo, ncomp = r, k = rep(k, r), nneg = FALSE, center = TRUE)
)
nspr_load <- nspr_res$rotation
cat("nsprcomp | NNZ:", colSums(abs(nspr_load) > 0),
		"| FVE:", round(fraction_variance_explained(pitprops, nspr_load), 4),
		"| Orth:", format(feasibility_violation_off(pitprops, nspr_load, 0),
											scientific = TRUE, digits = 3),
		"| Time:", round(t_nspr["elapsed"], 4), "s\n")

## --- nscumcomp (total nonzeros = r * k; uses pseudo-data) ---
## Note: nscumcomp() may fail on datasets with high inter-variable correlations
## ("Co-linear principal axes"). We use tryCatch() for graceful failure reporting.
max_attempts <- 20L
nscum_res <- NULL
last_err <- NULL

t_nscum <- system.time({
  for (attempt in seq_len(max_attempts)) {
    fit <- tryCatch(
      nsprcomp::nscumcomp(
        X_pseudo, ncomp = r, k = r * k, nneg = FALSE,
        center = TRUE
      ),
      error = function(e) e
    )

    if (!inherits(fit, "error")) {
      nscum_res <- fit
      break
    }

    last_err <- fit$message
  }
})

nscum_load <- if (!is.null(nscum_res)) nscum_res$rotation else NULL

if (is.null(nscum_load)) {
  cat("nscumcomp failed after", max_attempts, "attempts. Last error:", last_err, "\n")
} else {
  cat("nscumcomp | NNZ:", colSums(abs(nscum_load) > 0),
      "| FVE:", round(fraction_variance_explained(pitprops, nscum_load), 4),
      "| Orth:", format(feasibility_violation_off(pitprops, nscum_load, 0),
                        scientific = TRUE, digits = 3),
      "| Time:", round(t_nscum["elapsed"], 4), "s\n")
}

## --- Dense PCA ---
t_pca <- system.time(
	pca_res <- prcomp(t(pitprops), scale. = FALSE)  # treat as cov matrix
)
cat("Dense PCA (r=6) | FVE:",
		round(fraction_variance_explained(pitprops, pca_res$rotation[, 1:6]), 4),
		"| Time:", round(t_pca["elapsed"], 4), "s\n")

## --- Save results ---
results_pitprops <- data.frame(
	method  = c("msPCA", "elasticnet", "PMA", "sparsepca",
							"mixOmics", "nsprcomp", "nscumcomp", "Dense PCA"),
	nnz_pc1 = c(colSums(abs(mspca_res$x_best) > 0)[1],
							colSums(abs(enet_res$loadings)  > 0)[1],
							colSums(abs(pma_res$v)           > 0)[1],
							colSums(abs(spca_res$loadings)   > 0)[1],
							colSums(abs(mixo_load)           > 0)[1],
							colSums(abs(nspr_load)           > 0)[1],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[1] else NA,
							13L),
	nnz_pc2 = c(colSums(abs(mspca_res$x_best) > 0)[2],
							colSums(abs(enet_res$loadings)  > 0)[2],
							colSums(abs(pma_res$v)           > 0)[2],
							colSums(abs(spca_res$loadings)   > 0)[2],
							colSums(abs(mixo_load)           > 0)[2],
							colSums(abs(nspr_load)           > 0)[2],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[2] else NA,
							13L),
	nnz_pc3 = c(colSums(abs(mspca_res$x_best) > 0)[3],
							colSums(abs(enet_res$loadings)  > 0)[3],
							colSums(abs(pma_res$v)           > 0)[3],
							colSums(abs(spca_res$loadings)   > 0)[3],
							colSums(abs(mixo_load)           > 0)[3],
							colSums(abs(nspr_load)           > 0)[3],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[3] else NA,
							13L),
	nnz_pc4 = c(colSums(abs(mspca_res$x_best) > 0)[4],
							colSums(abs(enet_res$loadings)  > 0)[4],
							colSums(abs(pma_res$v)           > 0)[4],
							colSums(abs(spca_res$loadings)   > 0)[4],
							colSums(abs(mixo_load)           > 0)[4],
							colSums(abs(nspr_load)           > 0)[4],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[4] else NA,
							13L),
	nnz_pc5 = c(colSums(abs(mspca_res$x_best) > 0)[5],
							colSums(abs(enet_res$loadings)  > 0)[5],
							colSums(abs(pma_res$v)           > 0)[5],
							colSums(abs(spca_res$loadings)   > 0)[5],
							colSums(abs(mixo_load)           > 0)[5],
							colSums(abs(nspr_load)           > 0)[5],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[5] else NA,
							13L),
	nnz_pc6 = c(colSums(abs(mspca_res$x_best) > 0)[6],
							colSums(abs(enet_res$loadings)  > 0)[6],
							colSums(abs(pma_res$v)           > 0)[6],
							colSums(abs(spca_res$loadings)   > 0)[6],
							colSums(abs(mixo_load)           > 0)[6],
							colSums(abs(nspr_load)           > 0)[6],
							if (!is.null(nscum_load)) colSums(abs(nscum_load) > 0)[6] else NA,
							13L),
	fve     = round(c(fraction_variance_explained(pitprops, mspca_res$x_best),
							fraction_variance_explained(pitprops, enet_res$loadings),
							fraction_variance_explained(pitprops, pma_res$v),
							fraction_variance_explained(pitprops, spca_res$loadings),
							fraction_variance_explained(pitprops, mixo_load),
							fraction_variance_explained(pitprops, nspr_load),
							if (!is.null(nscum_load)) fraction_variance_explained(pitprops, nscum_load) else NA,
							fraction_variance_explained(pitprops, pca_res$rotation[, 1:6])), 3),
	orth    = c(feasibility_violation_off(pitprops, mspca_res$x_best, 0),
							feasibility_violation_off(pitprops, enet_res$loadings, 0),
							feasibility_violation_off(pitprops, pma_res$v, 0),
							feasibility_violation_off(pitprops, spca_res$loadings, 0),
							feasibility_violation_off(pitprops, mixo_load, 0),
							feasibility_violation_off(pitprops, nspr_load, 0),
							if (!is.null(nscum_load)) feasibility_violation_off(pitprops, nscum_load, 0) else NA,
							feasibility_violation_off(pitprops, pca_res$rotation[, 1:6], 0)),
	runtime = round(c(t_mspca["elapsed"], t_enet["elapsed"],
							t_pma["elapsed"],   t_spca["elapsed"],
							t_mixo["elapsed"],  t_nspr["elapsed"],
							t_nscum["elapsed"], t_pca["elapsed"]), 3)
)
write.csv(results_pitprops, "benchmarking/benchmarking_results_pitprops.csv",
					row.names = FALSE)
