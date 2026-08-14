############################################################
## In this notebook, we illustrate the use of the msPCA package on synthetic data.
## We generate a spiked-covariance matrix with 2 4-sparse spikes (or PCs).
##
## We generate datasets of increasing size consistent with this model, compute
## the empirical covariance matrix, and apply our msPCA algorithm on the data.
## We measure performance in terms of support recovery of the true PCs, fraction
## of variance explained by the PCs, and orthogonality violation.
##
## We compare with the spca method from the elasticnet package.
############################################################
required_packages <- c("mvtnorm", "msPCA", "elasticnet", "readr", "dplyr",
                       "ggplot2", "RColorBrewer")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(mvtnorm)
library(msPCA)
library(elasticnet)

set.seed(42)

## -----------------------------------------------------------
## Output destinations
##
## The two figures produced below are displayed in README.md and on the
## pkgdown home page, both of which read them from man/figures/ -- not
## from notebooks/. Writing only to notebooks/ is how the README came to
## be showing figures two months older than the results behind them, so
## every figure is now saved to both locations.
##
## The script may be sourced from the repository root or from
## replication/ (which is what run_all.R does), so locate the package
## root by looking for DESCRIPTION rather than assuming the working
## directory.
## -----------------------------------------------------------
find_pkg_root <- function() {
  for (cand in c(".", "..", "../..")) {
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(normalizePath(cand))
  }
  stop("Could not locate the package root: no DESCRIPTION found in '.', '..' ",
       "or '../..'. Run this script from the repository root or from ",
       "replication/.", call. = FALSE)
}

PKG_ROOT   <- find_pkg_root()
output_dir <- file.path(PKG_ROOT, "notebooks")
manfig_dir <- file.path(PKG_ROOT, "man", "figures")
for (d in c(output_dir, manfig_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

results_csv <- file.path(output_dir, "msPCA_synthetic_results.csv")
fig_ortho <- file.path(output_dir, "synthetic_orthogonality_violation.png")
fig_fve <- file.path(output_dir, "synthetic_variance_explained.png")

## Save a figure to notebooks/ (working copy) and to man/figures/ (the
## copy README.md and the pkgdown site actually display).
save_figure_both <- function(path, plot, ...) {
  ggplot2::ggsave(filename = path, plot = plot, ...)
  ggplot2::ggsave(filename = file.path(manfig_dir, basename(path)),
                  plot = plot, ...)
  cat("  wrote ", path, "\n    and ", file.path(manfig_dir, basename(path)),
      "\n", sep = "")
  invisible(NULL)
}

p = 50 #Dimension
r = 2 #Number of sparse PCs
k = 20 #Sparsity of each PC
b = 2 #Signal strength

#Construct random -1,0,+1 vectors with disjoint supports (use 1/sqrt(k) to ensure norm 1)
xtrue <- matrix(0,p,r)
xtrue[1:k,1] <- sign(runif(k) - .5)/sqrt(k)
xtrue[(k+1):(2*k),2] <- sign(runif(k) - .5)/sqrt(k)

shufflecoords <- sample(1:p) #Shuffle coordinates
xtrue <- xtrue[shufflecoords,]

#Compute support recover accuracy metric
accuracy <- function(xtrue, U){
  a1 <- sum( (abs(U) > 0)*(abs(xtrue) > 0) )
  a2 <- sum( (abs(U[,c(2,1)]) > 0)*(abs(xtrue) > 0) )
  max(a1,a2) / sum(abs(xtrue) > 0)
}

#True covariance matrix
Strue <- b*outer(xtrue[,1],xtrue[,1]) + b*outer(xtrue[,2],xtrue[,2]) +diag(p)


## We comare the algorithms on varying
resdf <- data.frame(n=numeric(), package=character(),
                   sparsity=numeric(), accuracy=numeric(),
                   ortho_viol=numeric(), varexplained=numeric())
icol <- 1
for (iter in 1:5){ #Number of replications

  Xfull <- mvtnorm::rmvnorm(2000, mean = numeric(p), sigma = Strue)

  for (n in c(25,50,75,100,150,200,250,500,1000,1500)){ #Size of the dataset used to recover the PCs
  #for (n in seq(25, 250, by = 25)){ #Size of the dataset used to recover the PCs
    print(n)

    X <- Xfull[1:n,]
    S <- cor(X)

    #msPCA algorithm
    set.seed(42)
    mspca_results <- msPCA::mspca(S, r, rep(k,r), verbose=FALSE)
    resdf[icol,] <- list(n,"msPCA",
                   sum(abs(mspca_results$x_best) > 0), #Sparsity level
                   accuracy(xtrue, mspca_results$x_best), #Support recovery accuracy
                   feasibility_violation_off(S, mspca_results$x_best, 0), #Orthogonality violation
                   fraction_variance_explained(Strue,mspca_results$x_best)) #Fraction of variance explained

    icol <- icol +1

    #msPCA algorithm working on the X matrix directly
    set.seed(42)
    mspca_X_results <- msPCA::mspca(X, 2, c(k,k), verbose=FALSE, type="X")

    resdf[icol,] <- list(n,"msPCA - X",
                         sum(abs(mspca_X_results$x_best) > 0), #Sparsity level
                         accuracy(xtrue, mspca_X_results$x_best), #Support recovery accuracy
                         feasibility_violation_off(S, mspca_X_results$x_best, 0), #Orthogonality violation
                         fraction_variance_explained(Strue,mspca_X_results$x_best)) #Fraction of variance explained

    icol <- icol +1

    #elasticnet algorithm
    enet_results <- elasticnet::spca(S, 2, sparse="varnum", para = c(k,k), type="Gram")
    resdf[icol,] <- list(n,"elasticnet",
                      sum(abs(enet_results$loadings) > 0),
                      accuracy(xtrue, enet_results$loadings),
                      feasibility_violation_off(S, enet_results$loadings, 0),
                      fraction_variance_explained(Strue,enet_results$loadings))
    icol <- icol +1

  }
 }

library(readr)
write_csv(resdf, results_csv) # For saving results in test/ when run from repo root

resdf <- read_csv(results_csv)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
std.error <- function(x){sd(x)/sqrt(length(x))}

sumdf <- resdf %>%
  group_by(n,package) %>%
  summarise(accuracy=mean(accuracy),
            ortho_viol_se = std.error(ortho_viol),
            ortho_viol_mean =mean(ortho_viol),
            varexplained_se = std.error(varexplained),
            varexplained_mean =mean(varexplained),
            nreplicates=n()
            ) %>%
  ungroup() %>%
  mutate(ortho_viol_min=ortho_viol_mean-2*ortho_viol_se,
         ortho_viol_max=ortho_viol_mean+2*ortho_viol_se,
         varexplained_min=varexplained_mean-2*varexplained_se,
         varexplained_max=varexplained_mean+2*varexplained_se
         )

p_ortho <- sumdf %>%
  #filter(n <= 500) %>%
  ggplot() + aes(x=n, y=ortho_viol_mean, group=package, color=package,shape=package) +
  geom_line(size=0.7) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=ortho_viol_min,ymax=ortho_viol_max), width=0.1) +
  theme_minimal() +
  scale_color_brewer(palette="Set1") +
  labs(x="Sample size n", y="Orthogonality violation", group="Package", shape="Package", color="Package") +
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

print(p_ortho)
save_figure_both(fig_ortho, p_ortho, width = 8, height = 5, dpi = 300)

p_fve <- sumdf %>%
  #filter(n <= 500) %>%
  ggplot() + aes(x=n, y=varexplained_mean, group=package, color=package, shape=package) +
  geom_line(size=0.7) + geom_point(size=3) +
  geom_errorbar(aes(ymin=varexplained_min,ymax=varexplained_max), width=0.1) +
  theme_minimal() +
  geom_hline(aes(yintercept=6/54, linetype="Information theoretic upper bound"), color="black") +
  scale_color_brewer(palette="Set1") +
  labs(x="Sample size n", y="Out-of-sample fraction of variance explained",
       linetype="", group="Package", shape="Package", color="Package") +
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"))

print(p_fve)
save_figure_both(fig_fve, p_fve, width = 8, height = 5, dpi = 300)

