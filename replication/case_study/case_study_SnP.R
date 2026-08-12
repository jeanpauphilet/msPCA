## ============================================================
## Section 5 Case Study: S&P 500 Factor Analysis
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Dataset: SnP_returns_cleaned.csv
##   423 stocks, Jan 2010 -- Dec 2019 (2,515 trading days)
##   Daily log-returns; market component removed via spectral
##   deflation before applying msPCA.
## ============================================================

library("msPCA")
library("RSpectra")
library("readr")
library("dplyr")
library("ggplot2")
library("nsprcomp")

## -----------------------------------------------------------
## 1. Load data and compute the correlation matrix
## -----------------------------------------------------------
df_returns <- read_csv("case_study/SnP_returns_cleaned.csv") |>
  filter(Date < "2020-01-01")
X  <- as.matrix(df_returns[, -1])  # drop Date column
Sigma  <- cor(X)

## -----------------------------------------------------------
## 2. Remove market component via spectral deflation
##
## The leading eigenvector v1 of S captures the "market factor"
## that loads positively on virtually every stock.  Projecting
## S onto the orthogonal complement of v1 suppresses this
## market-wide movement and exposes cross-sectional structure
## (sector / style effects).
## -----------------------------------------------------------
v1 <- eigs_sym(Sigma, k = 1, which = "LA")$vectors[, 1]
v1 <- v1 / sqrt(sum(v1^2))           # ensure unit norm
P  <- diag(length(v1)) - tcrossprod(v1)
SigmaR <- crossprod(P, Sigma) %*% P          # SigmaR = P Sigma P^T
rownames(SigmaR) <- colnames(SigmaR) <- colnames(X)  # preserve tickers

XR <- X %*% P
colnames(XR) <- colnames(SigmaR)

## -----------------------------------------------------------
## 3. Sweep over k for both constraint types
##    feasibilityConstraintType = 0  : orthogonality
##    feasibilityConstraintType = 1  : zero pairwise correlation
## -----------------------------------------------------------
ks_grid <- seq(5, 35, by = 5)
r       <- 4

run_sweep <- function(ctype) {
    do.call(rbind, lapply(ks_grid, function(k) {
      cat("  k =", k, "\n")
      set.seed(42)
      res <- mspca(SigmaR, r = r, ks = rep(k, r), maxIter = 100,
                   verbose = FALSE, 
                   feasibilityConstraintType = ctype
                   )
      data.frame(
        k              = k,
        constraint     = if (ctype == 0) "orthogonality"
                         else "zero-correlation",
        fve            = fraction_variance_explained(SigmaR, res$x_best),
        orth_violation = feasibility_violation_off(SigmaR, res$x_best, 0),
        total_pwcorr   = feasibility_violation_off(SigmaR, res$x_best, 1),
        nnz            = sum(abs(res$x_best) > 0)
      )
    }))
  }

cat("Running orthogonality sweep...\n")
res_orth <- run_sweep(0)
cat("Running zero-correlation sweep...\n")
res_corr <- run_sweep(1)


## Also run nsprcomp at same sparsity budgets for comparison
res_nsprcomp <- do.call(rbind, lapply(ks_grid, function(k) {
    cat("  nsprcomp k =", k, "\n")
    fit <- nsprcomp::nsprcomp(XR, ncomp = r, k = rep(k, r), nneg = FALSE, center = TRUE, scale. = TRUE)
        
    data.frame(
      k              = k,
      constraint     = "nsprcomp",
      fve            = fraction_variance_explained(SigmaR, fit$rotation),
      orth_violation = feasibility_violation_off(SigmaR, fit$rotation, 0),
      total_pwcorr   = feasibility_violation_off(SigmaR, fit$rotation, 1),
      nnz            = sum(abs(fit$rotation) > 0)
    )
  }))

results_df <- rbind(res_orth, res_corr, res_nsprcomp)

## -----------------------------------------------------------
## 3b. Save results table for paper
## -----------------------------------------------------------
write.csv(results_df, "case_study/snp_varyingk_results.csv", row.names = FALSE)


## -----------------------------------------------------------
## 4. Produce figures (saved to figures/)
## -----------------------------------------------------------

## Helper theme
theme_jss <- function() {
    theme_minimal(base_size = 10) +
      theme(plot.background  = element_rect(fill = "white",
                                            color = NA),
            panel.background = element_rect(fill = "white",
                                            color = NA),
            axis.text.x      = element_text(size = 11),
            axis.text.y      = element_text(size = 11),
            axis.title.x     = element_text(size = 11),
            axis.title.y     = element_text(size = 11),
            legend.position  = "bottom")
  }

pal   <- c("orthogonality"    = "#0173B2",
              "zero-correlation" = "#029E73",
              "nsprcomp"       = "#DE8F05")
ltype <- c("orthogonality"    = "solid",
              "zero-correlation" = "dashed",
              "nsprcomp"       = "dotted")
shp   <- c("orthogonality"    = 16,
              "zero-correlation" = 15,
              "nsprcomp"       = 17)

## Figure: FVE vs k
p_fve <- ggplot(results_df,
    aes(x = k, y = fve,
        color = constraint, linetype = constraint,
        shape = constraint)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = ltype) +
    scale_shape_manual(values = shp) +
    labs(x = "Sparsity budget k",
          y = "FVE",
         color = NULL, linetype = NULL, shape = NULL) +
    theme_jss()
ggsave("case_study/snp_fve.pdf", p_fve,
          width = 4, height = 3.0, units = "in")

## Figure: Orthogonality violation vs k
p_orth <- ggplot(results_df,
    aes(x = k, y = orth_violation,
        color = constraint, linetype = constraint,
        shape = constraint)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = ltype) +
    scale_shape_manual(values = shp) +
    labs(x = "Sparsity budget k",
         y = "Orthogonality violation ",
         color = NULL, linetype = NULL, shape = NULL) +
    theme_jss()
ggsave("case_study/snp_orth.pdf", p_orth,
          width = 4, height = 3.0, units = "in")

## Figure: Pairwise Sigma-correlation vs k
p_pwcorr <- ggplot(results_df,
    aes(x = k, y = total_pwcorr,
        color = constraint, linetype = constraint,
        shape = constraint)) +
    geom_line(linewidth = 1.1) + geom_point(size = 2.5) +
    scale_color_manual(values = pal) +
    scale_linetype_manual(values = ltype) +
    scale_shape_manual(values = shp) +
    labs(x = "Sparsity budget k",
         y = "Zero-correlation violation",
         color = NULL, linetype = NULL, shape = NULL) +
    theme_jss()
ggsave("case_study/snp_pwcorr.pdf", p_pwcorr,
          width = 4, height = 3.0, units = "in")

## -----------------------------------------------------------
## 5. Focal analysis at k = 10 (orthogonality constraint)
## -----------------------------------------------------------
set.seed(42)
res_k10 <- mspca(SigmaR, r = r, ks = rep(10, r),
                   verbose = TRUE,
                   feasibilityConstraintType = 0)
print_mspca(res_k10, SigmaR)
out <- capture.output(print_mspca(res_k10, SigmaR))
writeLines(out, "case_study/snp_k10_orthogonality_loadings.txt")

set.seed(42)
res_k10_corr <- mspca(SigmaR, r = r, ks = rep(10, r),
                   verbose = TRUE,
                   feasibilityConstraintType = 1)
print_mspca(res_k10_corr, SigmaR)
out <- capture.output(print_mspca(res_k10_corr, SigmaR))
writeLines(out, "case_study/snp_k10_correlation_loadings.txt")

## -----------------------------------------------------------
## 5b. Panel 1: Sparsity / composition heatmap at k = 10
## -----------------------------------------------------------

## Tidy helper: reshape a p x r loading matrix into long format
tidy_loadings <- function(res, label) {
  V <- res$x_best
  V[abs(V) < 1e-4] <- 0
  colnames(V) <- paste0("PC", seq_len(ncol(V)))
  rownames(V) <- rownames(SigmaR)
  as.data.frame(V) |>
    tibble::rownames_to_column("ticker") |>
    tidyr::pivot_longer(-ticker,
                        names_to  = "PC",
                        values_to = "loading") |>
    dplyr::mutate(constraint = label)
}

df_heat <- dplyr::bind_rows(
  tidy_loadings(res_k10,      "Orthogonal loadings"),
  tidy_loadings(res_k10_corr, "Uncorrelated PCs")
)

## Restrict to tickers that are non-zero in at least one
## method / PC combination
active_tickers <- df_heat |>
  dplyr::filter(abs(loading) > 0) |>
  dplyr::pull(ticker) |>
  unique()

df_heat_sub <- df_heat |>
  dplyr::filter(ticker %in% active_tickers) |>
  dplyr::mutate(PC = factor(PC, levels = paste0("PC", seq_len(r))))

## Order tickers: for each (method, PC) pair in sequence
##   orth-PC1, orth-PC2, orth-PC3, orth-PC4,
##   corr-PC1, corr-PC2, corr-PC3, corr-PC4
## add only stocks not yet listed, sorted by signed loading value
## (descending: most positive first, most negative last).
V_orth_full <- res_k10$x_best
V_corr_full <- res_k10_corr$x_best
rownames(V_orth_full) <- rownames(SigmaR)
rownames(V_corr_full) <- rownames(SigmaR)

ticker_order <- character(0)
for (V in list(V_orth_full, V_corr_full)) {
  for (k in seq_len(r)) {
    v_k       <- V[, k]
    new_ticks <- names(v_k)[abs(v_k) > 0 & !(names(v_k) %in% ticker_order)]
    new_ticks <- new_ticks[order(v_k[new_ticks], decreasing = TRUE)]
    ticker_order <- c(ticker_order, new_ticks)
  }
}

df_heat_sub <- df_heat_sub |>
  dplyr::mutate(ticker = factor(ticker, levels = rev(ticker_order)))

## Symmetric colour limits
lim <- max(abs(df_heat_sub$loading))

p_heat <- ggplot(df_heat_sub,
    aes(x = PC, y = ticker, fill = loading)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low      = "#0173B2",    # blue  → negative loading
    mid      = "white",
    high     = "#D55E00",    # red   → positive loading
    midpoint = 0,
    limits   = c(-lim, lim),
    name     = "Loading"
  ) +
  facet_wrap(~ constraint, ncol = 2) +
  labs(x = NULL, y = NULL) +
  theme_jss() +
  theme(
    axis.text.y      = element_text(size = 5.5, family = "mono"),
    axis.text.x      = element_text(size = 8),
    strip.text       = element_text(size = 9, face = "bold"),
    legend.key.width = unit(0.7, "cm"),
    panel.grid       = element_blank()
  )

ggsave("case_study/snp_heatmap.pdf", p_heat,
       width = 7.0, height = 6.0, units = "in")