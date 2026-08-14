## ============================================================
## Section 5 Case Study: S&P 500 Factor Analysis
## JSS replication script for "msPCA: Multiple Sparse PCA in R"
##
## Data: the deflated correlation matrix ships with the package as
##   data(snp500): 423 stocks, Jan 2010 -- Dec 2019 (2,515 trading
##   days), daily log-returns, market component removed via spectral
##   deflation.  Everything based on msPCA runs from this matrix
##   alone.
##
##   The raw returns (case_study/SnP_returns_cleaned.csv, ~36 MB) are
##   too large to ship and are only needed for the nsprcomp reference
##   curve, which requires a data matrix rather than a correlation
##   matrix.  Section 3 below detects whether that file is present
##   and falls back to cached nsprcomp results otherwise.
## ============================================================

library("msPCA")
library("ggplot2")
## nsprcomp and RSpectra are called with :: and only in the branch of
## section 3 that needs the raw returns; they are not required to run
## the rest of the script.

## -----------------------------------------------------------
## 0. Output destinations
##
## Everything this script produces is consumed in two places:
##   * replication/case_study/ -- PDFs and CSVs, for the paper;
##   * the package sources     -- PNGs under vignettes/figures/ and the
##     pre-computed sweep under inst/vignette-data/, for the vignette
##     and the pkgdown site.
##
## Writing only to the first is how the vignette came to be three days
## out of date against the replication outputs, still describing
## factors the re-run had replaced. save_figure() below and the copy at
## the end of section 3b keep the two in step, so a re-run of this
## script refreshes the website as well.
##
## Paths assume the script is sourced from the replication/ directory,
## which is what run_all.R and the README instruct.
## -----------------------------------------------------------
PKG_ROOT  <- ".."
FIG_DIR   <- file.path(PKG_ROOT, "vignettes", "figures")
VDATA_DIR <- file.path(PKG_ROOT, "inst", "vignette-data")

for (d in c("case_study", FIG_DIR, VDATA_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

## Save one figure twice: PDF for the paper, PNG for the vignette.
## 150 dpi matches the existing assets -- a 4 x 3 in panel becomes
## 600 x 450 px -- so figure sizing in the vignette is unaffected.
save_figure <- function(name, plot, width, height, dpi = 150) {
  ggsave(file.path("case_study", paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in")
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi)
  cat("  wrote case_study/", name, ".pdf and ", FIG_DIR, "/", name, ".png\n",
      sep = "")
  invisible(NULL)
}

## -----------------------------------------------------------
## 1. Load the market-deflated correlation matrix
##
## The leading eigenvector v1 of the raw correlation matrix Sigma
## captures the "market factor" that loads positively on virtually
## every stock.  Projecting Sigma onto the orthogonal complement of
## v1 suppresses this market-wide movement and exposes
## cross-sectional structure (sector / style effects):
##   SigmaR = P Sigma P^T,   P = I - v1 v1^T.
## This deflation is already applied in the shipped data object; see
## ?snp500 and data-raw/snp500.R for the derivation.
## -----------------------------------------------------------
data("snp500", package = "msPCA")
SigmaR <- snp500

## -----------------------------------------------------------
## 2. Sweep over k for both constraint types
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


## -----------------------------------------------------------
## 3. nsprcomp reference curve at the same sparsity budgets
##
## nsprcomp needs a data matrix, not a correlation matrix, so it
## cannot be run from snp500 alone.  If the raw returns file is
## available we rebuild the deflated data matrix XR and run the
## method, caching the results; otherwise we read the cached
## results back; if neither exists we warn and skip nsprcomp.
## -----------------------------------------------------------
returns_file  <- "case_study/SnP_returns_cleaned.csv"
nsprcomp_file <- "case_study/snp_nsprcomp_results.csv"

have_nsprcomp <- requireNamespace("nsprcomp", quietly = TRUE) &&
                 requireNamespace("RSpectra", quietly = TRUE)

if (file.exists(returns_file) && have_nsprcomp) {

  cat("Raw returns found: running nsprcomp...\n")

  df_returns <- utils::read.csv(returns_file, check.names = FALSE)
  df_returns <- df_returns[df_returns$Date < "2020-01-01", ]
  X     <- as.matrix(df_returns[, -1])   # drop Date column
  Sigma <- stats::cor(X)

  ## Reconstruct the deflation projector P = I - v1 v1^T used for snp500
  v1 <- RSpectra::eigs_sym(Sigma, k = 1, which = "LA")$vectors[, 1]
  v1 <- v1 / sqrt(sum(v1^2))             # ensure unit norm
  P  <- diag(length(v1)) - tcrossprod(v1)

  XR <- X %*% P
  colnames(XR) <- colnames(SigmaR)

  res_nsprcomp <- do.call(rbind, lapply(ks_grid, function(k) {
    cat("  nsprcomp k =", k, "\n")
    fit <- nsprcomp::nsprcomp(XR, ncomp = r, k = rep(k, r), nneg = FALSE,
                              center = TRUE, scale. = TRUE)

    data.frame(
      k              = k,
      constraint     = "nsprcomp",
      fve            = fraction_variance_explained(SigmaR, fit$rotation),
      orth_violation = feasibility_violation_off(SigmaR, fit$rotation, 0),
      total_pwcorr   = feasibility_violation_off(SigmaR, fit$rotation, 1),
      nnz            = sum(abs(fit$rotation) > 0)
    )
  }))

  ## Cache so the comparison remains reproducible without the raw data
  write.csv(res_nsprcomp, nsprcomp_file, row.names = FALSE)

} else if (file.exists(nsprcomp_file)) {

  cat("Raw returns not found: loading cached nsprcomp results from ",
      nsprcomp_file, "\n", sep = "")
  res_nsprcomp <- utils::read.csv(nsprcomp_file, stringsAsFactors = FALSE)

} else {

  warning("nsprcomp cannot be run (raw returns file '", returns_file,
          "' missing, or packages nsprcomp/RSpectra not installed) and no ",
          "cached results were found at '", nsprcomp_file, "'; the nsprcomp ",
          "comparison will be omitted from the results table and figures. ",
          "Rebuild the returns file from the source Kaggle dataset (see ",
          "section (B) of data-raw/snp500.R) to run it.", call. = FALSE)
  res_nsprcomp <- NULL

}

results_df <- rbind(res_orth, res_corr, res_nsprcomp)

## -----------------------------------------------------------
## 3b. Save results table for paper
## -----------------------------------------------------------
write.csv(results_df, "case_study/snp_varyingk_results.csv", row.names = FALSE)

## The vignette reads this table via system.file("vignette-data", ...), so the
## packaged copy has to be refreshed alongside the replication copy. Without
## this the vignette silently keeps rendering an older sweep.
file.copy("case_study/snp_varyingk_results.csv",
          file.path(VDATA_DIR, "snp_varyingk_results.csv"), overwrite = TRUE)
cat("  wrote case_study/snp_varyingk_results.csv and ",
    VDATA_DIR, "/snp_varyingk_results.csv\n", sep = "")


## -----------------------------------------------------------
## 4. Produce figures
##    Each is written twice by save_figure(): PDF to case_study/ for
##    the paper, PNG to vignettes/figures/ for the vignette and site.
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
save_figure("snp_fve", p_fve, width = 4, height = 3.0)

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
save_figure("snp_orth", p_orth, width = 4, height = 3.0)

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
save_figure("snp_pwcorr", p_pwcorr, width = 4, height = 3.0)

## -----------------------------------------------------------
## 5. Focal analysis at k = 10 (orthogonality constraint)
## -----------------------------------------------------------
set.seed(42)
res_k10 <- mspca(SigmaR, r = r, ks = rep(10, r),
                   verbose = TRUE,
                   feasibilityConstraintType = 0)
## print()/summary() take their figures from the fitted object; SigmaR does not
## need to be passed back in. summary() reports the violations under the
## constraint type that was enforced (orthogonality here) and labels them.
print(res_k10)
summary(res_k10)
out <- capture.output(print(res_k10))
writeLines(out, "case_study/snp_k10_orthogonality_loadings.txt")

set.seed(42)
res_k10_corr <- mspca(SigmaR, r = r, ks = rep(10, r),
                   verbose = TRUE,
                   feasibilityConstraintType = 1)
## Zero-correlation fit: summary() picks up feasibilityConstraintType = 1 from
## the object, so this reports pairwise correlations, not orthogonality.
print(res_k10_corr)
summary(res_k10_corr)
out <- capture.output(print(res_k10_corr))
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

save_figure("snp_heatmap", p_heat, width = 7.0, height = 6.0)