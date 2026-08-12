## data-raw/snp500.R -----------------------------------------------------------
##
## Builds data/snp500.rda: the market-deflated correlation matrix of daily
## log-returns for 423 S&P 500 constituents, January 2010 - December 2019.
##
## Run once from the package root:
##     source("data-raw/snp500.R")
##
## Two paths are provided:
##   (A) DEFAULT - rebuild from the derived matrix shipped in data-raw/.
##       Offline, deterministic, takes a second. Use this unless you want to
##       re-derive the numbers from the raw prices yourself.
##   (B) FULL PROVENANCE - rebuild from the raw Kaggle price file. Requires a
##       Kaggle account and reproduces (A) exactly for the stated date window.
##
## Source data: "S&P500 daily update dataset", CC0 1.0 (public domain),
##   https://www.kaggle.com/datasets/yash16jr/s-and-p500-daily-update-dataset
## Because the upstream dataset is updated daily, path (B) is only reproducible
## for the fixed window used here; path (A) is the archival record.
## ---------------------------------------------------------------------------

## ---- (A) Rebuild from the derived matrix (default) --------------------------

snp500 <- as.matrix(
  read.csv(xzfile("data-raw/snp500_corr.csv.xz"),
           row.names = 1, check.names = FALSE)
)
colnames(snp500) <- rownames(snp500)
storage.mode(snp500) <- "double"

## Enforce exact symmetry (guards against round-trip drift in the last digit).
snp500 <- (snp500 + t(snp500)) / 2

## ---- Sanity checks ----------------------------------------------------------

stopifnot(
  is.matrix(snp500),
  identical(dim(snp500), c(423L, 423L)),
  identical(rownames(snp500), colnames(snp500)),
  isSymmetric(snp500),
  all(is.finite(snp500))
)

ev <- eigen(snp500, symmetric = TRUE, only.values = TRUE)$values
stopifnot(
  abs(min(abs(ev))) < 1e-8,        # market factor removed => rank p - 1
  min(ev) > -1e-8                  # positive semidefinite up to rounding
)
message(sprintf("snp500: %d x %d, trace = %.4f, rank = %d",
                nrow(snp500), ncol(snp500), sum(diag(snp500)), sum(ev > 1e-8)))

## ---- Write package data -----------------------------------------------------

usethis::use_data(snp500, overwrite = TRUE, compress = "xz")
## Equivalent without usethis:
##   save(snp500, file = "data/snp500.rda", compress = "xz", version = 2)


## ---- (B) Full provenance: rebuild from raw Kaggle prices --------------------
## Not run. Download SnP_returns_cleaned.csv from the Kaggle link above (or
## regenerate it with code/case_study/, which converts adjusted closes to daily
## log-returns and drops tickers without a complete history over the window).

if (FALSE) {

  library("readr")
  library("dplyr")
  library("RSpectra")

  df_returns <- read_csv("SnP_returns_cleaned.csv") |>
    filter(Date < "2020-01-01")

  X <- as.matrix(df_returns[, -1])          # 2515 x 423 daily log-returns
  Sigma <- cor(X)

  ## S&P 500 returns are dominated by a market factor loading positively on
  ## nearly every stock. Remove it so the sparse components describe
  ## cross-sectional (sector and style) structure instead.
  v1 <- eigs_sym(Sigma, k = 1, which = "LA")$vectors[, 1]
  v1 <- v1 / sqrt(sum(v1^2))
  P  <- diag(length(v1)) - tcrossprod(v1)

  snp500 <- crossprod(P, Sigma) %*% P
  snp500 <- (snp500 + t(snp500)) / 2
  rownames(snp500) <- colnames(snp500) <- colnames(X)

  usethis::use_data(snp500, overwrite = TRUE, compress = "xz")
}
