## R CMD check results

0 errors | 0 warnings | 1 note

The check was run with `R CMD check --as-cran msPCA_0.5.1.tar.gz` on macOS
Sonoma 14.8.5 with R 4.6.0 (aarch64-apple-darwin23).

## Summary of changes

Version 0.5.1 adds the `snp500` dataset and two CRAN vignettes covering a
full S&P 500 case study and the algorithm and implementation details. It also
adds a website-only benchmarking article and replication scripts, which are
excluded from the CRAN source package.

The release also updates package metadata and documentation, including the
minimum R version, compressed lazy data, the GitHub issue tracker, and the
corresponding generated dataset documentation.

## Test environments

* macOS Sonoma 14.8.5, R 4.6.0 (local)

## Notes

The single NOTE concerns HTML validation:

* `Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.`

This is an environment-specific limitation on the local machine. The package
builds successfully, passes all tests and examples, and passes the remaining
CRAN checks.
