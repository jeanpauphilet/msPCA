## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a maintenance update (v0.4.1) focused on documentation consistency, dependency declarations for examples, and R-side performance improvements.
* The package implements the algorithm described in Cory-Wright and Pauphilet (2026) <doi:10.48550/arXiv.2209.14790>.

## Changes in this version

* Standardized function man page titles to consistent title style.
* Removed unnecessary `library(datasets)` calls from examples while keeping explicit `datasets::mtcars` usage.
* Added `datasets` to `Suggests` to align example dependencies with CRAN guidance.
* Improved efficiency in `feasibility_violation_off()` using `crossprod()` formulations.
* Optimized `variance_explained_perPC()` and `fraction_variance_explained()` with vectorized matrix formulations.
* Updated `print_mspca()` to use vectorized sparsity counting and added a `digits` argument for user-configurable print precision.

## Test environments

* local macOS Sonoma 14.8.5, R 4.4.3
* local macOS Sonoma 14.8.5, R CMD build + R CMD check --as-cran (source tarball)

## Notes

* `unable to verify current time` is an environment-specific check note on local macOS.
* `README.md`/`NEWS.md` cannot be checked without `pandoc`; this is due to local tooling and does not affect package functionality.
