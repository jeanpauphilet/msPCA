## R CMD check results

0 errors | 0 warnings | 1 note

## Summary of changes

Version 0.5.0 includes significant API improvements and additional S3 methods:

* Added dual input mode for `mspca()` and `tpm()`:
    - `type = "Sigma"` for covariance/correlation matrices (default, preserves backward compatibility)
    - `type = "X"` for raw data matrices with O(np) complexity
* Added raw-data preprocessing controls: `center`, `scale`, `divisor`, and comprehensive validation
* Added S3 classes and methods: `print.mspca()` and `summary.mspca()`
* Added `variance_explained` and `total_variance` to output
* Implemented covariance operator abstraction for scalability
* Updated citation with journal publication DOI: https://doi.org/10.1287/opre.2023.0598
* Added comprehensive vignette with worked example

## Test environments

* macOS Sonoma 14.8.5, R 4.6.0 (local)
* R CMD check --as-cran passed locally

## Notes

The single NOTE concerns HTML validation and math rendering tools:
- `Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.`
- `Skipping checking math rendering: package 'V8' unavailable`

These are environment-specific limitations on the local machine and do not indicate issues with the package itself. The package builds successfully with `devtools::build()` and passes all R CMD check tests.

## Technical terms in spell check

The spell check output includes technical domain-specific terms and abbreviations that are correct:
- Acronyms: TPM (Truncated Power Method), TPW (old acronym), PSD (Positive Semidefinite), FVE (Fraction of Variance Explained)
- Mathematical terms: decorrelation, semidefiniteness, uncorrelatedness, interpretability, scalability
- Publication references: doi, opre (Operations Research journal)
- Package/function names: mspca, roxygen, mtcars, Zhang (author)

All are appropriate in context.
