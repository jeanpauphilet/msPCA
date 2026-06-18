## R CMD check results

0 errors | 0 warnings | 2 notes

This update includes API and documentation improvements, including raw-data input support and S3 print/summary methods for msPCA objects.

## Changes in this version

* Added dual input mode for `mspca()` and `tpm()`:
	- `type = "Sigma"` for covariance/correlation matrices
	- `type = "X"` for raw data matrices
* Added raw-data preprocessing controls (`center`, `scale`, `divisor`) and validation checks.
* Added S3 methods `print.mspca()` and `summary.mspca()`; kept `print_mspca()` as a deprecated wrapper for backward compatibility.
* Added `variance_explained` and `total_variance` in `mspca()` outputs.
* Added a worked vignette and refreshed package documentation.

## Test environments

* local macOS Sonoma 14.8.5, R 4.6.0
* local macOS Sonoma 14.8.5, R CMD build + R CMD check --as-cran (source tarball)

## Notes

* CRAN incoming feasibility note:
	- Days since last update: 6
	- This is an expected resubmission timing note.
* HTML/manual notes are environment/tooling-related on the local machine:
	- HTML validation skipped because system `tidy` is not recent enough.
	- Math rendering check skipped because package `V8` is unavailable locally.
