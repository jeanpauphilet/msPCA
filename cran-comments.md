## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a minor update (v0.2.0) with new features and improvements.
* The package implements the algorithm described in Cory-Wright and Pauphilet (2022) <doi:10.48550/arXiv.2209.14790>.

## Changes in this version

* Added support for no-correlation constraints between PCs as well as orthogonality constraints. User chooses between orthogonality and uncorrelatedness constraints via the `feasibilityConstraintType` parameter to `msPCA()`.
* Renamed return field from `orthogonality_violation` to `feasibility_violation` to support both constraint types. 
* Renamed function `feasibility_violation()` as `orthogonality_violation()` to be more explicit
* Created function `pairwise_correlation()`
* Improved warning messages and renamed return fields for clarity

## Test environments

* local macOS Sonoma 14.7.4, R 4.5.2
* win-builder (R-devel, Windows Server 2022)
* rhub: macOS, Ubuntu, Windows
