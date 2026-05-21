## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a minor update (v0.4.0) focused on API clarity and documentation maintenance.
* The package implements the algorithm described in Cory-Wright and Pauphilet (2022) <doi:10.48550/arXiv.2209.14790>.

## Changes in this version

* Renamed hyperparameters controlling truncated power method restart budgets for clearer and more consistent API naming.
* Documentation polish across function docs and package materials.
* Removed `pairwise_correlation()` and `orthogonality_violation()` and replaced them with a unified `feasibility_violation_off()` helper.

## Test environments

* local macOS Sonoma 14.7.4, R 4.5.2
* win-builder (R-devel, Windows Server 2022)
* rhub: macOS, Ubuntu, Windows
