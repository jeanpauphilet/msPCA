## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a minor update (v0.3.0) focused on scalability and package maintenance.
* The package implements the algorithm described in Cory-Wright and Pauphilet (2022) <doi:10.48550/arXiv.2209.14790>.

## Changes in this version

* Improved scalability of `mspca()` and `tpw()` through algorithmic and implementation optimizations.
* Function `mspca()` now accepts a new hyper-parameter `restartsAfterFirstIter` that limits the number of random restarts for the truncated power method after the first outer iteration
* Improved scaling of the penalty parameters for the case of zero-correlation constraints

## Test environments

* local macOS Sonoma 14.7.4, R 4.5.2
* win-builder (R-devel, Windows Server 2022)
* rhub: macOS, Ubuntu, Windows
