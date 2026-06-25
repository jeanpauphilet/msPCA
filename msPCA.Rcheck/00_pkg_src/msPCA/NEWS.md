# msPCA 0.5.0

* `mspca()` and `tpm()` now take two possible inputs: the covariance/correlation matrix or the data matrix directly. In practice, the functions take single generic argument `M` together with a `type = c("Sigma", "X")` selector. `type = "Sigma"` (the default) treats `M` as a covariance/correlation matrix (p x p); `type = "X"` treats `M` as a raw data matrix (n observations x p variables). The `"Sigma"` default preserves the behaviour of existing matrix-based calls.
* The raw-data path applies the algorithm to the data directly: each product `Sigma %*% beta` is evaluated as `t(X) %*% (X %*% beta) / (n - 1)` at cost O(np), and the p x p matrix is never materialized. This substantially improves scalability when `n << p`. The covariance back-end was refactored behind a covariance-operator abstraction (`DenseOp` / `GramOp`) shared by both input modes.
* Added preprocessing controls for `type = "X"`: `center`, `scale` (covariance vs correlation), and `divisor` ("n-1" or "n").
* Added validation for both input modes: a `Sigma` input is checked for squareness, symmetry and positive semidefiniteness (`checkPSD`, `symTolerance`, `psdTolerance`); an `X` input is checked for finiteness, dimensions and (when scaling) zero-variance columns.
* `mspca()` results now include `variance_explained` (per-PC) and `total_variance`; `X`-mode results also record `inputType`, `center`, `scale`, `divisor`, `nObs` and `p`.
* `mspca()` and `tpm()` now return S3 objects of class `"mspca"` and `"tpm"` respectively, enabling use of standard R generics.
* Added `print.mspca()`: S3 print method displaying the sparse loading matrix restricted to the union of active variables, the percentage of variance explained per PC, and the number of non-zero loadings. Replaces the removed `print_mspca()`.
* Added `summary.mspca()`: produces a per-PC table of sparsity, variance explained, FVE, and cumulative FVE, followed by the full pairwise feasibility violation matrix.
* Updated citation

# msPCA 0.4.1

* Standardized function man page titles to consistent title style.
* Removed unnecessary `library(datasets)` calls from examples while keeping explicit `datasets::mtcars` usage, and added `datasets` to `Suggests` to align example dependencies with CRAN guidance.
* Improved efficiency and clarity of R code 
* Added a vignette

# msPCA 0.4.0

* Renamed hyperparameters controlling truncated power method restart budgets for clearer and more consistent API naming.
* Documentation polish across function docs and package materials.
* Removed `pairwise_correlation()` and `orthogonality_violation()` and replaced them with a unified `feasibility_violation_off()` helper for feasibility diagnostics across constraint types.

# msPCA 0.3.0

* Improved scalability of `mspca()` and `tpw()` through algorithmic and implementation optimizations.
* Function `mspca()` now accepts a new hyper-parameter `minRestartTPM` that limits the number of random restarts for the truncated power method after the first outer iteration
* Improved scaling of the penalty parameters for the case of zero-correlation constraints
* Fixed incorrect acronym for truncated power method (TPW <- TPM) 

# msPCA 0.2.0

* Added support for no-correlation constraints between PCs as well as orthogonality constraints. User chooses between orthogonality and uncorrelatedness constraints via the `feasibilityConstraintType` parameter to `msPCA()`.
* Renamed return field from `orthogonality_violation` to `feasibility_violation` to support both constraint types. 
* Renamed function `feasibility_violation()` as `orthogonality_violation()` to be more explicit
* Created function `pairwise_correlation()`
* Added warning message when no feasible solution is found


# msPCA 0.1.0

* Initial CRAN release
