# msPCA 0.5.1

## Bug fixes

* `summary.mspca()` no longer reports orthogonality violations for models fitted under uncorrelatedness constraints. The method defaulted `feasibilityConstraintType = 0` and `mspca()` did not record the type used to fit, so a user who fitted with `feasibilityConstraintType = 1` and called `summary()` without repeating the setting silently received the wrong diagnostic. The fitted type is now stored in the object and is the default for all diagnostics. Reported by a referee.
* `summary(object, feasibilityConstraintType = 1)` on a result with no `C` argument previously substituted the identity for the covariance matrix in the scalar violation, returning the orthogonality figure under the uncorrelatedness label, while the pairwise matrix errored on `NULL`. Both paths now use the values stored at fit time and need no `C`.

## New behaviour

* Uncorrelatedness violations are now normalized by the total variance `tr(Sigma)`, i.e. the pairwise term is `|u_t' Sigma u_s| / tr(Sigma)` rather than `|u_t' Sigma u_s|`. Because the loadings are unit-norm, the unnormalized quantity scaled linearly with `Sigma`, so `feasibilityTolerance` was effectively tighter on a covariance matrix with large variances than on the corresponding correlation matrix, and the same data expressed in different units gave different stopping behaviour. The normalized measure is invariant to a rescaling of `Sigma` and reads as a fraction of the total variance, on the same scale as the fraction of variance explained. This affects `mspca(..., feasibilityConstraintType = 1)` (the `feasibility_violation` field, the stopping rule, and the dual step size), `feasibility_violation_off()`, and the `uncorrelatedness` matrix in `nonredundancy`. Orthogonality violations are unchanged, as is the diagonal `|u_t' u_t - 1|` part of the solver's measure. Numerical results under `feasibilityConstraintType = 1` will differ from earlier versions unless the input has `tr(Sigma) = 1`; in particular a correlation matrix, for which `tr(Sigma) = p`, is also affected.
* The per-component penalty weights now follow a single rule for both constraint types: `weight_t = (u_t' Sigma u_t) / ||C u_t||^2`, fixed at the first outer iteration. The penalty a component can contribute is then bounded by `lambda * (u_t' Sigma u_t)` whatever `C` is, so `lambda` is a dimensionless penalty-to-objective ratio and the dual step-size schedule behaves identically under both constraints. Under orthogonality `C = I` and `||u_t|| = 1`, so this reduces exactly to the previous `weight_t = u_t' Sigma u_t` and orthogonality fits are bit-for-bit unchanged. Under uncorrelatedness it replaces the previous `weight_t = 1`, which — now that `C = Sigma / tr(Sigma)` — left the penalty weaker than the objective by a factor of order `p^2` and made the solver hit `maxIter` before reaching `feasibilityTolerance`. With the new weights, uncorrelatedness fits that previously terminated infeasible at the default `maxIter = 200` now converge; on `cor(mtcars)` with `r = 3, ks = 4` the violation goes from 4.5e-02 to 9.8e-05, and the `snp500` correlation matrix with `r = 3, ks = 10` reaches 2.6e-06 where neither the previous nor the pre-normalization behaviour reached tolerance at all.
* `mspca()` now records `feasibilityConstraintType` and `nonredundancy` in the returned object. `nonredundancy` holds two r x r matrices, `orthogonality` (\eqn{|u_t^\top u_s|}) and `uncorrelatedness` (\eqn{|u_t^\top \Sigma u_s| / \mathrm{tr}(\Sigma)}), computed after the columns have been sorted by explained variance, so indices always match the reported PC order. Both are stored regardless of which constraint was enforced, so inspecting the solution under the other definition requires neither the covariance matrix nor a refit.
* For `type = "X"` these are formed via `t(X) %*% (X %*% v)`, so the p x p matrix is still never materialized.
* `summary.mspca()` gains `feasibility_perPC` (per-PC maximum violation against any other PC), a `max_violation` column in the summary table, and `feasibilityConstraintType` / `fittedConstraintType` fields. Its printed output now names the constraint definition in use and points to the stored matrices for the other one.
* Passing `feasibilityConstraintType` explicitly to `summary.mspca()` still works, for inspecting a solution under the definition that was not enforced, but now warns when it differs from the fitted type.
* `mspca()` validates that `feasibilityConstraintType` is 0 or 1, and coerces it to integer.
* `type = "Sigma"` fits now carry the variable names on `x_best`, as `type = "X"` fits already did. `print()` therefore labels the loadings correctly without being passed `C`, which together with the stored variance figures and violation matrices means neither `print()` nor `summary()` needs the covariance matrix for any input type.

## Documentation

* Documented that `object$feasibility_violation` (from the solver) and the violation reported by `summary()` are different statistics: the former sums over all pairs `t >= s` and includes the normalization terms on the diagonal, the latter is the strictly off-diagonal sum of `feasibility_violation_off()`. The solver value is what `feasibilityTolerance` is compared against; the off-diagonal value is the redundancy diagnostic.
* The `C` argument of `summary.mspca()` is soft-deprecated. Every figure the method reports is now stored in the object; `C` is used only for objects fitted with msPCA 0.5.0 or earlier and will be removed in a future release.
* The "Algorithm and implementation notes" vignette now states the feasibility violation separately for each constraint type, including the `tr(C)` normalization; the introductory vignette shows `summary()` picking up the fitted constraint type and the stored `nonredundancy` matrices.
* Documentation fixes: `README.md` referred to the removed `print_mspca()` and to a `test/` directory that is now `notebooks/`; the `mtcars` vignette stated the wrong default for `scale` in `mspca()`.

## Data, vignettes and infrastructure

* Added the `snp500` dataset: the market-deflated correlation matrix of daily log-returns for 423 S&P 500 constituents, January 2010 - December 2019 (423 x 423, xz-compressed). Derived from a CC0-licensed Kaggle dataset by `data-raw/snp500.R`.
* Added the vignette "Case study: sparse factors in S&P 500 returns", a full application of `mspca()` to `snp500` comparing the two non-redundancy constraints.
* Added the vignette "Algorithm and implementation notes", documenting the optimization problem, both algorithms, the implicit matrix-vector implementation, computational complexity, and guidance on parameter choices.
* Added the website-only article "Benchmarking against other sparse PCA packages", comparing `mspca()` against seven competing implementations on four real datasets. It lives in `vignettes/articles/` and is not part of the CRAN build.
* Added `replication/`, the scripts reproducing the benchmarking and case-study results. Build-ignored.
* `_pkgdown.yml` now specifies the article ordering and a grouped reference index.
* `DESCRIPTION` gains `Depends: R (>= 3.5)`, `LazyData: true`, `LazyDataCompression: xz` and `BugReports`.

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
