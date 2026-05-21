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
