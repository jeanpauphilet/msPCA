# msPCA 0.2.0

* Added support for no-correlation constraints between PCs as well as orthogonality constraints. User chooses between orthogonality and uncorrelatedness constraints via the `feasibilityConstraintType` parameter to `msPCA()`.
* Renamed return field from `orthogonality_violation` to `feasibility_violation` to support both constraint types. 
* Renamed function `feasibility_violation()` as `orthogonality_violation()` to be more explicit
* Created function `pairwise_correlation()`
* Added warning message when no feasible solution is found


# msPCA 0.1.0

* Initial CRAN release
