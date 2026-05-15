# Scalability Improvement Plan — `msPCA` package

Target file: `code/msPCA/src/msPCA_R_CPP.cpp` (plus `Makevars` and minor changes in `ConstantArguments.h`).

The goal is to bring `mspca()` and `tpw()` runtime into the same order of magnitude as `nsprcomp::nsprcomp` on the benchmarking datasets in `code/benchmarking/` (mtcars, pitprops, breast, riboflavin). Each change below lists the *what*, the *why*, the *expected gain*, the *risk*, and the *validation* step.

The plan is organised in three waves: high-impact algorithmic changes first, then targeted micro-optimisations, then build-system / numerical hygiene. Within each wave, changes are ordered so that each one can be merged and benchmarked independently.

---

## Wave 1 — High-impact algorithmic changes

### 1.1 Remove the full `SelfAdjointEigenSolver` call inside the double loop

**Current behaviour.** For each outer iteration and each PC `t`, the code instantiates `Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(sigma_current)` (cost O(n³)) for two purposes: (i) shift the diagonal to make `sigma_current` PSD, (ii) obtain its leading eigenvector as a warm start for the truncated power method (TPW).

**Why this is unnecessary.**
- The PSD diagonal shift adds a constant `λ₀` to `xᵀ M x`, which does not change the argmax over the sparse unit sphere. The TPW iteration `x ← Truncate(Mx)` is invariant to a diagonal shift up to scale.
- After the first outer iteration, `x_current.col(t)` is an excellent warm start, strictly better than the leading eigenvector of the deflated matrix.

**Change.**
- Drop the diagonal shift entirely. Delete the `solver`/`lambda0` block.
- On iteration 1, build `beta0` with a handful of plain power iterations (≈ 20 steps) against the operator defined in §1.2.
- On iterations ≥ 2, set `beta0 = x_current.col(t)`.

**Expected gain.** This is the single largest win. For n=200, r=5, maxIter=200 this block alone costs ≈ 10⁹ flops per run. Empirically ≥ 5× speedup on medium-size problems.

**Risk.** Iteration 1 needs a deterministic, reasonable seed. Mitigation: keep a couple of random restarts (see §1.3) on iteration 1 only.

**Validation.** On mtcars / pitprops, the recovered support and explained-variance numbers should match the existing CSVs to within numerical tolerance.

---

### 1.2 Stop materialising `sigma_current`; keep the deflation factored

**Current behaviour.** Every PC update rebuilds the n×n deflated matrix as `Sigma - Σₛ λ wₛ wₛᵀ`, symmetrises it, and feeds it to TPW. Each rank-1 update is O(n²); for r PCs that's O(r·n²) per PC, plus an O(n²) symmetrisation.

**Change.** Replace the materialised matrix by a matvec functor used inside `iterativeTruncHeuristic`:

```cpp
// W: n x (r-1) matrix of the "w_s" vectors for s != t
// d: (r-1) vector of theLambda * weights[s]
auto applyM = [&](const Eigen::VectorXd& beta) -> Eigen::VectorXd {
  Eigen::VectorXd y = Sigma * beta;            // O(n^2)
  Eigen::VectorXd c = W.transpose() * beta;    // O(r n)
  y.noalias() -= W * (d.asDiagonal() * c);     // O(r n)
  return y;
};
```

Then `iterativeTruncHeuristic`, `singlePCHeuristic`, and `evaluate` operate via `applyM` rather than a dense matrix.

**Expected gain.** Eliminates O(r·n²) build cost per PC update and one O(n²) allocation + symmetrisation. Mat-vec cost inside TPW stays at O(n²) but with a smaller constant and fewer allocations.

**Risk.** API ripple: `singlePCHeuristic` and `iterativeTruncHeuristic` need to take a `std::function` (or templated functor) instead of an `Eigen::MatrixXd&`. Keep a thin overload that builds the functor from a matrix, so `truncatedPowerMethod` (single-PC) is unchanged.

**Validation.** Unit test: for a small Σ and known W, verify that `applyM(beta)` matches the explicit `sigma_current * beta` to 1e-12.

---

### 1.3 Gate the random restarts inside `singlePCHeuristic`

**Current behaviour.** Each call to `singlePCHeuristic` runs up to `maxIterTPM` (default 200) random restarts, each doing 100 truncated-power iterations. That's up to 20 000 matvecs per PC per outer iteration. nsprcomp uses ~50–200 iterations *total* per PC.

**Change.**
- Add an `int outerIter` parameter to `singlePCHeuristic`.
- Set `restartBudget = (outerIter == 1) ? maxIterTPM : k_small` where `k_small` is a small constant (default 2).
- Optional refinement: skip restarts entirely if the warm-started run already produced an objective ≥ the previous outer-iteration's per-PC objective.

**Expected gain.** 10–50× reduction in the per-PC inner work for outer iterations ≥ 2, which is where most of the wall time goes once the warm start is good.

**Risk.** Solution quality on hard instances may regress in the very early iterations. Mitigation: keep `k_small` exposed as a parameter (e.g. `restartsAfterFirstIter`) with sensible default, and keep the full restart budget for iteration 1.

**Validation.** On the pitprops and breast benchmarks, confirm that final variance-explained matches the existing results.

---

### 1.4 Add a convergence check to `iterativeTruncHeuristic`

**Current behaviour.** Hard-coded 100 iterations regardless of progress.

**Change.**
```cpp
for (int i = 0; i < 100; ++i) {
  Eigen::VectorXd prev = beta;
  beta = truncate(applyM(beta), k);   // truncate normalises in place
  if ((beta - prev).squaredNorm() < 1e-12) break;
}
```

A support-stability check (no index changed in the top-k for two consecutive iters) is even cheaper if profiling shows the squared-norm dominates.

**Expected gain.** 2–5× reduction in inner iterations on most instances (TPW typically converges in 5–20 iterations once warm-started).

**Risk.** None significant if the tolerance is conservative.

**Validation.** Compare iteration counts via a small debug counter on the benchmarking instances before and after.

---

## Wave 2 — Targeted micro-optimisations

### 2.1 Rewrite `truncateVector` / `selectperm2` to avoid R-side allocation

**Current behaviour.** `selectperm2` builds an `Rcpp::NumericVector` via `push_back` inside the hot path; `truncateVector` makes two full-vector copies.

**Change.**
- Make both functions purely C++ (no Rcpp types), returning a `std::vector<int>` or operating in place.
- Replace `std::partial_sort` (O(n log k)) with `std::nth_element` (O(n) average).
- Skip the second copy: write directly into the output vector.

Reference implementation:

```cpp
inline Eigen::VectorXd truncateVector(const Eigen::VectorXd& v, int k) {
  const int n = static_cast<int>(v.size());
  if (k >= n) { Eigen::VectorXd y = v; y.normalize(); return y; }
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::nth_element(idx.begin(), idx.begin() + k, idx.end(),
                   [&](int a, int b){ return std::abs(v[a]) > std::abs(v[b]); });
  Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < k; ++i) y[idx[i]] = v[idx[i]];
  y.normalize();
  return y;
}
```

**Expected gain.** Modest in absolute terms but called thousands of times per run; typically 1.5–3× speedup of the inner iteration when n is large.

**Risk.** None.

**Validation.** Unit test against the current implementation: identical output up to ordering inside the top-k.

---

### 2.2 Cache `Sigma * x_current` per outer iteration

**Current behaviour.** `ofv_overall = evaluate(x_current, Sigma)` and `fnviolation(..., feasibilityConstraintType==1)` both recompute `Σ x_current` (an n×r matvec).

**Change.** Compute `SX = Sigma * x_current` once at the bottom of each outer iteration; reuse for the trace objective and for the uncorrelatedness violation. When building the deflation operator for the next outer iteration, the columns `Σ uₛ` (needed when `feasibilityConstraintType == 1`) are also already in `SX`.

**Expected gain.** Small (one matrix product saved per outer iter) but trivial to implement.

**Risk.** None.

---

### 2.3 Replace `absoluteDouble` and clean up `evaluate`

**Change.**
- Delete `absoluteDouble`; use `std::abs` (it inlines and may auto-vectorise).
- Rewrite `evaluate(vec, A)` as `x.dot(A * x)` instead of `(xᵀ A x)[0]` — avoids the temporary 1×1 matrix.
- Rewrite the matrix overload of `evaluate` as `(A.selfadjointView<Eigen::Upper>() * x).cwiseProduct(x).sum()` if profiling shows the trace form dominates.

**Expected gain.** Minor, but it cleans up the hot path.

**Risk.** None.

---

### 2.4 Skip the explicit symmetrisation

`sigma_current = (sigma_current + sigma_current.transpose()) / 2;` is an O(n²) op preserving a property (symmetry) that is already exact up to FP rounding for symmetric rank-1 subtractions from a symmetric Σ. Once Wave 1 §1.2 is in, this disappears anyway. Track it here so it doesn't sneak back in.

---

## Wave 3 — Build-system and numerical hygiene

### 3.1 Add Eigen-friendly compile flags to `Makevars`

**Change.**
```
PKG_CXXFLAGS = -DEIGEN_NO_DEBUG -DNDEBUG
CXX_STD = CXX17
```

`-DEIGEN_NO_DEBUG` disables Eigen's per-access bounds checks, which compound inside the inner loops. `-DNDEBUG` removes `assert()` overhead in any third-party headers. Both are CRAN-acceptable.

**Do not** put `-O3` or `-march=native` in the package `Makevars` (CRAN policy). If you want to publish faster wheels, document the option for users to add to `~/.R/Makevars`.

**Expected gain.** 1.2–2× on tight Eigen loops.

**Risk.** None on supported platforms.

---

### 3.2 Reproducible RNG and `RNGScope`

Inside `singlePCHeuristic`, `R::rnorm` is called many times. Wrap the loop body with `Rcpp::RNGScope` (or `BEGIN_RCPP`/`END_RCPP` if Rcpp doesn't already wrap it) and consider drawing all random doubles in a single `rnorm(n*k_restarts)` call when `k_restarts` is known, then reshaping in place. Minor, but it makes runs reproducible under `set.seed()` from R and reduces per-scalar RNG overhead.

---

### 3.3 Optional: thin SVD path when Σ is supplied as `X` instead of covariance

`nsprcomp` takes the data matrix `X` (n_obs × p) and works against thin matvecs against `X` and `Xᵀ`. When `n_obs ≪ p` (e.g. riboflavin) this is dramatically cheaper than touching Σ at all. Add an optional code path that accepts `X` and uses `applyM(beta) = Xᵀ (X beta) / (n_obs - 1) - <deflation terms>`. This is a larger change and only worth it if benchmarks show Σ-based matvecs as the remaining bottleneck.

---

## Implementation order and milestones

Concretely, I propose to land the work as four PR-sized commits, each independently benchmarkable against the existing CSVs in `code/benchmarking/`:

**Commit A — Wave 1.1 + 1.2 + Wave 3.1.** Remove eigensolver, factor the deflation, add compile flags. Expected to deliver most of the speedup.

**Commit B — Wave 1.3 + 1.4.** Gate restarts and add convergence check inside `iterativeTruncHeuristic`. Expected to roughly halve the remaining inner-loop time.

**Commit C — Wave 2 (all four).** Truncation rewrite, caching, micro-cleanups. Polishes the hot path.

**Commit D — Optional: Wave 3.3.** Data-matrix path. Only if benchmarks still show a gap with nsprcomp on the wide datasets (riboflavin, breast).

## Validation protocol

Each commit must pass three checks before merge:

1. **Numerical regression.** Run `code/benchmarking/notebook_{mtcars,pitprops,breast,riboflavin}.R` and confirm that `objective_value`, `feasibility_violation`, and the support sets match the committed CSVs to within tolerance (`abs diff < 1e-4` on objective, identical supports up to label permutation).

2. **Unit tests.** Anything in `code/msPCA/test/` plus new tests for `applyM` correctness (§1.2) and `truncateVector` equivalence (§2.1).

3. **Wall-clock benchmark.** Time `mspca()` against `nsprcomp::nsprcomp()` on the four datasets and record results in a new `benchmarking_results_scalability.csv`. Target: within 2× of nsprcomp on pitprops and mtcars, within 3× on breast and riboflavin (the wide datasets where Wave 3.3 may be needed).

## Out of scope

- Algorithmic redesign of the outer iterative-deflation scheme. The penalty-based deflation is the core contribution of the paper; we keep its semantics intact and only change how each step is computed.
- Parallelism. Single-thread improvements first; OpenMP for the `r`-PC loop and for batched matvecs is a natural follow-up but is not part of this plan.
