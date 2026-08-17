############################################################
## bench_utils.R -- runtime, memory and solution-quality
## instrumentation for the msPCA benchmarking notebooks.
##
## ---------------------------------------------------------
## Memory
## ---------------------------------------------------------
## Why not gc() / bench::mark(mem_alloc)?
##   msPCA does its work in RcppEigen buffers allocated with C++
##   `new`, which R's garbage collector never accounts for. The same
##   applies to the embedded CPython heap used by the reticulate /
##   scikit-learn comparison. Any R-level memory metric therefore
##   under-reports msPCA and sklearn while fully charging pure-R
##   competitors such as elasticnet -- an unfair comparison.
##   We use getrusage(RUSAGE_SELF).ru_maxrss instead: the OS
##   high-water resident set size, which counts the R heap, Eigen,
##   BLAS and the Python interpreter alike.
##
## Why ONE repetition per process?
##   ru_maxrss is a monotone high-water mark for the whole process,
##   so several methods in one session, or even several repetitions
##   of one method, contaminate each other. Each repetition therefore
##   gets its own fresh R process via callr::r(). The algorithm is
##   invoked exactly as many times as before, so this costs only the
##   extra R startups, and in exchange every repetition yields an
##   INDEPENDENT peak-RSS sample.
##
##   That independence matters, because peak RSS is not deterministic
##   even for deterministic code: it depends on when the garbage
##   collector happens to run. Measured across two identical runs of
##   Pitprops, with identical seeds and bit-identical FVE and
##   orthogonality, nscumcomp's working set came out as 72.4 MB and
##   then 211.2 MB, and mixOmics as 0.6 MB and then 13.5 MB. A single
##   sample is not reportable. Memory is therefore summarised as a
##   MEDIAN over repetitions with min/max alongside, and the spread
##   should be quoted whenever it is wide.
##
## Why there is no warm-up
##   mem_delta_mb includes whatever a package allocates on its first
##   call -- lazy loading of namespaces, compiled code, S4 method
##   tables -- which is roughly 4 MB. The obvious fix is to call the
##   method once on a small problem before reading the baseline, so
##   that cost lands in the baseline instead. Do not do this. It was
##   tried and removed.
##
##   Calling anything before the baseline leaves R's garbage collector
##   with higher heap growth targets. gc(reset = TRUE, full = TRUE)
##   resets the max-used counters but NOT those targets, so the
##   measured run then collects less often and reaches a higher RSS
##   high-water mark -- and the increase in peak exceeds the increase
##   in baseline. Measured on Pitprops, warming up moved mixOmics from
##   0.6 to 11.6 MB, nsprcomp from 44.7 to 114.0 MB and nscumcomp from
##   72.4 to 275.8 MB, while msPCA (Sigma) went DOWN, 0.8 to 0.1 MB.
##
##   The inflation hits methods that allocate on R's heap and spares
##   those that do not, msPCA working in Eigen and sklearn in CPython.
##   It therefore biases the comparison in msPCA's favour, which is
##   exactly what this harness exists to prevent.
##
##   Carrying ~4 MB of loading cost in every row is the lesser evil,
##   and it is small next to the figures that carry the comparison
##   (272 MB for elasticnet at p = 500, 641 MB for msPCA on Sigma at
##   p = 4088). A single-repetition check on mtcars confirmed that the
##   ~21 MB elasticnet reports at p = 11 is genuine allocation, not
##   loading: warming up removed only about 4 MB of it.
##
## ---------------------------------------------------------
## Repetitions and seeds
## ---------------------------------------------------------
## Each method is run `reps` times, each in its own process, with a
## DIFFERENT seed per repetition (43, 44, 45, ... by default).
## set.seed() is called in the child immediately before the call, so
## benchmarked functions must NOT set their own seed. A function that
## needs the seed as an explicit argument -- scikit-learn's
## `random_state`, which does not read R's RNG -- declares a formal
## named `seed` and receives it.
##
## Consequences:
##   * runtime is the MEDIAN over repetitions;
##   * memory is the MEDIAN over repetitions, with min/max reported,
##     because peak RSS is not reproducible run to run (see above);
##   * FVE and orthogonality are medians, with min/max so that
##     solution stability is visible. Deterministic methods
##     (elasticnet, PMA, sparsepca, mixOmics, amanpg, prcomp) collapse
##     to a zero-width range, which doubles as a check that the
##     harness is behaving; nsprcomp and nscumcomp use random
##     initialisation and will not.
##   * nnz is taken from the representative repetition, defined as
##     the one whose FVE is closest to the median FVE.
##   * n_fail counts repetitions in which the method returned NULL
##     (nscumcomp's "Co-linear principal axes" failure).
##
## ---------------------------------------------------------
## Column glossary (written to each results CSV)
## ---------------------------------------------------------
##   runtime_s       -- median elapsed seconds over reps
##   runtime_min_s / runtime_max_s
##                   -- range of elapsed seconds over reps. Wall-clock
##                      timing is subject to machine conditions as well
##                      as to seed-dependent iteration counts, so quote
##                      the range whenever it is wide rather than
##                      presenting the median to three digits.
##   working_set_mb  -- input_mb + mem_delta_mb, median over reps; the
##                      headline memory figure. Total memory
##                      attributable to the method, excluding the
##                      interpreter and dependency-loading floor.
##   working_set_min_mb / working_set_max_mb
##                   -- range of that quantity over repetitions. Quote
##                      it whenever it is wide; a wide range means the
##                      point estimate should not be leaned on.
##   input_mb        -- object.size() of the matrix handed to the
##                      method, so the Sigma (p x p) versus X (n x p)
##                      difference is visible rather than hidden.
##   mem_delta_mb    -- median peak minus baseline: the solver's own
##                      working memory, on top of its inputs. Includes
##                      each package's first-call loading cost
##                      (roughly 4 MB); see "Why there is no warm-up".
##   peak_rss_mb     -- median absolute high-water RSS.
##   baseline_rss_mb -- median RSS after packages, setup() and inputs,
##                      but before the method is called.
##   gc_max_mb       -- R-level gc() maximum, median. Cross-check
##                      ONLY; not comparable across methods, for the
##                      reason given at the top of this file.
##
## Requires: callr, Rcpp (both already needed to build msPCA).
## Tested on macOS (Apple M2) and Linux. Not supported on Windows;
## there the memory columns come back NA and runtimes are unaffected.
############################################################

if (!requireNamespace("callr", quietly = TRUE))
  stop("bench_utils.R requires the 'callr' package. install.packages('callr')")

## Thread pinning, mirroring run_all.R, re-applied inside each child
## because callr::r() starts a clean process.
.BENCH_THREAD_ENV <- list(
  OPENBLAS_NUM_THREADS   = "1",
  OMP_NUM_THREADS        = "1",
  MKL_NUM_THREADS        = "1",
  BLIS_NUM_THREADS       = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

## Base seed; repetition i uses .BENCH_SEED0 + i - 1.
.BENCH_SEED0 <- 43L

## Store for the values chosen by tune_parameter(); emptied whenever this file
## is sourced, which each notebook does exactly once.
.BENCH_TUNING <- new.env(parent = emptyenv())

## ---------------------------------------------------------------
## Compiled peak-RSS probe.
##
## Compiled ONCE in the parent into a persistent cacheDir; each child
## then calls sourceCpp() on the same file with the same cacheDir and
## Rcpp reuses the cached shared object (dyn.load only, no compiler).
## Without this, every child would pay a full C++ compile -- and there
## are now reps x methods children per notebook.
## ---------------------------------------------------------------
.BENCH_CACHE_DIR <- file.path(tempdir(), "mspca_bench_cache")
.BENCH_SRC_FILE  <- file.path(.BENCH_CACHE_DIR, "peak_rss.cpp")

.bench_init_probe <- function() {
  dir.create(.BENCH_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(.BENCH_SRC_FILE)) {
    writeLines(c(
      '#include <Rcpp.h>',
      '#include <sys/resource.h>',
      '',
      '// [[Rcpp::export]]',
      'double peak_rss_bytes() {',
      '  struct rusage ru;',
      '  getrusage(RUSAGE_SELF, &ru);',
      '#ifdef __APPLE__',
      '  return (double) ru.ru_maxrss;          /* macOS reports bytes     */',
      '#else',
      '  return (double) ru.ru_maxrss * 1024.0; /* Linux reports kilobytes */',
      '#endif',
      '}'
    ), .BENCH_SRC_FILE)
  }
  tryCatch({
    Rcpp::sourceCpp(.BENCH_SRC_FILE, cacheDir = .BENCH_CACHE_DIR, env = new.env())
    TRUE
  }, error = function(e) {
    warning("Could not compile the peak-RSS probe; memory columns will be NA. ",
            conditionMessage(e), call. = FALSE)
    FALSE
  })
}

.BENCH_MEM_OK <- .bench_init_probe()


## ---------------------------------------------------------------
## The body of one repetition, executed in its own R process.
## Must be self-contained: callr does not ship the caller's globals.
## ---------------------------------------------------------------
.bench_one_run <- function(fun, inputs, globals, setup, packages, seed,
                           thr, srcfile, cachedir, mem_ok) {

  do.call(Sys.setenv, thr)
  for (p in packages) library(p, character.only = TRUE)

  ## Helpers visible to fun() and setup(), but not charged as inputs.
  env <- list2env(globals, parent = globalenv())
  environment(fun) <- env

  if (!is.null(setup)) {
    environment(setup) <- env
    sv <- setup()
    if (!is.null(sv)) inputs$setup <- sv
  }

  peak_rss_bytes <- function() NA_real_
  if (mem_ok) {
    probe <- new.env()
    try(Rcpp::sourceCpp(srcfile, cacheDir = cachedir, env = probe), silent = TRUE)
    if (!is.null(probe$peak_rss_bytes)) peak_rss_bytes <- probe$peak_rss_bytes
  }

  ## Input footprint excludes the setup handle (a module pointer).
  charged <- inputs[setdiff(names(inputs), "setup")]
  input_bytes <- sum(vapply(charged, function(z) as.numeric(object.size(z)), 0))

  if ("seed" %in% names(formals(fun))) inputs$seed <- seed

  gc(reset = TRUE, full = TRUE)
  baseline <- peak_rss_bytes()

  set.seed(seed)
  tt  <- system.time(res <- do.call(fun, inputs))

  peak  <- peak_rss_bytes()
  gcmax <- sum(gc()[, "max used"] * c(56, 8))  # Ncells, Vcells -> bytes

  list(time = tt[["elapsed"]], baseline = baseline, peak = peak,
       input_bytes = input_bytes, gc_max = gcmax,
       result = if (is.null(res)) NULL else as.matrix(res))
}


#' Benchmark one method, one repetition per isolated R process.
#'
#' @param fun      Function containing the call to benchmark. It receives the
#'                 elements of `inputs` as named arguments and MUST return a
#'                 plain p x r loadings matrix (already unit-normalised where
#'                 the package does not do so), or NULL if the method failed.
#'                 Returning the full model object would serialise it back
#'                 across the process boundary, which is wasteful and
#'                 impossible for Python objects.
#'                 `fun` must be self-contained: callr does not ship the
#'                 caller's global environment, so anything it needs arrives
#'                 via `inputs` or `globals`.
#'                 `fun` must NOT call set.seed(); the harness seeds each
#'                 repetition. If it declares a formal named `seed`, the
#'                 repetition's seed is passed in (for solvers such as
#'                 sklearn that take an explicit random_state).
#' @param inputs   Named list of data passed to `fun`. Materialised in the child
#'                 BEFORE the baseline reading and counted in `input_mb`.
#' @param globals  Named list of helper objects (e.g. unit_norm) injected into
#'                 the child. NOT counted in `input_mb`.
#' @param setup    Optional zero-argument function run in the child after the
#'                 packages are attached and before the baseline reading. Use it
#'                 for anything whose cost should be excluded from the method's
#'                 own memory delta -- notably initialising reticulate and
#'                 importing sklearn, which costs ~150 MB of NumPy/SciPy that
#'                 has nothing to do with the sparse-PCA solver. If it returns
#'                 a non-NULL value, that value is passed to `fun` as the
#'                 argument `setup`.
#' @param packages Character vector of packages to attach in the child.
#' @param reps     Number of repetitions; each runs in a separate process.
#' @param label    Method name as it should appear in the results table.
#'
#' @return list with runtime_s, runtime_all, seeds, the median memory figures
#'         and their ranges, `results` (per-repetition loadings), `result`
#'         (repetition 1, for convenience) and `error`.
bench_method <- function(fun, inputs = list(), globals = list(), setup = NULL,
                         packages = character(), reps = 3L,
                         label = "method") {

  reps  <- as.integer(reps)
  seeds <- .BENCH_SEED0 + seq_len(reps) - 1L
  cat("  [bench]", label, "")

  runs <- vector("list", reps)
  errs <- character(0)
  for (i in seq_len(reps)) {
    runs[[i]] <- tryCatch(
      callr::r(.bench_one_run,
               args = list(fun = fun, inputs = inputs, globals = globals,
                           setup = setup, packages = packages, seed = seeds[i],
                           thr = .BENCH_THREAD_ENV, srcfile = .BENCH_SRC_FILE,
                           cachedir = .BENCH_CACHE_DIR, mem_ok = .BENCH_MEM_OK),
               show = FALSE),
      error = function(e) { errs <<- c(errs, conditionMessage(e)); NULL })
    cat(if (is.null(runs[[i]])) "x" else ".")
  }

  ok <- !vapply(runs, is.null, logical(1))
  if (!any(ok)) {
    cat(" FAILED\n")
    warning(sprintf("%s failed in all %d repetitions: %s", label, reps, errs[1]),
            call. = FALSE)
    return(list(label = label, runtime_s = NA_real_, runtime_min_s = NA_real_,
                runtime_max_s = NA_real_, runtime_all = NA_real_,
                seeds = seeds, working_set_mb = NA_real_,
                working_set_min_mb = NA_real_, working_set_max_mb = NA_real_,
                input_mb = NA_real_, mem_delta_mb = NA_real_,
                peak_rss_mb = NA_real_, baseline_rss_mb = NA_real_,
                gc_max_mb = NA_real_, results = vector("list", reps),
                result = NULL, error = errs[1]))
  }
  cat(" done\n")

  good  <- runs[ok]
  pull  <- function(f) vapply(good, function(z) as.numeric(z[[f]]), 0)
  mb    <- function(x) x / 1024^2
  med   <- function(x) round(stats::median(x), 1)

  times    <- pull("time")
  deltas   <- mb(pull("peak") - pull("baseline"))
  input    <- mb(pull("input_bytes"))[1]
  working  <- input + deltas

  ## Per-repetition loadings, NULL for repetitions whose process died.
  results <- vector("list", reps)
  results[ok] <- lapply(good, `[[`, "result")

  list(label              = label,
       runtime_s          = round(stats::median(times), 3),
       runtime_min_s      = round(min(times), 3),
       runtime_max_s      = round(max(times), 3),
       runtime_all        = round(times, 3),
       seeds              = seeds,
       working_set_mb     = med(working),
       working_set_min_mb = round(min(working), 1),
       working_set_max_mb = round(max(working), 1),
       input_mb           = round(input, 1),
       mem_delta_mb       = med(deltas),
       peak_rss_mb        = med(mb(pull("peak"))),
       baseline_rss_mb    = med(mb(pull("baseline"))),
       gc_max_mb          = med(mb(pull("gc_max"))),
       results            = results,
       result             = results[[which(ok)[1]]],
       error              = if (length(errs)) errs[1] else NA_character_)
}


#' Choose a tuning parameter automatically, by sparsity.
#'
#' `msPCA`, `elasticnet`, `mixOmics` and `nsprcomp` accept an exact cardinality
#' per component and need no tuning. `PMA`, `sparsepca`, `amanpg` and
#' `scikit-learn` expose only a penalty magnitude or an l1 bound, so the value
#' that comes closest to the target sparsity has to be searched for. Doing that
#' by eye and pasting a constant into the script is how those constants go stale
#' -- correcting the input to PMA::SPC() silently invalidated a `sumabsv` that
#' had been tuned against the old input, and the mismatch survived a full re-run.
#' This function removes the manual step: the sweep and the selection happen in
#' the same place, every time the notebook runs.
#'
#' Selection minimises the total absolute deviation of the realised
#' per-component cardinality from `target`. Ties are broken on FVE when `C` is
#' supplied, so that among equally-sparse candidates the better solution wins.
#'
#' @param fit    Function of one grid value returning a p x r loadings matrix,
#'               or NULL / an error if that value is infeasible.
#' @param grid   Values to try. A numeric vector, or a list when the parameter
#'               is itself a vector; `grid[[i]]` is passed to `fit`.
#' @param target Integer vector of length r: the desired nonzeros per component.
#' @param C      Optional matrix for the FVE tie-break and for reporting.
#' @param label  Name used in the console output.
#' @param quiet  Suppress the per-value trace, keeping only the selection line.
#'
#' @return The selected element of `grid`.
tune_parameter <- function(fit, grid, target, C = NULL, label = "parameter",
                           quiet = FALSE) {
  n    <- length(grid)
  dev  <- rep(NA_real_, n)
  fve  <- rep(NA_real_, n)
  live <- rep(FALSE, n)          # did the candidate return r NON-EMPTY components?
  nnzs <- vector("list", n)

  cat("  tuning ", label, " over ", n, " values...\n", sep = "")
  for (i in seq_len(n)) {
    L <- tryCatch(fit(grid[[i]]), error = function(e) NULL)
    if (is.null(L) || !is.matrix(L)) next
    nnz     <- colSums(abs(L) > 0)
    nnzs[[i]] <- nnz
    dev[i]  <- sum(abs(nnz - target))
    live[i] <- all(nnz > 0)
    if (!is.null(C)) fve[i] <- fraction_variance_explained(C, L)
    if (!quiet)
      cat("    ", format(grid[[i]]), " | NNZ: ", paste(nnz, collapse = " "),
          " | deviation: ", dev[i],
          if (!is.null(C)) paste0(" | FVE: ", round(fve[i], 4)) else "",
          if (!live[i]) "  [rejected: empty component]" else "", "\n", sep = "")
  }

  if (all(is.na(dev)))
    stop("tune_parameter(): every value in the grid for ", label, " failed.",
         call. = FALSE)

  ## A penalty large enough to empty a component has not produced a sparse
  ## solution -- it has produced FEWER THAN r COMPONENTS, which is a different
  ## and worse thing. Scoring by |nnz - target| alone rewards exactly that: an
  ## empty component costs only k, while destroying the variance that component
  ## would have explained. On mtcars this let sparsepca select (4,4,0) over
  ## (7,4,1), dropping FVE from 0.585 to 0.427, and amanpg select (5,0,0) over
  ## (7,4,2), dropping it from 0.649 to 0.311. Such candidates are excluded.
  usable <- which(live)
  if (length(usable) == 0L) {
    warning("Every candidate for ", label, " left at least one component empty; ",
            "falling back to the smallest deviation overall. The method cannot ",
            "return ", length(target), " non-empty components anywhere on this ",
            "grid -- widen it, or report the method as unable to meet the ",
            "budget.", call. = FALSE)
    usable <- which(!is.na(dev))
  }

  best <- usable[dev[usable] == min(dev[usable], na.rm = TRUE)]
  if (length(best) > 1 && !all(is.na(fve[best])))
    best <- best[which.max(fve[best])]
  best <- best[1]

  n_rejected <- sum(!live & !is.na(dev))
  cat("  -> ", label, " = ", format(grid[[best]]),
      " | NNZ: ", paste(nnzs[[best]], collapse = " "),
      " | deviation ", dev[best], " from target ", paste(target, collapse = " "),
      if (!is.na(fve[best])) paste0(" | FVE ", round(fve[best], 4)) else "",
      if (n_rejected) paste0("  (", n_rejected, " of ", n,
                             " candidates rejected for empty components)") else "",
      "\n", sep = "")

  ## A selection sitting on the edge of the grid usually means the grid is too
  ## narrow and the real optimum lies outside it. Say so rather than silently
  ## returning a boundary value.
  if (dev[best] > 0 && (best == min(usable) || best == max(usable)))
    warning(label, " was selected at the ",
            if (best == min(usable)) "lower" else "upper",
            " end of its usable grid with a non-zero deviation (", dev[best],
            "). The grid may be too narrow -- widen it and re-run.",
            call. = FALSE)

  ## Record the choice so tuning_table() can write it out; this is what makes
  ## the article's "Tuning:" lines readable off a file instead of transcribed.
  assign(label, list(parameter = label, value = grid[[best]],
                     nnz = nnzs[[best]], deviation = dev[best],
                     fve = fve[best], violation = NA_real_,
                     at_grid_edge = FALSE, n_rejected = n_rejected),
         envir = .BENCH_TUNING)

  grid[[best]]
}


#' Tune a VECTOR-valued penalty, one coordinate at a time.
#'
#' For methods whose penalty is a vector with one entry per component --
#' `amanpg::spca.amanpg()`'s `lambda1` is the only one here -- a single shared
#' value cannot deliver balanced sparsity. With `lambda2 = Inf` the amanpg
#' formulation soft-thresholds the columns of Sigma %*% A at lambda1[j], and
#' those columns scale with the eigenvalue of component j. A threshold that
#' leaves k nonzeros on the leading component therefore annihilates the trailing
#' ones as the spectrum decays: tuning a shared value produced (5,0,0) on mtcars
#' and (8,4,2,0,0,0) on Pitprops.
#'
#' Since the package exposes per-component penalties, we tune them. Methods
#' whose interface offers per-component control -- elasticnet, mixOmics,
#' nsprcomp, msPCA -- all receive it, so this is parity rather than preferential
#' treatment.
#'
#' The search is coordinate descent: start from the best shared value, then
#' sweep each coordinate in turn against the FULL objective (total deviation
#' across all components), since changing one entry moves the whole solution.
#' Two passes are enough in practice; the second rarely moves anything.
#'
#' @inheritParams tune_parameter
#' @param passes Number of coordinate-descent sweeps.
#' @return Numeric vector of length `length(target)`.
tune_parameter_vector <- function(fit, grid, target, C = NULL,
                                  label = "parameter", passes = 2L) {
  r <- length(target)

  score <- function(vec) {
    L <- tryCatch(fit(vec), error = function(e) NULL)
    if (is.null(L) || !is.matrix(L)) return(list(dev = Inf, nnz = NULL, fve = NA_real_))
    nnz <- colSums(abs(L) > 0)
    list(dev = if (any(nnz == 0)) Inf else sum(abs(nnz - target)),
         nnz = nnz,
         fve = if (!is.null(C)) fraction_variance_explained(C, L) else NA_real_)
  }

  cat("  tuning ", label, " (vector of ", r, ") by coordinate descent\n", sep = "")

  ## Start from the best shared value, so the search begins somewhere sensible.
  shared <- lapply(grid, function(g) score(rep(g, r)))
  devs   <- vapply(shared, function(z) z$dev, 0)
  cur    <- rep(grid[[which.min(devs)]], r)
  best   <- score(cur)
  if (!is.finite(best$dev)) {
    ## No shared value keeps every component alive; start from the smallest
    ## penalty in the grid, which is the least likely to empty anything.
    cur  <- rep(min(unlist(grid)), r)
    best <- score(cur)
  }
  cat("    start: ", paste(format(cur), collapse = " "),
      " | NNZ: ", paste(best$nnz, collapse = " "),
      " | deviation: ", best$dev, "\n", sep = "")

  for (pass in seq_len(passes)) {
    moved <- FALSE
    for (j in seq_len(r)) {
      for (g in grid) {
        cand <- cur; cand[j] <- g
        s <- score(cand)
        better <- s$dev < best$dev ||
          (s$dev == best$dev && !is.na(s$fve) && !is.na(best$fve) && s$fve > best$fve)
        if (better) { cur <- cand; best <- s; moved <- TRUE }
      }
    }
    cat("    pass ", pass, ": ", paste(format(cur), collapse = " "),
        " | NNZ: ", paste(best$nnz, collapse = " "),
        " | deviation: ", best$dev,
        if (!is.na(best$fve)) paste0(" | FVE: ", round(best$fve, 4)) else "", "\n", sep = "")
    if (!moved) break
  }

  if (!is.finite(best$dev))
    warning(label, ": no combination on this grid keeps all ", r,
            " components non-empty. The method cannot meet the budget here.",
            call. = FALSE)

  cat("  -> ", label, " = ", paste(format(cur), collapse = " "),
      " | NNZ: ", paste(best$nnz, collapse = " "),
      " | deviation ", best$dev, " from target ", paste(target, collapse = " "),
      if (!is.na(best$fve)) paste0(" | FVE ", round(best$fve, 4)) else "", "\n", sep = "")

  assign(label, list(parameter = label, value = cur, nnz = best$nnz,
                     deviation = best$dev, fve = best$fve, violation = NA_real_,
                     at_grid_edge = FALSE, n_rejected = NA_integer_),
         envir = .BENCH_TUNING)
  cur
}


#' Tune a parameter that trades variance explained against FEASIBILITY.
#'
#' `nsprcomp::nscumcomp()` is the only competing function with a knob on
#' non-redundancy: `gamma` penalises divergence from orthonormality of the
#' loadings. Its default is 0, i.e. no penalty at all -- the worst possible
#' setting for the orthogonality column, and leaving it there would mean
#' reporting that `nscumcomp` is infeasible while never having asked it to be
#' feasible. That is not a comparison worth publishing.
#'
#' Sparsity is not at stake here: `nscumcomp` takes its cardinality budget
#' separately and `gamma` does not change it. What `gamma` buys is feasibility,
#' at a cost in FVE. The selection rule is therefore: among values that bring
#' the orthogonality violation within `tol` -- the same tolerance `msPCA`
#' enforces by default -- keep the one with the highest FVE. This asks "at
#' equal feasibility, how much variance does each method explain?", which is
#' the like-for-like question. If no value on the grid reaches `tol`, the one
#' with the smallest violation is kept and a warning is issued, since the
#' honest report is then that the method cannot reach the tolerance at all.
#'
#' @param fit   Function of one grid value returning a p x r loadings matrix.
#' @param grid  Values to try.
#' @param C     Matrix against which FVE and the violation are evaluated.
#' @param ctype Constraint type passed to feasibility_violation_off().
#' @param tol   Feasibility target; defaults to msPCA's own 1e-4.
#' @param label Name used in the console output and the tuning table.
#'
#' @return The selected element of `grid`.
tune_for_feasibility <- function(fit, grid, C, ctype = 0, tol = 1e-4,
                                 label = "parameter", quiet = FALSE) {
  n    <- length(grid)
  viol <- rep(NA_real_, n)
  fve  <- rep(NA_real_, n)
  nnzs <- vector("list", n)

  cat("  tuning ", label, " for feasibility (target <= ", format(tol),
      ") over ", n, " values...\n", sep = "")
  for (i in seq_len(n)) {
    L <- tryCatch(fit(grid[[i]]), error = function(e) NULL)
    if (is.null(L) || !is.matrix(L)) next
    nnzs[[i]] <- colSums(abs(L) > 0)
    viol[i]   <- feasibility_violation_off(C, L, ctype)
    fve[i]    <- fraction_variance_explained(C, L)
    if (!quiet)
      cat("    ", format(grid[[i]]), " | violation: ",
          format(viol[i], scientific = TRUE, digits = 3),
          " | FVE: ", round(fve[i], 4), "\n", sep = "")
  }

  if (all(is.na(viol)))
    stop("tune_for_feasibility(): every value in the grid for ", label,
         " failed.", call. = FALSE)

  feasible <- which(viol <= tol)
  if (length(feasible)) {
    best <- feasible[which.max(fve[feasible])]
  } else {
    best <- which.min(viol)
    warning(label, ": no value on the grid brings the violation within ",
            format(tol), ". The smallest achievable is ",
            format(viol[best], scientific = TRUE, digits = 3),
            " -- report the method as unable to reach the tolerance rather ",
            "than as merely untuned.", call. = FALSE)
  }

  cat("  -> ", label, " = ", format(grid[[best]]),
      " | violation ", format(viol[best], scientific = TRUE, digits = 3),
      " | FVE ", round(fve[best], 4),
      if (!length(feasible)) "  (tolerance NOT reached)" else "", "\n", sep = "")

  assign(label, list(parameter = label, value = grid[[best]],
                     nnz = nnzs[[best]], deviation = NA_real_,
                     fve = fve[best], violation = viol[best],
                     at_grid_edge = (best == 1L || best == n) && !length(feasible),
                     n_rejected = NA_integer_),
         envir = .BENCH_TUNING)

  grid[[best]]
}


#' Selected tuning values for this notebook, as a data frame.
#'
#' Populated by tune_parameter(); reset each time bench_utils.R is sourced, so
#' one notebook's selections cannot leak into another's.
tuning_table <- function() {
  keys <- ls(.BENCH_TUNING)
  if (!length(keys)) return(NULL)
  rows <- mget(keys, envir = .BENCH_TUNING)
  do.call(rbind, lapply(rows, function(z) data.frame(
    parameter    = z$parameter,
    value        = paste(format(z$value), collapse = " "),
    selected_nnz = paste(z$nnz, collapse = " "),
    deviation    = z$deviation,
    fve          = round(z$fve, 4),
    violation    = z$violation,
    at_grid_edge = z$at_grid_edge,
    n_rejected   = z$n_rejected,
    stringsAsFactors = FALSE)))
}


#' Summarise solution quality across the repetitions of one benchmark.
#'
#' @return list(nnz, fve, fve_min, fve_max, orth, orth_max, n_fail, rep_used)
#'         The representative repetition (`rep_used`) is the one whose FVE is
#'         closest to the median FVE; nnz is read from it.
bench_quality <- function(b, S, r) {
  ok <- !vapply(b$results, is.null, logical(1))
  n_fail <- sum(!ok)

  if (!any(ok)) {
    return(list(nnz = rep(NA_real_, r), fve = NA_real_, fve_min = NA_real_,
                fve_max = NA_real_, orth = NA_real_, orth_max = NA_real_,
                n_fail = n_fail, rep_used = NA_integer_))
  }

  idx   <- which(ok)
  fves  <- vapply(b$results[idx], function(L) fraction_variance_explained(S, L), 0)
  orths <- vapply(b$results[idx], function(L) feasibility_violation_off(S, L, 0), 0)

  ## Representative repetition: the one closest to the median FVE. With an even
  ## number of repetitions the median falls between two runs, so "closest"
  ## rather than "equal to" keeps nnz sourced from an actual run.
  med  <- stats::median(fves)
  pick <- idx[which.min(abs(fves - med))]

  list(nnz      = colSums(abs(b$results[[pick]]) > 0)[seq_len(r)],
       fve      = med,
       fve_min  = min(fves),
       fve_max  = max(fves),
       orth     = stats::median(orths),
       orth_max = max(orths),
       n_fail   = n_fail,
       rep_used = pick)
}


#' One-line console report, matching the cat() style of the notebooks.
report_bench <- function(b, S) {
  r <- if (is.null(b$result)) 0L else ncol(b$result)
  if (r == 0L) { cat(b$label, "| FAILED in all repetitions\n"); return(invisible(b)) }
  q <- bench_quality(b, S, r)
  cat(b$label,
      "| NNZ:", q$nnz,
      "| FVE:", round(q$fve, 4),
      if (q$fve_max - q$fve_min > 0)
        paste0("[", round(q$fve_min, 4), ", ", round(q$fve_max, 4), "]") else "(stable)",
      "| Orth:", format(q$orth, scientific = TRUE, digits = 3),
      "| Time:", b$runtime_s, "s",
      paste0("[", b$runtime_min_s, ", ", b$runtime_max_s, "]"),
      "| Working set:", b$working_set_mb, "MB",
      paste0("[", b$working_set_min_mb, ", ", b$working_set_max_mb, "]"),
      if (q$n_fail > 0) paste0("| FAILED ", q$n_fail, "/", length(b$results)) else "",
      "\n")
  invisible(b)
}


#' Assemble timing/memory columns for the results data frame.
#'
#' working_set_mb = input_mb + mem_delta_mb is the headline memory figure:
#' total memory attributable to the method, with the R interpreter and the
#' package/Python loading floor excluded. It is the column that carries the
#' Sigma (p x p) versus X (n x p) comparison, because the input matrix is
#' where most of that difference lives. Do not quote input_mb and
#' mem_delta_mb separately as well -- that double-counts. Always check
#' working_set_min_mb / working_set_max_mb before leaning on the median:
#' peak RSS is not reproducible run to run for methods that allocate heavily
#' on R's heap.
bench_table <- function(benches) {
  do.call(rbind, lapply(benches, function(b) data.frame(
    method             = b$label,
    runtime_s          = b$runtime_s,
    runtime_min_s      = b$runtime_min_s,
    runtime_max_s      = b$runtime_max_s,
    working_set_mb     = b$working_set_mb,
    working_set_min_mb = b$working_set_min_mb,
    working_set_max_mb = b$working_set_max_mb,
    input_mb           = b$input_mb,
    mem_delta_mb       = b$mem_delta_mb,
    peak_rss_mb        = b$peak_rss_mb,
    baseline_rss_mb    = b$baseline_rss_mb,
    gc_max_mb          = b$gc_max_mb,
    stringsAsFactors = FALSE)))
}


#' Sparsity / quality columns, summarised across repetitions and NA-safe.
#'
#' fve and orth are medians over repetitions; fve_min/fve_max and orth_max give
#' the spread. A zero-width FVE range means the method is deterministic under
#' reseeding. n_fail counts repetitions that returned NULL.
#'
#' @param benches list of bench_method() results
#' @param S       matrix against which FVE and orthogonality are evaluated
#'                (note: for Pitprops this is the published correlation
#'                matrix, not the covariance of the pseudo-data)
#' @param r       number of components
quality_table <- function(benches, S, r) {
  qs <- lapply(benches, bench_quality, S = S, r = r)

  nnz <- t(vapply(qs, function(q) as.numeric(q$nnz), numeric(r)))
  colnames(nnz) <- paste0("nnz_pc", seq_len(r))

  data.frame(
    nnz,
    fve      = round(vapply(qs, function(q) q$fve, 0), 3),
    fve_min  = round(vapply(qs, function(q) q$fve_min, 0), 3),
    fve_max  = round(vapply(qs, function(q) q$fve_max, 0), 3),
    orth     = vapply(qs, function(q) q$orth, 0),
    orth_max = vapply(qs, function(q) q$orth_max, 0),
    n_fail   = vapply(qs, function(q) q$n_fail, 0L)
  )
}
