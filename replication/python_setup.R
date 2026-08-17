## ============================================================
## python_setup.R
##
## Single place where the Python side of the replication is configured.
## Source this file BEFORE any Python is touched (see "Timing" below).
##
## WHY THIS FILE EXISTS
##
## Up to reticulate 1.40 the package searched the machine for an interpreter:
## RETICULATE_PYTHON, then `python3` on PATH, then any virtualenv/conda
## environment it could find. If that interpreter had scikit-learn installed,
## the comparison ran.
##
## From reticulate 1.41 onwards the default is the opposite. Unless an
## interpreter is named explicitly, reticulate provisions its OWN ephemeral
## environment through `uv` (on macOS under
## ~/Library/Caches/org.R-project.R/R/reticulate/uv/), containing nothing but
## what the session has declared via py_require(). A machine whose system
## python3 has scikit-learn will therefore still report
##
##   scikit-learn was not found for the Python interpreter at
##   '.../reticulate/uv/cache/archive-v0/.../bin/python3'
##
## because that managed interpreter -- not python3 -- is what reticulate
## imports from.
##
## Declaring the requirement here makes the comparison self-provisioning and
## pinned to the version quoted in the benchmarking article: uv downloads
## scikit-learn on first use and caches it, on any machine, with no manual pip
## step and no dependence on whatever Python the reader happens to have.
##
## TIMING
##
## py_require() must run before Python is initialised (the first import() or
## py_module_available() call). Afterwards packages can only be added, and the
## Python version can no longer be pinned. run_all.R sources this file first;
## each script that touches Python also sources it, so the scripts remain
## runnable one at a time.
##
## ESCAPE HATCHES (both honoured without editing this file)
##
##   RETICULATE_PYTHON=/path/to/python
##       Use that interpreter instead. It must already have scikit-learn.
##   RETICULATE_USE_MANAGED_VENV=no
##       Disable managed environments entirely and restore the pre-1.41
##       behaviour of discovering `python3` on PATH.
##
## A missing scikit-learn is NOT fatal anywhere in this replication: the
## sklearn rows are dropped and every other method still runs.
## ============================================================

SKLEARN_SPEC   <- "scikit-learn==1.6.1"  # version quoted in the article
PY_VERSION     <- "3.11"                 # Python version quoted in the article

## Declare the requirement for reticulate's managed environment. Skipped when
## the user has named their own interpreter, when managed environments are
## switched off, or when Python has already been initialised (in which case
## py_require() could no longer pin versions anyway).
.py_declare_requirements <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
  if (utils::packageVersion("reticulate") < "1.41.0") return(invisible(FALSE))
  if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) return(invisible(FALSE))
  if (identical(tolower(Sys.getenv("RETICULATE_USE_MANAGED_VENV")), "no"))
    return(invisible(FALSE))
  if (isTRUE(reticulate::py_available())) return(invisible(FALSE))
  ok <- tryCatch({
    reticulate::py_require(SKLEARN_SPEC, python_version = PY_VERSION)
    TRUE
  }, error = function(e) FALSE)
  invisible(ok)
}

.py_declare_requirements()

## Initialise Python and report what was found. Returns a list with
##   available : is scikit-learn importable?
##   python    : the interpreter actually in use (empty string if none)
##   version   : scikit-learn version, or NA
## The `python` field is the resolved interpreter from py_config(), not the
## RETICULATE_PYTHON environment variable, so it can be handed to the
## benchmark subprocesses and point them at the same environment.
sklearn_setup <- function() {
  out <- list(available = FALSE, python = "", version = NA_character_)
  if (!requireNamespace("reticulate", quietly = TRUE)) return(out)

  user_python <- Sys.getenv("RETICULATE_PYTHON")
  if (nzchar(user_python))
    try(reticulate::use_python(user_python, required = TRUE), silent = TRUE)

  out$available <- isTRUE(tryCatch(reticulate::py_module_available("sklearn"),
                                   error = function(e) FALSE))
  out$python <- tryCatch({
    p <- reticulate::py_config()$python
    if (is.null(p)) user_python else p
  }, error = function(e) user_python)

  if (out$available)
    out$version <- tryCatch(
      as.character(reticulate::import("sklearn")$`__version__`),
      error = function(e) NA_character_)
  out
}

## Explanatory text for the scripts to emit when scikit-learn is unavailable.
sklearn_missing_message <- function(py) {
  paste0(
    "scikit-learn was not found for the Python interpreter at '", py$python,
    "'. The sklearn.decomposition.SparsePCA comparison will be omitted from ",
    "this script's results; every other method still runs.\n",
    "  * reticulate >= 1.41 uses its own uv-managed environment rather than ",
    "the system python3, so installing scikit-learn into python3 will not ",
    "help. python_setup.R declares ", SKLEARN_SPEC, " via py_require(); if ",
    "this warning appeared, that declaration did not take effect -- most ",
    "often because Python had already been initialised earlier in the ",
    "session (restart R and source run_all.R from the start), or because uv ",
    "could not reach the network to download it.\n",
    "  * To use your own interpreter instead, set RETICULATE_PYTHON to one ",
    "that already has scikit-learn, or set RETICULATE_USE_MANAGED_VENV=no to ",
    "restore discovery of python3 on PATH.")
}
