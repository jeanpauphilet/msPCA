# Every package whose version the benchmarking article quotes must appear here,
# otherwise the article's Environment section cannot be checked against the
# recorded environment. amanpg and reticulate were previously quoted in the
# article but not recorded; callr drives the benchmarking harness.
compared_pkgs <- c(
  "msPCA", "elasticnet", "PMA", "sparsepca",
  "mixOmics", "nsprcomp", "amanpg", "reticulate", "callr"
)

# Add dataset/helper packages used in a given script if relevant:
# compared_pkgs <- c(compared_pkgs, "MASS", "hdi", "datasets")

pkg_versions <- data.frame(
  package = compared_pkgs,
  installed = vapply(compared_pkgs, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(
    compared_pkgs,
    function(p) {
      if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

r_platform <- data.frame(
  R_version = R.version$version.string,
  platform = R.version$platform,
  os = R.version$os,
  arch = R.version$arch,
  system = R.version$system,
  machine = Sys.info()[["machine"]],
  release = Sys.info()[["release"]],
  sysname = Sys.info()[["sysname"]],
  stringsAsFactors = FALSE
)

# scikit-learn is quoted in the article too, but lives on the Python side;
# record it here so the whole Environment section is verifiable from one file.
# python_setup.R declares the scikit-learn requirement for reticulate's
# managed environment and must be sourced before Python is initialised -- this
# script is the first thing run_all.R runs, so this is where Python first comes
# up in a full replication run.
if (file.exists("python_setup.R")) source("python_setup.R")

sklearn_version <- tryCatch({
  if (exists("sklearn_setup", mode = "function")) {
    sklearn_setup()$version
  } else if (requireNamespace("reticulate", quietly = TRUE)) {
    as.character(reticulate::import("sklearn")$`__version__`)
  } else NA_character_
}, error = function(e) NA_character_)

pkg_versions <- rbind(
  pkg_versions,
  data.frame(package = "scikit-learn (Python)",
             installed = !is.na(sklearn_version),
             version = sklearn_version, stringsAsFactors = FALSE)
)

# Save machine-readable tables
write.csv(pkg_versions, "package_versions.csv", row.names = FALSE)
write.csv(r_platform, "r_platform_info.csv", row.names = FALSE)

# Save complete session details (includes loaded package versions + locale + BLAS/LAPACK etc.)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")