compared_pkgs <- c(
  "msPCA", "elasticnet", "PMA", "sparsepca",
  "mixOmics", "nsprcomp"
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

# Save machine-readable tables
write.csv(pkg_versions, "package_versions.csv", row.names = FALSE)
write.csv(r_platform, "r_platform_info.csv", row.names = FALSE)

# Save complete session details (includes loaded package versions + locale + BLAS/LAPACK etc.)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")