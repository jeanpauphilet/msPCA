# Contributing to msPCA

Thank you for your interest in contributing to **msPCA**. This document covers
everything you need to get the development environment running and submit
high-quality changes.

---

## Table of contents

1. [Development setup](#development-setup)
2. [Running tests](#running-tests)
3. [Code style](#code-style)
4. [Documentation](#documentation)
5. [Submitting changes](#submitting-changes)
6. [Reporting issues](#reporting-issues)

---

## Development setup

### Prerequisites

- R ≥ 4.1 (tested on release and devel)
- A C++14-capable compiler (GCC ≥ 7, Clang ≥ 5, or MSVC 2017+)
- The following R packages:

```r
install.packages(c("Rcpp", "RcppEigen", "devtools", "testthat", "covr"))
```

### Clone and install

```bash
git clone https://github.com/jeanpauphilet/msPCA.git
cd msPCA
```

```r
# From the R console, inside the package directory:
devtools::install(dependencies = TRUE)
```

### Compile C++ code after edits

```r
devtools::load_all()   # recompiles src/ and sources R/ in one step
```

---

## Running tests

```r
devtools::test()
# or equivalently:
testthat::test_local()
```

All tests live in `tests/testthat/`. The test files are:

| File | What it covers |
|------|---------------|
| `test-mspca.R` | `mspca()` output structure, Sigma and X paths, real-data smoke tests |
| `test-tpm.R` | `tpm()` output structure, both input paths |
| `test-validators.R` | Error paths for invalid inputs (non-square, non-PSD, non-finite, …) |
| `test-s3methods.R` | `print.mspca`, `summary.mspca`, and helper functions |

### Coverage report

```r
covr::report()   # opens an HTML coverage report in the browser
```

The project targets ≥ 90% line coverage.

---

## Code style

### R

- Follow [tidyverse style](https://style.tidyverse.org/) for R code.
- Exported functions and helpers must have roxygen2 documentation blocks.
- All roxygen directives use Markdown (enabled via `Roxygen: list(markdown = TRUE)` in `DESCRIPTION`).
- Run `devtools::document()` after any change to roxygen comments.

### C++

- Standard: C++14.
- Match the existing brace and spacing style in `src/`.
- Every source file (`*.cpp`, `*.h`) must carry the SPDX header:

```cpp
// SPDX-License-Identifier: MIT
// Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet
```

### R source files

```r
# SPDX-License-Identifier: MIT
# Copyright (C) 2025 Ryan Cory-Wright, Jean Pauphilet
```

> Auto-generated files (`R/RcppExports.R`, `src/RcppExports.cpp`) are excluded
> from the header requirement — do not edit them by hand.

---

## Documentation

Rebuild Rd files after changing roxygen comments:

```r
devtools::document()
```

Rebuild the pkgdown site (requires `pkgdown`):

```r
pkgdown::build_site()
```

The `docs/` folder is committed to `main` and served as GitHub Pages — do not
delete it.

---

## Submitting changes

1. Fork the repository and create a branch from `main`:
   ```bash
   git checkout -b fix/my-bugfix
   ```
2. Make your changes and add or update tests as appropriate.
3. Ensure `devtools::check()` passes with no new errors, warnings, or notes
   (the CI matrix runs `R CMD check --as-cran` on Linux, macOS, and Windows).
4. Push your branch and open a pull request against `main`.
5. Fill in the pull request description, referencing any related issues.

---

## Reporting issues

Please use the GitHub issue templates:

- **Bug report** — for reproducible errors or unexpected behaviour.
- **Feature request** — for suggestions or ideas.

For security vulnerabilities, contact the maintainer directly by e-mail rather
than opening a public issue.
