# Contributing to kamila

Thank you for your interest in contributing to **`kamila`**! This document provides development guidelines, testing standards, and operational instructions for all contributors.

---

## 1. Development Workflow

1. **Create a Feature Branch**: Never commit directly to the `master` branch. Always branch off `master` with a descriptive name (e.g., `feature/new-metric`, `fix/edge-case-nan`).
2. **Make Changes**:
   - Write clean, documented R code following tidyverse/CRAN style.
   - For C++ changes via Rcpp, ensure memory safety and run `Rcpp::compileAttributes()` to regenerate bindings.
   - Document all exported functions using `roxygen2` comments (`#'`) in `R/*.R`. Do **not** manually edit `.Rd` files in `man/`.
3. **Log Changes to News**:
   - Every enhancement, bug fix, breaking change, or dependency update **must** be logged in [`NEWS.md`](NEWS.md) following CRAN package news standards.
4. **Run Pre-Commit Verification**: Ensure all local checks pass before committing (see below).
5. **Submit a Pull Request**: Provide a clear description of the changes, referencing any relevant issue numbers.

---

## 2. Pre-Commit Verification & Quality Standards

Before committing or submitting a pull request, you **must** execute and pass the following checks:

### 1. Unit Tests
All unit tests must pass with 0 errors and 0 failures:
```r
devtools::test()
```

### 2. Code Coverage (Strict 100% Policy)
Code coverage **must not decrease below 100%** from any commit:
```r
covr::package_coverage()
```
- Every new function, control flow branch, edge case, and error condition must have corresponding unit tests in `tests/testthat/`.

### 3. Linting & Code Style
The codebase strictly adheres to linter standards. There must be 0 linter warnings or errors:
```r
lintr::lint_package()
```

### 4. Full CRAN Package Check
Verify that the package builds cleanly with 0 errors, 0 warnings, and 0 notes:
```r
devtools::check(cran = TRUE)
```

### 5. Multi-Platform & Sanitizer Checks (R-hub v2)
Prior to submitting releases or making low-level C++ changes, run pre-flight checks across multiple architectures and memory sanitizers (e.g., ASAN/UBSAN, Valgrind, Windows UCRT, Linux containers) using R-hub v2:
```r
# Verify setup and GitHub PAT
rhub::rhub_doctor()

# Run interactive or targeted checks
rhub::rhub_check()
```

---

## 3. Key Development Commands

| Action | Command |
| :--- | :--- |
| **Load Package** | `Rscript -e "devtools::load_all()"` |
| **Recompile C++ (Rcpp)** | `Rscript -e "Rcpp::compileAttributes(); devtools::clean_dll()"` |
| **Generate Documentation** | `Rscript -e "devtools::document()"` |
| **Run Unit Tests** | `Rscript -e "devtools::test()"` |
| **Check Code Coverage** | `Rscript -e "covr::package_coverage()"` |
| **Lint Codebase** | `Rscript -e "lintr::lint_package()"` |
| **Run Full Package Check** | `Rscript -e "devtools::check(cran = TRUE)"` |
| **Validate R-hub Setup** | `Rscript -e "rhub::rhub_doctor()"` |
| **Run R-hub v2 Checks** | `Rscript -e "rhub::rhub_check()"` |

---

## 4. Git Commit Guidelines

Commit messages should follow [Conventional Commits](https://www.conventionalcommits.org/):
- `feat:` New features or parameter additions
- `fix:` Bug fixes or numerical error corrections
- `docs:` Documentation updates (roxygen2 comments, README, NEWS.md)
- `test:` Adding or modifying unit tests
- `refactor:` Code refactoring without behavioral changes
- `ci:` CI/CD pipeline modifications
