# Changelog

All notable changes to the **kamila** R package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.3] - 2026-09-07

### Added
- **Standardized Multi-Start Benchmark Runs**: Standardized all compared clustering methods in the interactive Shiny benchmark dashboard (`inst/shiny/horserace/app.R`) to use 5 random initialization starts (`numInit = 5`, `nstart = 5`, `nrep = 5`, `nbKeep = 5`) for fair computational and performance comparisons.
- **Fast C++ Prediction Strength**: Re-implemented prediction strength evaluation in C++ (`calcPsCpp` via Rcpp) with $O(N + K^2)$ algorithmic complexity, replacing previous $O(N^2)$ pairwise operations (#12, #42).
- **Parallel Prediction Strength**: Added multi-core parallel processing support for prediction strength cross-validation runs using `numCores` parameter in `kamila()` (#46).
- **Non-Mixed Data Validation**: Added explicit error validation and informative messaging when non-mixed data (continuous-only or categorical-only) are passed to mixed-data clustering functions (#19).
- **Prediction Strength Documentation**: Comprehensive vignettes, usage examples, and documentation on interpreting prediction strength for cluster number selection (#15, #40).
- **Modern CI/CD Pipelines**: Replaced AppVeyor with GitHub Actions workflows for multi-OS R CMD check, Codecov coverage tracking, and automated linting.
- **100% Test Coverage & Regression Tests**: Added unit tests and snapshot regression tests achieving 100% test coverage across R and C++ codebases (#26, #32, #42, #46).

### Fixed
- **Unseen Categorical Levels in Prediction**: `classifyKamila()` now checks test categorical factors against training levels, throwing an informative error identifying column names and unseen levels (#16, #36).
- **Factor Level Ordering Alignment**: `classifyKamila()` re-aligns test factor levels to match training factor ordering to prevent incorrect index mapping in conditional probability tables (#16, #36).
- **Radial KDE Density at Origin**: Resolved mathematical division-by-zero boundary issue where distance 0 produced infinite log-likelihood (`-Inf`) (#9, #43).
- **NA / NaN Input Validation**: Added strict input validation for missing and invalid values across all exported functions (`kamila()`, `gmsClust()`, `classifyKamila()`, `dummyCodeFactorDf()`, `genMixedData()`) (#3, #39).
- **Modha-Spangler Degenerate Cases**: Added explicit error handling in `gmsClust()` for degenerate cluster assignments and division-by-zero objective calculations (#7, #38).
- **Single-Column Subsetting**: Added `drop = FALSE` in prediction strength subsetting routines to avoid dimensionality loss on single-column factor dataframes (#14, #30).
- **Linter Compliance**: Resolved all 651+ lintr issues and enforced clean styling.

---

## [0.1.2] - 2020-03-10

### Fixed
- Fixed dimension dropping during data frame column subsetting by adding `drop = FALSE` (#20).

### Changed
- Updated package metadata, maintainer information, and paper citations in `DESCRIPTION`.

---

## [0.1.1.4] - 2020-03-10

### Changed
- Updated data frame construction for R 4.0.0 compatibility (handling default `stringsAsFactors = FALSE`).

---

## [0.1.1.3] - 2019-03-14

### Changed
- Updated random number generator test seeds for consistency and reproducibility with R 3.6.0.

---

## [0.1.1.2] - 2018-02-17

### Added
- Added Journal of Statistical Software (JSS) citation and DOI (<doi:10.18637/jss.v083.i13>).

### Fixed
- Fixed namespace importation for `quantile` from `stats`.
- Removed stale build artifacts (`src/symbols.rds`).

---

## [0.1.1.1] - 2016-08-18

### Added
- Initial CRAN release candidate preparation.
- Core KAMILA (KA-means for MIXed LArge data sets) algorithm implementation (`kamila()`).
- Modha-Spangler clustering algorithm (`gmsClust()`).
- High-performance C++ backend routines with Rcpp (`dptmCpp()`, `wkmeans()`).
- Mixed-type synthetic dataset generator (`genMixedData()`).
- Cluster classification for new observations (`classifyKamila()`).
- Prediction strength cluster validation method.
- Full roxygen2 documentation and unit tests.
- GPL-3 license.

---

## [0.1.0] - 2015-10-06

### Added
- Initial project release.
