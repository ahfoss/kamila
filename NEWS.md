# kamila 0.1.3

## New Features

* **End-to-End C++ Iteration Engine**: Core KAMILA while-convergence loop implemented in C++ (`kamilaLoopCpp`) with pre-allocated scratch buffers, achieving near-zero heap memory allocations during iteration.
* **Native C++ Linear Binning & Gaussian Convolution**: Fast $O(N)$ histogram accumulation and discrete Gaussian convolution directly in C++, removing the dependency on calling R's `KernSmooth::bkde` inside the loop.
* **$O(N)$ Quantile Selection & Cache-Coherent Streaming**: Replaced $O(N \log N)$ sorting with `std::nth_element` for exact type-7 quantile bandwidth determination, inlined distance/minDist calculations, and implemented contiguous column-streaming for categorical lookups.
* **Fast C++ Prediction Strength**: Re-implemented prediction strength evaluation in C++ (`calcPsCpp` via Rcpp) with $O(N + K^2)$ algorithmic complexity, replacing previous $O(N^2)$ pairwise operations.
* **Parallel Prediction Strength**: Added multi-core parallel processing support for prediction strength cross-validation runs using `numCores` parameter in `kamila()`.
* **Non-Mixed Data Validation**: Added explicit error validation and informative messaging when non-mixed data (continuous-only or categorical-only) are passed to mixed-data clustering functions.
* **Documentation & Vignettes**: Added comprehensive documentation and usage examples on interpreting prediction strength for cluster number selection.
* **Modern CI/CD Pipelines**: Replaced AppVeyor with GitHub Actions workflows for multi-OS R CMD check, Codecov coverage tracking, and automated linting.
* **100% Test Coverage & Regression Tests**: Added unit tests and snapshot regression tests achieving 100% test coverage across R and C++ codebases.

## Bug Fixes

* **Unseen Categorical Levels in Prediction**: `classifyKamila()` now checks test categorical factors against training levels, throwing an informative error identifying column names and unseen levels.
* **Factor Level Ordering Alignment**: `classifyKamila()` re-aligns test factor levels to match training factor ordering to prevent incorrect index mapping in conditional probability tables.
* **Radial KDE Density at Origin**: Resolved mathematical division-by-zero boundary issue where distance 0 produced infinite log-likelihood (`-Inf`).
* **NA / NaN Input Validation**: Added strict input validation for missing and invalid values across all exported functions (`kamila()`, `gmsClust()`, `classifyKamila()`, `dummyCodeFactorDf()`, `genMixedData()`).
* **Modha-Spangler Degenerate Cases**: Added explicit error handling in `gmsClust()` for degenerate cluster assignments and division-by-zero objective calculations.
* **Single-Column Subsetting**: Added `drop = FALSE` in prediction strength subsetting routines to avoid dimensionality loss on single-column factor dataframes.
* **Linter Compliance**: Resolved all lintr issues and enforced clean styling across R code.

# kamila 0.1.2

## Bug Fixes

* Fixed dimension dropping during data frame column subsetting by adding `drop = FALSE`.

## Other Changes

* Updated package metadata, maintainer information, and paper citations in `DESCRIPTION`.

# kamila 0.1.1.4

## Minor Changes

* Updated data frame construction for R 4.0.0 compatibility (handling default `stringsAsFactors = FALSE`).

# kamila 0.1.1.3

## Minor Changes

* Updated random number generator test seeds for consistency and reproducibility with R 3.6.0.

# kamila 0.1.1.2

## New Features

* Added Journal of Statistical Software (JSS) citation and DOI (<doi:10.18637/jss.v083.i13>).

## Bug Fixes

* Fixed namespace importation for `quantile` from `stats`.
* Removed stale build artifacts (`src/symbols.rds`).

# kamila 0.1.1.1

## New Features

* Initial CRAN release candidate preparation.
* Core KAMILA (KA-means for MIXed LArge data sets) algorithm implementation (`kamila()`).
* Modha-Spangler clustering algorithm (`gmsClust()`).
* High-performance C++ backend routines with Rcpp (`dptmCpp()`, `wkmeans()`).
* Mixed-type synthetic dataset generator (`genMixedData()`).
* Cluster classification for new observations (`classifyKamila()`).
* Prediction strength cluster validation method.
* Full roxygen2 documentation and unit tests.
* GPL-3 license.

# kamila 0.1.0

## New Features

* Initial project release.
