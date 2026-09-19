# kamila 0.2.0

* **CRAN Milestone Release**: Version 0.2.0 is the official CRAN release consolidating all major enhancements, C++ engine optimizations, algorithmic speedups, bug fixes, and infrastructure improvements developed across versions 0.1.3 and 0.1.4 since the last CRAN release (0.1.2).

# kamila 0.1.4

## Added
* **CRAN Release News (`NEWS.md`)**: Added structured `NEWS.md` file adhering to CRAN release notes standards and `newsmd` guidelines.
* **End-to-End C++ Iteration Engine (`kamilaLoopCpp`)**: Implemented the core while convergence loop in C++ with pre-allocated scratch buffers (`distMat`, `minDist`, `catLogLiks`, `allLogLiks`, `membOld`, `membNew`), achieving near-zero heap memory allocations in the iteration loop (#49, #61).
* **Native C++ Linear Binning & Gaussian Convolution**: Integrated fast $O(N)$ histogram accumulation and discrete Gaussian convolution directly in C++, removing the dependency on calling R's `KernSmooth::bkde` inside the loop (#49, #61).
* **$O(N)$ Quantile Selection & Cache-Coherent Streaming**: Replaced $O(N \log N)$ sorting with `std::nth_element` for exact type-7 quantile bandwidth determination, inlined distance/minDist calculations, and implemented contiguous column-streaming for categorical lookups, demonstrating statistically significant performance superiority across Small, Medium, and Large datasets (#49, #61).
* **Interactive Horse-Race Web App (WebR / Shinylive)**: Built and deployed a zero-install interactive benchmark web app running client-side via Shinylive and WebAssembly (#62, #64).
* **Statistical Superiority Testing Framework**: Added automated benchmarking and performance regression testing suite (`inst/benchmarks/run_superiority_benchmark.R`) with CI verification workflow (#63).

## Changed
* **CRAN Compliance & Cleanup**: Updated `.Rbuildignore` to ignore non-package root files and fixed `inherits(fac, "factor")` in `R/misc_functions.R`.
* **Documentation & Badges**: Added CRAN version badge, GPL-3 license badge, quick installation guide, and reproducible quick-start clustering example to `README.md`.

# kamila 0.1.3

## Added
* **Restored Progress Indicator with Non-Dismissible Notification**: Restored the original floating `withProgress` / `incProgress` notification indicator for benchmark runs, configured with hidden close button styling (`.shiny-notification-close { display: none !important; }`) to maintain full progress visibility during execution.
* **Standardized Multi-Start Benchmark Runs**: Standardized all compared clustering methods in the interactive Shiny benchmark dashboard (`inst/shiny/horserace/app.R`) to use 5 random initialization starts (`numInit = 5`, `nstart = 5`, `nrep = 5`, `nbKeep = 5`) for fair computational and performance comparisons.
* **Fast C++ Prediction Strength**: Re-implemented prediction strength evaluation in C++ (`calcPsCpp` via Rcpp) with $O(N + K^2)$ algorithmic complexity, replacing previous $O(N^2)$ pairwise operations (#12, #42).
* **Parallel Prediction Strength**: Added multi-core parallel processing support for prediction strength cross-validation runs using `numCores` parameter in `kamila()` (#46).
* **Non-Mixed Data Validation**: Added explicit error validation and informative messaging when non-mixed data (continuous-only or categorical-only) are passed to mixed-data clustering functions (#19).
* **Prediction Strength Documentation**: Comprehensive vignettes, usage examples, and documentation on interpreting prediction strength for cluster number selection (#15, #40).
* **Modern CI/CD Pipelines**: Replaced AppVeyor with GitHub Actions workflows for multi-OS R CMD check, Codecov coverage tracking, and automated linting.
* **100% Test Coverage & Regression Tests**: Added unit tests and snapshot regression tests achieving 100% test coverage across R and C++ codebases (#26, #32, #42, #46).

## Fixed
* **Unseen Categorical Levels in Prediction**: `classifyKamila()` now checks test categorical factors against training levels, throwing an informative error identifying column names and unseen levels (#16, #36).
* **Factor Level Ordering Alignment**: `classifyKamila()` re-aligns test factor levels to match training factor ordering to prevent incorrect index mapping in conditional probability tables (#16, #36).
* **Radial KDE Density at Origin**: Resolved mathematical division-by-zero boundary issue where distance 0 produced infinite log-likelihood (`-Inf`) (#9, #43).
* **NA / NaN Input Validation**: Added strict input validation for missing and invalid values across all exported functions (`kamila()`, `gmsClust()`, `classifyKamila()`, `dummyCodeFactorDf()`, `genMixedData()`) (#3, #39).
* **Modha-Spangler Degenerate Cases**: Added explicit error handling in `gmsClust()` for degenerate cluster assignments and division-by-zero objective calculations (#7, #38).
* **Single-Column Subsetting**: Added `drop = FALSE` in prediction strength subsetting routines to avoid dimensionality loss on single-column factor dataframes (#14, #30).
* **Linter Compliance**: Resolved all 651+ lintr issues and enforced clean styling.

# kamila 0.1.2

## Fixed
* Fixed dimension dropping during data frame column subsetting by adding `drop = FALSE` (#20).

## Changed
* Updated package metadata, maintainer information, and paper citations in `DESCRIPTION`.

# kamila 0.1.1.4

## Changed
* Updated data frame construction for R 4.0.0 compatibility (handling default `stringsAsFactors = FALSE`).

# kamila 0.1.1.3

## Changed
* Updated random number generator test seeds for consistency and reproducibility with R 3.6.0.

# kamila 0.1.1.2

## Added
* Added Journal of Statistical Software (JSS) citation and DOI (<doi:10.18637/jss.v083.i13>).

## Fixed
* Fixed namespace importation for `quantile` from `stats`.
* Removed stale build artifacts (`src/symbols.rds`).

# kamila 0.1.1.1

## Added
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

## Added
* Initial project release.
