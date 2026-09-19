
# kamila [![CRAN Version](https://www.r-pkg.org/badges/version/kamila)](https://cran.r-project.org/package=kamila) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![R-CMD-check](https://github.com/ahfoss/kamila/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ahfoss/kamila/actions/workflows/R-CMD-check.yaml) [![test-coverage](https://github.com/ahfoss/kamila/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/ahfoss/kamila/actions/workflows/test-coverage.yaml) [![Codecov test coverage](https://codecov.io/gh/ahfoss/kamila/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ahfoss/kamila) [![lint](https://github.com/ahfoss/kamila/actions/workflows/lint.yaml/badge.svg)](https://github.com/ahfoss/kamila/actions/workflows/lint.yaml) [![Shinylive App](https://img.shields.io/badge/Shinylive-Interactive_Demo-blue?logo=r)](https://ahfoss.github.io/kamila/) [![CRAN_Status_Badge](https://cranlogs.r-pkg.org/badges/grand-total/kamila)](https://cran.r-project.org/package=kamila) [![CRAN_Status_Badge](https://cranlogs.r-pkg.org/badges/kamila)](https://cran.r-project.org/package=kamila)

`kamila` implements methods for clustering mixed-type data (continuous and nominal categorical variables), specifically **KAMILA** (KA-means for MIXed LArge datasets) and **Modha-Spangler clustering**. Special attention is paid to the problem of equitably balancing the contribution of continuous and categorical variables without requiring artificial dummy coding.

## Installation

### Stable CRAN Release
```r
install.packages("kamila")
```

### Development Version
```r
# install.packages("remotes")
remotes::install_github("ahfoss/kamila")
```

## Quick Start

```r
library(kamila)

# Generate synthetic mixed-type data (continuous + categorical)
set.seed(123)
dat <- genMixedData(
  sampSize = 200,
  nCon = 3,
  nCat = 2,
  nClust = 3,
  nIndepConCat = 0
)

# Run KAMILA clustering (e.g., K = 3 clusters with 5 random initializations)
kamRes <- kamila(
  conVar = dat$conVars,
  catFactor = dat$catVars,
  numClust = 3,
  numInit = 5
)

# Inspect cluster assignments
table(True = dat$trueID, Predicted = kamRes$finalMemb)
```

## Scientific Publications
For an in-depth discussion of the challenges involved in clustering mixed-type data, please see:
* [Foss, Markatou, Ray, and Heching (2016). A semiparametric method for clustering mixed data. **Machine Learning**, 105(3), 419-458. DOI: 10.1007/s10994-016-5575-7](https://link.springer.com/article/10.1007/s10994-016-5575-7)
* [Foss and Markatou (2018). kamila: Clustering Mixed-Type Data in R and Hadoop. **Journal of Statistical Software**, 83(13). DOI: 10.18637/jss.v083.i13](https://www.jstatsoft.org/article/view/v083i13)
* [Foss, Markatou, and Ray (2018). Distance Metrics and Clustering Methods for Mixed-Type Data. **International Statistical Review**. DOI: 10.1111/insr.12274.](https://onlinelibrary.wiley.com/doi/abs/10.1111/insr.12274)

## Updates & Changelog
For release notes and version history, see CHANGELOG.md.

## [Under Construction] Interactive Clustering Horse-Race App

Explore and benchmark KAMILA against competing mixed-type clustering algorithms (FlexMix, PAM + Gower's distance, k-prototypes, VarSelLCM) interactively:

### [🚀 Launch Live Web App (Client-Side WebAssembly)](https://ahfoss.github.io/kamila/)
> **No installation required:** Runs 100% client-side in your web browser via **WebR & Shinylive**.

---

### Run Locally in R
To run the dashboard locally:

```r
shiny::runGitHub(
  repo = "kamila",
  username = "ahfoss",
  ref = "master",
  subdir = "inst/shiny/horserace"
)
```

## Performance Benchmarking & Superiority Testing

For contributors proposing algorithmic speedups, `kamila` includes a statistical superiority testing framework. For details on how to run local benchmarks or trigger CI via PR labels, see SUPERIORITY_TESTING.md.
