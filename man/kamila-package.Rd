\name{kamila-package}
\alias{kamila-package}
\docType{package}
\title{
Clustering for mixed continuous and categorical data sets
}
\description{
A collection of methods for clustering mixed type data, including KAMILA (KAy-means for MIxed LArge data) and a flexible implementation of Modha-Spangler clustering
}

\author{
Alex Foss and Marianthi Markatou

Maintainer: Alex Foss <alexanderhfoss@gmail.com>
}
\references{
AH Foss, M Markatou, B Ray, and A Heching (in press). A semiparametric method for clustering mixed data. Machine Learning, DOI: 10.1007/s10994-016-5575-7.

DS Modha and S Spangler (2003). Feature weighting in k-means clustering. Machine Learning 52(3), 217-237.
}

\examples{
\dontrun{
# Generate a synthetic mixed-type dataset
set.seed(123)
dat <- genMixedData(100, 2, 2,
  nCatLevels = 4, nConWithErr = 1, nCatWithErr = 1,
  popProportions = c(0.3, 0.7), conErrLev = 0.3, catErrLev = 0.2
)
conVars <- data.frame(scale(dat$conVars))
catVarsFac <- data.frame(lapply(as.data.frame(dat$catVars), factor))
catVarsDum <- dummyCodeFactorDf(catVarsFac)

# Modha-Spangler clustering with default Hartigan-Wong algorithm
gmsResHw <- gmsClust(conVars, catVarsDum, nclust = 2)

# Modha-Spangler clustering with Lloyd algorithm
# NOTE searchDensity should be >= 10 for optimal performance:
# this is just a syntax demo
gmsResLloyd <- gmsClust(conVars, catVarsDum, nclust = 2,
  algorithm = "Lloyd", searchDensity = 3)

# KAMILA clustering with fixed number of clusters
kamRes <- kamila(conVars, catVarsFac, numClust = 2, numInit = 10)

# KAMILA clustering with prediction strength estimation for number of clusters
kamPsRes <- kamila(conVars, catVarsFac, numClust = 2:3, numInit = 10,
  calcNumClust = "ps", numPredStrCvRun = 5, predStrThresh = 0.5)
kamPsRes$nClust$bestNClust

# Plot results using base R graphics
plot(conVars[, 1], conVars[, 2], col = kamRes$finalMemb, pch = 19,
  main = "KAMILA Clusters", xlab = "Continuous Variable 1",
  ylab = "Continuous Variable 2")
}
}
