######################
# Helper functions implementing observation-level sample weights
# (e.g. survey weights) for the KAMILA algorithm.
#
# The weighted estimators below reduce exactly to their unweighted
# counterparts in kamila.R / cppfunctions.cpp when all weights equal 1.
#
# References:
#  Foss A, Markatou M, Ray B, Heching A; A semiparametric method for
#    clustering mixed data. Machine Learning, 105(3), 419-458. 2016.
#  Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and Hadoop.
#    Journal of Statistical Software, 83(13). 2018.

######################
# Weighted quantile via linear interpolation of the weighted empirical CDF.
#' @importFrom stats approx
weightedQuantile <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  if (length(x) == 1) return(rep(x, length(probs)))
  stats::approx(x = cw, y = x, xout = probs, rule = 2, ties = "ordered")$y
}

######################
# Weighted analogue of stats::bw.nrd0 (Silverman's rule of thumb).
# Uses the Kish effective sample size neff = (sum w)^2 / sum(w^2)
# in place of n, and weighted sd / weighted IQR.
weightedBwNrd0 <- function(x, w) {
  if (length(x) < 2L) stop("need at least 2 data points")
  wsum <- sum(w)
  neff <- wsum^2 / sum(w^2)
  mu <- sum(w * x) / wsum
  sdw <- sqrt(sum(w * (x - mu)^2) / wsum)
  iqrw <- diff(weightedQuantile(x, w, c(0.25, 0.75)))
  lo <- min(sdw, iqrw / 1.349)
  # fallback chain mirrors stats::bw.nrd0
  if (lo == 0) lo <- sdw
  if (lo == 0) lo <- abs(x[1L])
  if (lo == 0) lo <- 1
  0.9 * lo * neff^(-0.2)
}

######################
# Weighted analogue of the Rcpp function aggregateMeans:
# cluster centroids as weighted means. Empty clusters get 0 rows,
# matching aggregateMeans behavior.
weightedAggregateMeans <- function(conVar, membNew, kk, obsWeights) {
  conVar <- as.matrix(conVar)
  pp <- ncol(conVar)
  outMat <- matrix(0, nrow = kk, ncol = pp)
  wsums <- rowsum(obsWeights, group = membNew)
  wx <- rowsum(conVar * obsWeights, group = membNew)
  presentClusts <- as.integer(rownames(wsums))
  outMat[presentClusts, ] <- wx / as.vector(wsums)
  outMat
}

######################
# R replication of the Rcpp function smooth2dTable (categorical kernel
# smoothing of a two-way table), operating on weighted counts.
smooth2dTableWeighted <- function(inputTab, catBw) {
  dim1 <- nrow(inputTab)
  dim2 <- ncol(inputTab)
  colS <- colSums(inputTab)
  offCounts1 <- sweep(-inputTab, 2, colS, "+")
  midMat <- (1 - catBw) * inputTab + catBw / (dim1 - 1) * offCounts1
  rowS <- rowSums(midMat)
  offCounts2 <- sweep(-midMat, 1, rowS, "+")
  (1 - catBw) * midMat + catBw / (dim2 - 1) * offCounts2
}

######################
# Weighted analogue of the Rcpp function jointTabSmoothedList:
# joint (cluster x level) tables of summed observation weights,
# smoothed with the categorical kernel.
weightedJointTabSmoothedList <- function(catFactorNum, membNew, numLev, catBw, kk, obsWeights) {
  qq <- ncol(catFactorNum)
  outList <- vector("list", qq)
  membF <- factor(membNew, levels = 1:kk)
  for (q in 1:qq) {
    levF <- factor(catFactorNum[, q], levels = 1:numLev[q])
    tab <- tapply(obsWeights, list(membF, levF), sum, default = 0)
    tab <- matrix(
      tab
      ,nrow = kk
      ,ncol = numLev[q]
      ,dimnames = list(clust = 1:kk, level = 1:numLev[q])
    )
    if (catBw != 0) {
      tab <- smooth2dTableWeighted(tab, catBw)
    }
    outList[[q]] <- tab
  }
  outList
}
