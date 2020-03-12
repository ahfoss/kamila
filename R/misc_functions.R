
# For some reason standard importing isn't working for plyr.
# Workaround:
#' @import plyr

#################
# Dummy coding of a single factor variable.
# Clustering dummy coding, where each variable level gets its own
# dummy variable to ensure that distances between each level are
# consistent; as opposed to regression coding with an intercept
# and a dropped level etc.
dummyCodeOneVar <- function(fac) {
  if (class(fac) != 'factor') fac <- factor(fac)
  lev <- levels(fac)
  mapply(
    lev,
    FUN=function(ll) {as.numeric(fac==ll)}
  )
}


#' Dummy coding of a data frame of factor variables
#'
#' Given a data frame of factor variables, this function returns a numeric
#' matrix of 0--1 dummy-coded variables.
#'
#' @export
#' @param dat A data frame of factor variables
#' @return A numeric matrix of 0--1 dummy coded variables
#' @examples
#' dd <- data.frame(a=factor(1:8), b=factor(letters[1:8]), stringsAsFactors = TRUE, stringsAsFactors = TRUE)
#' dummyCodeFactorDf(dd)
dummyCodeFactorDf <- function(dat) {
  catTypes <- sapply(dat,class)
  if (!all(catTypes=='factor')) {
    stop('Input data frame must have only factor variables.')
  }
  outMat <- Reduce(cbind,lapply(dat,dummyCodeOneVar))
  levNames <- lapply(dat,levels)
  colnames(outMat) <- paste(
    rep(colnames(dat), times=lapply(levNames,length)),
    unlist(levNames),
    sep = '_'
  )
  return(outMat)
}


#################
# Squared Euclidean distance between two rows of a data frame
squaredEuc <- function(v1,v2) {
  sum( (v1-v2)^2 )
}


#################
# Distance between a data frame and a single centroid using a user-specified
# distance function.
distFromData2Centroid <- function(dat,centroid,distFun) {
  numRows <- nrow(dat)
  outDists <- rep(Inf,numRows)
  for (i in 1:numRows) {
    outDists[i] <- distFun(dat[i,], centroid)
  }
  return(outDists)
# can't use apply since it converts from factor to character
#  as.vector(apply(
#    X = dat,
#    MARGIN = 1,
#    FUN = function(vec) distFun(vec,centroid)
#  ))
}


#################
# Cluster-wise distances between points of a dataset and their respective
# cluster centroids. 
# Centroids must be a k X p data frame, where k is # clusters and p is
# the number of variables in dat.
# The cluster memberships must be integers giving the corresponding row of
# centroids.
withinClusterDist <- function(dat,centroids,distFun,memberships) {
  if (is.null(centroids)) {
    stop('Formal parameter centroids cannot be NULL.')
  } else if (nrow(centroids) < 1) {
    stop('Must include at least one centroid.')
  }
  centroids <- as.data.frame(centroids)
  clusterDists <- ddply(
    .data = cbind(as.data.frame(dat),clust_uNiQuE=factor(memberships)),
    ~ clust_uNiQuE,
    function(dd) data.frame(dist=sum(distFromData2Centroid(
      dat=data.frame(dd[,-ncol(dd)]),
      centroid=data.frame(centroids[dd$clust_uNiQuE[1],]),
      distFun=distFun
    )), stringsAsFactors = TRUE)
  )
  return(sum(clusterDists$dist))
}




#################
## Cosine distance between (X) an nxp matrix of px1 vectors and (Y)
## a single px1 vector
## Note vectors must already be normalized to length 1
#cosDist <- function(xx,yy) {
#  2 * (1 - xx %*% yy)
#}

#################
## orthant mean
## input is n x p matrix
## rows must already be normed
#orthMean <- function(xx) {
#  cs <- colSums(xx)
#  cs / sqrt(sum(cs^2))
#}

#################
## distortion sum, euclidean
#sumDistEuc <- function(xx) {
#  mn <- colMeans(xx)
#  cent <- xx - rep(1,nrow(xx)) %o% mn
#  sum(cent^2)
#}

#################
## distortion sum, spherical
#sumDistSph<- function(xx) {
#  om <- orthMean(xx)
#  sum(cosDist(xx,om))
#}


#################
# Weighted k-means for mixed-type data.
# conData and catData must be coercible to data frames.
# catData can be a data frame of factor variables or dummy-coded categorical
# variables. If factor variables, catData is dummy coded.
# conWeight must be an element of [0,1], and gives the weight applied to the
# continuous variables. The categorical variables are assigned a weight of
# 1 - conWeight.
# ... optional arguments passed to kmeans. An input argument of 'centers' is
# ignored.

#' Weighted k-means for mixed-type data
#'
#' Weighted k-means for mixed continuous and categorical variables. A 
#' user-specified weight \code{conWeight} controls the relative contribution of the
#' variable types to the cluster solution.
#'
#' A simple adaptation of \code{stats::kmeans} to mixed-type data.  Continuous
#' variables are multiplied by the input parameter \code{conWeight}, and categorical
#' variables are multipled by \code{1-conWeight}. If factor variables are input to
#' \code{catData}, they are transformed to 0-1 dummy coded variables with the function
#' \code{dummyCodeFactorDf}.
#'
#' @export
#' @importFrom stats kmeans
#' @param conData The continuous variables. Must be coercible to a data frame.
#' @param catData The categorical variables, either as factors or dummy-coded variables. Must be coercible to a data frame.
#' @param conWeight The continuous weight; must be between 0 and 1. The categorical weight is \code{1-conWeight}.
#' @param nclust The number of clusters.
#' @param ... Optional arguments passed to \code{kmeans}.
#' @return A stats::kmeans results object, with additional slots \code{conCenters} and \code{catCenters} giving the actual centers adjusted for the weighting process.
#' @seealso \code{\link{dummyCodeFactorDf}}
#' @seealso \code{\link[stats]{kmeans}}
#' @examples
#' # Generate toy data set with poor quality categorical variables and good
#' # quality continuous variables.
#' set.seed(1)
#' dat <- genMixedData(200, nConVar=2, nCatVar=2, nCatLevels=4, nConWithErr=2,
#'   nCatWithErr=2, popProportions=c(.5,.5), conErrLev=0.3, catErrLev=0.8)
#' catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
#' conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)
#'
#' # A clustering that emphasizes the continuous variables
#' r1 <- with(dat,wkmeans(conDf, catDf, 0.9, 2))
#' table(r1$cluster, dat$trueID)
#'
#' # A clustering that emphasizes the categorical variables; note argument
#' # passed to the underlying stats::kmeans function
#' r2 <- with(dat,wkmeans(conDf, catDf, 0.1, 2, nstart=4))
#' table(r2$cluster, dat$trueID)
wkmeans <- function(
  conData,
  catData,
  conWeight,
  nclust,
  ...
  ) {
  conData <- as.data.frame(conData)
  catData <- as.data.frame(catData)
  catTypes <- sapply(catData,class)
  if (!all(catTypes=='factor') && !all(catTypes %in% c('integer','numeric'))) {
    stop('Argument catData must be a data frame with all factor variables or all numeric variables.')
  }
  if (!all(sapply(conData,class) %in% c('integer','numeric'))) {
    stop('Argument conData must be data frame with all integer/numeric types.')
  }
  if ( !is.numeric(conWeight) || conWeight > 1 || conWeight < 0 ) {
    stop('Argument conWeight must be numeric and in [0,1].')
  }
  nclust <- as.integer(nclust)
  if (nclust <= 0) {
    stop('Argument nclust must be a positive integer')
  }
  if (catTypes[1] == 'factor') {
    catData <- dummyCodeFactorDf(catData)
  }
  dotArgs <- list(...)
  dotArgs$centers <- NULL
  clustRes <- do.call(
    what = function(...) {
      kmeans(
        x = cbind(conData * conWeight,catData * (1-conWeight)),
        centers=nclust,
        ...
      )
    },
    args = dotArgs
  )
  conVarInds <- 1:ncol(conData)
  # "Reconstitute" means from their scaled versions
  clustRes$conCenters <- clustRes$centers[,conVarInds] / conWeight
  clustRes$catCenters <- clustRes$centers[,-conVarInds] / (1-conWeight)
  return(clustRes)
}


#################
# define a function that maps clusters
# to true class (i.e. as computed in purity)
# Input: two factor variables
clust2class <- function(clust, trueClass) {
  if (length(clust) != length(trueClass)) stop('input vectors must be same length.')
  tab <- table(clust, trueClass)
  maxInds <- apply(tab,1,which.max)
  newVar <- factor(maxInds[clust])
  return(newVar)
}

#################
# hand written function to calculate purity;
# replicates NMF::purity().
# Note that purity is equivalent to calculating
# micro-averaged precision (which is equivalent
# to micro-averaged recall in this context), if
# using the same cluster-class assignment rule as
# in purity calculations.
myPurity <- function(clust, trueClass) {
  predClass <- clust2class(clust, trueClass)
  tab <- table(predClass, trueClass)
  return(sum(diag(tab)) / sum(tab))
}

#################
# Calculates macro-precision from cluster assignments
# and true class variable; use assigment rule from
# Modha-Spangler paper (same rule used in purity 
# calculations).
macroPrecRec <- function(clust, trueClass) {
  predClass <- clust2class(clust, trueClass)
  tab <- table(predClass, trueClass)
  nLev <- nrow(tab)
  if (nLev != ncol(tab)) stop('predClass and trueClass have different number of levels')
  precisionSum <- 0
  recallSum <- 0
  for (i in 1:nLev) {
    precisionSum <- precisionSum + tab[i,i] / sum(tab[i,])
    recallSum <- recallSum + tab[i,i] / sum(tab[,i])
  }
  macroP <- precisionSum / nLev
  macroR <- recallSum / nLev
  return(list(macroP=macroP, macroR=macroR))
}
