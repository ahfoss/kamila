
#################
# dummy coding of a single factor variable
dummyCodeOneVar <- function(fac) {
  fac <- factor(fac)
  lev <- levels(fac)
  mapply(
    lev,
    FUN=function(ll) {as.numeric(fac==ll)}
  )
}


#################
# Dummy coding of a data frame of factor variables
dummyCodeFactorDf <- function(dat) {
  catTypes <- sapply(dat,class)
  if (!all(catTypes=='factor')) {
    stop('Input data frame must have only factor variables.')
  }
  Reduce(cbind,lapply(dat,dummyCodeOneVar2))
}


#################
# Squared Euclidean distance between two rows of a data frame
squaredEuc <- function(v1,v2) {
  dist(rbind(v1,v2))^2
}


#################
# Distance between a data frame and a single centroid using a user-specified
# distance function.
distFromData2Centroid <- function(dat,centroid,distFun) {
  as.vector(apply(
    X = dat,
    MARGIN = 1,
    FUN = function(vec) distFun(vec,centroid)
  ))
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
    )))
  )
  return(sum(clusterDists$dist))
}




#################
# Cosine distance between (X) an nxp matrix of px1 vectors and (Y)
# a single px1 vector
# Note vectors must already be normalized to length 1
#cosDist <- function(xx,yy) {
#  2 * (1 - xx %*% yy)
#}

#################
# orthant mean
# input is n x p matrix
# rows must already be normed
#orthMean <- function(xx) {
#  cs <- colSums(xx)
#  cs / sqrt(sum(cs^2))
#}

#################
# distortion sum, euclidean
sumDistEuc <- function(xx) {
  mn <- colMeans(xx)
  cent <- xx - rep(1,nrow(xx)) %o% mn
  sum(cent^2)
}

#################
# distortion sum, spherical
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
#' @export

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
  # nclust scalar convertible to integer
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
  clustRes$conCenters <- clustRes$centers[,conVarInds]
  clustRes$catCenters <- clustRes$centers[,-conVarInds]
  return(clustRes)
}


#' A general implementation of Modha-Spangler clustering for mixed-type data.
#'
#' Modha-Spangler clustering estimates the optimal weighting for continuous
#' vs categorical variables using a brute-force search strategy.
#'
#' Modha-Spangler clustering uses a brute-force search strategy to estimate
#' the optimal weighting for continuous vs categorical variables. This
#' implementation admits an arbitrary clustering function and arbitrary
#' objective functions for continuous and categorical variables.
#' 
#' The input parameter clustFun must be a function accepting inputs 
#' (conData, catData, conWeight, nclust, ...) and returning a list containing
#' (at least) the elements cluster, conCenters, and catCenters. The list element
#' "cluster" contains cluster memberships denoted by the integers 1:nclust. The
#' list elements "conCenters" and "catCenters" must be data frames whose rows
#' denote cluster centroids. The function clustFun must allow nclust = 1, in
#' which case $centers returns a data frame with a single row.
#' Input parameters conDist and catDist are functions that must each take two
#' data frame rows as input and return a scalar distance measure.
#' @export
#' @param conData A data frame of continuous variables.
#' @param catData A data frame of categorical variables; the allowable variable types depend on the specific clustering function used.
#' @param nclust An integer specifying the number of clusters.
#' @param searchDensity An integer determining the number of distinct cluster weightings evaluated in the brute-force search.
#' @param clustFun The clustering function to be applied.
#' @param conDist The continuous distance function used to construct the objective function.
#' @param catDist The categorical distance function used to construct the objective function.

gmsClust <- function(
  conData,
  catData,
  nclust,
  searchDensity = 10,
  clustFun = wkmeans,
  conDist = squaredEuc,
  catDist = squaredEuc,
  ...
  ) {
  # variable tests
  conData <- as.data.frame(conData)
  catData <- as.data.frame(catData)
  
  # initializations
  bestObj <- Inf
  weights <- seq(
    from = 1/(1+searchDensity),
    to = 1-1/(1+searchDensity),
    by = 1/(1+searchDensity)
  )
  objFun <- rep(NaN,length(weights))
  Qcon <- rep(NaN,length(weights))
  Qcat <- rep(NaN,length(weights))
  nullClustering <- clustFun(
    conData=conData,
    catData=catData,
    conWeight=1/2,
    nclust=1,
    ...
  )
  totalConDist <- sum(distFromData2Centroid(
    dat=conData,
    centroid=nullClustering$conCenters,
    distFun = conDist
  ))
  totalCatDist <- sum(distFromData2Centroid(
    dat=catData,
    centroid=nullClustering$catCenters,
    distFun = catDist
  ))
  for (i in 1:length(weights)) {
    currentClustering <- clustFun(
      conData=conData,
      catData=catData,
      conWeight=weights[i],
      nclust=nclust,
      ...
    )
    withinConDist <- withinClusterDist(
      conData,
      centroids = currentClustering$conCenters,
      distFun = conDist,
      memberships = currentClustering$cluster
    )
    withinCatDist <- withinClusterDist(
      catData,
      centroids = currentClustering$catCenters,
      distFun = catDist,
      memberships = currentClustering$cluster
    )
    Qcon[i] <- withinConDist / (totalConDist - withinConDist)
    Qcat[i] <- withinCatDist / (totalCatDist - withinCatDist)
    # If w/in cluster distortion is worse than total distortion,
    # the clustering is "infinitely" bad.
    if (Qcon[i] < 0) Qcon[i] <- Inf
    if (Qcat[i] < 0) Qcat[i] <- Inf
    objFun[i] <- Qcon[i] * Qcat[i]
    if (i==1 || objFun[i] < bestObj) {
      bestInd <- i
      bestObj <- objFun[i]
      bestRes <- currentClustering
    }
  }
  if (any(objFun==0)) {
    warning('At least one entry of zero in the objective function; 
    is nclust >= the number of categorical variable 
    level combinations?')
  }
  return(list(
    results = bestRes,
    objFun = objFun,
    Qcon = Qcon,
    Qcat = Qcat,
    bestInd = bestInd,
    weights = weights
  ))
}


#################
## Modha & Spangler optimal weight clustering
## Note that categorical data should be factors
#owClust <- function(conData,categFact,samplingInt=0.1,centers,iter.max,nstart) {
#  numUnique <- unique(cbind(conData,categFact))
#  if (nrow(numUnique) < centers) {
#    stop('more cluster centers than distinct data points.')
#  }
#  ncon <- ncol(conData)
#  ncat <- ncol(categFact)
#  conData <- scale(conData)
#  conDataDf <- as.data.frame(conData)
#  catDum1 <- Reduce(cbind,lapply(categFact,dummyCodeOneVar))
#
#  #print(head(catDum1))
#  #print(dist(head(catDum1))^2)
#
#  catDum2 <- catDum1 / sqrt(ncat/2)
#  ndum <- ncol(catDum2)
#
#  alphas <- seq(0,1,by=samplingInt)
#  alphas <- alphas[-c(1,length(alphas))]
#  objFun <- rep(NaN,length(alphas))
#  Qcon <- rep(NaN,length(alphas))
#  Qcat <- rep(NaN,length(alphas))
#
#  totConDist <- sumDistEuc(conData)
#  totCatDist <- sumDistSph(catDum2)
#
#  for (i in 1:length(alphas)) {
#    # run ith kmeans
#    currentKmRes <- kmeans(
#      x = cbind(
#        alphas[i]*conData
#       ,(1-alphas[i])*catDum1
#      )
#     ,centers = centers
#     ,iter.max = iter.max
#     ,nstart = nstart
#     ,algorithm = 'Hartigan-Wong'
#    )
#    # calculate continuous distortion
#    wnConDistDf <- ddply(
#      .data=cbind(as.data.frame(conDataDf),ind=factor(currentKmRes$cluster))
#     ,~ind
#     ,function(dd) data.frame(dist=sumDistEuc(as.matrix(dd[-(1+ncon)])))
#    )
#    wnConDist <- sum(wnConDistDf$dist)
#    bwConDist <- totConDist - wnConDist
#
#    # calculate categorical distortion
#    wnCatDistDf <- ddply(
#      .data=cbind(as.data.frame(catDum2),ind=factor(currentKmRes$cluster))
#     ,~ind
#     ,function(dd) data.frame(dist=sumDistSph(as.matrix(dd[-(1+ndum)])))
#    )
#    wnCatDist <- sum(wnCatDistDf$dist)
#    bwCatDist <- totCatDist - wnCatDist
#    Qcon[i] <- wnConDist/bwConDist
#    Qcat[i] <- wnCatDist/bwCatDist
#    objFun[i] <- Qcon[i] * Qcat[i]
#    if (objFun[i] < 0) objFun[i] <- Inf
#
#    if (i==1 || objFun[i] < bestObj) {
#      bestInd <- i
#      bestObj <- objFun[i]
#      bestRes <- currentKmRes
#    }
#  }
#  return(list(res=bestRes,objFun=objFun,Qcon=Qcon,Qcat=Qcat,bestInd=bestInd,alphas=alphas))
#}

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
