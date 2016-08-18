
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
    to = 1 - 1/(1+searchDensity),
    by = 1/(1+searchDensity)
  )
  objFun <- rep(NaN,length(weights))
  Qcon <- rep(NaN,length(weights))
  Qcat <- rep(NaN,length(weights))
  nullConClustering <- clustFun(
    conData = conData,
    catData = rep(catData[1,1], nrow(catData)),
    conWeight = 1,
    nclust = 1,
    ...
  )
  nullCatClustering <- clustFun(
    conData = rep(conData[1,1], nrow(conData)),
    catData = catData,
    conWeight = 0,
    nclust = 1,
    ...
  )
  totalConDist <- sum(distFromData2Centroid(
    dat=conData,
    centroid=nullConClustering$conCenters,
    distFun = conDist
  ))
  totalCatDist <- sum(distFromData2Centroid(
    dat=catData,
    centroid=nullCatClustering$catCenters,
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

