
# categorical coding for k-means
# Scaled such that euclidean distance between any two levels is one
dummyCodeOneVar <- function(fac) {
  fac <- factor(fac)
  cbind(
    as.numeric(fac==levels(fac)[1])
   ,model.matrix( ~ fac)[,-1]
  )/sqrt(2)
}

#################
# Cosine distance between (X) an nxp matrix of px1 vectors and (Y)
# a single px1 vector
# Note vectors must already be normalized to length 1
cosDist <- function(xx,yy) {
  2 * (1 - xx %*% yy)
}

#################
# orthant mean
# input is n x p matrix
# rows must already be normed
orthMean <- function(xx) {
  cs <- colSums(xx)
  cs / sqrt(sum(cs^2))
}

#################
# distortion sum, euclidean
sumDistEuc <- function(xx) {
  mn <- apply(xx,2,mean)
  cent <- xx - rep(1,nrow(xx)) %o% mn
  sum(cent^2)
}

#################
# distortion sum, spherical
sumDistSph<- function(xx) {
  om <- orthMean(xx)
  sum(cosDist(xx,om))
}

#################
# Modha & Spangler clustering
# Note that categorical data should be factors
#' @export
msClust <- function(conData,categFact,samplingInt=0.1,centers,iter.max,nstart) {
  numUnique <- unique(cbind(conData,categFact))
  if (nrow(numUnique) < centers) {
    stop('more cluster centers than distinct data points.')
  }
  ncon <- ncol(conData)
  ncat <- ncol(categFact)
  conData <- scale(conData)
  conDataDf <- as.data.frame(conData)
  catDum1 <- Reduce(cbind,lapply(categFact,dummyCodeOneVar))

  #print(head(catDum1))
  #print(dist(head(catDum1))^2)

  catDum2 <- catDum1 / sqrt(ncat/2)
  ndum <- ncol(catDum2)

  alphas <- seq(0,1,by=samplingInt)
  alphas <- alphas[-c(1,length(alphas))]
  objFun <- rep(NaN,length(alphas))
  Qcon <- rep(NaN,length(alphas))
  Qcat <- rep(NaN,length(alphas))

  totConDist <- sumDistEuc(conData)
  totCatDist <- sumDistSph(catDum2)

  for (i in 1:length(alphas)) {
    # run ith kmeans
    currentKmRes <- kmeans(
      x = cbind(
        alphas[i]*conData
       ,(1-alphas[i])*catDum1
      )
     ,centers = centers
     ,iter.max = iter.max
     ,nstart = nstart
     ,algorithm = 'Hartigan-Wong'
    )
    # calculate continuous distortion
    wnConDistDf <- ddply(
      .data=cbind(as.data.frame(conDataDf),ind=factor(currentKmRes$cluster))
     ,~ind
     ,function(dd) data.frame(dist=sumDistEuc(as.matrix(dd[-(1+ncon)])))
    )
    wnConDist <- sum(wnConDistDf$dist)
    bwConDist <- totConDist - wnConDist

    # calculate categorical distortion
    wnCatDistDf <- ddply(
      .data=cbind(as.data.frame(catDum2),ind=factor(currentKmRes$cluster))
     ,~ind
     ,function(dd) data.frame(dist=sumDistSph(as.matrix(dd[-(1+ndum)])))
    )
    wnCatDist <- sum(wnCatDistDf$dist)
    bwCatDist <- totCatDist - wnCatDist
    Qcon[i] <- wnConDist/bwConDist
    Qcat[i] <- wnCatDist/bwCatDist
    objFun[i] <- Qcon[i] * Qcat[i]
    if (objFun[i] < 0) objFun[i] <- Inf

    if (i==1 || objFun[i] < bestObj) {
      bestInd <- i
      bestObj <- objFun[i]
      bestRes <- currentKmRes
    }
  }
  return(list(res=bestRes,objFun=objFun,Qcon=Qcon,Qcat=Qcat,bestInd=bestInd,alphas=alphas))
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
