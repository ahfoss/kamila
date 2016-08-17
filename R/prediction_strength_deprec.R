
# Deprecated functions

# Uses function classifyKamila in kamila.r
# Uses function medeaWgts in medea.r

# create function for calculating prediction strength for a given 
# test clustering size k and train/test classification size k.
# Note: Denote cluster test clusters A_k1, A_k2, ..., A_kk
prediction.strength <- function(testClustering, classifyTest) {
  # Count # pairs in each A_kj classified into the same group using classfication rule
  tab <- table(testClustering,classifyTest)
  nPairsAgree = choose(tab,2)

  # for the jth test cluster, sum over nPairsAgree:
  rs = rowSums(nPairsAgree)

  # divide by # pairs in each test cluster:
  meanAgree = rs / choose(table(testClustering),2)
  meanAgree[is.nan(meanAgree)] <- 0

  # calculate minimum
  return(min(meanAgree))
}



# Function to classify new observations given kmeans results object 
classify.kmeans <- function(res,newData) {
  centers <- res$centers
  kk <- nrow(centers)
  pp <- ncol(centers)

  dists <- mapply(
    function(xx) {
      mahalanobis(
        x=newData
       ,center=centers[xx,]
       ,cov=diag(pp)
      )^(1/2)
    }
   ,1:kk
  )

  membership <- apply(dists,1,which.min)
  return(membership)
}

# simple wrapper for kamila function for passing into nclust1x
# Note dat is a 2-element list with continuous in 1st and categorical in 2nd
# Note that return object must have memberships listed in field "cluster"
kamilaMethod <- function(dat,k) {
  res <- kamila(
    conVar=dat[[1]]
   ,catFactor=dat[[2]]
   ,numClust=k
   ,numInit=10
  )
  res$cluster <- res$finalMemb
  return(res)
}



# Function designed for cyclical variables, e.g. day of week.
# Recodes to equidistant points on the unit circle.
cyclicalCoding <- function(invar) {
  minDetected <- min(invar)
  maxDetected <- max(invar)
  tmp <- 2*invar*pi/(maxDetected+1-minDetected)
  return(cbind(cos(tmp),sin(tmp)))
}


# wrapper for kamila + medea functions for passing into nclust1x
# Note that \code{dat} is a 3-element list with continuous variables as first
# element, categorical as 2nd, and optional cyclical variables as 3rd.
# Cyclical variables must be coded as factors, unit circle recoding
# is done automatically.
# All variable components must be data frames.
# If no cyclical variables, then leave third position empty.
# Note that return object must have memberships listed in field "cluster"
kamilaMedeaMethod <- function(dat,k) {
  numConVar <- ncol(dat[[1]])
  numCatVar <- ncol(dat[[2]])
  if (length(dat) == 3) {
    cyclicVarPresent <- TRUE
    numCycVar <- ncol(dat[[3]])
  } else if (length(dat) == 2) {
    cyclicVarPresent <- FALSE
    numCycVar <- NA
  } else {
    stop('Error in function kamilaMedeaMethod: data must be list of length 2 or 3')
  }

  if (cyclicVarPresent) {
    # convert cyc vars from int to factor for MEDEA
    cycFactors <- dat[[3]]
    for (i in 1:numCycVar) cycFactors[,i] <- factor(dat[[3]][,i])
    medeaVar <- as.data.frame(cbind(dat[[1]],dat[[2]],cycFactors))
    # convert cyc vars from int to numerically coded for KAMILA
    cyclicRecoded <- as.data.frame(lapply(dat[[3]],cyclicalCoding))
    kamilaConVar <- as.data.frame(cbind(dat[[1]],cyclicRecoded))
    # match up indices of medea weights with variables as input into kamila function
    kamilaConWeightInds <- c(
      1:numConVar
     ,rep(
        (1+numConVar+numCatVar):(numConVar+numCatVar+numCycVar)
       ,each=2
      )
    )
  } else {
    medeaVar <- as.data.frame(cbind(dat[[1]],dat[[2]]))
    kamilaConVar <- as.data.frame(dat[[1]])
    kamilaConWeightInds <- 1:numConVar
  }


  medeaTmp <- medeaWgts(
    dat = medeaVar
   ,associationMethod = 'famd'
   ,verbosity=FALSE
  )

  res <- kamila(
    conVar = kamilaConVar
   ,catFactor = dat[[2]]
   ,numClust = k
   ,numInit = 10
   ,conWeights = medeaTmp$wgts[kamilaConWeightInds]
   ,catWeights = medeaTmp$wgts[(numConVar+1):(numConVar+numCatVar)]
  )

  res$cluster <- res$finalMemb
  return(res)
}

# function for subsetting matrix of continuous variables
#
getSubsetCon <- function(dat,inds) dat[inds,]

# function for subsetting matrix of mixed variables
# Data is list, 1st element continuous, 2nd element categorical,
# and possible 3rd element cyclical
getSubsetMix <- function(dat,inds) {
  if (length(dat)==3) {
    return(list(
      dat[[1]][inds,]
     ,data.frame(dat[[2]][inds,])
     ,data.frame(dat[[3]][inds,])
    ))
  } else if (length(dat)==2) {
    return(list(
      dat[[1]][inds,]
     ,data.frame(dat[[2]][inds,])
    ))
  } else {
    stop('Error in function getSubsetMix: data must be list of length 2 or 3')
  }

}

# function for selecting number of clusters using prediction strength
# Default nfold 2 as recommended in tibshirani 2005
# This executes one iteration.
nclust1x <- function(
  inData
 ,clustMethod = function(dat,k) kmeans(x=dat,centers=k)
 ,classMethod = classify.kmeans
 ,subsetMethod = getSubsetCon
 ,kmax = 4
 ,psThresh = 0.8
 ,verbose = FALSE
# ,nfold = 2 # currently only two-fold
)
{

  if (is.list(inData)) {
    nn <- nrow(inData[[1]])
  } else {
    nn <- nrow(inData)
  }

  # 2-fold cross-validation:
  nfold <- 2
  #rawInd <- rep(1:nfold,ceiling(nn / nfold))
  #trainInd <- sample(rawInd,size=nn,replace=FALSE)
  trainInd <- sample(1:nn,size=ceiling(nn/nfold),replace=FALSE)

  sizeVec <- 2:kmax
  psMatrix <- matrix(nrow=length(sizeVec),ncol=2,dimnames=list(nclust=sizeVec,fold=1:2))

  # for k = 2, 3, ..., K_max number of clusters
  for (i in 1:length(sizeVec)) {

    if (verbose) cat('\n Now starting k =',sizeVec[i],'\n')

    # Cluster training set into k clusters, obtain classification rule
    trainingDataTmp <- subsetMethod(inData,trainInd)
    trainRes <- clustMethod(dat=trainingDataTmp,k=sizeVec[i])

    # Cluster test set individually obtaining clusters A_k1, A_k2, ..., A_kk
    testRes <- clustMethod(subsetMethod(inData,-trainInd),sizeVec[i])

    # Classify test set using classification rule
    testClass <- classMethod(trainRes,newData=subsetMethod(inData,-trainInd))

    # Obtain prediction strenght
    psMatrix[i,1] <- prediction.strength(testRes$cluster,testClass)

    # FOLD 2
    trainRes.f2 <- clustMethod(subsetMethod(inData,-trainInd),sizeVec[i])
    testRes.f2 <- clustMethod(subsetMethod(inData,trainInd),sizeVec[i])
    testClass.f2 <- classMethod(trainRes.f2,newData=subsetMethod(inData,trainInd))
    psMatrix[i,2] <- prediction.strength(testRes.f2$cluster,testClass.f2)

    
  }

  # for each fold, average
  psVec <- apply(psMatrix,1,mean)

  return(list(psVec=psVec,psMatrix=psMatrix))
}

nclustFull <- function(
  inData
 ,clustMethod = function(dat,k) kmeans(x=dat,centers=k)
 ,classMethod = classify.kmeans
 ,subsetMethod = getSubsetCon
 ,kmax = 4
 ,psThresh = 0.8
 ,nrep = 10
# ,nfold = 2 # currently only two-fold
)
{
  allPsMatrix <- matrix(
    nrow=kmax-1
   ,ncol=nrep
   ,dimnames=list(nclust=2:kmax,rep=1:nrep)
  )

  for (iter in 1:nrep) {

    ithRes <- nclust1x(
      inData = inData
     ,clustMethod = clustMethod
     ,classMethod = classMethod
     ,subsetMethod = subsetMethod
     ,kmax = kmax
     ,psThresh = psThresh
     # ,nfold = 2 # currently only two-fold
    )
    allPsMatrix[,iter] <- ithRes$psVec
  }

  psVec <- apply(allPsMatrix,1,mean)
  overThresh <- psVec > psThresh
  if (any(overThresh)) {
    bestK <- names(psVec)[max(which(overThresh))]
  } else {
    bestK <- names(psVec)[which.max(psVec)]
  }

  return(list(
    k=as.numeric(bestK)
   ,psVec=psVec
   ,scores=allPsMatrix
   ,thresh=psThresh
   ,kmax=kmax
   ,nrep=nrep
  ))
}

# Uses Hmisc for errbar
plot_nclustFull <- function(obj) {
  if ('Hmisc' %in% installed.packages()) {
    library(Hmisc)
    ses <- apply(obj$scores,1,sd) #/sqrt(obj$nrep)
    with(obj,errbar(x=2:kmax,y=psVec,yplus=psVec+ses,yminus=psVec-ses,ylim=c(0,1)))
    abline(h=obj$thresh,lty=2)
    text(x=obj$k,y=1,labels='*',cex=3)
  } else {
    stop('plot.nclustFull requires package Hmisc to run')
  }
}

