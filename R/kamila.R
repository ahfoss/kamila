

#library(abind)
#library(KernSmooth) # for bkde
#library(gtools) # for rdirichlet


#############################
# compile rcpp functions, including
# - dptm, used in dptmCpp for calculating weighted Euc distances
# - rowMin, used to take row mins of a matrix
# - rowMax, used to take row maxs of a matrix
# - rowMaxInds, gives indices of row maxs of a matrix
# - sumMatList, sums all matrices in a list of matrices
# - getIndividualLogProbs, matches factor level log probs to observed data
# - aggregateMeans, implement restricted but speedy version of aggregate()
# - jointTabSmoothedList and two helper functions
#library(Rcpp)

#' @useDynLib kamila
#' @importFrom Rcpp sourceCpp

#tryCatch(
#  sourceCpp("./src/cppfunctions.cpp")
#  ,error = function(e){}
#)

# kdat must be a data frame with factor rows, levels 1:L;
# factors must have the correct levels specified to ensure
# missing cells are not lost
# This version deprecated; use my Rcpp implementation
myCatKern <- function(kdat,bw,tabOnly=TRUE) {
  tt <- table(kdat)
  dims <- dim(tt)
  ndim <- length(dims)
  nn <- nrow(kdat)
  if (length(bw) == 1) bw <- rep(bw,ndim)
  for (dd in 1:ndim) {
    dimCounts <- apply(X=tt,MARGIN=(1:ndim)[-dd],FUN=sum)
    l1 = list()
    for (i in 1:dims[dd]) {
      dimInds <- dims
      dimInds[dd] <- 1
      rotatedMat <- array(dimCounts,dim=dimInds)
      l1[[i]] <- rotatedMat
    }
    l1$along <- dd
    offCounts <- do.call(abind,l1)
    offCounts  <- offCounts - tt
    tt <- (1-bw[dd])*tt + bw[dd]/(dims[dd]-1)*offCounts
  }
  if (tabOnly) return(tt)
  preds <- mapply(
    FUN = function(ii) {
      l2 <- as.list(c(0,as.numeric(kdat[ii,])))
      l2[[1]] <- tt/nn
      do.call('[',l2)
    }
    ,ii=1:nn
  )
  return(list(preds=preds,tab=tt))
}

############
# function for initializing means
initMeans <- function(conVar,method=conInitMethod,numClust) {
  if (method=='sample') {
    return(
      sapply(conVar,function(xx) sample(xx,size=numClust,replace=TRUE))
    )
  } else if (method=='runif') {
    ranges <- sapply(conVar,range)
    return(
      apply(ranges,2,function(xx) runif(numClust,min=xx[1],max=xx[2]))
    )
  } else {
    stop('Unrecognized mean initialization method')
  }
}


######################
# remove bumps (i.e. make sure xx is nondecreasing)
# Version 1 starts from right side moving left, replaces any decrease
# with the previous value
rmBump1 <- function(xx) {
  for (i in (length(xx)-1):1) {
    if (xx[i] < xx[i+1]) xx[i] <- xx[i+1]
  }
  return(xx)
}
# Version 2 starts from the left side moving right, replaces any increase
# with the earlier value
rmBump2 <- function(xx) {
  for (i in 2:length(xx)) {
    if (xx[i] > xx[i-1]) xx[i] <- xx[i-1]
  }
  return(xx)
}

######################
# Calculate weighted Euclidean distances from a set of n points
# in p-space to a set of k means
distPointsToMeans <- function(pts,myMeans,wgts) {
  ppDim <- ncol(pts)
  if (ppDim != ncol(myMeans)) stop('Dimensionality of pts and myMeans must be equal')
  if (ppDim != length(wgts)) stop('Dimensionality of pts must equal # weights')
  kkMean <- nrow(myMeans)
  
  sapply(
    1:kkMean
    ,function(kk) {
      diff_k = sapply(
        1:ppDim
        ,function(jj) {wgts[jj]*(pts[,jj] - myMeans[kk,jj])}
      )
      apply(diff_k,1,function(ii) sqrt(sum(ii^2)))
    }
  )
}

##########################
# Rcpp implementation of the above function
dptmCpp <- function(pts,myMeans,wgts) {
  pts <- as.matrix(pts)
  ppDim <- ncol(pts)
  if (ppDim != ncol(myMeans)) stop('Dimensionality of pts and myMeans must be equal')
  if (ppDim != length(wgts)) stop('Dimensionality of pts must equal number of weights')
  kkMean <- nrow(myMeans)
  nn <- nrow(pts)
  
  dptm(pts,myMeans,wgts,ppDim,kkMean,nn)
}


######################
# psort: poor man's approximation to sample quantile
# Underestimates the quantile, but this bias decreases
# as sample size increases.
# Used in radialKDE function below, where we don't need
# an exact quantile anyway.
# Not used since gain in computation 
# time is minimal.
psort <- function(xx,pp) {
  nn = length(xx)
  jj = floor(nn*pp)
  if (jj <= 1) return(min(xx))
  return(sort(xx,partial=jj)[jj])
}

######################
# Estimate density of radii
# Takes a vector of distances to mean, estimates radial KDE
# then evaluates at given vector of points
# pdim is the number of continuous variables used
# returnFun causes a resampling function to be returned
radialKDE <- function(radii,evalPoints,pdim,returnFun=FALSE) {
  MAXDENS <- 1
  # Note using a chosen constant for bw reduces time by about 7%
  radialBW <- bw.nrd0(radii)
  radKDE <- KernSmooth::bkde(
    x = radii
    ,kernel = "normal"
    ,bandwidth = radialBW
    #   ,range.x = c(0, max(radii))
    ,range.x = c(0,max(evalPoints))
  )
  
  # remove any zero and negative density estimates
  newY <- radKDE$y
  nonnegTF <- newY > 0
  if (any(!nonnegTF)) {
    minPos <- min(newY[nonnegTF])
    newY[!nonnegTF] <- minPos/100
  }
  
  # at bottom 5th percentile, replace with line through (0,0) and (q05,f(q05)).
  # This removes substantial variability in output near zero
  quant05 <- quantile(x = radKDE$x, probs = 0.05)
  #quant05 <- psort(xx = radKDE$x, pp = 0.05)
  coordsLtQ05 <- radKDE$x < quant05
  maxPt <- max(which(coordsLtQ05))
  newY[coordsLtQ05] <- radKDE$x[coordsLtQ05] * (newY[maxPt]/radKDE$x[maxPt]) # y = x * sl
  
  # radial Jacobian transformation; up to proportionality constant
  radY <- c(0,newY[-1] / radKDE$x[-1]^(pdim-1))
  
  # replace densities over MAXDENS with MAXDENS
  overMax <- radY > MAXDENS
  radY[overMax] <- MAXDENS
  
  # remove bumps (i.e. make sure nondecreasing)
  # (eventually may want to replace with better strategy)
  #  radY <- rmBump1(radY)
  
  # normalize to area 1
  binwidthX <- diff(radKDE$x[1:2])
  densR <- radY/(binwidthX * sum(radY))
  
  # now create resampling function
  resampler <- approxfun(x=radKDE$x, y=densR,rule=1:2, method='linear')
  kdes <- resampler(evalPoints)
  if (!returnFun) {
    resampler <- NULL
  }
  
  #return(list(kdes=resampler(evalPoints),resampler=resampler))
  return(list(kdes=kdes,resampler=resampler))
}

# Main kamila function
#
# Note conVar and catFactor must be dataframes
# Categorical initialization using random draws from a dirichlet(alpha=rep(1,nlev))
# 
# Weighting scheme:
#    If no weights are desired, set all weights to 1 (the default setting).
#    Let a_1, ..., a_p denote the weights for p continuous variables.
#    Let b_1, ..., b_q denote the weights for q categorical variables.
#    Currently, continuous weights are applied during the calculation
#    of Euclidean distance, as:
#        sqrt( a_1^2(x1-y1)^2 + ... + a_p^2(xp-yp)^2 )
#    Categorical weights are applied to the log-likelihoods obtained
#    by the level probabilities given cluster membership as:
#        logLikCat_k = b_1*log(P(lev_1|clust_k)) + ... + b_q*log(P(lev_q|clust_k)) for each k
#    Total log likelihood for the kth cluster is obtained by weighting
#    the single continuous log-likelihood by the mean of all continuous
#    weights plus logLikCat_k:
#        total_k = mean(a_1 + ... + a_p)*logLikCon_k + logLikCat_k
#    Note that weights between 0 and 1 are admissible; weights equal to zero
#    completely remove a variable's influence on the clustering; weights equal
#    to 1 leave a variable's contribution unchanged. Weights between 0 and 1
#    may not be comparable across continuous and categorical variables.

#' KAMILA clustering of mixed-type data.
#' 
#' KAMILA is an iterative clustering method that equitably balances the
#' contribution of continuous and categorical variables.
#' 
#' KAMILA (KAy-means for MIxed LArge data sets) is an iterative clustering
#' method that equitably balances the contribution of the continuous and
#' categorical variables. It uses a kernel density estimation technique to
#' flexibly model spherical clusters in the continuous domain, and uses a
#' multinomial model in the categorical domain.
#' 
#' Weighting scheme: If no weights are desired, set all weights to 1 (the
#' default setting). Let a_1, ..., a_p denote the weights for p continuous 
#' variables. Let b_1, ..., b_q denote the weights for q categorical variables.
#'  Currently, continuous weights are applied during the calculation of 
#' Euclidean distance, as:
#        sqrt( a_1^2(x1-y1)^2 + ... + a_p^2(xp-yp)^2 )
#' Categorical weights are applied to the log-likelihoods obtained by the 
#' level probabilities given cluster membership as:
#        logLikCat_k = b_1*log(P(lev_1|clust_k)) + ... + b_q*log(P(lev_q|clust_k)) for each k
#' Total log likelihood for the kth cluster is obtained by weighting the
#' single continuous log-likelihood by the mean of all continuous weights
#' plus logLikCat_k:
#        total_k = mean(a_1 + ... + a_p)*logLikCon_k + logLikCat_k
#' Note that weights between 0 and 1 are admissible; weights equal to zero
#' completely remove a variable's influence on the clustering; weights equal
#' to 1 leave a variable's contribution unchanged. Weights between 0 and 1
#' may not be comparable across continuous and categorical variables.
#' @export
#' @param conVar A data frame of continuous variables.
#' @param catFactor A data frame of factors.
#' @param numClust The number of clusters returned by the algorithm.
#' @param numInit The number of initializations used.
#' @param conWeights A vector of continuous weights for the continuous variables.
#' @param catWeights A vector of continuous weights for the categorical variables.
#' @param maxIter The maximum number of iterations in each run.
#' @param conInitMethod Character: he method used to initialize each run.
#' @param catBw The bandwidth used for the categorical kernel.
#' @param verbose Logical: whether detailed results should be printed and returned.
#' @return A list of the following results objects:
#' \item{finalMemb}{A numeric vector with cluster assignment indicated by integer.}
#' \item{numIter}{}
#' \item{finalLogLik}{The pseudo log-likelihood of the returned clustering.}
#' \item{finalObj}{}
#' \item{finalCenters}{}
#' \item{finalProbs}{}
#' \item{input}{Object with the given input parameter values.}
#' \item{verbose}{An optionally returned object with more detailed information.}

kamila <- function(
  conVar
  ,catFactor
  ,numClust
  ,numInit
  ,conWeights = rep(1,ncol(conVar))
  ,catWeights = rep(1,ncol(catFactor))
  ,maxIter = 25
  ,conInitMethod = 'runif'
  ,catBw = 0.025
  ,verbose = FALSE
) {
  
  
  # Deprecated option
  returnResampler = FALSE
  
  # (0) extract data characteristics, checks
  if (!is.data.frame(conVar) | !is.data.frame(catFactor)) {
    stop("Input datasets must be dataframes")
  }
  if (max(c(conWeights,catWeights)) > 1 | min(c(conWeights,catWeights)) < 0) stop('Weights must be in [0,1]')
  
  numObs <- nrow(conVar)
  numConVar <- ncol(conVar)
  if (nrow(catFactor) != numObs) {
    stop("Number of observations in con and cat vars don't match")
  }
  numCatVar <- ncol(catFactor)
  
  numLev <- sapply(catFactor,function(xx) length(levels(xx)))
  
  # for later use, convert to numeric matrix by level codes
  catFactorNumeric <- sapply(catFactor,as.numeric,simplify=TRUE)
  
  numIterVect <- rep(NaN,numInit)
  totalLogLikVect <- rep(NaN,numInit)
  catLogLikVect <- rep(NaN,numInit)
  winDistVect <- rep(NaN,numInit)
  totalDist <- sum(dptmCpp(
    pts=conVar
    ,myMeans=matrix(colMeans(conVar),nrow=1)
    ,wgts=conWeights #rep(1,numConVar)
  ))
  objectiveVect <- rep(NaN,numInit)
  
  # for verbose output, list of memberships for each init, iter
  if (verbose) membLongList = rep(list(list()),numInit)
  
  
  # Apply continuous weighting vector directly to variables
  #for (colInd in 1:ncol(conVar)) conVar[,colInd] <- conVar[,colInd] * conWeights[colInd]
  
  # (1) loop over each initialization
  for (init in 1:numInit) {
    # initialize means; matrix size numClust X numConVar
    means_i <- initMeans(conVar=conVar,method=conInitMethod,numClust=numClust)
    #print(means_i)
    
    # initialize probabilities; list length numCatVar of numClust X numLev matrices
    # generate each level prob from dirichlet(alpha=rep(1,nlev))
    # These are P( level | clust )
    logProbsCond_i <- lapply(
      numLev
      ,function(xx) {
        matrix(
          data=log(gtools::rdirichlet(n=numClust,alpha=rep(1,xx)))
          ,nrow=numClust
          ,dimnames=list(clust=1:numClust,level=1:xx)
        )
      }
    )
    
    
    
    #print('logProbsCond_i')
    #print(logProbsCond_i)
    
    # initialize structures for iterative procedure
    membOld <- membNew <- rep(0,numObs)
    numIter <- 0
    degenerateSoln <- FALSE
    
    # Loop until convergence
    # loop for minimum two iterations, while all memberships are NOT UNchanged, while still under max # iter
    while(
      ( (numIter < 3) | !all(membOld == membNew))
      & (numIter < maxIter)
    ) { 
      
      
      numIter <- numIter + 1
      #print(paste('Iteration',numIter))
      
      #print('conVar')
      #str(conVar)
      #print('means_i')
      #str(means_i)
      #print('conWeights')
      #str(conWeights)
      #print('initMeans Structure')
      #str(initMeans(conVar,conInitMethod,2))
      #str(initMeans(cbind(conVar,conVar),conInitMethod,2))
      
      #2 calc weighted euclidean distances to means
      # result is numObs X numClust matrix
      #str(means_i)
      #str(conVar)
      
      dist_i <- dptmCpp(pts=conVar,myMeans=means_i,wgts=conWeights)
      #dist_i <- dptmCpp(pts=conVar,myMeans=means_i,wgts=rep(1,numConVar)) # no weights
      #print('Distance X cluster matrix')
      #print(head(dist_i,n=15))
      
      #3 extract min distances
      #minDist_i <- apply(dist_i,1,min)
      minDist_i <- rowMin(dist_i)
      
      #4 RKDE of min distances
      #5 Evaluate all distances with kde
      logDistRadDens_vec <- log(
        radialKDE(
          radii=minDist_i
          ,evalPoints=c(dist_i)
          ,pdim=numConVar
          ,returnFun = returnResampler
        )$kdes
      )
      logDistRadDens_i <- matrix(
        logDistRadDens_vec
        ,nrow=numObs
        ,ncol=numClust
      )
      #print('log Density X cluster')
      #print(head(logDistRadDens_i))
      if (FALSE) { #(verbose & numIter == 1) {
        hist(minDist_i)
        plot(c(dist_i),c(exp(logDistRadDens_i)))
        myXs <- seq(0,10,by=0.01)
        plot(myXs,log(radialKDE(radii=minDist_i,evalPoints=myXs,pdim=numConVar)$kdes),main='Radial Density')
      }
      
      
      #6 categorical likelihoods
      # this step takes the log probs for each level of each
      # categorical variable (list length Q, elements k X l_q)
      # and creates a list of length Q, each element is n X k
      # giving the log prob for nth observed level for variable q, cluster k
      
      # old method
      #catLikWeight <- function(ind) catWeights[ind]*t(logProbsCond_i[[ind]][,as.numeric(catFactor[,ind])])
      #individualLogProbs = lapply(
      #  X = 1:numCatVar
      # ,FUN = catLikWeight
      #)
      
      # New Rcpp method
      individualLogProbs <- getIndividualLogProbs(
        catFactorNum = catFactorNumeric
        ,catWeights = catWeights
        ,logProbsCond_i = logProbsCond_i
      )
      
      #print('exp(individualLogProbs)')
      #print(lapply(individualLogProbs,function(xx) head(exp(xx))))
      
      #7 log likelihood eval for each point, for each of k clusters (WEIGHTED)
      # output is n X k matrix of log likelihoods
      catLogLiks <- Reduce(f='+',x=individualLogProbs)
      
      # sumMatList is cpp function that replaces the above; not clear if faster than Reduce
      #catLogLiks <- sumMatList(individualLogProbs)
      
      # scaling relative con to cat likelihood by weights
      #allLogLiks <- mean(conWeights)*logDistRadDens_i + catLogLiks
      
      # NOT scaling relative to weights
      allLogLiks <- logDistRadDens_i + catLogLiks
      #print('allLogLiks')
      #print(head(allLogLiks))
      
      
      #8 partition data into clusters
      membOld <- membNew
      
      #membNew <- apply(allLogLiks,1,which.max)
      membNew <- rowMaxInds(allLogLiks)
      
      #print('membNew')
      #str(membNew)
      
      #9 calculate new means: k X p matrix
      #means_i <- as.matrix(aggregate(x=conVar,by=list(membNew),FUN=mean)[,-1])
      means_i <- aggregateMeans(
        conVar = as.matrix(conVar)
        ,membNew = membNew
        ,kk = numClust
      )
      
      #10.1 new joint probability table, possibly with kernel estimator
      #10.2 Also calculate conditional probabilities: P(lev | clust)
      
      # New Rcpp implementation
      jointProbsList <- jointTabSmoothedList(catFactorNumeric,membNew,numLev,catBw,kk=numClust)
      logProbsCond_i <- lapply(jointProbsList,FUN=function(xx) log(xx/rowSums(xx)))
      
      ### Old approach:
      #logProbsCond_i <- lapply(
      #  X = as.list(catFactor)
      # ,FUN = function(xx) {
      #    jointTab <- myCatKern(
      #      kdat=data.frame(clust=membNew,level=xx)
      #     ,bw=catBw
      #     ,tabOnly=TRUE
      #    )
      #    #jointTab <- table(membNew,xx) # without kernel
      #    return( log(jointTab / rowSums(jointTab)) )
      #  }
      #)
      
      
      #print('exp(logProbsCond_i)')
      #print(lapply(logProbsCond_i,function(xx) head(exp(xx))))
      
      # print current plot and metrics, if requested
      if (FALSE) { #(verbose) {
        #catLikTmp <- sum(apply(catLogLiks,1,max))
        catLikTmp <- sum(rowMax(catLogLiks))
        winDistTmp <- sum(minDist_i)
        winToBetRatTmp <- winDistTmp / (totalDist - winDistTmp)
        if (winToBetRatTmp < 0) winToBetRatTmp <- 100
        objectiveTmp <- winToBetRatTmp * catLikTmp 
        
        plot(conVar[,1],conVar[,2],col=membNew,main=paste('Init',init,'; Iter',numIter) )
        points(means_i,pch=18,col='blue',cex=2)
        legend('topleft',legend=c(
          paste('catLik =',round(catLikTmp ,2))
          ,paste('w/b dist =',round(winDistTmp / (totalDist - winDistTmp),2))
          ,paste('objective =',round(objectiveTmp,2))
        ))
      }
      
      # store every solution if verbose=true
      if (verbose) {
        membLongList[[init]][[numIter]] <- membOld
      }
      
      #11 test for degenerate solution, calculate total log-likelihood (WEIGHTED)
      if(length(unique(membNew)) < numClust) {
        degenerateSoln <- TRUE 
        break
      }
    }
    
    
    # Store log likelihood for each initialization
    if (degenerateSoln) {
      totalLogLikVect[init] <- -Inf
    } else {
      #totalLogLikVect[init] <- sum(apply(allLogLiks,1,max))
      totalLogLikVect[init] <- sum(rowMax(allLogLiks))
    }
    
    # store num init
    numIterVect[init] <- numIter
    
    # other useful internal measures of cluster quality
    catLogLikVect[init] <- sum(rowMax(catLogLiks))
    winDistVect[init] <- sum(minDist_i)
    winToBetRat <- winDistVect[init] / (totalDist - winDistVect[init])
    if (winToBetRat < 0) winToBetRat <- 100
    objectiveVect[init] <- winToBetRat * catLogLikVect[init]
    
    # Store current solution if objective beats all others
    if (
      (init == 1)
      | (init > 1 & objectiveVect[init] > max(objectiveVect[1:(init-1)]))
    ) {
      finalLogLik <- totalLogLikVect[init]
      finalObj <- objectiveVect[init]
      finalMemb <- membNew
      finalCenters <- matrix(
        data=as.matrix(means_i)
        ,nrow=nrow(means_i)
        ,ncol=numConVar
        ,dimnames=list(
          cluster=paste('Clust',1:nrow(means_i))
          ,variableMean=paste('Mean',1:numConVar))
      )
      # rescale by weights
      #    finalCenters <- finalCenters * matrix(rep(1/conWeights,times=numClust),byrow=TRUE,nrow=numClust)
      finalProbs <- lapply(logProbsCond_i,exp)
      names(finalProbs) <- paste('Categorical Variable',1:numCatVar)
      finalClustSize <- table(membNew)
    }
    
    # store last membership also
    if (verbose) membLongList[[init]][[numIter+1]] <- membNew
    
  }
  
  
  #12 Prepare output data structure
  if (verbose) {
    optionalOutput <- list(
      totalLogLikVect =totalLogLikVect
      ,catLogLikVect = catLogLiks
      ,winDistVect = winDistVect
      ,totalDist = totalDist
      ,objectiveVect = objectiveVect
      ,membLongList = membLongList
    )
  } else {
    optionalOutput <- list()
  }
  
  inputList <- list(
    conVar = conVar
    ,catFactor = catFactor
    ,numClust = numClust
    ,numInit = numInit
    ,conWeights = conWeights
    ,catWeights = catWeights
    ,maxIter = maxIter
    ,conInitMethod = conInitMethod
    ,catBw = catBw
    ,verbose = verbose
  )
  
  
  return(
    list(
      finalMemb=as.numeric(finalMemb)
      ,numIter=numIterVect
      ,finalLogLik=finalLogLik
      ,finalObj=finalObj
      ,finalCenters=finalCenters
      ,finalProbs=finalProbs
      ,input=inputList
      ,verbose=optionalOutput
    )
  )
}



# Function designed for cyclical variables, e.g. day of week.
# Recodes to equidistant points on the unit circle.
#############################################
### duplicated in predictionStrength.r ######
#### correct this redundancy ################
#############################################
cyclicalCoding <- function(invar) {
  minDetected <- min(invar)
  maxDetected <- max(invar)
  tmp <- 2*invar*pi/(maxDetected+1-minDetected)
  return(cbind(cos(tmp),sin(tmp)))
}

# Function to take kamila results object and classify new points
# Note newData is a list with two elements, continuous data and dataframe of 
# factors of categorial varables.
classifyKamila <- function(obj,newData) {
  if (length(newData) == 3) {
    cyclicRecoded <- as.data.frame(lapply(newData[[3]],cyclicalCoding))
    newCon <- as.data.frame(cbind(newData[[1]],cyclicRecoded))
  } else if (length(newData) == 2) {
    newCon <- newData[[1]]
  } else {
    stop('Error in function classifyKamila: newData must be list of length 2 or 3')
  }
  newCatFactor <- newData[[2]]
  
  ########################################
  # 1) reconstruct classification rule
  ########################################
  # calculate (weighted) distance to means
  distances <- with(obj,dptmCpp(
    pts=input$conVar
    ,myMeans=finalCenters
    ,wgts=input$conWeights
  ))
  minDistances <- apply(distances,1,min)
  #print(minDistances)
  
  ########################################
  # 2) classify new points
  ########################################
  
  # calculate distances to each mean from each new point
  newDistances <- with(obj,dptmCpp(
    pts = newCon
    ,myMeans=finalCenters
    ,wgts=input$conWeights
  ))
  #print(newDistances)
  
  # calculate continuous log density KD estimates
  logRadDens <- matrix(
    log(radialKDE(radii=minDistances,evalPoints=c(newDistances),pdim=ncol(newCon))$kdes)
    ,nrow=nrow(newCon)
    ,ncol=nrow(obj$finalCenters)
  )
  #print('logRadDens')
  #print(logRadDens)
  
  # calculate categorical probabilities
  logClustProbs <- lapply(obj$finalProbs,log)
  logCatKProbs = with(obj,lapply(
    X = 1:ncol(newCatFactor)
    ,FUN = function(ind) {
      input$catWeights[ind]*t(logClustProbs[[ind]][,as.numeric(newCatFactor[,ind])])
    }
  )) # this is a list of q matrices, each n x k; with elements log-likelihood
  #print('logCatKProbs')
  #print(logCatKProbs)
  
  catLogLiks <- Reduce(f='+',x=logCatKProbs)
  #print('catLogLiks')
  #print(catLogLiks)
  
  # maximum "likelihood" classification
  combinedLogLik <- logRadDens + catLogLiks
  #print('combinedLogLik')
  #print(combinedLogLik)
  
  membership <- as.numeric(apply(combinedLogLik,1,which.max))
  #print('membership')
  #print(membership)
  
  #pairs(cbind(obj$verbose$catLogLikVect,catLogLiks))
  return(membership)
}
