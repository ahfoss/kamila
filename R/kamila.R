

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
#' @importFrom KernSmooth bkde
#' @importFrom gtools rdirichlet
#' @importFrom abind abind

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
    offCounts <- do.call(abind::abind,l1)
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
initMeans <- function(conVar,method,numClust) {
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
# with the previous value.
# Not used; doesn't appear to be necessary.
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
##########################
#' Calculate distances from a set of points to a set of centroids
#'
#' A function that calculates a NxM matrix of distances between a NxP set of
#' points and a MxP set of points.
#'
#' @export
#' @param pts A matrix of points
#' @param myMeans A matrix of centroids, must have same ncol as pts
#' @param wgts A Px1 vector of variable weights
#' @return A MxP matrix of distances
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
#' @importFrom stats bw.nrd0 approxfun quantile
radialKDE <- function(radii,evalPoints,pdim,returnFun=FALSE) {
  MAXDENS <- 1
  # Note using a chosen constant for bw reduces time by about 7%
  radialBW <- bw.nrd0(radii)
  radKDE <- bkde(
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
  # not used currently
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
#' Estimating the number of clusters: Default is no estimation method. Setting
#' calcNumClust to 'ps' uses the prediction strength method of Tibshirani & 
#' Walther (J. of Comp. and Graphical Stats. 14(3), 2005). There is no perfect
#' method for estimating the number of clusters; PS tends to give a smaller
#' number than, say, BIC based methods for large sample sizes. The user must
#' specify the number of cross-validation runs and the threshold for
#' determining the number of clusters. The smaller the threshold, the larger
#' the number of clusters selected.
#' @export
#' @importFrom stats runif sd setNames
#' @param conVar A data frame of continuous variables.
#' @param catFactor A data frame of factors.
#' @param numClust The number of clusters returned by the algorithm.
#' @param numInit The number of initializations used.
#' @param conWeights A vector of continuous weights for the continuous variables.
#' @param catWeights A vector of continuous weights for the categorical variables.
#' @param maxIter The maximum number of iterations in each run.
#' @param conInitMethod Character: The method used to initialize each run.
#' @param catBw The bandwidth used for the categorical kernel.
#' @param verbose Logical: Whether detailed results should be printed and returned.
#' @param calcNumClust Character: Method for selecting the number of clusters.
#' @param numPredStrCvRun Numeric: Number of CV runs for prediction strength method. Ignored unless calcNumClust == 'ps'
#' @param predStrThresh Numeric: Threshold for prediction strength method. Ignored unless calcNumClust == 'ps'
#' @return A list with the following results objects:
#' \item{finalMemb}{A numeric vector with cluster assignment indicated by integer.}
#' \item{numIter}{}
#' \item{finalLogLik}{The pseudo log-likelihood of the returned clustering.}
#' \item{finalObj}{}
#' \item{finalCenters}{}
#' \item{finalProbs}{}
#' \item{input}{Object with the given input parameter values.}
#' \item{nClust}{An object describing the results of selecting the number of clusters, empty if calcNumClust == 'none'.}
#' \item{verbose}{An optionally returned object with more detailed information.}
#' @examples
#' # Generate toy data set with poor quality categorical variables and good
#' # quality continuous variables.
#' set.seed(1)
#' dat <- genMixedData(200, nConVar = 2, nCatVar = 2, nCatLevels = 4,
#'   nConWithErr = 2, nCatWithErr = 2, popProportions = c(.5, .5),
#'   conErrLev = 0.3, catErrLev = 0.8)
#' catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
#' conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)
#'
#' kamRes <- kamila(conDf, catDf, numClust = 2, numInit = 10)
#'
#' table(kamRes$finalMemb, dat$trueID)
#' @references Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and Hadoop. Journal of Statistical Software, 83(13). 2018. doi: 10.18637/jss.v083.i13

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
  ,calcNumClust = 'none'
  ,numPredStrCvRun = 10
  ,predStrThresh = 0.8
) {

  if (calcNumClust == 'none') {
    if (length(numClust) != 1) {
      stop('Input parameter numClust must be length 1 if calcNumClust == "none"')
    }
    # main function
  
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
        
        # New Rcpp method
        individualLogProbs <- getIndividualLogProbs(
          catFactorNum = catFactorNumeric
          ,catWeights = catWeights
          ,logProbsCond_i = logProbsCond_i
        )
        
        #7 log likelihood eval for each point, for each of k clusters (WEIGHTED)
        # output is n X k matrix of log likelihoods
        catLogLiks <- Reduce(f='+',x=individualLogProbs)
        
        # sumMatList is cpp function that replaces the above; not clear if faster than Reduce
        #catLogLiks <- sumMatList(individualLogProbs)
        
        # scaling relative con to cat likelihood by weights
        #allLogLiks <- mean(conWeights)*logDistRadDens_i + catLogLiks
        
        # NOT scaling relative to weights
        allLogLiks <- logDistRadDens_i + catLogLiks
        
        #8 partition data into clusters
        membOld <- membNew
        
        #membNew <- apply(allLogLiks,1,which.max)
        membNew <- rowMaxInds(allLogLiks)
        
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
        
        # print current plot and metrics, if requested
        if (FALSE) { #(verbose) {
          catLikTmp <- sum(rowMax(catLogLiks))
          winDistTmp <- sum(dist_i[cbind(1:numObs,membNew)])
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
        totalLogLikVect[init] <- sum(rowMax(allLogLiks))
      }
      
      # store num init
      numIterVect[init] <- numIter
      
      # other useful internal measures of cluster quality
      catLogLikVect[init] <- sum(rowMax(catLogLiks))
      winDistVect[init] <- sum(dist_i[cbind(1:numObs, membNew)])
      winToBetRat <- winDistVect[init] / (totalDist - winDistVect[init])
      if (winToBetRat < 0) winToBetRat <- 100
      # Note catLogLik is negative, larger is better
      # Note WSS/BSS is positive, with smaller better
      # Thus, we maximize their product.
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
        ,nClust=list()
      )
    )
  } else if (calcNumClust == 'ps') {
    # Recursive call to kamila implementing prediction strength method.

    # Test that numClust is an integer vector.
    allEqInteger <- all(numClust == as.integer(numClust))
    if ( is.na(allEqInteger) ||
         !allEqInteger ||
         !all(sapply(numClust, is.numeric)) ||
         length(numClust) != length(unique(numClust)) ) {
      stop('Input parameter numClust must be a vector of unique integers
         if input parameter calcNumClust == "ps"')
    }
    if (length(numClust) == 1) {
      warning('Input parameter numClust is a scalar; the prediction strength 
        method is probably not appropriate or desired')
    }
    if (any(numClust > floor(nrow(conVar)/2))) {
      stop('The number of clusters in input parameter numClust cannot exceed
        one-half of the sample size.')
    }

    # Test that predStrThresh is within (0,1).
    if ( length(predStrThresh) != 1 ||
         is.na(predStrThresh) ||
         predStrThresh <= 0 ||
         predStrThresh >= 1 ) {
      stop('Input parameter predStrThresh must be scalar in (0,1)')
    }

    # Test that numPredStrCvRun is a valid number of cv runs.
    if ( length(numPredStrCvRun) != 1 ||
         is.na(numPredStrCvRun) || 
         numPredStrCvRun != as.integer(numPredStrCvRun) ||
         numPredStrCvRun < 1 ) {
      stop('Input parameter numPredStrCvRun must be a positive integer.')
    }

    psCvRes <- matrix(
      NaN,
      nrow = length(numClust),
      ncol = numPredStrCvRun,
      dimnames = list(
        nClust = numClust,
        CVRun = paste('Run', 1:numPredStrCvRun)
      )
    )
    
    # Implement CV procedure
    numObs <- nrow(conVar)
    numInTest <- floor(numObs/2)
    for (cvRun in 1:numPredStrCvRun) {
      for (ithNcInd in 1:length(numClust)) {
        # generate cv indices
        testInd <- sample(numObs, size=numInTest, replace=FALSE)

        # cluster test data
        testClust <- kamila(
           conVar = conVar[testInd,,drop=FALSE],
           catFactor = catFactor[testInd,,drop=FALSE],
           numClust = numClust[ithNcInd],
           numInit = numInit,
           conWeights = conWeights,
           catWeights = catWeights,
           maxIter = maxIter,
           conInitMethod = conInitMethod,
           catBw = catBw,
           verbose = FALSE
        )
        
        # cluster training data
        trainClust <- kamila(
          conVar = conVar[-testInd,,drop=FALSE],
          catFactor = catFactor[-testInd,,drop=FALSE],
          numClust = numClust[ithNcInd],
          numInit = numInit,
          conWeights = conWeights,
          catWeights = catWeights,
          maxIter = maxIter,
          conInitMethod = conInitMethod,
          catBw = catBw,
          verbose = FALSE
        )
        
        # Generate a list of indices for each cluster within test data.
        # Note there are two index systems floating around:
	# Indices of the full data set slicing out test set (testInd)
        # and indices of cluster membership within test data (testIndList).
	# If the length of the outputs are equal, mapply defaults to returning
	# a matrix unless default of SIMPLIFY parameter is changed to FALSE.
        testIndList <- mapply(
          x = 1:numClust[ithNcInd],
          function(x) which(testClust$finalMemb == x),
	  SIMPLIFY = FALSE
        )
        
        # Allocate test data based on training clusters.
        teIntoTr <- classifyKamila(
          trainClust,
          list(conVar[testInd,,drop=FALSE],catFactor[testInd,,drop=FALSE])
        )
        
        # Initialize D matrix.
        dMat <- matrix(
          NaN,
          nrow = numInTest,
          ncol = numInTest
        )
        
        # Calculate D matrix.
        # replace with RCpp
        for (i in 1:(numInTest-1)) {
          for (j in (i+1):numInTest) {
            dMat[i,j] <- teIntoTr[i] == teIntoTr[j]
          }
        }
        
        # Initialize proportions to zero.
        psProps <- rep(0,numClust[ithNcInd])

        # Calculate prediction strength proportions.
        # replace with RCpp
        for (cl in 1:numClust[ithNcInd]) {
          clustN <- length(testIndList[[cl]])
	  if (clustN > 1) {
##################################
# Construct dMat separately within each cluster
# initialize dMat_cl
# dMat_cl[i,j] <- teIntoTr[testIndList[[cl]][i]] == teIntoTr[testIndList[[cl]][j]]
##################################
            for (i in 1:(clustN-1)) {
              for (j in (i+1):clustN) {
                psProps[cl] <- psProps[cl] + dMat[testIndList[[cl]][i], testIndList[[cl]][j]]
              }
            }
	  }
	  if (clustN < 2) {
	    # if cluster size is 1 or zero, not applicable
	    psProps[cl] <- NA
	  } else {
            # * 2 since only upper triangle of dMat is used.
            psProps[cl] <- psProps[cl] / (clustN*(clustN-1)) * 2
	  }
        }

        # Calculate and update prediction strength results.
	# Remove NaN/NAs for empty clusters.
	# min(...,na.rm=T) returns Inf if entire vector is NA.
	# Workaround is simple but WEIRD behavior.
        psCvRes[ithNcInd, cvRun] <- ifelse(
	  test = all(is.na(psProps)),
	  yes = NA,
	  no = min(psProps,na.rm=TRUE)
	)

      } # end clusters
    } # end cv runs

    # Calculate CV estimate of prediction strength for each cluster size.
    avgPredStr    <- apply(psCvRes,1,mean,na.rm=TRUE)
    stdErrPredStr <- apply(psCvRes,1,sd,na.rm=TRUE) / sqrt(numPredStrCvRun)

    # Calculate final number of clusters: largest # clust such that avg+sd
    # score is above the threshold.
    psValues <- avgPredStr + stdErrPredStr
    clustAboveThresh <- psValues > predStrThresh

    if (all(!clustAboveThresh)) {
      warning('No cluster size is above prediction strength threshold.
        Consider lowering the ps threshold; returning the cluster size
	corresponding to the highest ps value.')
      kfinal <- numClust[which.max(psValues)]
    } else {
      kfinal <- max(numClust[clustAboveThresh])
    }

    # Do the final clustering.
    outRes <- kamila(
      conVar = conVar,
      catFactor = catFactor,
      numClust=kfinal,
      numInit = numInit,
      conWeights = conWeights,
      catWeights = catWeights,
      maxIter = maxIter,
      conInitMethod = conInitMethod,
      catBw = catBw,
      verbose = FALSE
    )
    outRes$nClust <- list(
      bestNClust = kfinal,
      psValues = setNames(psValues, numClust),
      avgPredStr = avgPredStr,
      stdErrPredStr = stdErrPredStr,
      psCvRes = psCvRes
    )

    return(outRes)

  } else {
    stop('Currently calcNumClust must be either "none" or "ps"')
  }
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

#' Classify new data into existing KAMILA clusters
#'
#' A function that classifies a new data set into existing KAMILA clusters
#' using the output object from the kamila function.
#'
#' A function that takes obj, the output from the kamila function, and newData,
#' a list of length 2, where the first element is a data frame of continuous
#' variables, and the second element is a data frame of categorical factors.
#' Both data frames must have the same format as the original data used
#' to construct the kamila clustering.
#' @export
#' @param obj An output object from the kamila function.
#' @param newData A list of length 2, with first element a data frame of continuous variables, and second element a data frame of categorical factors.
#' @return An integer vector denoting cluster assignments of the new data points.
#' @examples
#' # Generate toy data set
#' set.seed(1234)
#' dat1 <- genMixedData(400, nConVar = 2, nCatVar = 2, nCatLevels = 4,
#'   nConWithErr = 2, nCatWithErr = 2, popProportions = c(.5,.5),
#'   conErrLev = 0.2, catErrLev = 0.2)
#' # Partition the data into training/test set
#' trainingIds <- sample(nrow(dat1$conVars), size = 300, replace = FALSE)
#' catTrain <- data.frame(apply(dat1$catVars[trainingIds,], 2, factor), stringsAsFactors = TRUE)
#' conTrain <- data.frame(scale(dat1$conVars)[trainingIds,], stringsAsFactors = TRUE)
#' catTest <- data.frame(apply(dat1$catVars[-trainingIds,], 2, factor), stringsAsFactors = TRUE)
#' conTest <- data.frame(scale(dat1$conVars)[-trainingIds,], stringsAsFactors = TRUE)
#' # Run the kamila clustering procedure on the training set
#' kamilaObj <- kamila(conTrain, catTrain, numClust = 2, numInit = 10)
#' table(dat1$trueID[trainingIds], kamilaObj$finalMemb)
#' # Predict membership in the test data set
#' kamilaPred <- classifyKamila(kamilaObj, list(conTest, catTest))
#' table(dat1$trueID[-trainingIds], kamilaPred)
#' @references Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and Hadoop. Journal of Statistical Software, 83(13). 2018. doi: 10.18637/jss.v083.i13
classifyKamila <- function(obj, newData) {
  #if (length(newData) == 3) {
  #  cyclicRecoded <- as.data.frame(lapply(newData[[3]],cyclicalCoding))
  #  newCon <- as.data.frame(cbind(newData[[1]],cyclicRecoded))
  #} else
  if (length(newData) == 2) {
    newCon <- newData[[1]]
  } else {
    stop('Error in function classifyKamila: newData must be list of length 2')
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
  
  ########################################
  # 2) classify new points
  ########################################
  
  # calculate distances to each mean from each new point
  newDistances <- with(obj,dptmCpp(
    pts = newCon
    ,myMeans=finalCenters
    ,wgts=input$conWeights
  ))
  
  # calculate continuous log density KD estimates
  logRadDens <- matrix(
    log(radialKDE(radii=minDistances,evalPoints=c(newDistances),pdim=ncol(newCon))$kdes)
    ,nrow=nrow(newCon)
    ,ncol=nrow(obj$finalCenters)
  )
  
  # calculate categorical probabilities
  logClustProbs <- lapply(obj$finalProbs,log)
  logCatKProbs = with(obj,lapply(
    X = 1:ncol(newCatFactor)
    ,FUN = function(ind) {
      input$catWeights[ind]*t(logClustProbs[[ind]][,as.numeric(newCatFactor[,ind])])
    }
  )) # this is a list of q matrices, each n x k; with elements log-likelihood
  
  catLogLiks <- Reduce(f='+',x=logCatKProbs)
  
  # maximum "likelihood" classification
  combinedLogLik <- logRadDens + catLogLiks
  
  membership <- as.numeric(apply(combinedLogLik,1,which.max))
  
  return(membership)
}
