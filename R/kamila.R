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
# library(Rcpp)

#' @useDynLib kamila
#' @importFrom Rcpp sourceCpp
#' @importFrom KernSmooth bkde
#' @importFrom gtools rdirichlet
#' @importFrom abind abind

# tryCatch(
#  sourceCpp("./src/cppfunctions.cpp")
#  ,error = function(e){}
# )

# kdat must be a data frame with factor rows, levels 1:L;
# factors must have the correct levels specified to ensure
# missing cells are not lost
# This version deprecated; use my Rcpp implementation
myCatKern <- function(kdat, bw, tabOnly = TRUE) {
  tt <- table(kdat)
  dims <- dim(tt)
  ndim <- length(dims)
  nn <- nrow(kdat)
  if (length(bw) == 1) bw <- rep(bw, ndim)
  for (dd in 1:ndim) {
    dimCounts <- apply(X = tt, MARGIN = (1:ndim)[-dd], FUN = sum)
    l1 <- list()
    for (i in 1:dims[dd]) {
      dimInds <- dims
      dimInds[dd] <- 1
      rotatedMat <- array(dimCounts, dim = dimInds)
      l1[[i]] <- rotatedMat
    }
    l1$along <- dd
    offCounts <- do.call(abind::abind, l1)
    offCounts <- offCounts - tt
    tt <- (1 - bw[dd]) * tt + bw[dd] / (dims[dd] - 1) * offCounts
  }
  if (tabOnly) {
    return(tt)
  }
  preds <- mapply(
    FUN = function(ii) {
      l2 <- as.list(c(0, as.numeric(kdat[ii, ])))
      l2[[1]] <- tt / nn
      do.call("[", l2)
    },
    ii = 1:nn
  )
  return(list(preds = preds, tab = tt))
}

############
# function for initializing means
initMeans <- function(conVar, method, numClust) {
  if (method == "sample") {
    return(
      sapply(conVar, function(xx) sample(xx, size = numClust, replace = TRUE))
    )
  } else if (method == "runif") {
    ranges <- sapply(conVar, range)
    return(
      apply(ranges, 2, function(xx) runif(numClust, min = xx[1], max = xx[2]))
    )
  } else {
    stop("Unrecognized mean initialization method")
  }
}

# Setup worker nodes in a cluster with library paths and package namespace
setupClusterWorkers <- function(cl) {
  lp <- .libPaths()
  parallel::clusterCall(cl, function(p) .libPaths(p), lp)
  pkgPath <- tryCatch(
    getNamespaceInfo(asNamespace("kamila"), "path"),
    error = function(e) ""
  )
  isDev <- tryCatch(
    requireNamespace("pkgload", quietly = TRUE) && pkgload::is_dev_package("kamila"),
    error = function(e) FALSE
  )
  parallel::clusterCall(cl, function(is_dev, path) {
    if (is_dev && nzchar(path)) {
      # #nocov start
      pkgload::load_all(path, quiet = TRUE)
      # #nocov end
    } else {
      library(kamila)
    }
  }, isDev, pkgPath)
}

# Execute one prediction strength cross-validation run
calcSinglePsCvRun <- function(
  cvRun,
  numObs,
  numInTest,
  hasCon,
  hasCat,
  conVar,
  catFactor,
  numClust,
  numInit,
  conWeights,
  catWeights,
  maxIter,
  conInitMethod,
  catBw
) {
  resCol <- numeric(length(numClust))
  for (ithNcInd in seq_along(numClust)) {
    # generate cv indices
    testInd <- sample(numObs, size = numInTest, replace = FALSE)

    testCon <- if (hasCon) conVar[testInd, , drop = FALSE] else NULL
    testCat <- if (hasCat) catFactor[testInd, , drop = FALSE] else NULL
    trainCon <- if (hasCon) conVar[-testInd, , drop = FALSE] else NULL
    trainCat <- if (hasCat) catFactor[-testInd, , drop = FALSE] else NULL

    # cluster test data
    testClust <- kamila(
      conVar = testCon,
      catFactor = testCat,
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
      conVar = trainCon,
      catFactor = trainCat,
      numClust = numClust[ithNcInd],
      numInit = numInit,
      conWeights = conWeights,
      catWeights = catWeights,
      maxIter = maxIter,
      conInitMethod = conInitMethod,
      catBw = catBw,
      verbose = FALSE
    )

    # Allocate test data based on training clusters.
    testDataClassify <- if (hasCon && hasCat) {
      list(testCon, testCat)
    } else if (hasCon) {
      testCon
    } else {
      testCat
    }

    teIntoTr <- classifyKamila(
      trainClust,
      testDataClassify
    )

    # Calculate prediction strength proportions using Rcpp function.
    # Uses exact combinatorial identity via contingency table counts C_{k,m}:
    # psProps[k] = sum_m [C_{k,m} * (C_{k,m} - 1)] / [n_k * (n_k - 1)],
    # evaluating pair co-membership in O(N + K^2) time instead of O(N^2) loops.
    psProps <- calcPsCpp(
      testClust$finalMemb,
      teIntoTr,
      numClust[ithNcInd]
    )

    # Calculate and update prediction strength results.
    resCol[ithNcInd] <- ifelse(
      test = all(is.na(psProps)),
      yes = NA,
      no = min(psProps, na.rm = TRUE)
    )
  }
  return(resCol)
}


######################
# remove bumps (i.e. make sure xx is nondecreasing)
# Version 1 starts from right side moving left, replaces any decrease
# with the previous value.
# Not used; doesn't appear to be necessary.
# rmBump1 <- function(xx) {
#   for (i in (length(xx) - 1):1) {
#     if (xx[i] < xx[i + 1]) xx[i] <- xx[i + 1]
#   }
#   return(xx)
# }
# Version 2 starts from the left side moving right, replaces any increase
# with the earlier value
# rmBump2 <- function(xx) {
#   for (i in 2:length(xx)) {
#     if (xx[i] > xx[i - 1]) xx[i] <- xx[i - 1]
#   }
#   return(xx)
# }

######################
# Calculate weighted Euclidean distances from a set of n points
# in p-space to a set of k means
# distPointsToMeans <- function(pts, myMeans, wgts) {
#   ppDim <- ncol(pts)
#   if (ppDim != ncol(myMeans)) stop("Dimensionality of pts and myMeans must be equal")
#   if (ppDim != length(wgts)) stop("Dimensionality of pts must equal # weights")
#   kkMean <- nrow(myMeans)
#
#   sapply(
#     1:kkMean,
#     function(kk) {
#       diff_k <- sapply(
#         1:ppDim,
#         function(jj) {
#           wgts[jj] * (pts[, jj] - myMeans[kk, jj])
#         }
#       )
#       apply(diff_k, 1, function(ii) sqrt(sum(ii^2)))
#     }
#   )
# }

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
dptmCpp <- function(pts, myMeans, wgts) {
  if (!is.matrix(pts)) pts <- as.matrix(pts)
  ppDim <- ncol(pts)
  if (ppDim != ncol(myMeans)) stop("Dimensionality of pts and myMeans must be equal")
  if (ppDim != length(wgts)) stop("Dimensionality of pts must equal number of weights")
  kkMean <- nrow(myMeans)
  nn <- nrow(pts)

  dptm(pts, myMeans, wgts, ppDim, kkMean, nn)
}


######################
# psort: poor man's approximation to sample quantile
# Underestimates the quantile, but this bias decreases
# as sample size increases.
# Used in radialKDE function below, where we don't need
# an exact quantile anyway.
# Not used since gain in computation
# time is minimal.
# psort <- function(xx, pp) {
#   nn <- length(xx)
#   jj <- floor(nn * pp)
#   if (jj <= 1) {
#     return(min(xx))
#   }
#   return(sort(xx, partial = jj)[jj])
# }

######################
# Estimate density of radii
# Takes a vector of distances to mean, estimates radial KDE
# then evaluates at given vector of points
# pdim is the number of continuous variables used
# returnFun causes a resampling function to be returned
#' @importFrom stats bw.nrd0 approxfun quantile
radialKDE <- function(radii, evalPoints, pdim, returnFun = FALSE) {
  MAXDENS <- 1
  # Note using a chosen constant for bw reduces time by about 7%
  radialBW <- bw.nrd0(radii)
  maxEval <- max(evalPoints)
  radKDE <- bkde(
    x = radii,
    kernel = "normal",
    bandwidth = radialBW
    #   ,range.x = c(0, max(radii))
    , range.x = c(0, maxEval)
  )

  if (!returnFun) {
    kdes <- interpRadialKde(
      y = radKDE$y,
      maxEval = maxEval,
      pdim = pdim,
      evalPoints = evalPoints
    )
    return(list(kdes = kdes, resampler = NULL))
  }

  # remove any zero and negative density estimates
  newY <- radKDE$y
  nonnegTF <- newY > 0
  if (any(!nonnegTF)) {
    minPos <- min(newY[nonnegTF])
    newY[!nonnegTF] <- minPos / 100
  }

  # at bottom 5th percentile, replace with line through (0,0) and (q05,f(q05)).
  # This removes substantial variability in output near zero
  quant05 <- quantile(x = radKDE$x, probs = 0.05)
  # quant05 <- psort(xx = radKDE$x, pp = 0.05)
  coordsLtQ05 <- radKDE$x < quant05
  maxPt <- max(which(coordsLtQ05))
  newY[coordsLtQ05] <- radKDE$x[coordsLtQ05] * (newY[maxPt] / radKDE$x[maxPt]) # y = x * sl

  # radial Jacobian transformation; up to proportionality constant
  # Issue #9 fix: avoid zero density at r = 0 by right-continuous extension
  radY_rest <- newY[-1] / radKDE$x[-1]^(pdim - 1)
  radY <- c(radY_rest[1], radY_rest)

  # replace densities over MAXDENS with MAXDENS
  overMax <- radY > MAXDENS
  radY[overMax] <- MAXDENS

  # remove bumps (i.e. make sure nondecreasing)
  # not used currently
  #  radY <- rmBump1(radY)

  # normalize to area 1
  binwidthX <- diff(radKDE$x[1:2])
  densR <- radY / (binwidthX * sum(radY))

  # now create resampling function
  resampler <- approxfun(x = radKDE$x, y = densR, rule = 1:2, method = "linear")
  kdes <- resampler(evalPoints)
  kdes <- pmax(kdes, min(densR))

  # return(list(kdes=resampler(evalPoints),resampler=resampler))
  return(list(kdes = kdes, resampler = resampler))
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

#' KAMILA clustering of mixed-type, continuous-only, or categorical-only data.
#'
#' KAMILA is an iterative clustering method that equitably balances the
#' contribution of continuous and categorical variables. It also natively
#' supports continuous-only and categorical-only data.
#'
#' KAMILA (KAy-means for MIxed LArge data sets) is an iterative clustering
#' method that equitably balances the contribution of the continuous and
#' categorical variables. It uses a kernel density estimation technique to
#' flexibly model spherical clusters in the continuous domain, and uses a
#' multinomial model in the categorical domain.
#'
#' In addition to mixed-type data, KAMILA natively supports single-modality
#' data:
#' \itemize{
#'   \item \strong{Continuous-only data}: When \code{catFactor = NULL}, KAMILA
#'   clusters continuous data using semi-parametric radial kernel density
#'   estimation (RKDE), selecting the best initialization via continuous
#'   pseudo log-likelihood maximization.
#'   \item \strong{Categorical-only data}: When \code{conVar = NULL}, KAMILA
#'   clusters categorical data using fast C++ smoothed multinomial estimation,
#'   selecting the best initialization via categorical log-likelihood
#'   maximization.
#' }
#'
#' Weighting scheme: If no weights are desired, set all weights to 1 (the
#' default setting). Let a_1, ..., a_p denote the weights for p continuous
#' variables. Let b_1, ..., b_q denote the weights for q categorical variables.
#' Currently, continuous weights are applied during the calculation of
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
#' the number of clusters selected. Cluster pair co-membership evaluation is
#' accelerated via an Rcpp routine (calcPsCpp) that calculates exact agreement
#' proportions via contingency table partitioning in O(N + K^2) time and O(K^2)
#' memory, which is mathematically and numerically identical to the pairwise
#' formulation without creating an O(N^2) distance matrix.
#'
#' Prediction strength scores range from 0 to 1, reflecting cluster stability and
#' reproducibility across cross-validation folds. Higher scores indicate greater
#' cluster reproducibility. Relative differences in prediction strength scores across
#' candidate values of \code{numClust} provide insight into cluster structure: a sharp
#' drop in prediction strength when transitioning from \emph{k} to \emph{k+1} clusters
#' suggests that \emph{k+1} creates unstable or arbitrary partitions. Following Tibshirani
#' & Walther (2005), the largest candidate \emph{k} whose score (mean plus standard error)
#' meets or exceeds \code{predStrThresh} is selected.
#' @export
#' @importFrom stats runif sd setNames
#' @importFrom parallel makeCluster stopCluster parLapply clusterCall clusterExport clusterSetRNGStream
#' @param conVar An optional data frame of continuous variables. At least one of
#'   \code{conVar} or \code{catFactor} must be specified.
#' @param catFactor An optional data frame of factors. At least one of
#'   \code{conVar} or \code{catFactor} must be specified.
#' @param numClust The number of clusters returned by the algorithm.
#' @param numInit The number of initializations used.
#' @param conWeights A vector of continuous weights for the continuous variables.
#' @param catWeights A vector of weights for the categorical variables.
#' @param maxIter The maximum number of iterations in each run.
#' @param conInitMethod Character: The method used to initialize each run.
#' @param catBw The bandwidth used for the categorical kernel.
#' @param verbose Logical: Whether detailed results should be printed and returned.
#' @param calcNumClust Character: Method for selecting the number of clusters.
#' @param numPredStrCvRun Numeric: Number of CV runs for prediction strength method. Ignored unless calcNumClust == 'ps'
#' @param predStrThresh Numeric: Threshold for prediction strength method. Ignored unless calcNumClust == 'ps'
#' @param numCores Numeric or cluster: Number of CPU cores to use for parallel
#'   execution, or a cluster object created by \code{parallel::makeCluster}.
#'   Defaults to 1 (sequential execution). Ignored unless
#'   \code{calcNumClust == 'ps'}.
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
#' dat <- genMixedData(200,
#'   nConVar = 2, nCatVar = 2, nCatLevels = 4,
#'   nConWithErr = 2, nCatWithErr = 2, popProportions = c(.5, .5),
#'   conErrLev = 0.3, catErrLev = 0.8
#' )
#' catDf <- data.frame(apply(dat$catVars, 2, factor), stringsAsFactors = TRUE)
#' conDf <- data.frame(scale(dat$conVars), stringsAsFactors = TRUE)
#'
#' # Standard KAMILA clustering with a specified number of clusters
#' kamRes <- kamila(conDf, catDf, numClust = 2, numInit = 10)
#' table(kamRes$finalMemb, dat$trueID)
#'
#' \dontrun{
#' # KAMILA clustering with prediction strength estimation for number of clusters
#' kamPsRes <- kamila(
#'   conDf, catDf,
#'   numClust = 2:4, numInit = 10,
#'   calcNumClust = "ps", numPredStrCvRun = 10, predStrThresh = 0.8
#' )
#' kamPsRes$nClust$bestNClust
#' table(kamPsRes$finalMemb, dat$trueID)
#' }
#' @references Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and
#'   Hadoop. Journal of Statistical Software, 83(13). 2018.
#'   doi: 10.18637/jss.v083.i13
kamila <- function(
  conVar = NULL,
  catFactor = NULL,
  numClust,
  numInit,
  conWeights = NULL,
  catWeights = NULL,
  maxIter = 25,
  conInitMethod = "runif",
  catBw = 0.025,
  verbose = FALSE,
  calcNumClust = "none",
  numPredStrCvRun = 10,
  predStrThresh = 0.8,
  numCores = 1
) {
  hasCon <- !is.null(conVar)
  hasCat <- !is.null(catFactor)

  if (!hasCon && !hasCat) {
    stop("At least one of conVar or catFactor must be specified.")
  }

  if (hasCon) {
    if (!is.data.frame(conVar)) {
      stop("Input dataset conVar must be a dataframe.")
    }
    if (ncol(conVar) < 1) {
      stop("Input dataset conVar must have at least 1 column.")
    }
    if (anyNA(conVar)) {
      stop("Input dataset conVar contains missing values (NA). Missing values are not supported.")
    }
    numConVar <- ncol(conVar)
    if (is.null(conWeights)) {
      conWeights <- rep(1, numConVar)
    }
    if (length(conWeights) != numConVar) {
      stop("Length of conWeights must equal number of continuous variables")
    }
    if (max(conWeights) > 1 || min(conWeights) < 0) {
      stop("Weights must be in [0,1]")
    }
  } else {
    numConVar <- 0
    conWeights <- numeric(0)
  }

  if (hasCat) {
    if (!is.data.frame(catFactor)) {
      stop("Input dataset catFactor must be a dataframe.")
    }
    if (ncol(catFactor) < 1) {
      stop("Input dataset catFactor must have at least 1 column.")
    }
    if (anyNA(catFactor)) {
      stop("Input dataset catFactor contains missing values (NA). Missing values are not supported.")
    }
    numCatVar <- ncol(catFactor)
    if (is.null(catWeights)) {
      catWeights <- rep(1, numCatVar)
    }
    if (length(catWeights) != numCatVar) {
      stop("Length of catWeights must equal number of categorical variables")
    }
    if (max(catWeights) > 1 || min(catWeights) < 0) {
      stop("Weights must be in [0,1]")
    }
  } else {
    numCatVar <- 0
    catWeights <- numeric(0)
  }

  if (hasCon && hasCat) {
    if (nrow(catFactor) != nrow(conVar)) {
      stop("Number of observations in con and cat vars don't match")
    }
    numObs <- nrow(conVar)
  } else if (hasCon) {
    numObs <- nrow(conVar)
  } else {
    numObs <- nrow(catFactor)
  }

  if (calcNumClust == "none") {
    if (length(numClust) != 1) {
      stop('Input parameter numClust must be length 1 if calcNumClust == "none"')
    }
    # main function

    # Deprecated option
    returnResampler <- FALSE

    if (hasCon) {
      conVarMat <- as.matrix(conVar)
    }

    if (hasCat) {
      numLev <- as.integer(sapply(catFactor, function(xx) length(levels(xx))))
      catFactorNumeric <- matrix(
        as.integer(sapply(catFactor, as.integer, simplify = TRUE)),
        nrow = numObs,
        ncol = numCatVar
      )
    }

    numIterVect <- rep(NaN, numInit)
    totalLogLikVect <- rep(NaN, numInit)
    catLogLikVect <- rep(NaN, numInit)
    winDistVect <- rep(NaN, numInit)
    if (hasCon) {
      totalDist <- sum(dptm(
        pts = conVarMat,
        myMeans = matrix(colMeans(conVarMat), nrow = 1),
        wgts = conWeights,
        ppDim = numConVar,
        kkMean = 1,
        nn = numObs
      ))
    } else {
      totalDist <- NaN
    }
    objectiveVect <- rep(NaN, numInit)

    # for verbose output, list of memberships for each init, iter
    if (verbose) membLongList <- rep(list(list()), numInit)

    # (1) loop over each initialization
    for (init in 1:numInit) {
      if (hasCon) {
        means_i <- initMeans(conVar = conVar, method = conInitMethod, numClust = numClust)
      } else {
        means_i <- NULL
      }

      if (hasCat) {
        logProbsCond_i <- lapply(
          numLev,
          function(xx) {
            matrix(
              data = log(gtools::rdirichlet(n = numClust, alpha = rep(1, xx))),
              nrow = numClust,
              dimnames = list(clust = 1:numClust, level = 1:xx)
            )
          }
        )
      } else {
        logProbsCond_i <- list()
      }

      # initialize structures for iterative procedure
      membOld <- membNew <- rep(0, numObs)
      numIter <- 0
      degenerateSoln <- FALSE

      # Loop until convergence
      while (
        ((numIter < 3) || !all(membOld == membNew)) &&
          (numIter < maxIter)
      ) {
        numIter <- numIter + 1

        if (hasCon) {
          dist_i <- dptm(
            pts = conVarMat,
            myMeans = means_i,
            wgts = conWeights,
            ppDim = numConVar,
            kkMean = numClust,
            nn = numObs
          )
          minDist_i <- rowMin(dist_i)

          logDistRadDens_vec <- log(
            radialKDE(
              radii = minDist_i,
              evalPoints = c(dist_i),
              pdim = numConVar,
              returnFun = returnResampler
            )$kdes
          )
          logDistRadDens_i <- matrix(
            logDistRadDens_vec,
            nrow = numObs,
            ncol = numClust
          )
        }

        if (hasCat) {
          catLogLiks <- calcCatLogLiks(
            catFactorNum = catFactorNumeric,
            catWeights = catWeights,
            logProbsCond_i = logProbsCond_i
          )
        }

        if (hasCon && hasCat) {
          allLogLiks <- logDistRadDens_i + catLogLiks
        } else if (hasCon) {
          allLogLiks <- logDistRadDens_i
        } else {
          allLogLiks <- catLogLiks
        }

        # partition data into clusters
        membOld <- membNew
        membNew <- as.integer(rowMaxInds(allLogLiks))

        # calculate new means / probabilities
        if (hasCon) {
          means_i <- aggregateMeans(
            conVar = conVarMat,
            membNew = membNew,
            kk = numClust
          )
        }

        if (hasCat) {
          logProbsCond_i <- updateLogProbs(
            catFactorNum = catFactorNumeric,
            membNew = membNew,
            numLev = numLev,
            catBw = catBw,
            kk = numClust
          )
        }

        if (verbose) {
          membLongList[[init]][[numIter]] <- membOld
        }

        if (any(tabulate(membNew, numClust) == 0L)) {
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

      numIterVect[init] <- numIter

      # other useful internal measures of cluster quality
      if (hasCat) catLogLikVect[init] <- sum(rowMax(catLogLiks))
      if (hasCon) {
        winDistVect[init] <- sum(dist_i[1:numObs + (membNew - 1L) * numObs])
        winToBetRat <- winDistVect[init] / (totalDist - winDistVect[init])
        if (winToBetRat < 0) winToBetRat <- 100
      }

      if (hasCon && hasCat) {
        objectiveVect[init] <- winToBetRat * catLogLikVect[init]
      } else if (hasCon) {
        objectiveVect[init] <- totalLogLikVect[init]
      } else {
        objectiveVect[init] <- if (degenerateSoln) -Inf else catLogLikVect[init]
      }

      # Store current solution if objective beats all others
      if (
        (init == 1) ||
          (init > 1 && objectiveVect[init] > max(objectiveVect[1:(init - 1)]))
      ) {
        finalLogLik <- totalLogLikVect[init]
        finalObj <- objectiveVect[init]
        finalMemb <- membNew
        if (hasCon) {
          finalCenters <- matrix(
            data = as.matrix(means_i),
            nrow = nrow(means_i),
            ncol = numConVar,
            dimnames = list(
              cluster = paste("Clust", seq_len(nrow(means_i))),
              variableMean = paste("Mean", seq_len(numConVar))
            )
          )
        } else {
          finalCenters <- NULL
        }

        if (hasCat) {
          finalProbs <- lapply(logProbsCond_i, exp)
          names(finalProbs) <- paste("Categorical Variable", seq_len(numCatVar))
        } else {
          finalProbs <- list()
        }
        finalClustSize <- table(membNew)
      }

      if (verbose) membLongList[[init]][[numIter + 1]] <- membNew
    }

    # 12 Prepare output data structure
    if (verbose) {
      optionalOutput <- list(
        totalLogLikVect = totalLogLikVect,
        catLogLikVect = if (hasCat) catLogLiks else NULL,
        winDistVect = if (hasCon) winDistVect else NULL,
        totalDist = if (hasCon) totalDist else NULL,
        objectiveVect = objectiveVect,
        membLongList = membLongList
      )
    } else {
      optionalOutput <- list()
    }

    inputList <- list(
      conVar = conVar,
      catFactor = catFactor,
      numClust = numClust,
      numInit = numInit,
      conWeights = conWeights,
      catWeights = catWeights,
      maxIter = maxIter,
      conInitMethod = conInitMethod,
      catBw = catBw,
      verbose = verbose
    )

    return(
      list(
        finalMemb = as.numeric(finalMemb),
        numIter = numIterVect,
        finalLogLik = finalLogLik,
        finalObj = finalObj,
        finalCenters = finalCenters,
        finalProbs = finalProbs,
        input = inputList,
        verbose = optionalOutput,
        nClust = list()
      )
    )
  } else if (calcNumClust == "ps") {
    # Recursive call to kamila implementing prediction strength method.

    # Test that numClust is an integer vector.
    allEqInteger <- all(numClust == as.integer(numClust))
    if (
      is.na(allEqInteger) ||
        !allEqInteger ||
        !all(sapply(numClust, is.numeric)) ||
        length(numClust) != length(unique(numClust))
    ) {
      stop('Input parameter numClust must be a vector of unique integers
         if input parameter calcNumClust == "ps"')
    }
    if (length(numClust) == 1) {
      warning("Input parameter numClust is a scalar; the prediction strength
        method is probably not appropriate or desired")
    }
    if (any(numClust > floor(numObs / 2))) {
      stop("The number of clusters in input parameter numClust cannot exceed
        one-half of the sample size.")
    }

    # Test that predStrThresh is within (0,1).
    if (
      length(predStrThresh) != 1 ||
        is.na(predStrThresh) ||
        predStrThresh <= 0 ||
        predStrThresh >= 1
    ) {
      stop("Input parameter predStrThresh must be scalar in (0,1)")
    }

    # Test that numPredStrCvRun is a valid number of cv runs.
    if (
      length(numPredStrCvRun) != 1 ||
        is.na(numPredStrCvRun) ||
        numPredStrCvRun != as.integer(numPredStrCvRun) ||
        numPredStrCvRun < 1
    ) {
      stop("Input parameter numPredStrCvRun must be a positive integer.")
    }

    # Test that numCores is a valid number of cores or a cluster object.
    if (!inherits(numCores, "cluster")) {
      if (
        length(numCores) != 1 ||
          is.na(numCores) ||
          !is.numeric(numCores) ||
          numCores != as.integer(numCores) ||
          numCores < 1
      ) {
        stop("Input parameter numCores must be a positive integer or a cluster object.")
      }
    }

    psCvRes <- matrix(
      NaN,
      nrow = length(numClust),
      ncol = numPredStrCvRun,
      dimnames = list(
        nClust = numClust,
        CVRun = paste("Run", 1:numPredStrCvRun)
      )
    )

    # Implement CV procedure
    numInTest <- floor(numObs / 2)
    isClustObj <- inherits(numCores, "cluster")
    effectiveCores <- if (isClustObj) length(numCores) else min(as.integer(numCores), numPredStrCvRun)

    if (!isClustObj && (numCores == 1 || effectiveCores == 1)) {
      for (cvRun in 1:numPredStrCvRun) {
        psCvRes[, cvRun] <- calcSinglePsCvRun(
          cvRun = cvRun,
          numObs = numObs,
          numInTest = numInTest,
          hasCon = hasCon,
          hasCat = hasCat,
          conVar = conVar,
          catFactor = catFactor,
          numClust = numClust,
          numInit = numInit,
          conWeights = conWeights,
          catWeights = catWeights,
          maxIter = maxIter,
          conInitMethod = conInitMethod,
          catBw = catBw
        )
      }
    } else {
      if (isClustObj) {
        cl <- numCores
      } else {
        cl <- parallel::makeCluster(effectiveCores)
        on.exit(parallel::stopCluster(cl), add = TRUE)
      }

      setupClusterWorkers(cl)
      parallel::clusterSetRNGStream(cl)

      parallel::clusterExport(
        cl = cl,
        varlist = c(
          "calcSinglePsCvRun", "numObs", "numInTest", "hasCon", "hasCat",
          "conVar", "catFactor", "numClust", "numInit", "conWeights",
          "catWeights", "maxIter", "conInitMethod", "catBw"
        ),
        envir = environment()
      )

      psResList <- parallel::parLapply(
        cl = cl,
        X = 1:numPredStrCvRun,
        fun = function(cvRun) {
          calcSinglePsCvRun(
            cvRun = cvRun,
            numObs = numObs,
            numInTest = numInTest,
            hasCon = hasCon,
            hasCat = hasCat,
            conVar = conVar,
            catFactor = catFactor,
            numClust = numClust,
            numInit = numInit,
            conWeights = conWeights,
            catWeights = catWeights,
            maxIter = maxIter,
            conInitMethod = conInitMethod,
            catBw = catBw
          )
        }
      )

      for (cvRun in 1:numPredStrCvRun) {
        psCvRes[, cvRun] <- psResList[[cvRun]]
      }
    }

    # Calculate CV estimate of prediction strength for each cluster size.
    avgPredStr <- apply(psCvRes, 1, mean, na.rm = TRUE)
    stdErrPredStr <- if (numPredStrCvRun > 1) {
      apply(psCvRes, 1, sd, na.rm = TRUE) / sqrt(numPredStrCvRun)
    } else {
      setNames(rep(0, length(numClust)), numClust)
    }

    # Calculate final number of clusters: largest # clust such that avg+sd
    # score is above the threshold.
    psValues <- avgPredStr + stdErrPredStr
    clustAboveThresh <- psValues > predStrThresh

    if (all(!clustAboveThresh)) {
      warning("No cluster size is above prediction strength threshold.
        Consider lowering the ps threshold; returning the cluster size
	corresponding to the highest ps value.")
      kfinal <- numClust[which.max(psValues)]
    } else {
      kfinal <- max(numClust[clustAboveThresh])
    }

    # Do the final clustering.
    outRes <- kamila(
      conVar = conVar,
      catFactor = catFactor,
      numClust = kfinal,
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
  tmp <- 2 * invar * pi / (maxDetected + 1 - minDetected)
  return(cbind(cos(tmp), sin(tmp)))
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
#' which contains the new observations to classify. For mixed-type models,
#' \code{newData} must be a list of length 2 (continuous data frame and
#' categorical factor data frame). For single-modality models (continuous-only
#' or categorical-only), \code{newData} can be supplied either as a data frame
#' or as a list containing the single data frame.
#' @export
#' @param obj An output object from the kamila function.
#' @param newData For mixed-type data, a list of length 2, with first element
#'   a data frame of continuous variables, and second element a data frame of
#'   categorical factors. For continuous-only or categorical-only models, either
#'   a single data frame or a list containing the single data frame.
#' @return An integer vector denoting cluster assignments of the new data points.
#' @examples
#' # Generate toy data set
#' set.seed(1234)
#' dat1 <- genMixedData(400,
#'   nConVar = 2, nCatVar = 2, nCatLevels = 4,
#'   nConWithErr = 2, nCatWithErr = 2, popProportions = c(.5, .5),
#'   conErrLev = 0.2, catErrLev = 0.2
#' )
#' # Partition the data into training/test set
#' trainingIds <- sample(nrow(dat1$conVars), size = 300, replace = FALSE)
#' catTrain <- data.frame(apply(dat1$catVars[trainingIds, ], 2, factor), stringsAsFactors = TRUE)
#' conTrain <- data.frame(scale(dat1$conVars)[trainingIds, ], stringsAsFactors = TRUE)
#' catTest <- data.frame(apply(dat1$catVars[-trainingIds, ], 2, factor), stringsAsFactors = TRUE)
#' conTest <- data.frame(scale(dat1$conVars)[-trainingIds, ], stringsAsFactors = TRUE)
#' # Run the kamila clustering procedure on the training set
#' kamilaObj <- kamila(conTrain, catTrain, numClust = 2, numInit = 10)
#' table(dat1$trueID[trainingIds], kamilaObj$finalMemb)
#' # Predict membership in the test data set
#' kamilaPred <- classifyKamila(kamilaObj, list(conTest, catTest))
#' table(dat1$trueID[-trainingIds], kamilaPred)
#' @references Foss A, Markatou M; kamila: Clustering Mixed-Type Data in R and
#'   Hadoop. Journal of Statistical Software, 83(13). 2018.
#'   doi: 10.18637/jss.v083.i13
classifyKamila <- function(obj, newData) {
  hasCon <- is.list(obj) && !is.null(obj$input$conVar) && ncol(obj$input$conVar) > 0
  hasCat <- is.list(obj) && !is.null(obj$input$catFactor) && ncol(obj$input$catFactor) > 0

  if (
    !is.list(obj) || is.null(obj$input) ||
      (!hasCon && !hasCat) ||
      (hasCon && is.null(obj$finalCenters)) ||
      (hasCat && is.null(obj$finalProbs))
  ) {
    stop("Error in function classifyKamila: obj must be a valid kamila object")
  }

  if (is.list(newData) && !is.data.frame(newData)) {
    if (length(newData) == 2) {
      newCon <- newData[[1]]
      newCatFactor <- newData[[2]]
    } else if (length(newData) == 1) {
      if (hasCon && !hasCat) {
        newCon <- newData[[1]]
        newCatFactor <- NULL
      } else if (!hasCon && hasCat) {
        newCon <- NULL
        newCatFactor <- newData[[1]]
      } else {
        stop("Error in function classifyKamila: newData must be list of length 2 for mixed data")
      }
    } else {
      stop("Error in function classifyKamila: newData list must have length 1 or 2")
    }
  } else if (is.data.frame(newData) || (hasCon && !hasCat && is.matrix(newData))) {
    if (hasCon && !hasCat) {
      newCon <- newData
      newCatFactor <- NULL
    } else if (!hasCon && hasCat) {
      newCon <- NULL
      newCatFactor <- newData
    } else {
      stop("Error in function classifyKamila: newData must be list of length 2 for mixed data")
    }
  } else {
    stop("Error in function classifyKamila: newData must be a data frame or list of data frames")
  }

  if (hasCon) {
    if (!is.data.frame(newCon) && !is.matrix(newCon)) {
      stop("Error in function classifyKamila: newData continuous element must be a data frame or matrix")
    }
    newCon <- as.data.frame(newCon)
    if (ncol(newCon) < 1) {
      stop("Error in function classifyKamila: data frames in newData must have at least 1 column")
    }
    if (ncol(newCon) != ncol(obj$input$conVar)) {
      stop(paste(
        "Error in function classifyKamila: number of continuous columns in newData",
        "does not match model"
      ))
    }
  }

  if (hasCat) {
    if (!is.data.frame(newCatFactor)) {
      stop("Error in function classifyKamila: elements of newData must be data frames")
    }
    if (ncol(newCatFactor) < 1) {
      stop("Error in function classifyKamila: data frames in newData must have at least 1 column")
    }
    if (ncol(newCatFactor) != ncol(obj$input$catFactor)) {
      stop(paste(
        "Error in function classifyKamila: number of categorical columns in newData",
        "does not match model"
      ))
    }

    numCatVar <- ncol(newCatFactor)
    for (ind in seq_len(numCatVar)) {
      trLevels <- levels(obj$input$catFactor[[ind]])
      colName <- colnames(newCatFactor)[ind]
      varDesc <- if (!is.null(colName) && nchar(colName) > 0) {
        paste0("'", colName, "' (column ", ind, ")")
      } else {
        paste0("column ", ind)
      }

      valChar <- as.character(newCatFactor[[ind]])
      uniqVals <- unique(valChar[!is.na(valChar)])
      unseenLevels <- setdiff(uniqVals, trLevels)

      if (length(unseenLevels) > 0) {
        formattedLevels <- paste(paste0("'", unseenLevels, "'"), collapse = ", ")
        stop(
          "Error in function classifyKamila: Categorical variable ",
          varDesc,
          " contains level(s) not present in training data: ",
          formattedLevels
        )
      }

      newCatFactor[[ind]] <- factor(newCatFactor[[ind]], levels = trLevels)
    }
  }

  if (hasCon && hasCat) {
    if (nrow(newCon) != nrow(newCatFactor)) {
      stop("Error in function classifyKamila: number of observations in con and cat vars don't match")
    }
  }

  if (hasCon) {
    conVarMat <- as.matrix(obj$input$conVar)
    distances <- dptm(
      pts = conVarMat,
      myMeans = obj$finalCenters,
      wgts = obj$input$conWeights,
      ppDim = ncol(conVarMat),
      kkMean = nrow(obj$finalCenters),
      nn = nrow(conVarMat)
    )
    minDistances <- rowMin(distances)

    newConMat <- as.matrix(newCon)
    newDistances <- dptm(
      pts = newConMat,
      myMeans = obj$finalCenters,
      wgts = obj$input$conWeights,
      ppDim = ncol(newConMat),
      kkMean = nrow(obj$finalCenters),
      nn = nrow(newConMat)
    )

    logRadDens <- matrix(
      log(radialKDE(radii = minDistances, evalPoints = c(newDistances), pdim = ncol(newConMat))$kdes),
      nrow = nrow(newConMat),
      ncol = nrow(obj$finalCenters)
    )
  }

  if (hasCat) {
    numCatVar <- ncol(newCatFactor)
    logClustProbs <- lapply(obj$finalProbs, log)
    newCatFactorNum <- matrix(
      as.integer(unlist(lapply(newCatFactor, as.integer), use.names = FALSE)),
      nrow = nrow(newCatFactor),
      ncol = numCatVar
    )
    catLogLiks <- calcCatLogLiks(
      catFactorNum = newCatFactorNum,
      catWeights = obj$input$catWeights,
      logProbsCond_i = logClustProbs
    )
  }

  if (hasCon && hasCat) {
    combinedLogLik <- logRadDens + catLogLiks
  } else if (hasCon) {
    combinedLogLik <- logRadDens
  } else {
    combinedLogLik <- catLogLiks
  }

  membership <- as.numeric(rowMaxInds(combinedLogLik))

  return(membership)
}
