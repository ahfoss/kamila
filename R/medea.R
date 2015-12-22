
clustOneVar <- function(vec,k,nstart,iter.max) {
  nDistinctPoints <- length(unique(vec))
  if (nDistinctPoints < k) {
    if (nDistinctPoints == 1) {
      stop('Error in clustOneVar: variables must have more than one value.')
    }
    #print(length(unique(vec)))
    k <- nDistinctPoints
  }
  return(kmeans(
    vec,
    centers=k,
    nstart=nstart,
    iter.max=iter.max
  )$cluster)
}


# set variable using input args from list, with defaults
# Returns the value from input (list) named "name", with the
# given default if no element in input is named "name".
setArgs <- function(input, name, default) {
  output <- default
  nameBool <- (names(input) == name)
  if (sum(nameBool) > 1) {
    stop('multiple matches to name in function setArgs')
  }
  if (any(nameBool)) {
    output <- input[[which(nameBool)]] # could replace with syntax input[[name]] but test first
  }
  return(output)
}


# helper function for medea with famd option
propVar <- function(X1,X2) {
  if (is.numeric(X1) && is.numeric(X2)) {
    depVar <- scale(X1)
    indVar <- scale(X2)
    return(summary(lm(depVar ~ indVar))$r.squared)
  }
  if (is.numeric(X1) && is.factor(X2)) {
    depVar <- scale(X1)
    indVar <- model.matrix(~ 0 + X2)
    return(summary(lm(depVar ~ indVar))$r.squared)
  }
  if (is.factor(X1) && is.numeric(X2)) {
    depVar <- scale(X2)
    indVar <- model.matrix(~ 0 + X1)
    return(summary(lm(depVar ~ indVar))$r.squared)
  }
  if (is.factor(X1) && is.factor(X2)) {
    X1d <- model.matrix(~ 0 + X1)
    X2d <- model.matrix(~ 0 + X2)
    if (ncol(X1d) < 2 || ncol(X2d) < 2) stop('Factor must have more than one level')

    # original method
    #H1 <- X1d %*% solve(t(X1d) %*% X1d) %*% t(X1d)
    #H2 <- X2d %*% solve(t(X2d) %*% X2d) %*% t(X2d)
    #nn <- nrow(X1d)
    #ninv <- 1/nn
    #Jmat <- matrix(ninv, nrow=nn, ncol=nn)
    #ratio1 <- sum(diag(t(X1d) %*% (H2-Jmat) %*% X1d)) / sum(diag(t(X1d) %*% (diag(nn)-Jmat) %*% X1d))
    #ratio2 <- sum(diag(t(X2d) %*% (H1-Jmat) %*% X2d)) / sum(diag(t(X2d) %*% (diag(nn)-Jmat) %*% X2d))

    # faster method
    ssr1 <- 0
    sst1 <- 0
    for (i in 1:ncol(X1d)) {
      res1 <- lm(X1d[,i] ~ X2d)
      aov1 <- anova(res1)
      sss1 <- aov1$`Sum Sq`
      ssr1 <- ssr1 + sum(sss1[-length(sss1)])
      sst1 <- sst1 + sum(sss1)
    }
    rSq1 <- ssr1 / sst1

    ssr2 <- 0
    sst2 <- 0
    for (i in 1:ncol(X2d)) {
      res2 <- lm(X2d[,i] ~ X1d)
      aov2 <- anova(res2)
      sss2 <- aov2$`Sum Sq`
      ssr2 <- ssr2 + sum(sss2[-length(sss2)])
      sst2 <- sst2 + sum(sss2)
    }
    rSq2 <- ssr2 / sst2

    return(rSq1/2 + rSq2/2)
  } 
  stop('Unrecognized input types.')
}

# helper function for medea with famd option
propVarAssocMat <- function(dat) {
  pp <- ncol(dat)
  outMat <- matrix(0,nrow=pp,ncol=pp)
  for (i in 1:(pp-1)) {
    for (j in (i+1):pp) {
      outMat[i,j] <- propVar(dat[,i],dat[,j])
    }
  }
  return(outMat + t(outMat))
}


# Options: ari+kmeans, ari+quantile, or famd

#' MEDEA weight calculation
#' 
#' @export
#' @param dat A data frame of numeric and factor variables.
#' @param discretizationMethod The discretization method used.
#' @param associationMethod The method for calculating the association matrix.
#' @param verbosity Whether to include extended results objects.
#' @param ... Optional input to the discretization function.
#' @return A list containing a vector \code{wgts} of variable-specific weights with
#' length equal to \code{ncol(dat)}.
medeaWgts <- function(
  dat,
  discretizationMethod='kmeans',
  associationMethod='ari',
  verbosity=FALSE,
  ...
) {
  varargs <- list(...)
  DISC_METHODS <- c('kmeans','quantile','')
  ASSOC_METHODS <- c('ari','famd')
  whichCon <- which(!(sapply(dat,class) == 'factor'))
  
  if (associationMethod == 'ari') {
    # discretize variables
    if (discretizationMethod == 'kmeans') {
      # define discFun in terms of clustOneVar
      discFun <- function(x) {
        centers <- setArgs(varargs, 'centers', 4)
        nstart <- setArgs(varargs, 'nstart', 10)
        iter.max <- setArgs(varargs, 'iter.max', 20)
        return(clustOneVar(
          vec=x,
          k=centers,
          nstart=nstart,
          iter.max=iter.max))
      }
    } else if (discretizationMethod == 'quantile') {
      # define discFun in terms of quantile
      discFun <- function(x) {
        probs <- setArgs(varargs, 'probs', c(1:3)/4)
        qs <- unique(quantile(x, probs))
        return(as.numeric(cut(x, breaks=c(-Inf, qs, Inf))))
      }
    } else {
      stop('unknown discretizationMethod specified in function medeaWgts')
    }
    # sapply with discFun to generate categorical data frame
    if (length(whichCon) > 1) {
      clustered <- sapply(
        dat[,whichCon]
        ,function(xx) discFun(xx)
      )
    } else if (length(whichCon) == 1) {
      clustered <- discFun(dat[,whichCon])
    }
    dat[,whichCon] <- clustered
    dat <- as.data.frame(lapply(dat,factor))
    # use ari to generate association matrix
    nvar <- ncol(dat)
    assocMat <- matrix(nrow=nvar,ncol=nvar)
    diag(assocMat) <- 1
    for (i in 1:(nvar-1)) {
      for (j in (i+1):nvar) {
        assocMat[i,j] <- assocMat[j,i] <- mclust::adjustedRandIndex(dat[,i],dat[,j])
      }
    }
    # analyze matrix
    eigs <- eigen(assocMat - diag(nvar))
    wgts <- abs(eigs$vectors[,1])
    results <- list(wgts=wgts)
    if (verbosity & discretizationMethod %in% DISC_METHODS) {
      results$verbose <- list(ARI=assocMat, eigs=eigs, cluster=clustered)
    }
    return(results)
  } else if (associationMethod == 'famd') {
    assocMat <- propVarAssocMat(dat)
    # note the association matrix already has zeros
    # along the main diagonal as required by medea
    eigs <- eigen(assocMat)
    wgts <- abs(eigs$vectors[,1])
    results <- list(wgts=wgts)
    return(results)
  }
}

