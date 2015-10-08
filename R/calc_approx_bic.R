
# function for calculating approximate BIC
calcApproxBIC <- function(resObj, BIC=TRUE) {
  nObs <- length(resObj$finalMemb)
  nConP <- prod(dim(resObj$finalCenters))
  nCatP <- sum(sapply(resObj$finalProbs,function(xx) nrow(xx)*( ncol(xx)-1 ) )) # subtract one param for each cat var
  nParmTotal <- nConP + nCatP + 1 # change if different density technique is used
  if (BIC) {
    thisCrit <- log(nObs)*nParmTotal - 2*resObj$finalLogLik
  } else {
    thisCrit <- 2*nParmTotal - 2*resObj$finalLogLik
  }
  return(list(
    criteria = thisCrit
   ,k = nParmTotal
   ,logLik = resObj$finalLogLik
   ,nConParm = nConP
   ,nCatParm = nCatP
  ))
}
