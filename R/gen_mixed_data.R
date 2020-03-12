
#' Generate simulated mixed-type data with cluster structure.
#'
#' This function simulates mixed-type data sets with a latent cluster
#' structure, with continuous and nominal variables.
#'
#' This function simulates mixed-type data sets with a latent cluster
#' structure. Continuous variables follow a normal mixture model, and
#' categorical variables follow a multinomial mixture model. Overlap of the
#' continuous and categorical variables (i.e. how clear the cluster structure
#' is) can be manipulated by the user. Overlap between two clusters is the area
#' of the overlapping region defined by their densities (or, for categorical
#' variables, the summed height of overlapping segments defined by their point
#' masses). The default overlap level is 0.01 (i.e. almost perfect separation).
#' A user-specified number of continuous and categorical variables can be
#' specified to be "error variables" with arbitrary overlap within 0.01 and 1.00
#' (where 1.00 corresponds to complete overlap).
#' NOTE: Currently, only two populations (clusters) are supported. While exact
#' control of overlap between two clusters is straightforward, controlling the
#' overlap between the K choose 2 pairwise combinations of clusters is a more
#' difficult task.
#'
#' @export
#' @importFrom stats qnorm rnorm
#' @param sampSize Integer: Size of the simulated data set.
#' @param nConVar The number of continuous variables.
#' @param nCatVar The number of categorical variables.
#' @param nCatLevels Integer: The number of categories per categorical variables. Currently must be a multiple of the number of populations specified in popProportions.
#' @param nConWithErr Integer: The number of continuous variables with error.
#' @param nCatWithErr Integer: The number of categorical variables with error.
#' @param popProportions A vector of scalars that sums to one. The length gives the number of populations (clusters), with values denoting the prior probability of observing a member of the corresponding population. NOTE: currently only two populations are supported.
#' @param conErrLev A scalar between 0.01 and 1 denoting the univariate overlap between clusters on the continuous variables specified to have error.
#' @param catErrLev Univariate overlap level for the categorical variables with error.
#' @return A list with the following elements:
#' \item{trueID}{Integer vector giving population (cluster) membership of each observation}
#' \item{trueMus}{Mean parameters used for population (cluster) centers in the continuous variables}
#' \item{conVars}{The continuous variables}
#' \item{errVariance}{Variance parameter used for continuous error distribution}
#' \item{popProbsNoErr}{Multinomial probability vectors for categorical variables without measurement error}
#' \item{popProbsWithErr}{Multinomial probability vectors for categorical variables with measurement error}
#' \item{catVars}{The categorical variables}
#' @examples
#' dat <- genMixedData(100, 2, 2, nCatLevels=4, nConWithErr=1, nCatWithErr=1,
#'   popProportions=c(0.3,0.7), conErrLev=0.3, catErrLev=0.2)
#' with(dat,plot(conVars,col=trueID))
#' with(dat,table(data.frame(catVars[,1:2],trueID, stringsAsFactors = TRUE)))
genMixedData = function(
  sampSize
 ,nConVar
 ,nCatVar
 ,nCatLevels
 ,nConWithErr
 ,nCatWithErr
 ,popProportions
 ,conErrLev
 ,catErrLev
) {

#####################################
# ARgument checking
#####################################

# Only two populations currently supported
if (length(popProportions) > 2) {
  # MixSim will need to be used
  stop('Error in genMixedData: More than two populations not currently supported')
}

# check number of continuous error variables specified
if (is.logical(nConWithErr)) {
  if (length(nConWithErr) != nConVar) {
    stop('Error in genMixedData: length of nConWithErr must equal nConVar')
  }
} else if ( length(nConWithErr) != 1 | !is.numeric(nConWithErr) ) {
  # else if nConWithErr is not a scalar of length 1
  stop('Error in genMixedData: improper input to nConWithErr')
  } else if (nConWithErr > nConVar) {
    stop('Error in genMixedData: nConWithErr must be less than nConVar')
} else {
  # Replace scalar nConWithErr with boolean; default to first with error
  nConWithErr = rep(c(T,F), c(nConWithErr,nConVar - nConWithErr))
}


# Check that nCatLevels is a multiple of number of populations
if (nCatLevels %% length(popProportions) != 0) {
  stop('Error in genMixedData: number of categorical levels must
        be a multiple of number of populations')
}

# Set variables of general interest
OVERLAP_DEFAULT = 0.01 # overlap

#####################################
# Generate true population membership
#####################################
trueID = sample(
  x = 1:length(popProportions)
 ,size = sampSize
 ,replace = T
 ,prob = popProportions
)

#####################################
# Continuous variables
#####################################

# Note that mean1 = 0, mean2 = \pm 2 * qnorm(x/2)
# guarantees overlap of 100*x% when sd=1
trueMus = c(0,-2*qnorm(OVERLAP_DEFAULT/2))
muVec = trueMus[trueID]

# Generate continuous variables
conVars = matrix(
  rnorm(sampSize*nConVar,mean=muVec,sd=1)
 ,nrow=sampSize
 ,ncol=nConVar
)

# Use formula in my IBM presentation with closed-form
# solution for a sigma achieving target overlap:
#  sig = (mu1 - mu2) / (2*qnorm(overlap/2))
# where mu1 < mu2; 
# we then subtract one from sig since variables already
# have variance 1
errVariance = (abs(diff(trueMus)) / (2 * qnorm(conErrLev/2)))^2 - 1

# Code to check the previous theoretical overlap:
#with(res,pnorm(abs(diff(trueMus))/2,mean=trueMus[2],sd=sqrt(errVariance+1)))
# Code to check that desired variance is attained:
#tapply(res$conVars[,1],res$trueID,var)


# Add continuous error
for (i in 1:nConVar) {
  if (nConWithErr[i]) {
    ithError = rnorm(sampSize,mean=0,sd=sqrt(errVariance))
    conVars[,i] = conVars[,i] + ithError
  }
}


#####################################
# Categorical variables
#####################################

# Generate error-free prob vecs
rightCatProb = (1-OVERLAP_DEFAULT/2) / (nCatLevels/2)
wrongCatProb = (OVERLAP_DEFAULT/2) / (nCatLevels/2)
popProb1 = rep(c(rightCatProb,wrongCatProb),each=nCatLevels/2)
popProb2 = rep(c(wrongCatProb,rightCatProb),each=nCatLevels/2)

# Number observed pop1
numInPop1 = sum(trueID==1)

if (nCatWithErr < nCatVar) {
  # Simulate population 1 no error
  pop1CatNoErr = sample(
    x = 1:nCatLevels
   ,size = numInPop1*(nCatVar - nCatWithErr)
   ,replace = T
   ,prob = popProb1
  )

  # simulate population 2 no error
  pop2CatNoErr = sample(
    x = 1:nCatLevels
   ,size = (sampSize - numInPop1)*(nCatVar - nCatWithErr)
   ,replace = T
   ,prob = popProb2
  )

  # combine no error data
  catNoErr = rep(NA,sampSize*(nCatVar - nCatWithErr))
  vectIdNoErr = rep(trueID,nCatVar - nCatWithErr)
  catNoErr[vectIdNoErr==1] = pop1CatNoErr
  catNoErr[vectIdNoErr==2] = pop2CatNoErr

  catNoErr = matrix(
    catNoErr
   ,ncol = nCatVar - nCatWithErr
  )
} else {
  # Empty if no variables error-free
  catNoErr = c()
}

if (nCatWithErr > 0) {
  # Generate error prone prob vecs
  rightCatProbErr = (1-catErrLev/2) / (nCatLevels/2)
  wrongCatProbErr = (catErrLev/2) / (nCatLevels/2)
  popProb1Err = rep(c(rightCatProbErr,wrongCatProbErr),each=nCatLevels/2)
  popProb2Err = rep(c(wrongCatProbErr,rightCatProbErr),each=nCatLevels/2)

  # Simulate pops 1 and 2 with error
  pop1CatWithErr = sample(
    x = 1:nCatLevels
   ,size = numInPop1*nCatWithErr
   ,replace = T
   ,prob = popProb1Err
  )
  pop2CatWithErr = sample(
    x = 1:nCatLevels
   ,size = (sampSize - numInPop1)*nCatWithErr
   ,replace = T
   ,prob = popProb2Err
  )

  # combine error data pops 1 and 2
  catWithErr = rep(NA,sampSize*nCatWithErr)
  vectIdWithErr = rep(trueID,nCatWithErr)
  catWithErr[vectIdWithErr==1] = pop1CatWithErr
  catWithErr[vectIdWithErr==2] = pop2CatWithErr

  catWithErr = matrix(
    catWithErr
   ,ncol = nCatWithErr
  )
  
} else {
  # Empty if no variables with error
  catWithErr = c()
}



# Combine no error with error
catVars = cbind(catNoErr,catWithErr)

return(
  list(
    trueID = trueID
   ,trueMus = trueMus
   ,conVars = conVars
   ,errVariance = errVariance
   ,popProbsNoErr = data.frame(popProb1=popProb1,popProb2=popProb2, stringsAsFactors = TRUE)
   ,popProbsWithErr = data.frame(popProb1Err=popProb1Err,popProb2Err=popProb2Err, stringsAsFactors = TRUE)
   ,catVars = catVars
  )
)

}


