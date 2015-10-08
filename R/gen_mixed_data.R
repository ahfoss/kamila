
###########################################################
# genMixedData.r
#
# Generate a dataset of mixed variables with underlying 
# cluster structure.
# ------------------------------------------------------
#
# Distribution of true data is currently independent normal.
# Separation is currently fixed at 1%.
#
# Error distribution is currently independent normal.
#
# Error distribution is currently equal among all continuous
# variables measured with error, and equal among all
# categorical variables measured with error.
#
# ------------------------------------------------------
# Input args:
#
# sampSize = {integer}
#    'N', that is, number of draws from the mixture distribution.
#
# popProportions = {vector of scalars that sums to 1}
#    Simultaneously specify number and frequency of populations.
#    length(popProportions) gives number of underlying populations,
#    with each element corresponding to the prior probability
#    of an observation being drawn from that population.
#    NOTE THAT CURRENTLY ONLY 2 POPULATIONS ARE SUPPORTED
#
# conErrDist = 'mvn'
#    Continuous error distribution. Currently only independent
#    Normal supported.
# 
# nConVar = {integer}
#    Number of continuous variables
#
# nConWithErr = {integer <= nConVar} or {boolean vector}
#    Number of continuous variables with error, or vector
#    specifying which categorical variables have measurement error.
#
# conErrLev = {numeric \in [0.01,1]}
#    Amount of overlap introduced by measurement error of
#    continuous variables.
#
# nCatVar = {integer}
#    Number of categorical variables
#
# nCatLevels = {integer}
#    A scalar with the number of categories per categorical
#    Variable. Must be a multiple of 2.
#    In the future, could have a vector of level numbers,
#    but assigning populations to categorical levels is
#    ambiguous when nCatLevels is not a multiple of the number
#    of populations, so more detailed input or more extensive
#    assumptions would be needed.
#
# nCatWithErr = {integer <= nCatVar}
#    Number of categorical variables with error
#
# catErrType = c('none','all','adjacent','early')
#    Type of measurement error in the categorical data.
#    Option 'none' bypasses confusion matrices entirely.
#
# catErrParams
#    List of parameters corresponding to catErrType. Option 'all'
#    doesn't require any parameters, option 'adjacent' and
#    'early' require a tuple c(epsilon,k), where epsilon
#    is \in [0,1] and k is an integer. See genConfMat documentation
#    for details.
#


# TO BE MOVED TO SEPARATE FUNCTION 
# Now always output as matrix; should be converted later as needed
# since it doesn't save time to do it here.

# output = c('matrix','dfFactor')
#    Should data be output as matrices (numeric for continuous,
#    character for categorical with integer values); or output
#    as data frames (numeric vars for continuous, factors for
#    categorical).
###########################################################


genMixedData = function(
  sampSize
 ,nConVar
 ,nCatVar
 ,nCatLevels
 ,nConWithErr
 ,nCatWithErr
 ,popProportions
 ,conErrDist = 'mvn'
 ,conErrLev
 ,catErrType
 ,catErrParams
) {

# Function that generates confusion matrices
source('genConfMat.r')

#####################################
# ARgument checking
#####################################

# Only two populations currently supported
if (length(popProportions) > 2) {
  # MixSim will need to be used
  stop('Error in genMixedData: More than two populations not currently supported')
}

# Input continuous distribution
if (conErrDist != 'mvn') {
  stop('Error in genMixedData: conErrDist must be mvn')
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
OVERLAP = 0.01 # overlap

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
trueMus = c(0,-2*qnorm(OVERLAP/2))
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
rightCatProb = (1-OVERLAP/2) / (nCatLevels/2)
wrongCatProb = (OVERLAP/2) / (nCatLevels/2)
popProb1 = rep(c(rightCatProb,wrongCatProb),each=nCatLevels/2)
popProb2 = rep(c(wrongCatProb,rightCatProb),each=nCatLevels/2)
# NOTE: Currently fixed at 1% overlap

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
  # use genConfMat
  confMat = genConfMat(
    catErrType = catErrType
   ,nLev = nCatLevels
   ,epsilon = catErrParams$epsilon
   ,k = catErrParams$k
  )

  # Note currently assuming uniform probability
  # for each categorical level within a population
  popProbErr1 = c(popProb1 %*% confMat)
  popProbErr2 = c(popProb2 %*% confMat)
  
  # Simulate pops 1 and 2 with error
  pop1CatWithErr = sample(
    x = 1:nCatLevels
   ,size = numInPop1*nCatWithErr
   ,replace = T
   ,prob = popProbErr1
  )
  pop2CatWithErr = sample(
    x = 1:nCatLevels
   ,size = (sampSize - numInPop1)*nCatWithErr
   ,replace = T
   ,prob = popProbErr2
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


#### MOVE to SEPARATE FUNCTiON
#### Whether to output as factor or character
# Separate function for converting to dummy variables
# use data.frame(...,stringsAsFactors)

return(
  list(
    trueID = trueID
   ,trueMus = trueMus
   ,conVars = conVars
   ,errVariance = errVariance
   ,popProbsNoErr = data.frame(popProb1=popProb1,popProb2=popProb2)
   ,popProbsWithErr = data.frame(popProbErr1=popProbErr1,popProbErr2=popProbErr2)
   ,catVars = catVars
   ,confMat = confMat
  )
)

}


