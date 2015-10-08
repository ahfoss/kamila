
# A function that generates confusion matrices
# for several different categorical measurement
# error scenarios. A helper function for genMixedData.r

# Input arguments:
#
# nLev: number of categorical levels
# epsilon: for catErrType 'diag','adjacent','early';
#    given a particular categorical level, the probability
#    of an incorrect answer (approximately)
# k: 
#  if catErrType = 'adjacent' the number of categories
#    on either side of the correct with nonzero conditional
#    probability.
#  if catErrType = 'early', the number of choices with
#    nonzero conditional probability of incorrect responses.
#
# catErrType = c('all','none','diag','adjacent','early')
#    Type of measurement error in the categorical data
#
# Error types:
#  all: each entry equal to 1/nLev
#  none: identity matrix
#  diag: Main diagonal equals (1-epsilon), all other entries
#        are equal to epsilon/(nLev-1)
#  adjacent: epsilon is spread out k categories to the left
#        and right of the true category.
#  early: first k categories are classified correctly with
#         probability 1. All remaining categories are mistaken
#         for the first k categories with probability epsilon/k,
#         correct is chosen with prob 1-epsilon. 
#  
#
# For condition 'adjacent':
# Note that if epsilon > 2k/(1+2k), then the conditional probability
# of the correct level is smaller than any conditional prob
# of other levels (except for first few and last few cases, possibly).

genConfMat = function(
  catErrType
 ,nLev
 ,epsilon
 ,k=NA
) {

  if (!is.na(k)) {
    if (k < 1) stop('Error in genConfMat: k must be >= 1')
  }

  if (catErrType == 'all') {
    return(matrix(rep(1/nLev,nLev^2),nrow=nLev))


  } else if (catErrType == 'none') {
    return(diag(nLev))


  } else if (catErrType == 'diag') {
    mat = epsilon/(nLev-1) * matrix(1,nLev,nLev)
    diag(mat) = 1-epsilon
    return(mat)


  } else if (catErrType == 'early') {
    mat = cbind(
      rbind(
        diag(k)
       ,matrix(epsilon/k,nrow=nLev-k,ncol=k)
      )
     ,rbind(
        matrix(0,nrow=k,ncol=nLev-k)
       ,diag(nLev-k)*(1-epsilon)
      )
    )
    return(mat)


  } else if (catErrType == 'adjacent') {
    if (k >= nLev) {
      stop('Error in genConfMat: k must be < nLev for
            adjacent condition')
    }
    if (epsilon > 2*k/(1+2*k)) {
      warning('genConfMat: strange behavior when
               epsilon > 2k/(1+2k) for adjacent condition')
    }
    # build off-diagonals
    offDiagErr = epsilon/(2*k)
    mat = matrix(0,nrow=nLev,ncol=nLev)
    for (i in 1:k){
      mat[seq(i*nLev+1,nLev^2-i,nLev+1)] = offDiagErr
      mat[seq(i+1,(nLev^2-nLev*i),nLev+1)] = offDiagErr
    }

    # build diagonals: dropped terms on the left
    tmp1 = c(k:1,rep(0,max(0,nLev-k)))
    tmp2 = tmp1[1:nLev]

    diag(mat) = (tmp2+rev(tmp2)) * offDiagErr + 1 - epsilon

    return(mat)

  } else {
    stop(
      paste('Unknown error distribution specified in catErrType: ',catErrType)
    )
  }
}

# Debugging code:
#mm = genConfMat(catErrType='adjacent',nLev=4,k=3,epsilon=.4); mm; image(mm); rowSums(mm)
