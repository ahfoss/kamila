#include <Rcpp.h>
using namespace Rcpp;

// https://github.com/RcppCore/Rcpp/issues/636
void R_init_kamila(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

// [[Rcpp::export]]
NumericMatrix dptm(
  NumericMatrix pts
 ,NumericMatrix myMeans
 ,NumericVector wgts
 ,int ppDim
 ,int kkMean
 ,int nn
)
{
  NumericMatrix outMat(nn,kkMean);

  for (int i=0; i < nn; i++){
    for (int j=0; j<kkMean; j++){
      double distij2 = 0;
      for (int p=0; p<ppDim; p++){
        distij2 += pow(
          wgts[p] * (pts(i,p) - myMeans(j,p))
         ,2
        );
      }
      outMat(i,j) = sqrt(distij2);
    }
  }
  return(outMat);
}

// [[Rcpp::export]]
NumericVector rowMax( NumericMatrix inMat )
{
  int nn = inMat.nrow(), pp = inMat.ncol();
  NumericVector outVec(nn);

  for (int i=0; i<nn; i++) {
    outVec[i] = inMat(i,0);
    // note j starts at 1 not zero
    for (int j=1; j<pp; j++) {
      if (inMat(i,j) > outVec[i]) {
        outVec[i] = inMat(i,j);
      }
    }
  }
  return(outVec);
}

// [[Rcpp::export]]
NumericVector rowMin( NumericMatrix inMat )
{
  int nn = inMat.nrow(), pp = inMat.ncol();
  NumericVector outVec(nn);

  for (int i=0; i<nn; i++) {
    outVec[i] = inMat(i,0);
    // note j starts at 1 not zero
    for (int j=1; j<pp; j++) {
      if (inMat(i,j) < outVec[i]) {
        outVec[i] = inMat(i,j);
      }
    }
  }
  return(outVec);
}

// [[Rcpp::export]]
NumericVector rowMaxInds( NumericMatrix inMat )
{
  int nn = inMat.nrow(), pp = inMat.ncol();
  NumericVector outVec(nn);

  for (int i=0; i<nn; i++) {
    outVec[i] = 1;
    // note j starts at 1 not zero
    for (int j=1; j<pp; j++) {
      if (inMat(i,j) > inMat(i,outVec[i]-1)) {
        outVec[i] = j+1;
      }
    }
  }
  return(outVec);
}



// [[Rcpp::export]]
NumericMatrix sumMatList( List x )
{
  int qq = x.size();
  NumericMatrix mat0 = x[0];
  int nn = mat0.nrow();
  int kk = mat0.ncol();

  // Note Rcpp initializes with zeroes
  NumericMatrix outMatrix(nn,kk);

  for (int q=0; q<qq; q++) {
    NumericMatrix qthMat = x[q];
    for (int n=0; n<nn; n++) {
      for (int k=0; k<kk; k++) {
        outMatrix(n,k) += qthMat(n,k);
        //outMatrix(n,k) = outMatrix(n,k) + qthMat(n,k);
      }
    }
  }

  return(outMatrix);
}

/* Pseudocode: getIndividualLogProbs
 * For q=0 where q < numCatVar
 *   ithVarCodes <- catFactor[,q]
 *   ithLogLiks <- logProbsCond_i[q]
 *   for (n in 1:nn-1) {
 *     for (k in 0:kk-1) {
 *       ithOutMat[n,k] <- catWeights[q] * ithLogLiks(k,ithVarCodes[n]-1)
 *     }
 *   }
 *   outList[q] <- ithOutMat
 * }
 */

// [[Rcpp::export]]
List getIndividualLogProbs(
  NumericMatrix catFactorNum
 ,NumericVector catWeights
 ,List logProbsCond_i
)
{
  int qq = catWeights.size();
  int nn = catFactorNum.nrow();
  NumericMatrix logProbs0 = logProbsCond_i(0);
  int kk = logProbs0.nrow();
  List outList(qq);

  for (int q=0; q<qq; q++) {
    NumericMatrix ithOutMat(nn,kk);
    NumericMatrix::Column ithVarCodes = catFactorNum(_,q);
    NumericMatrix ithLogLiks = logProbsCond_i[q];
    for (int n=0; n<nn; n++){
      for (int k=0; k<kk; k++) {
        ithOutMat(n,k) = catWeights[q] * ithLogLiks(k,ithVarCodes[n]-1);
      }
    }
    outList[q] = ithOutMat;
  }
  return(outList);
}



/* Pseudocode: aggregateMeans
 *
 * Input: conVar; n X p numeric matrix
 *        membNew; n X 1 numeric vector
 *
 * Procedure:
 * initialize outMat(kk,pp)
 * initialize countVec(kk)
 * for (n in 0:(nn-1)
 *   for (p in 0:(pp-1))
 *     outMat(membNew[n]-1,p) += conVar(n,p)
 *     countVec[membNew[n]-1] += 1
 * for (k in 0:(kk-1))
 *   for (p in 0:(pp-1))
 *     outMat(k,p) /= countVec[k]
 *
 * Output: k X p numeric matrix of means
 */

// [[Rcpp::export]]
NumericMatrix aggregateMeans(
  NumericMatrix conVar
 ,IntegerVector membNew
 ,int kk
)
{
  int pp = conVar.ncol(), nn = conVar.nrow();
  NumericVector countVec(kk);
  NumericMatrix outMat(kk,pp);

  for (int n=0; n<nn; n++) {
    countVec[membNew[n]-1] += 1;
    for (int p=0; p<pp; p++) {
      outMat(membNew[n]-1,p) += conVar(n,p);
    }
  }
  for (int k=0; k<kk; k++) {
    if (countVec[k] != 0) {
      for (int p=0; p<pp; p++) {
        outMat(k,p) /= countVec[k];
      }
    }
  }

  return(outMat);
}



/* Helper function for jointTabSmoothedList
 * implements table function in Rcpp
 * for two integer vectors with known number
 * of categories. Both must be coded 1:nc, where
 * nc is number of categories. Both must be same length.
 */
IntegerMatrix tabulateTwoIntVec(
  IntegerVector vec1
 ,IntegerVector vec2
 ,int nc1
 ,int nc2
 ,int nn
)
{
  IntegerMatrix outMat(nc1,nc2);
  for (int i=0; i<nn; i++) {
    outMat(vec1[i]-1,vec2[i]-1) += 1;
  }
  return(outMat);
}

/* Helper function for jointTabSmoothedList
 * function for categorical kernels
 * Replicates the results of np::npudens, except it's
 * much faster, and limited to two-dimensional tables.
 * See np package in R for more details.
 * Also, note that bandwidth is a scalar, and is replicated
 * for both dimensions of the table.
 */
NumericMatrix smooth2dTable(
  IntegerMatrix inputTab
 ,double catBw
 ,int nn
)
{
  int dim1 = inputTab.nrow(), dim2 = inputTab.ncol();
  NumericMatrix midMat(dim1,dim2);
  NumericMatrix outMat(dim1,dim2);

  // get colsums of original matrix
  IntegerVector colSums(dim2);
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      colSums[j] += inputTab(i,j);
    }
  }

  // Begin construction of output matrix
  // weighted sum of original matrix and offcounts
  int offCounts1;
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      offCounts1 = colSums[j] - inputTab(i,j);
      midMat(i,j) = (1-catBw)*inputTab(i,j) + catBw/(dim1-1)*offCounts1;
    }
  }

  // Get rowsums of new matrix
  NumericVector rowSums(dim1);
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      rowSums[i] += midMat(i,j);
    }
  }

  // Construct final matrix
  // Weighted sum of current outMat and offCounts2
  double offCounts2;
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      offCounts2 =  rowSums[i] - midMat(i,j);
      outMat(i,j) = (1-catBw)*midMat(i,j) + catBw/(dim2-1)*offCounts2;
    }
  }

  return(outMat);
}

/* Pseudocode: jointTabSmoothedList
 * Input: catFactorNum, integer coded matrix n X q factor variables
 *        membNew, integer coded nX1 vector of cluster memberships
 *        numLev, qX1 vector giving number of level codes for each variable
 *                (needed in case of missing levels)
 *        catBw, bandwidth for categorical kernel
 *        kk, number of clusters (needed in case of empty clusters)
 * Output: List of length Q, each element k X l_q matrix of smoothed table
 *
 * Procedure: 
 *   for q categorical variables:
 *     get crosstab of membNew and qth categorical variable
 *     smooth table using input bandwidth
 */

// [[Rcpp::export]]
List jointTabSmoothedList(
  IntegerMatrix catFactorNum
 ,IntegerVector membNew
 ,IntegerVector numLev
 ,double catBw
 ,int kk
)
{
  int qq = catFactorNum.ncol(), nn = catFactorNum.nrow();
  List outList(qq);
  for (int q=0; q<qq; q++) {
    IntegerMatrix::Column qthVar = catFactorNum(_,q);
    IntegerMatrix qthTabRaw = tabulateTwoIntVec(
      membNew,qthVar,kk,numLev[q],nn
    );
    if (catBw != 0) {
      outList(q) = smooth2dTable(qthTabRaw,catBw,nn);
    }
  }

  return(outList);
}

