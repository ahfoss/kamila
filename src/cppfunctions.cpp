#include <Rcpp.h>
using namespace Rcpp;

// https://github.com/RcppCore/Rcpp/issues/636
// #nocov start
void R_init_kamila(DllInfo* info) {
  R_registerRoutines(info, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}
// #nocov end

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
  NumericMatrix outMat(nn, kkMean);
  const double* p_pts = pts.begin();
  const double* p_means = myMeans.begin();
  const double* p_wgts = wgts.begin();
  double* p_out = outMat.begin();

  for (int j = 0; j < kkMean; ++j) {
    double* out_col = p_out + j * nn;
    std::fill(out_col, out_col + nn, 0.0);

    for (int p = 0; p < ppDim; ++p) {
      double w = p_wgts[p];
      if (w == 0.0) continue;
      double m = p_means[j + p * kkMean];
      const double* pts_col = p_pts + p * nn;

      if (w == 1.0) {
        for (int i = 0; i < nn; ++i) {
          double diff = pts_col[i] - m;
          out_col[i] += diff * diff;
        }
      } else {
        for (int i = 0; i < nn; ++i) {
          double diff = w * (pts_col[i] - m);
          out_col[i] += diff * diff;
        }
      }
    }

    for (int i = 0; i < nn; ++i) {
      out_col[i] = std::sqrt(out_col[i]);
    }
  }

  return outMat;
}

// [[Rcpp::export]]
NumericVector rowMax( NumericMatrix inMat )
{
  int nn = inMat.nrow(), pp = inMat.ncol();
  NumericVector outVec(nn);
  const double* p_in = inMat.begin();
  double* p_out = outVec.begin();

  for (int i=0; i<nn; i++) {
    p_out[i] = p_in[i];
  }
  for (int j=1; j<pp; j++) {
    const double* col_ptr = p_in + j * nn;
    for (int i=0; i<nn; i++) {
      if (col_ptr[i] > p_out[i]) {
        p_out[i] = col_ptr[i];
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
  const double* p_in = inMat.begin();
  double* p_out = outVec.begin();

  for (int i=0; i<nn; i++) {
    p_out[i] = p_in[i];
  }
  for (int j=1; j<pp; j++) {
    const double* col_ptr = p_in + j * nn;
    for (int i=0; i<nn; i++) {
      if (col_ptr[i] < p_out[i]) {
        p_out[i] = col_ptr[i];
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
  std::vector<double> maxVals(nn);
  const double* p_in = inMat.begin();
  double* p_out = outVec.begin();

  for (int i=0; i<nn; i++) {
    p_out[i] = 1.0;
    maxVals[i] = p_in[i];
  }
  for (int j=1; j<pp; j++) {
    const double* col_ptr = p_in + j * nn;
    double colIdx = j + 1.0;
    for (int i=0; i<nn; i++) {
      if (col_ptr[i] > maxVals[i]) {
        maxVals[i] = col_ptr[i];
        p_out[i] = colIdx;
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
NumericMatrix calcCatLogLiks(
  IntegerMatrix catFactorNum
 ,NumericVector catWeights
 ,List logProbsCond_i
)
{
  int qq = catWeights.size();
  int nn = catFactorNum.nrow();
  NumericMatrix logProbs0 = logProbsCond_i[0];
  int kk = logProbs0.nrow();

  NumericMatrix outMat(nn, kk);
  double* p_out = outMat.begin();
  const int* p_cat = catFactorNum.begin();

  // Pre-calculate weighted lookups for each q: table of size [nlev * kk]
  std::vector<std::vector<double>> weightedTabs(qq);
  for (int q = 0; q < qq; ++q) {
    double w = catWeights[q];
    NumericMatrix mat = logProbsCond_i[q];
    int nlev = mat.ncol();
    weightedTabs[q].resize(kk * nlev);
    const double* lp = mat.begin();
    for (int lev = 0; lev < nlev; ++lev) {
      for (int cl = 0; cl < kk; ++cl) {
        weightedTabs[q][lev * kk + cl] = w * lp[cl + lev * kk];
      }
    }
  }

  std::vector<double*> out_cols(kk);
  for (int cl = 0; cl < kk; ++cl) {
    out_cols[cl] = p_out + cl * nn;
  }

  for (int q = 0; q < qq; ++q) {
    if (catWeights[q] == 0.0) continue;
    const int* q_col = p_cat + q * nn;
    const double* w_tab = weightedTabs[q].data();

    if (kk == 4) {
      double* c0 = out_cols[0];
      double* c1 = out_cols[1];
      double* c2 = out_cols[2];
      double* c3 = out_cols[3];
      for (int i = 0; i < nn; ++i) {
        const double* w_row = w_tab + (q_col[i] - 1) * 4;
        c0[i] += w_row[0];
        c1[i] += w_row[1];
        c2[i] += w_row[2];
        c3[i] += w_row[3];
      }
    } else if (kk == 2) {
      double* c0 = out_cols[0];
      double* c1 = out_cols[1];
      for (int i = 0; i < nn; ++i) {
        const double* w_row = w_tab + (q_col[i] - 1) * 2;
        c0[i] += w_row[0];
        c1[i] += w_row[1];
      }
    } else if (kk == 3) {
      double* c0 = out_cols[0];
      double* c1 = out_cols[1];
      double* c2 = out_cols[2];
      for (int i = 0; i < nn; ++i) {
        const double* w_row = w_tab + (q_col[i] - 1) * 3;
        c0[i] += w_row[0];
        c1[i] += w_row[1];
        c2[i] += w_row[2];
      }
    } else {
      for (int i = 0; i < nn; ++i) {
        const double* w_row = w_tab + (q_col[i] - 1) * kk;
        for (int cl = 0; cl < kk; ++cl) {
          out_cols[cl][i] += w_row[cl];
        }
      }
    }
  }

  return outMat;
}

// [[Rcpp::export]]
NumericMatrix aggregateMeans(
  NumericMatrix conVar
 ,IntegerVector membNew
 ,int kk
)
{
  int pp = conVar.ncol(), nn = conVar.nrow();
  std::vector<double> countVec(kk, 0.0);
  const int* p_memb = membNew.begin();
  for (int n = 0; n < nn; ++n) {
    countVec[p_memb[n] - 1] += 1.0;
  }

  NumericMatrix outMat(kk, pp);
  double* p_out = outMat.begin();
  const double* p_con = conVar.begin();

  for (int p = 0; p < pp; ++p) {
    const double* col_ptr = p_con + p * nn;
    double* out_col = p_out + p * kk;
    for (int n = 0; n < nn; ++n) {
      out_col[p_memb[n] - 1] += col_ptr[n];
    }
    for (int k = 0; k < kk; ++k) {
      if (countVec[k] != 0.0) {
        out_col[k] /= countVec[k];
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
    } else {
      NumericMatrix numTab(kk, numLev[q]);
      for (int i = 0; i < kk; ++i) {
        for (int j = 0; j < numLev[q]; ++j) {
          numTab(i, j) = qthTabRaw(i, j);
        }
      }
      outList(q) = numTab;
    }
  }

  return(outList);
}

// [[Rcpp::export]]
List updateLogProbs(
  IntegerMatrix catFactorNum
 ,IntegerVector membNew
 ,IntegerVector numLev
 ,double catBw
 ,int kk
)
{
  int qq = catFactorNum.ncol();
  int nn = catFactorNum.nrow();
  const int* p_memb = membNew.begin();
  const int* p_cat = catFactorNum.begin();

  List outList(qq);

  for (int q = 0; q < qq; ++q) {
    int nlev = numLev[q];
    const int* q_col = p_cat + q * nn;

    // 1. Tabulate
    std::vector<int> rawTab(kk * nlev, 0);
    for (int i = 0; i < nn; ++i) {
      int row = p_memb[i] - 1;
      int col = q_col[i] - 1;
      rawTab[row + col * kk] += 1;
    }

    // 2. Smooth 2D table
    std::vector<double> outMat(kk * nlev, 0.0);
    if (catBw != 0.0) {
      std::vector<double> midMat(kk * nlev, 0.0);
      std::vector<int> colSums(nlev, 0);
      for (int j = 0; j < nlev; ++j) {
        for (int i = 0; i < kk; ++i) {
          colSums[j] += rawTab[i + j * kk];
        }
      }

      double bw_div_k = (kk > 1) ? (catBw / (kk - 1.0)) : 0.0;
      double one_minus_bw = 1.0 - catBw;
      for (int j = 0; j < nlev; ++j) {
        for (int i = 0; i < kk; ++i) {
          int offCounts1 = colSums[j] - rawTab[i + j * kk];
          midMat[i + j * kk] = one_minus_bw * rawTab[i + j * kk] + bw_div_k * offCounts1;
        }
      }

      std::vector<double> rowSums(kk, 0.0);
      for (int j = 0; j < nlev; ++j) {
        for (int i = 0; i < kk; ++i) {
          rowSums[i] += midMat[i + j * kk];
        }
      }

      double bw_div_lev = (nlev > 1) ? (catBw / (nlev - 1.0)) : 0.0;
      for (int j = 0; j < nlev; ++j) {
        for (int i = 0; i < kk; ++i) {
          double offCounts2 = rowSums[i] - midMat[i + j * kk];
          outMat[i + j * kk] = one_minus_bw * midMat[i + j * kk] + bw_div_lev * offCounts2;
        }
      }
    } else {
      for (size_t idx = 0; idx < rawTab.size(); ++idx) {
        outMat[idx] = static_cast<double>(rawTab[idx]);
      }
    }

    // 3. Normalize by rowSums and take log
    NumericMatrix logProbMat(kk, nlev);
    std::vector<double> finalRowSums(kk, 0.0);
    for (int j = 0; j < nlev; ++j) {
      for (int i = 0; i < kk; ++i) {
        finalRowSums[i] += outMat[i + j * kk];
      }
    }

    for (int j = 0; j < nlev; ++j) {
      for (int i = 0; i < kk; ++i) {
        double denom = finalRowSums[i];
        logProbMat(i, j) = (denom > 0.0) ? std::log(outMat[i + j * kk] / denom) : R_NegInf;
      }
    }

    outList[q] = logProbMat;
  }

  return(outList);
}

// [[Rcpp::export]]
NumericVector interpRadialKde(
  NumericVector y,
  double maxEval,
  int pdim,
  NumericVector evalPoints,
  bool takeLog = false
)
{
  int m = 401;
  double h = (maxEval > 0.0) ? (maxEval / (m - 1.0)) : 1.0;
  std::vector<double> x(m);
  for (int i = 0; i < m; ++i) {
    x[i] = i * h;
  }

  // 1. remove any zero and negative density estimates
  std::vector<double> newY(m);
  double minPos = 1e300;
  for (int i = 0; i < m; ++i) {
    if (y[i] > 0.0 && y[i] < minPos) minPos = y[i];
  }
  for (int i = 0; i < m; ++i) {
    newY[i] = (y[i] > 0.0) ? y[i] : (minPos / 100.0);
  }

  // 2. at bottom 5th percentile, replace with line through (0,0) and (q05, f(q05))
  // For 401 equally spaced points from 0, index 21 (0-based 20) is 5th percentile
  // coordsLtQ05 is 0..19, maxPt is 19
  double slope = (x[19] > 0.0) ? (newY[19] / x[19]) : 0.0;
  for (int i = 0; i < 20; ++i) {
    newY[i] = x[i] * slope;
  }

  // 3. radial Jacobian transformation; up to proportionality constant
  std::vector<double> radY(m);
  for (int i = 1; i < m; ++i) {
    radY[i] = newY[i] / std::pow(x[i], pdim - 1);
  }
  radY[0] = radY[1];

  // 4. replace densities over MAXDENS with MAXDENS (MAXDENS = 1.0)
  double sumRadY = 0.0;
  for (int i = 0; i < m; ++i) {
    if (radY[i] > 1.0) radY[i] = 1.0;
    sumRadY += radY[i];
  }

  // 5. normalize to area 1
  double minDensR = 1e300;
  std::vector<double> densR(m);
  double normFactor = h * sumRadY;
  for (int i = 0; i < m; ++i) {
    densR[i] = (normFactor > 0.0) ? (radY[i] / normFactor) : 0.0;
    if (densR[i] < minDensR) minDensR = densR[i];
  }

  std::vector<double> diffDensR(m - 1);
  for (int i = 0; i < m - 1; ++i) {
    diffDensR[i] = densR[i + 1] - densR[i];
  }

  // 6. linear interpolation at evalPoints with rule 1:2 and pmax(..., min(densR))
  int nEval = evalPoints.size();
  NumericVector kdes(nEval);
  const double* p_eval = evalPoints.begin();
  double* p_kdes = kdes.begin();

  double inv_h = (h > 0.0) ? (1.0 / h) : 0.0;
  if (takeLog) {
    double logVal0 = std::log((densR[0] > minDensR) ? densR[0] : minDensR);
    double logValMax = std::log((densR[m - 1] > minDensR) ? densR[m - 1] : minDensR);
    for (int i = 0; i < nEval; ++i) {
      double u = p_eval[i];
      if (u <= 0.0) {
        p_kdes[i] = logVal0;
      } else if (u >= maxEval) {
        p_kdes[i] = logValMax;
      } else {
        double pos = u * inv_h;
        int idx = static_cast<int>(pos);
        if (idx >= m - 1) idx = m - 2;
        double frac = pos - idx;
        double val = densR[idx] + frac * diffDensR[idx];
        p_kdes[i] = std::log((val > minDensR) ? val : minDensR);
      }
    }
  } else {
    double val0 = (densR[0] > minDensR) ? densR[0] : minDensR;
    double valMax = (densR[m - 1] > minDensR) ? densR[m - 1] : minDensR;
    for (int i = 0; i < nEval; ++i) {
      double u = p_eval[i];
      if (u <= 0.0) {
        p_kdes[i] = val0;
      } else if (u >= maxEval) {
        p_kdes[i] = valMax;
      } else {
        double pos = u * inv_h;
        int idx = static_cast<int>(pos);
        if (idx >= m - 1) idx = m - 2;
        double frac = pos - idx;
        double val = densR[idx] + frac * diffDensR[idx];
        p_kdes[i] = (val > minDensR) ? val : minDensR;
      }
    }
  }

  if (evalPoints.hasAttribute("dim")) {
    kdes.attr("dim") = evalPoints.attr("dim");
  }

  return kdes;
}

/* Pseudocode & Note: calcPsCpp
 * Input: testMemb, integer coded vector (1-indexed) of test cluster assignments
 *        teIntoTr, integer coded vector (1-indexed) of predicted test cluster assignments
 *        numClust, number of clusters
 * Output: NumericVector of prediction strength proportions psProps for each cluster.
 *
 * Mathematical Note:
 * Evaluates prediction strength (Tibshirani & Walther, 2005) for each cluster k:
 *   ps(k) = (1 / choose(n_k, 2)) * sum_{i < j in A_k} I(hat(y)_i == hat(y)_j)
 * By partitioning cluster A_k into predicted class counts C_{k, m}, the number of
 * agreeing pairs in A_k is sum_{m} choose(C_{k, m}, 2) = sum_{m} C_{k, m}(C_{k, m} - 1) / 2.
 * Thus:
 *   ps(k) = sum_{m} [ C_{k, m} * (C_{k, m} - 1) ] / [ n_k * (n_k - 1) ]
 * This computes the exact combinatorial identity in O(N + K^2) time and O(K^2) memory,
 * producing mathematically and numerically identical results to the legacy O(N^2) pairwise
 * distance matrix approach without allocating an N x N matrix or performing quadratic loops.
 */

// [[Rcpp::export]]
NumericVector calcPsCpp(
  IntegerVector testMemb,
  IntegerVector teIntoTr,
  int numClust
)
{
  int n = testMemb.size();
  int maxPred = 0;
  for (int i = 0; i < n; ++i) {
    if (teIntoTr[i] > maxPred) {
      maxPred = teIntoTr[i];
    }
  }
  if (maxPred == 0 || numClust <= 0) {
    return NumericVector(numClust, NA_REAL);
  }

  std::vector<double> counts(numClust * maxPred, 0.0);
  std::vector<double> clustSize(numClust, 0.0);

  for (int i = 0; i < n; ++i) {
    int tm = testMemb[i] - 1;
    int tr = teIntoTr[i] - 1;
    if (tm >= 0 && tm < numClust && tr >= 0 && tr < maxPred) {
      counts[tm * maxPred + tr] += 1.0;
      clustSize[tm] += 1.0;
    }
  }

  NumericVector psProps(numClust);
  for (int cl = 0; cl < numClust; ++cl) {
    double n_cl = clustSize[cl];
    if (n_cl < 2.0) {
      psProps[cl] = NA_REAL;
    } else {
      double nPairsAgree = 0.0;
      for (int m = 0; m < maxPred; ++m) {
        double cnt = counts[cl * maxPred + m];
        if (cnt >= 2.0) {
          nPairsAgree += cnt * (cnt - 1.0) / 2.0;
        }
      }
      double totalPairs = n_cl * (n_cl - 1.0) / 2.0;
      psProps[cl] = nPairsAgree / totalPairs;
    }
  }

  return psProps;
}

// [[Rcpp::export]]
List kamilaLoopCpp(
  Nullable<NumericMatrix> conVarMat_,
  Nullable<IntegerMatrix> catFactorNum_,
  NumericVector conWeights,
  NumericVector catWeights,
  Nullable<NumericMatrix> initMeans_,
  Nullable<List> initLogProbs_,
  IntegerVector numLev,
  double catBw,
  int numClust,
  int maxIter,
  bool verbose,
  bool hasCon,
  bool hasCat
)
{
  int nn = 0;
  int pp = 0;
  int qq = 0;

  NumericMatrix conVarMat;
  if (hasCon) {
    conVarMat = NumericMatrix(conVarMat_.get());
    nn = conVarMat.nrow();
    pp = conVarMat.ncol();
  }

  IntegerMatrix catFactorNum;
  if (hasCat) {
    catFactorNum = IntegerMatrix(catFactorNum_.get());
    nn = catFactorNum.nrow();
    qq = catFactorNum.ncol();
  }

  const double* p_con = hasCon ? conVarMat.begin() : nullptr;
  const int* p_cat = hasCat ? catFactorNum.begin() : nullptr;
  const double* p_con_wgts = hasCon ? conWeights.begin() : nullptr;

  // Pre-allocated scratch buffers (reused across all iterations)
  std::vector<double> distMat(hasCon ? (nn * numClust) : 0, 0.0);
  std::vector<double> minDist(hasCon ? nn : 0, 0.0);
  std::vector<double> catLogLiks(hasCat ? (nn * numClust) : 0, 0.0);
  std::vector<double> allLogLiks(nn * numClust, 0.0);
  std::vector<int> membOld(nn, 0);
  std::vector<int> membNew(nn, 0);
  std::vector<double> countVec(numClust, 0.0);

  std::vector<double> currentMeans(hasCon ? (numClust * pp) : 0, 0.0);
  if (hasCon) {
    NumericMatrix initM(initMeans_.get());
    const double* p_init_m = initM.begin();
    std::copy(p_init_m, p_init_m + numClust * pp, currentMeans.begin());
  }

  std::vector<std::vector<double>> logProbs(qq);
  if (hasCat) {
    List initL(initLogProbs_.get());
    for (int q = 0; q < qq; ++q) {
      NumericMatrix mat = initL[q];
      int nlev = numLev[q];
      logProbs[q].resize(numClust * nlev);
      const double* p_lp = mat.begin();
      std::copy(p_lp, p_lp + numClust * nlev, logProbs[q].begin());
    }
  }

  // Pre-allocated KDE workspace (m_grid = 401)
  const int m_grid = 401;
  std::vector<double> r_sorted(hasCon ? nn : 0);
  std::vector<double> gcounts(m_grid);
  std::vector<double> kappa(m_grid + 1);
  std::vector<double> y_kde(m_grid);
  std::vector<double> newY(m_grid);
  std::vector<double> radY(m_grid);
  std::vector<double> densR(m_grid);
  std::vector<double> diffDensR(m_grid - 1);

  // Pre-allocated Categorical update workspace
  std::vector<std::vector<double>> weightedTabs(qq);
  if (hasCat) {
    for (int q = 0; q < qq; ++q) {
      weightedTabs[q].resize(numLev[q] * numClust);
    }
  }
  int max_nlev = 0;
  for (int q = 0; q < qq; ++q) {
    if (numLev[q] > max_nlev) max_nlev = numLev[q];
  }
  std::vector<int> rawTab(numClust * max_nlev);
  std::vector<double> outMat(numClust * max_nlev);
  std::vector<double> midMat(numClust * max_nlev);
  std::vector<int> colSums(max_nlev);
  std::vector<double> rowSums(numClust);
  std::vector<double> finalRowSums(numClust);

  std::vector<IntegerVector> membHistory;

  int numIter = 0;
  bool degenerateSoln = false;

  while (true) {
    bool canStop = (numIter >= 3);
    if (canStop) {
      bool allEqual = true;
      for (int i = 0; i < nn; ++i) {
        if (membOld[i] != membNew[i]) {
          allEqual = false;
          break;
        }
      }
      if (allEqual) break;
    }
    if (numIter >= maxIter) break;

    numIter++;

    // 1. Continuous calculations: dptm, rowMin, radialKDE
    if (hasCon) {
      // (a) dptm
      const double* p_means = currentMeans.data();
      for (int j = 0; j < numClust; ++j) {
        double* out_col = distMat.data() + j * nn;
        std::fill(out_col, out_col + nn, 0.0);

        for (int p = 0; p < pp; ++p) {
          double w = p_con_wgts[p];
          if (w == 0.0) continue;
          double m_val = p_means[j + p * numClust];
          const double* pts_col = p_con + p * nn;

          if (w == 1.0) {
            for (int i = 0; i < nn; ++i) {
              double diff = pts_col[i] - m_val;
              out_col[i] += diff * diff;
            }
          } else {
            for (int i = 0; i < nn; ++i) {
              double diff = w * (pts_col[i] - m_val);
              out_col[i] += diff * diff;
            }
          }
        }

        for (int i = 0; i < nn; ++i) {
          out_col[i] = std::sqrt(out_col[i]);
        }
      }

      // (b) rowMin & maxEval
      double maxEval = distMat[0];
      for (int i = 0; i < nn; ++i) {
        minDist[i] = distMat[i];
      }
      for (int j = 0; j < numClust; ++j) {
        const double* col_ptr = distMat.data() + j * nn;
        for (int i = 0; i < nn; ++i) {
          double d = col_ptr[i];
          if (j > 0 && d < minDist[i]) minDist[i] = d;
          if (d > maxEval) maxEval = d;
        }
      }

      // (c) radial KDE: bw.nrd0
      double mean_r = 0.0;
      for (int i = 0; i < nn; ++i) {
        r_sorted[i] = minDist[i];
        mean_r += minDist[i];
      }
      mean_r /= nn;

      double var_r = 0.0;
      for (int i = 0; i < nn; ++i) {
        double diff = minDist[i] - mean_r;
        var_r += diff * diff;
      }
      double hi = (nn > 1) ? std::sqrt(var_r / (nn - 1.0)) : 0.0;

      std::sort(r_sorted.begin(), r_sorted.end());

      auto get_quantile_7 = [&](double p_prob) {
        double index = 1.0 + (nn - 1.0) * p_prob;
        int lo = static_cast<int>(std::floor(index));
        int hi_idx = static_cast<int>(std::ceil(index));
        double g = index - lo;
        if (lo == hi_idx || lo >= nn) {
          return r_sorted[std::min(lo - 1, nn - 1)];
        } else {
          return (1.0 - g) * r_sorted[lo - 1] + g * r_sorted[hi_idx - 1];
        }
      };

      double q25 = get_quantile_7(0.25);
      double q75 = get_quantile_7(0.75);
      double iqr = q75 - q25;

      double lo = std::min(hi, iqr / 1.34);
      if (lo <= 0.0 || std::isnan(lo)) {
        lo = hi;
        if (lo <= 0.0 || std::isnan(lo)) {
          lo = std::abs(minDist[0]);
          if (lo <= 0.0 || std::isnan(lo)) {
            lo = 1.0;
          }
        }
      }
      double h_bw = 0.9 * lo * std::pow(static_cast<double>(nn), -0.2);

      // (d) 1D linear binning directly into gcounts
      double delta_grid = (maxEval > 0.0) ? (maxEval / (m_grid - 1.0)) : 1.0;
      double inv_delta = (delta_grid > 0.0) ? (1.0 / delta_grid) : 0.0;
      std::fill(gcounts.begin(), gcounts.end(), 0.0);

      for (int i = 0; i < nn; ++i) {
        double r = minDist[i];
        if (r >= 0.0 && r < maxEval) {
          double pos = r * inv_delta;
          int l = static_cast<int>(pos);
          double rem = pos - l;
          if (l >= 0 && l < m_grid - 1) {
            gcounts[l] += (1.0 - rem);
            gcounts[l + 1] += rem;
          }
        }
      }

      // (e) Discrete Gaussian convolution
      double delta = delta_grid / h_bw;
      int L = (delta > 0.0) ? static_cast<int>(std::floor(4.0 / delta)) : m_grid;
      if (L > m_grid) L = m_grid;
      if (L < 0) L = 0;

      double sum_kappa = 0.0;
      for (int l = 0; l <= L; ++l) {
        double z = l * delta;
        kappa[l] = std::exp(-0.5 * z * z) * M_1_SQRT_2PI / (nn * h_bw);
        sum_kappa += (l == 0) ? kappa[l] : (2.0 * kappa[l]);
      }
      double tot = sum_kappa * delta_grid * nn;
      double inv_tot = (tot > 0.0) ? (1.0 / tot) : 0.0;
      for (int l = 0; l <= L; ++l) {
        kappa[l] *= inv_tot;
      }

      for (int i = 0; i < m_grid; ++i) {
        int j_min = std::max(0, i - L);
        int j_max = std::min(m_grid - 1, i + L);
        double val = 0.0;
        for (int j = j_min; j <= j_max; ++j) {
          val += gcounts[j] * kappa[std::abs(i - j)];
        }
        y_kde[i] = val;
      }

      // (f) Radial KDE post-processing
      double minPos = 1e300;
      for (int i = 0; i < m_grid; ++i) {
        if (y_kde[i] > 0.0 && y_kde[i] < minPos) minPos = y_kde[i];
      }
      for (int i = 0; i < m_grid; ++i) {
        newY[i] = (y_kde[i] > 0.0) ? y_kde[i] : (minPos / 100.0);
      }

      double x19 = 19.0 * delta_grid;
      double slope = (x19 > 0.0) ? (newY[19] / x19) : 0.0;
      for (int i = 0; i < 20; ++i) {
        newY[i] = (i * delta_grid) * slope;
      }

      for (int i = 1; i < m_grid; ++i) {
        double xi = i * delta_grid;
        radY[i] = newY[i] / std::pow(xi, pp - 1);
      }
      radY[0] = radY[1];

      double sumRadY = 0.0;
      for (int i = 0; i < m_grid; ++i) {
        if (radY[i] > 1.0) radY[i] = 1.0;
        sumRadY += radY[i];
      }

      double minDensR = 1e300;
      double normFactor = delta_grid * sumRadY;
      for (int i = 0; i < m_grid; ++i) {
        densR[i] = (normFactor > 0.0) ? (radY[i] / normFactor) : 0.0;
        if (densR[i] < minDensR) minDensR = densR[i];
      }

      for (int i = 0; i < m_grid - 1; ++i) {
        diffDensR[i] = densR[i + 1] - densR[i];
      }

      // (g) Write continuous log-densities directly into allLogLiks
      double inv_h = inv_delta;
      double logVal0 = std::log((densR[0] > minDensR) ? densR[0] : minDensR);
      double logValMax = std::log((densR[m_grid - 1] > minDensR) ? densR[m_grid - 1] : minDensR);

      int totalPoints = nn * numClust;
      for (int i = 0; i < totalPoints; ++i) {
        double u = distMat[i];
        if (u <= 0.0) {
          allLogLiks[i] = logVal0;
        } else if (u >= maxEval) {
          allLogLiks[i] = logValMax;
        } else {
          double pos = u * inv_h;
          int idx = static_cast<int>(pos);
          if (idx >= m_grid - 1) idx = m_grid - 2;
          double frac = pos - idx;
          double val = densR[idx] + frac * diffDensR[idx];
          allLogLiks[i] = std::log((val > minDensR) ? val : minDensR);
        }
      }
    }

    // 2. Categorical calculations: calcCatLogLiks
    if (hasCat) {
      std::fill(catLogLiks.begin(), catLogLiks.end(), 0.0);

      // Pre-calculate weighted lookups
      for (int q = 0; q < qq; ++q) {
        double w = catWeights[q];
        int nlev = numLev[q];
        const double* lp = logProbs[q].data();
        double* w_tab = weightedTabs[q].data();
        for (int lev = 0; lev < nlev; ++lev) {
          for (int cl = 0; cl < numClust; ++cl) {
            w_tab[lev * numClust + cl] = w * lp[cl + lev * numClust];
          }
        }
      }

      for (int q = 0; q < qq; ++q) {
        if (catWeights[q] == 0.0) continue;
        const int* q_col = p_cat + q * nn;
        const double* w_tab = weightedTabs[q].data();

        for (int i = 0; i < nn; ++i) {
          const double* w_row = w_tab + (q_col[i] - 1) * numClust;
          for (int cl = 0; cl < numClust; ++cl) {
            catLogLiks[i + cl * nn] += w_row[cl];
          }
        }
      }

      // Combine into allLogLiks
      if (hasCon) {
        for (int i = 0; i < nn * numClust; ++i) {
          allLogLiks[i] += catLogLiks[i];
        }
      } else {
        std::copy(catLogLiks.begin(), catLogLiks.end(), allLogLiks.begin());
      }
    }

    // 3. Partition data into clusters: membOld <- membNew, membNew <- rowMaxInds(allLogLiks)
    std::copy(membNew.begin(), membNew.end(), membOld.begin());

    for (int i = 0; i < nn; ++i) {
      double maxVal = allLogLiks[i];
      int maxIdx = 1;
      for (int cl = 1; cl < numClust; ++cl) {
        double val = allLogLiks[i + cl * nn];
        if (val > maxVal) {
          maxVal = val;
          maxIdx = cl + 1;
        }
      }
      membNew[i] = maxIdx;
    }

    // 4. Update means: aggregateMeans
    std::fill(countVec.begin(), countVec.end(), 0.0);
    for (int n = 0; n < nn; ++n) {
      countVec[membNew[n] - 1] += 1.0;
    }

    if (hasCon) {
      std::fill(currentMeans.begin(), currentMeans.end(), 0.0);
      for (int p = 0; p < pp; ++p) {
        const double* col_ptr = p_con + p * nn;
        double* out_col = currentMeans.data() + p * numClust;
        for (int n = 0; n < nn; ++n) {
          out_col[membNew[n] - 1] += col_ptr[n];
        }
        for (int k = 0; k < numClust; ++k) {
          if (countVec[k] != 0.0) {
            out_col[k] /= countVec[k];
          }
        }
      }
    }

    // 5. Update categorical probabilities: updateLogProbs
    if (hasCat) {
      for (int q = 0; q < qq; ++q) {
        int nlev = numLev[q];
        const int* q_col = p_cat + q * nn;

        // (a) Tabulate
        std::fill(rawTab.begin(), rawTab.begin() + numClust * nlev, 0);
        for (int i = 0; i < nn; ++i) {
          int row = membNew[i] - 1;
          int col = q_col[i] - 1;
          rawTab[row + col * numClust] += 1;
        }

        // (b) Smooth
        if (catBw != 0.0) {
          std::fill(colSums.begin(), colSums.begin() + nlev, 0);
          for (int j = 0; j < nlev; ++j) {
            for (int i = 0; i < numClust; ++i) {
              colSums[j] += rawTab[i + j * numClust];
            }
          }

          double bw_div_k = (numClust > 1) ? (catBw / (numClust - 1.0)) : 0.0;
          double one_minus_bw = 1.0 - catBw;
          for (int j = 0; j < nlev; ++j) {
            for (int i = 0; i < numClust; ++i) {
              int offCounts1 = colSums[j] - rawTab[i + j * numClust];
              midMat[i + j * numClust] = one_minus_bw * rawTab[i + j * numClust] + bw_div_k * offCounts1;
            }
          }

          std::fill(rowSums.begin(), rowSums.begin() + numClust, 0.0);
          for (int j = 0; j < nlev; ++j) {
            for (int i = 0; i < numClust; ++i) {
              rowSums[i] += midMat[i + j * numClust];
            }
          }

          double bw_div_lev = (nlev > 1) ? (catBw / (nlev - 1.0)) : 0.0;
          for (int j = 0; j < nlev; ++j) {
            for (int i = 0; i < numClust; ++i) {
              double offCounts2 = rowSums[i] - midMat[i + j * numClust];
              outMat[i + j * numClust] = one_minus_bw * midMat[i + j * numClust] + bw_div_lev * offCounts2;
            }
          }
        } else {
          for (int idx = 0; idx < numClust * nlev; ++idx) {
            outMat[idx] = static_cast<double>(rawTab[idx]);
          }
        }

        // (c) Normalize and take log
        std::fill(finalRowSums.begin(), finalRowSums.begin() + numClust, 0.0);
        for (int j = 0; j < nlev; ++j) {
          for (int i = 0; i < numClust; ++i) {
            finalRowSums[i] += outMat[i + j * numClust];
          }
        }

        double* lp_dest = logProbs[q].data();
        for (int j = 0; j < nlev; ++j) {
          for (int i = 0; i < numClust; ++i) {
            double denom = finalRowSums[i];
            lp_dest[i + j * numClust] = (denom > 0.0) ? std::log(outMat[i + j * numClust] / denom) : R_NegInf;
          }
        }
      }
    }

    if (verbose) {
      IntegerVector curOld(nn);
      std::copy(membOld.begin(), membOld.end(), curOld.begin());
      membHistory.push_back(curOld);
    }

    // Check degenerate
    for (int cl = 0; cl < numClust; ++cl) {
      if (countVec[cl] == 0.0) {
        degenerateSoln = true;
        break;
      }
    }
    if (degenerateSoln) break;
  }

  // After loop: compute summary measures
  double totalLogLik = 0.0;
  if (degenerateSoln) {
    totalLogLik = R_NegInf;
  } else {
    for (int i = 0; i < nn; ++i) {
      double maxVal = allLogLiks[i];
      for (int cl = 1; cl < numClust; ++cl) {
        double val = allLogLiks[i + cl * nn];
        if (val > maxVal) maxVal = val;
      }
      totalLogLik += maxVal;
    }
  }

  double catLogLik = 0.0;
  if (hasCat) {
    for (int i = 0; i < nn; ++i) {
      double maxVal = catLogLiks[i];
      for (int cl = 1; cl < numClust; ++cl) {
        double val = catLogLiks[i + cl * nn];
        if (val > maxVal) maxVal = val;
      }
      catLogLik += maxVal;
    }
  }

  double winDist = 0.0;
  if (hasCon) {
    for (int i = 0; i < nn; ++i) {
      int cl = membNew[i] - 1;
      winDist += distMat[i + cl * nn];
    }
  }

  IntegerVector finalMemb(nn);
  std::copy(membNew.begin(), membNew.end(), finalMemb.begin());

  RObject finalMeansObj = R_NilValue;
  if (hasCon) {
    NumericMatrix meansMat(numClust, pp);
    std::copy(currentMeans.begin(), currentMeans.end(), meansMat.begin());
    finalMeansObj = meansMat;
  }

  List finalLogProbsList(qq);
  if (hasCat) {
    for (int q = 0; q < qq; ++q) {
      int nlev = numLev[q];
      NumericMatrix m_q(numClust, nlev);
      std::copy(logProbs[q].begin(), logProbs[q].end(), m_q.begin());
      finalLogProbsList[q] = m_q;
    }
  }

  RObject catLogLiksMatObj = R_NilValue;
  if (verbose && hasCat) {
    NumericMatrix clm(nn, numClust);
    std::copy(catLogLiks.begin(), catLogLiks.end(), clm.begin());
    catLogLiksMatObj = clm;
  }

  List membHistList(membHistory.size());
  if (verbose) {
    for (size_t i = 0; i < membHistory.size(); ++i) {
      membHistList[i] = membHistory[i];
    }
  }

  return List::create(
    Named("numIter") = numIter,
    Named("degenerateSoln") = degenerateSoln,
    Named("finalMemb") = finalMemb,
    Named("finalMeans") = finalMeansObj,
    Named("finalLogProbs") = finalLogProbsList,
    Named("totalLogLik") = totalLogLik,
    Named("catLogLik") = catLogLik,
    Named("winDist") = winDist,
    Named("catLogLiks") = catLogLiksMatObj,
    Named("membHistory") = membHistList
  );
}

