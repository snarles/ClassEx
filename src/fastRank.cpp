#include <Rcpp.h>
using namespace Rcpp;

//' @title Computing the rank of the correct label given classification scores.
//' @param pmat Matrix of margins m_y(x) for test instances, rows are observations and columns are labels.
//' @param Zs True labels.
//' @export
// [[Rcpp::export]]
IntegerVector fastRank(NumericMatrix pmat, IntegerVector Zs) {
  int m = pmat.nrow(), 
    k = pmat.ncol();
  IntegerVector counts(m);
  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = 0; j1 < k; j1++) {
      if (pmat(i1, j1) >= pmat(i1, Zs[i1] - 1)) {
        counts[i1]++;
      }
    }
  }
  return counts; 
}

