#include <Rcpp.h>
using namespace Rcpp;

//' @title Computing the rank of the Euclidean distance between observation and correct centroid vs all centroids.
//' @param Xm Matrix of centroids, rows are centroids and columns are dimensions.
//' @param Ym Data matrix, rows are observations and columns are dimensions.
//' @param Zs Assignments of observations to centroids.
//' @export
// [[Rcpp::export]]
NumericVector fastRank1nn(NumericMatrix Xm, NumericMatrix Ym, IntegerVector Zs) {
  int n = Xm.nrow(), 
    k = Xm.ncol(),
    m = Ym.nrow();
  NumericVector counts(m);
  NumericVector rSq(m);
  for (int i1 = 0; i1 < m; i1++) {
    //compute distance between Y[i1, ] and X[i1, ]
    double sum = 0;
    for (int d1 = 0; d1 < k; d1++) {
      sum += pow(Xm(Zs[i1] - 1, d1) - Ym(i1, d1), 2.0);
    }
    rSq[i1] = sum;
  }
  for (int i1 = 0; i1 < m; i1++) {
    for (int j1 = 0; j1 < n; j1++) {
      if (j1 != Zs[i1] - 1) {
        //compute distance between Y[i1, ] and X[j1, ]
        double sum = 0;
        for (int d1 = 0; d1 < k; d1++) {
          sum += pow(Xm(j1, d1) - Ym(i1, d1), 2.0);
        }
        if (sum < rSq[i1]) {
          counts[i1]++;
        }
      }
    }
    counts[i1]++; // to convert to a rank
  }
  return counts; 
}

