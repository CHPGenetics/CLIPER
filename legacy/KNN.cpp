#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector determineOverlapCpp(Rcpp::IntegerMatrix m, int overlapCut) {
  const int nrow = m.nrow();
  const int ncol = m.ncol();
  const int twiceK = 2 * ncol;
  
  Rcpp::IntegerVector keepFlag(nrow);   
  Rcpp::NumericVector maxOverlapPerRow(nrow);  
  Rcpp::IntegerVector overlaps(nrow); 
  
  for (int i = 1; i < nrow; ++i) {
    if (i % 500 == 0) {
      Rcpp::Rcout << "Completed Computing KNN Overlap " << i << " of " << nrow << std::endl;
    }
    
    for (int j = 0; j < i; ++j) {
      if (keepFlag[j] != 0) {
        overlaps[j] = 0;
        continue;
      }
      Rcpp::IntegerVector row_i = m(i, Rcpp::_);
      Rcpp::IntegerVector row_j = m(j, Rcpp::_);
      overlaps[j] = twiceK - Rcpp::union_(row_i, row_j).size();
    }
    
    int currentMax = 0;
    for (int j = 0; j < i; ++j) {
      if (overlaps[j] > currentMax) currentMax = overlaps[j];
    }
    maxOverlapPerRow[i] = currentMax;
    
    if (currentMax > overlapCut) {
      keepFlag[i] = -1;
    }
  }
  
  return keepFlag;
}
