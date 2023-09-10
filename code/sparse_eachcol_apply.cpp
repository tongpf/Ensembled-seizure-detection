#include <RcppArmadillo.h>
// #include <RcppEigen.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;
// using namespace Eigen;

//[[Rcpp::export]]
NumericVector sparse_eachcol_apply(arma::sp_mat x, NumericVector y){
  const int n = x.n_cols;
  NumericVector f(n);
  for(int i=0;i<n;++i){
    sp_mat::const_iterator it     = x.begin_col(i);
    sp_mat::const_iterator it_end = x.end_col(i);
    
    for(; it != it_end; ++it)
    {
      //cout << "val: " << (*it)    << endl;
      //cout << "row: " << it.row() << endl;
      //cout << "col: " << it.col() << endl;
      f[i]=f[i]+(*it)*y[it.row()];
    }
  }
  return f;
}