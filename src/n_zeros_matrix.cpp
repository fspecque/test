#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
int64_t n_zeros_mat(const arma::mat &mat) {
  int64_t count = arma::accu(mat == 0);
  return count;
}

