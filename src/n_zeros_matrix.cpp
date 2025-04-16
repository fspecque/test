#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
int64_t n_zeros_dense_mat(const arma::mat &mat) {
  int64_t count = arma::accu(mat == 0);
  return count;
}

// [[Rcpp::export]]
int64_t n_zeros_sparse_mat(const arma::sp_mat &mat) {
  int64_t count = mat.n_elem - mat.n_nonzero;
  return count;
}

