#include <RcppArmadillo.h>
#include <queue>
#include <unordered_set>

using namespace Rcpp ;

typedef std::pair<double, int> e;

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
List dijkstra_cpp(const arma::sp_mat& m, const int k) {

  const int n_rows = m.n_rows;
  arma::uvec seq_idx = arma::regspace<arma::uvec>(0, 1, n_rows - 1);
  arma::umat knn_idx = arma::repmat(seq_idx, 1, k);
  arma::dmat knn_dist(n_rows, k, arma::fill::zeros);

  #ifdef _OPENMP
  #pragma omp parallel num_threads(omp_get_max_threads())
  #pragma omp for
  #endif
  for (int root = 0; root < n_rows; root++) {
    arma::sp_mat::const_iterator it;
    std::priority_queue<e, std::vector<e>, std::greater<e> > pq;
    std::unordered_set<int> reached;
    reached.insert(root);
    arma::sp_mat col(m.col(root));
    for (it = col.begin(); it != col.end(); ++it){
      pq.push(std::make_pair((*it), it.row()));
    }

    int i = 1;
    while(i < k) {
      e current = pq.top();
      int source = current.second;
      pq.pop();
      auto found = reached.find(source);
      if (found == reached.end()) {
        knn_idx(root,  i) = source;
        knn_dist(root, i) = current.first;
        i++;
        reached.insert(source);

        col = m.col(source);
        for (it = col.begin(); it != col.end(); ++it){
          pq.push(std::make_pair((*it) + current.first, it.row()));
        }
      } else if (pq.size() == 0) {
        break;
      }
    }
  }

  return(List::create(Named("knn.idx") = knn_idx + 1,
                      Named("knn.dist") = knn_dist));
}
