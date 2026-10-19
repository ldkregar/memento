// [[Rcpp::depends(RcppEigen)]]

// OpenMP is optional: when the toolchain supports it, src/Makevars passes
// $(SHLIB_OPENMP_CXXFLAGS) and Eigen parallelises the product. Without it,
// Eigen::setNbThreads() is a no-op and the product runs single-threaded.
#ifdef _OPENMP
#include <omp.h>
#endif
#include <RcppEigen.h>

// [[Rcpp::export]]
SEXP eigenMapMatMult2(const Eigen::Map<Eigen::MatrixXd> A,
                      Eigen::Map<Eigen::MatrixXd> B,
                      int n_cores = 8){

  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A * B;
  return Rcpp::wrap(C);
}
