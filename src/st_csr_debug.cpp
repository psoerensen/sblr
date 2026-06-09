// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_csr_common.h"

// [[Rcpp::export]]
Rcpp::DataFrame debug_stldcsr_xij_pairs(
  std::string ld_prefix,
  std::vector<double> xx,
  Rcpp::IntegerMatrix pairs
) {
 const int m = static_cast<int>(xx.size());

 if (m <= 0) {
  throw std::runtime_error("debug_stldcsr_xij_pairs: xx must be non-empty.");
 }

 if (pairs.ncol() != 2) {
  throw std::runtime_error("debug_stldcsr_xij_pairs: pairs must have two columns.");
 }

 STLDCSR ld = read_and_build_st_ld_csr(
  ld_prefix,
  m,
  xx
 );

 const int npairs = pairs.nrow();

 Rcpp::IntegerVector out_i(npairs);
 Rcpp::IntegerVector out_j(npairs);
 Rcpp::LogicalVector out_has_ij(npairs);
 Rcpp::LogicalVector out_has_ji(npairs);
 Rcpp::NumericVector out_xij_ij(npairs);
 Rcpp::NumericVector out_xij_ji(npairs);

 auto find_xij = [&](int i, int j, double& value) -> bool {
  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int jj = ld.idx[static_cast<std::size_t>(p)];

   if (jj == j) {
    value = static_cast<double>(ld.xij[static_cast<std::size_t>(p)]);
    return true;
   }
  }

  value = NA_REAL;
  return false;
 };

 for (int k = 0; k < npairs; ++k) {
  if (pairs(k, 0) == NA_INTEGER || pairs(k, 1) == NA_INTEGER) {
   throw std::runtime_error("debug_stldcsr_xij_pairs: pairs contains NA.");
  }

  const int i = pairs(k, 0) - 1;
  const int j = pairs(k, 1) - 1;

  if (i < 0 || i >= m || j < 0 || j >= m) {
   throw std::runtime_error("debug_stldcsr_xij_pairs: marker index out of range.");
  }

  double xij_ij = NA_REAL;
  double xij_ji = NA_REAL;

  const bool has_ij = find_xij(i, j, xij_ij);
  const bool has_ji = find_xij(j, i, xij_ji);

  out_i[k] = i + 1;
  out_j[k] = j + 1;
  out_has_ij[k] = has_ij;
  out_has_ji[k] = has_ji;
  out_xij_ij[k] = xij_ij;
  out_xij_ji[k] = xij_ji;
 }

 return Rcpp::DataFrame::create(
  Rcpp::Named("i") = out_i,
  Rcpp::Named("j") = out_j,
  Rcpp::Named("has_ij") = out_has_ij,
  Rcpp::Named("has_ji") = out_has_ji,
  Rcpp::Named("xij_ij") = out_xij_ij,
  Rcpp::Named("xij_ji") = out_xij_ji
 );
}
