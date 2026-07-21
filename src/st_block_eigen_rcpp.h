#ifndef SBLR_ST_BLOCK_EIGEN_RCPP_H
#define SBLR_ST_BLOCK_EIGEN_RCPP_H

#include <Rcpp.h>

#include <stdexcept>
#include <string>
#include <vector>

#include "st_block_eigen.h"

inline EigenFilterMode parse_block_eigen_filter_mode(const std::string& mode) {
  if (mode == "hard_truncate") return EigenFilterMode::hard_truncate;
  if (mode == "ridge_fixed") return EigenFilterMode::ridge_fixed;
  if (mode == "ridge_lw") return EigenFilterMode::ridge_lw;
  throw std::runtime_error(
    "eigen_filter must be one of 'hard_truncate', 'ridge_fixed', or 'ridge_lw'."
  );
}

inline Rcpp::DataFrame block_eigen_diagnostics_to_data_frame(
    const std::vector<BlockEigenDiag>& diagnostics) {
  const int count = static_cast<int>(diagnostics.size());
  Rcpp::IntegerVector start(count), size(count), n_kept(count);
  Rcpp::NumericVector mu_min(count), shrink(count);
  for (int i = 0; i < count; ++i) {
    const BlockEigenDiag& item = diagnostics[static_cast<std::size_t>(i)];
    start[i] = item.start;
    size[i] = item.size;
    n_kept[i] = item.n_kept;
    mu_min[i] = item.mu_min;
    shrink[i] = item.shrink;
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("start") = start,
    Rcpp::Named("size") = size,
    Rcpp::Named("n_kept") = n_kept,
    Rcpp::Named("mu_min") = mu_min,
    Rcpp::Named("shrink") = shrink
  );
}

#endif
