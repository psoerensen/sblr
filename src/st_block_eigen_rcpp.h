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

inline Rcpp::DataFrame block_low_rank_diagnostics_to_data_frame(
    const std::vector<BlockLowRankDiag>& diagnostics) {
  const int count = static_cast<int>(diagnostics.size());
  Rcpp::IntegerVector start(count), size(count), positive_rank(count), retained_rank(count),
    negative_count(count);
  Rcpp::NumericVector retained_rank_fraction(count), positive_mass(count), retained_mass(count),
    retained_mass_fraction(count), minimum_retained(count), maximum_omitted(count), negative_mass(count),
    tolerance(count, sblr::core::block_low_rank_eigenvalue_tolerance);
  for (int i = 0; i < count; ++i) {
    const BlockLowRankDiag& item = diagnostics[static_cast<std::size_t>(i)];
    start[i] = item.start; size[i] = item.size;
    positive_rank[i] = item.positive_rank; retained_rank[i] = item.retained_rank;
    retained_rank_fraction[i] = static_cast<double>(item.retained_rank) / item.size;
    positive_mass[i] = item.positive_mass; retained_mass[i] = item.retained_mass;
    retained_mass_fraction[i] = item.retained_mass / item.positive_mass;
    minimum_retained[i] = item.minimum_retained_eigenvalue;
    maximum_omitted[i] = item.maximum_omitted_eigenvalue;
    negative_count[i] = item.negative_eigenvalue_count;
    negative_mass[i] = item.negative_eigenvalue_mass;
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("block_start") = start, Rcpp::Named("block_size") = size,
    Rcpp::Named("positive_rank") = positive_rank, Rcpp::Named("retained_rank") = retained_rank,
    Rcpp::Named("retained_rank_fraction") = retained_rank_fraction,
    Rcpp::Named("positive_eigenvalue_mass") = positive_mass,
    Rcpp::Named("retained_eigenvalue_mass") = retained_mass,
    Rcpp::Named("retained_mass_fraction") = retained_mass_fraction,
    Rcpp::Named("minimum_retained_eigenvalue") = minimum_retained,
    Rcpp::Named("maximum_omitted_eigenvalue") = maximum_omitted,
    Rcpp::Named("negative_eigenvalue_count") = negative_count,
    Rcpp::Named("negative_eigenvalue_mass") = negative_mass,
    Rcpp::Named("eigenvalue_tolerance") = tolerance
  );
}

inline Rcpp::List block_low_rank_build_metadata(const BlockLowRankOperator& op) {
  return Rcpp::List::create(
    Rcpp::Named("block_count") = static_cast<int>(op.blocks.size()),
    Rcpp::Named("operator_storage_bytes") = op.operator_storage_bytes,
    Rcpp::Named("chain_residual_storage_bytes") =
      static_cast<double>(sizeof(double)) * op.reduced_dimension,
    Rcpp::Named("construction_workspace_bytes") = op.construction_workspace_bytes,
    Rcpp::Named("construction_time") = op.construction_seconds,
    Rcpp::Named("cross_product_time") = op.cross_product_seconds,
    Rcpp::Named("eigendecomposition_time") = op.eigendecomposition_seconds,
    Rcpp::Named("transformation_time") = op.transformation_seconds
  );
}

#endif
