#ifndef SBLR_ST_BLOCK_EIGEN_EXECUTION_H
#define SBLR_ST_BLOCK_EIGEN_EXECUTION_H

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "st_block_eigen.h"

// Resolved, binding-neutral retained-block construction input. Rcpp adapters
// remain responsible for validating and converting R objects exactly once.
struct BlockEigenExecutionInput {
  std::vector<std::string> bed_files;
  int n_bed = 0;
  std::vector<std::vector<int>> cls;
  std::vector<int> rows0;
  std::vector<double> af;
  std::vector<int> block_start;
  EigenFilterMode filter_mode = EigenFilterMode::hard_truncate;
  double eigen_tau = 0.01;
  double eigen_eta = 0.0;
  bool low_rank = false;
  double eigen_prop = 0.995;
  int ncores = 1;
};

struct PreparedBlockEigenOperator {
  BlockEigenDispatchOperator op;
  std::vector<BlockEigenDiag> dense_diagnostics;
  std::vector<BlockLowRankDiag> low_rank_diagnostics;
};

inline PreparedBlockEigenOperator prepare_block_eigen_operator(
  const BlockEigenExecutionInput& input,
  int marker_count,
  arma::mat& wy
) {
  PackedBedMatrix packed = read_bedfiles_to_packed_matrix(
    input.bed_files,
    input.n_bed,
    input.rows0.empty() ? nullptr : input.rows0.data(),
    static_cast<int>(input.rows0.size()),
    input.cls
  );
  if (packed.m != marker_count) {
    throw std::runtime_error("BED marker count does not match m.");
  }

  PreparedBlockEigenOperator prepared;
  prepared.op.low_rank = input.low_rank;
  if (input.low_rank) {
    prepared.op.retained = build_block_low_rank(
      packed, input.af, input.block_start, input.eigen_prop, wy,
      input.ncores, &prepared.low_rank_diagnostics
    );
  } else {
    prepared.op.dense = build_block_eigen(
      packed, input.af, input.block_start, input.filter_mode,
      input.eigen_tau, input.eigen_eta, wy, input.ncores,
      &prepared.dense_diagnostics
    );
  }
  return prepared;
}

#endif
