#ifndef ST_BLOCK_EIGEN_H
#define ST_BLOCK_EIGEN_H

#include <vector>

#include <RcppArmadillo.h>

#include "packed_bed.h"
#include "st_ld_operator.h"

enum class EigenFilterMode {
  hard_truncate = 0,
  ridge_fixed = 1,
  ridge_lw = 2
};

struct BlockEigenDiag {
  int start = 0;
  int size = 0;
  int n_kept = 0;
  double mu_min = 0.0;
  double shrink = 0.0;
};

BlockEigenOperator build_block_eigen(
    const PackedBedMatrix& G,
    const std::vector<double>& af,
    const std::vector<int>& block_start,
    EigenFilterMode mode,
    double tau,
    double eta,
    arma::mat& wy_mat,
    int nthreads,
    std::vector<BlockEigenDiag>* diag_out = nullptr);

#endif
