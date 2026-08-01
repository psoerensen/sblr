#ifndef ST_BLOCK_EIGEN_H
#define ST_BLOCK_EIGEN_H

#include <vector>

#include <armadillo>

#include "packed_bed.h"
#include "st_ld_operator.h"

using EigenFilterMode = sblr::core::BlockEigenFilterMode;
using BlockEigenDiag = sblr::core::BlockEigenBlockDiagnostics;
using BlockLowRankDiag = sblr::core::BlockLowRankDiagnostics;

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

BlockLowRankOperator build_block_low_rank(
    const PackedBedMatrix& G,
    const std::vector<double>& af,
    const std::vector<int>& block_start,
    double eigen_prop,
    arma::mat& wy_mat,
    int nthreads,
    std::vector<BlockLowRankDiag>* diag_out = nullptr);

#endif
