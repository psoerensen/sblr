#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <RcppArmadillo.h>

#include "st_bed_decode.h"
#include "st_block_eigen.h"

BlockLowRankOperator build_block_low_rank(
    const PackedBedMatrix& G,
    const std::vector<double>& af,
    const std::vector<int>& block_start,
    double eigen_prop,
    arma::mat& wy_mat,
    int nthreads,
    std::vector<BlockLowRankDiag>* diag_out) {
  const auto construction_start = std::chrono::steady_clock::now();
  const int n = G.n;
  const int m = G.m;
  const int nt = static_cast<int>(wy_mat.n_rows);
  if (!(eigen_prop > 0.0 && eigen_prop < 1.0) || !std::isfinite(eigen_prop))
    throw std::runtime_error("build_block_low_rank: eigen_prop must be finite and in (0, 1).");
  if (m <= 0 || n <= 0 || block_start.empty() || block_start.front() != 0)
    throw std::runtime_error("build_block_low_rank: invalid dimensions or block coverage.");
  if (static_cast<int>(af.size()) != m || static_cast<int>(wy_mat.n_cols) != m || nt <= 0)
    throw std::runtime_error("build_block_low_rank: score or allele-frequency dimensions are invalid.");
  for (int i = 0; i < m; ++i) {
    if (!std::isfinite(af[static_cast<std::size_t>(i)]) ||
        af[static_cast<std::size_t>(i)] <= 0.0 || af[static_cast<std::size_t>(i)] >= 1.0)
      throw std::runtime_error("build_block_low_rank: allele frequencies must be finite and in (0, 1).");
  }
  for (std::size_t block = 1; block < block_start.size(); ++block) {
    if (block_start[block] <= block_start[block - 1] || block_start[block] >= m)
      throw std::runtime_error("build_block_low_rank: block starts must be strictly increasing and cover markers.");
  }

  BlockLowRankOperator op;
  op.marker_count = static_cast<std::size_t>(m);
  op.trait_count = nt;
  op.blocks.reserve(block_start.size());
  op.block_of.assign(static_cast<std::size_t>(m), -1);
  op.local_of.assign(static_cast<std::size_t>(m), -1);
  op.diagonal.zeros(static_cast<arma::uword>(m));
  op.transformed_score_norm_squared.zeros(static_cast<arma::uword>(nt));
  if (diag_out != nullptr) diag_out->clear();

  std::vector<float> decoded;
  for (std::size_t group = 0; group < block_start.size(); ++group) {
    const int start = block_start[group];
    const int end = group + 1 < block_start.size() ? block_start[group + 1] : m;
    const int size = end - start;
    if (size <= 0) throw std::runtime_error("build_block_low_rank: empty block.");

    decoded.assign(static_cast<std::size_t>(size) * static_cast<std::size_t>(n), 0.0f);
    decode_packed_block_float(G, start, size, af.data(), decoded.data(), nthreads);
    op.construction_workspace_bytes = std::max(
      op.construction_workspace_bytes,
      static_cast<double>(sizeof(float)) * size * n +
        static_cast<double>(sizeof(double)) * (3.0 * size * size + 5.0 * size)
    );
    const auto cross_start = std::chrono::steady_clock::now();
    arma::mat cross(static_cast<arma::uword>(size), static_cast<arma::uword>(size),
                    arma::fill::zeros);
    for (int a = 0; a < size; ++a) {
      const float* za = decoded.data() + static_cast<std::size_t>(a) * n;
      for (int b = a; b < size; ++b) {
        const float* zb = decoded.data() + static_cast<std::size_t>(b) * n;
        double value = 0.0;
        for (int row = 0; row < n; ++row)
          value += static_cast<double>(za[row]) * static_cast<double>(zb[row]);
        cross(static_cast<arma::uword>(a), static_cast<arma::uword>(b)) = value;
        cross(static_cast<arma::uword>(b), static_cast<arma::uword>(a)) = value;
      }
    }
    op.cross_product_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - cross_start).count();

    const arma::vec diagonal = cross.diag();
    arma::vec sqrt_diagonal(static_cast<arma::uword>(size));
    arma::vec inverse_sqrt_diagonal(static_cast<arma::uword>(size));
    for (int local = 0; local < size; ++local) {
      const double value = diagonal(static_cast<arma::uword>(local));
      if (!std::isfinite(value) || value <= 0.0)
        throw std::runtime_error("build_block_low_rank: non-positive standardized diagonal.");
      sqrt_diagonal(static_cast<arma::uword>(local)) = std::sqrt(value);
      inverse_sqrt_diagonal(static_cast<arma::uword>(local)) = 1.0 / std::sqrt(value);
    }
    arma::mat correlation = cross;
    correlation.each_col() %= inverse_sqrt_diagonal;
    correlation.each_row() %= inverse_sqrt_diagonal.t();
    correlation = 0.5 * (correlation + correlation.t());
    if (!correlation.is_finite())
      throw std::runtime_error("build_block_low_rank: non-finite correlation matrix.");

    const auto eigen_start = std::chrono::steady_clock::now();
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    if (!arma::eig_sym(eigenvalues, eigenvectors, correlation))
      throw std::runtime_error("build_block_low_rank: eig_sym failed.");
    op.eigendecomposition_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - eigen_start).count();

    BlockLowRankDiag diagnostic;
    diagnostic.start = start;
    diagnostic.size = size;
    for (arma::uword idx = 0; idx < eigenvalues.n_elem; ++idx) {
      const double value = eigenvalues(idx);
      if (value > sblr::core::block_low_rank_eigenvalue_tolerance) {
        ++diagnostic.positive_rank;
        diagnostic.positive_mass += value;
      } else if (value < -sblr::core::block_low_rank_eigenvalue_tolerance) {
        ++diagnostic.negative_eigenvalue_count;
        diagnostic.negative_eigenvalue_mass += -value;
      }
    }
    if (diagnostic.positive_rank == 0 || !(diagnostic.positive_mass > 0.0))
      throw std::runtime_error("build_block_low_rank: block has no positive eigenvalues.");

    std::vector<arma::uword> retained;
    retained.reserve(static_cast<std::size_t>(diagnostic.positive_rank));
    for (arma::sword idx = static_cast<arma::sword>(eigenvalues.n_elem) - 1; idx >= 0; --idx) {
      const double value = eigenvalues(static_cast<arma::uword>(idx));
      if (value <= sblr::core::block_low_rank_eigenvalue_tolerance) continue;
      retained.push_back(static_cast<arma::uword>(idx));
      diagnostic.retained_mass += value;
      if (diagnostic.retained_mass / diagnostic.positive_mass > eigen_prop) break;
    }
    diagnostic.retained_rank = static_cast<int>(retained.size());
    diagnostic.minimum_retained_eigenvalue =
      eigenvalues(retained.back());
    if (diagnostic.retained_rank < diagnostic.positive_rank) {
      const arma::sword omitted = static_cast<arma::sword>(retained.back()) - 1;
      if (omitted >= 0)
        diagnostic.maximum_omitted_eigenvalue = eigenvalues(static_cast<arma::uword>(omitted));
    }

    const auto transformation_start = std::chrono::steady_clock::now();
    sblr::core::BlockLowRankBlock block;
    block.start = start;
    block.size = size;
    block.rank = diagnostic.retained_rank;
    block.residual_offset = op.reduced_dimension;
    block.factor.resize(static_cast<std::size_t>(size) * retained.size());
    block.transformed_score.zeros(static_cast<arma::uword>(nt), retained.size());
    for (int k = 0; k < block.rank; ++k) {
      const arma::uword eig = retained[static_cast<std::size_t>(k)];
      const double root = std::sqrt(eigenvalues(eig));
      for (int local = 0; local < size; ++local) {
        block.factor[static_cast<std::size_t>(local) * retained.size() +
                     static_cast<std::size_t>(k)] = static_cast<float>(
          root * eigenvectors(static_cast<arma::uword>(local), eig) *
          sqrt_diagonal(static_cast<arma::uword>(local)));
      }
      for (int trait = 0; trait < nt; ++trait) {
        double transformed = 0.0;
        for (int local = 0; local < size; ++local) {
          transformed += eigenvectors(static_cast<arma::uword>(local), eig) *
            inverse_sqrt_diagonal(static_cast<arma::uword>(local)) *
            wy_mat(static_cast<arma::uword>(trait), static_cast<arma::uword>(start + local));
        }
        block.transformed_score(static_cast<arma::uword>(trait), static_cast<arma::uword>(k)) =
          transformed / root;
        op.transformed_score_norm_squared(static_cast<arma::uword>(trait)) +=
          (transformed / root) * (transformed / root);
      }
    }

    for (int local = 0; local < size; ++local) {
      const int marker = start + local;
      op.block_of[static_cast<std::size_t>(marker)] = static_cast<int>(group);
      op.local_of[static_cast<std::size_t>(marker)] = local;
      double runtime_diagonal = 0.0;
      for (int k = 0; k < block.rank; ++k) {
        const double value = block.q(k, local);
        runtime_diagonal += value * value;
      }
      if (!std::isfinite(runtime_diagonal) || runtime_diagonal <= 0.0)
        throw std::runtime_error("build_block_low_rank: non-positive runtime diagonal.");
      op.diagonal(static_cast<arma::uword>(marker)) = runtime_diagonal;
      for (int trait = 0; trait < nt; ++trait) {
        double projected = 0.0;
        for (int k = 0; k < block.rank; ++k)
          projected += block.q(k, local) * block.transformed_score(
            static_cast<arma::uword>(trait), static_cast<arma::uword>(k));
        wy_mat(static_cast<arma::uword>(trait), static_cast<arma::uword>(marker)) = projected;
      }
    }
    op.reduced_dimension += static_cast<arma::uword>(block.rank);
    op.transformation_seconds += std::chrono::duration<double>(
      std::chrono::steady_clock::now() - transformation_start).count();
    op.blocks.push_back(std::move(block));
    if (diag_out != nullptr) diag_out->push_back(diagnostic);
  }
  op.operator_storage_bytes =
    static_cast<double>(sizeof(int)) * (op.block_of.size() + op.local_of.size()) +
    static_cast<double>(sizeof(double)) *
      (op.diagonal.n_elem + op.transformed_score_norm_squared.n_elem);
  for (const sblr::core::BlockLowRankBlock& block : op.blocks) {
    op.operator_storage_bytes += static_cast<double>(sizeof(float)) * block.factor.size() +
      static_cast<double>(sizeof(double)) * block.transformed_score.n_elem;
  }
  op.construction_seconds = std::chrono::duration<double>(
    std::chrono::steady_clock::now() - construction_start).count();
  sblr::core::validate_block_low_rank_storage(op);
  return op;
}
