#ifndef SBLR_BLR_BLOCK_LOW_RANK_H
#define SBLR_BLR_BLOCK_LOW_RANK_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sblr {
namespace core {

constexpr double block_low_rank_eigenvalue_tolerance = 1e-10;

struct BlockLowRankDiagnostics {
  int start = 0;
  int size = 0;
  int positive_rank = 0;
  int retained_rank = 0;
  double positive_mass = 0.0;
  double retained_mass = 0.0;
  double minimum_retained_eigenvalue = 0.0;
  double maximum_omitted_eigenvalue = 0.0;
  int negative_eigenvalue_count = 0;
  double negative_eigenvalue_mass = 0.0;
};

struct BlockLowRankBlock {
  int start = 0;
  int size = 0;
  int rank = 0;
  arma::uword residual_offset = 0;
  // Column-major rank by marker storage: q(local marker) is contiguous.
  std::vector<float> factor;
  // One transformed score row per input trait.
  arma::mat transformed_score;

  inline double q(int eigen, int local) const {
    return static_cast<double>(factor[
      static_cast<std::size_t>(local) * static_cast<std::size_t>(rank) +
      static_cast<std::size_t>(eigen)
    ]);
  }
};

struct BlockLowRankChainState {
  arma::rowvec residual;
};

struct BlockLowRankStorage {
  std::size_t marker_count = 0;
  int trait_count = 0;
  arma::uword reduced_dimension = 0;
  std::vector<BlockLowRankBlock> blocks;
  std::vector<int> block_of;
  std::vector<int> local_of;
  arma::rowvec diagonal;
  arma::rowvec transformed_score_norm_squared;
  double construction_seconds = 0.0;
  double cross_product_seconds = 0.0;
  double eigendecomposition_seconds = 0.0;
  double transformation_seconds = 0.0;
  double construction_workspace_bytes = 0.0;
  double operator_storage_bytes = 0.0;

  inline const arma::rowvec& diag() const { return diagonal; }
  inline arma::uword residual_size() const { return reduced_dimension; }

  inline double corrected_rhs(int marker, double beta_old,
                              const arma::rowvec& residual) const {
    const int group = block_of.data()[marker];
    const int local = local_of.data()[marker];
    const BlockLowRankBlock& block = blocks.data()[group];
    const float* factor_ptr = block.factor.data() +
      static_cast<std::size_t>(local) * static_cast<std::size_t>(block.rank);
    const double* residual_ptr = residual.memptr() + block.residual_offset;
    double value = diagonal.memptr()[marker] * beta_old;
#ifdef _OPENMP
#pragma omp simd reduction(+:value)
#endif
    for (int k = 0; k < block.rank; ++k) {
      value += static_cast<double>(factor_ptr[k]) * residual_ptr[k];
    }
    return value;
  }

  inline void apply_difference(int marker, double difference,
                               arma::rowvec& residual) const {
    const int group = block_of.data()[marker];
    const int local = local_of.data()[marker];
    const BlockLowRankBlock& block = blocks.data()[group];
    const float* factor_ptr = block.factor.data() +
      static_cast<std::size_t>(local) * static_cast<std::size_t>(block.rank);
    double* residual_ptr = residual.memptr() + block.residual_offset;
#ifdef _OPENMP
#pragma omp simd
#endif
    for (int k = 0; k < block.rank; ++k) {
      residual_ptr[k] -= static_cast<double>(factor_ptr[k]) * difference;
    }
  }

  inline void rebuild(int trait, const arma::rowvec&, const arma::rowvec& effects,
                      arma::rowvec& residual) const {
    if (trait < 0 || trait >= trait_count)
      throw std::out_of_range("low-rank trait index is out of range.");
    residual.set_size(reduced_dimension);
    for (const BlockLowRankBlock& block : blocks) {
      for (int k = 0; k < block.rank; ++k) {
        double value = block.transformed_score(
          static_cast<arma::uword>(trait), static_cast<arma::uword>(k));
        for (int local = 0; local < block.size; ++local) {
          value -= block.q(k, local) * effects(
            static_cast<arma::uword>(block.start + local));
        }
        residual(block.residual_offset + static_cast<arma::uword>(k)) = value;
      }
    }
  }

  inline double projected_score_dot(int, const arma::rowvec& effects,
                                    const arma::rowvec& projected_score) const {
    return arma::dot(effects, projected_score);
  }

  inline double quadratic_form(const arma::rowvec& effects) const {
    // Validation/reference operation: ordinary retained-low-rank MCMC
    // iterations must use fitted_quadratic() instead of traversing Q columns.
    double result = 0.0;
    for (const BlockLowRankBlock& block : blocks) {
      for (int k = 0; k < block.rank; ++k) {
        double fitted = 0.0;
        for (int local = 0; local < block.size; ++local) {
          fitted += block.q(k, local) * effects(
            static_cast<arma::uword>(block.start + local));
        }
        result += fitted * fitted;
      }
    }
    return result;
  }

  inline double fitted_quadratic(int trait, const arma::rowvec&,
                                 const arma::rowvec&,
                                 const arma::rowvec& residual) const {
    if (trait < 0 || trait >= trait_count)
      throw std::out_of_range("low-rank trait index is out of range.");
    double result = 0.0;
    const double* residual_ptr = residual.memptr();
    for (const BlockLowRankBlock& block : blocks) {
      const double* transformed_ptr = block.transformed_score.memptr();
      for (int k = 0; k < block.rank; ++k) {
        const double fitted = transformed_ptr[
          static_cast<std::size_t>(trait) +
          static_cast<std::size_t>(k) * static_cast<std::size_t>(trait_count)
        ] - residual_ptr[block.residual_offset + static_cast<arma::uword>(k)];
        result += fitted * fitted;
      }
    }
    return result;
  }

  inline double residual_sse(int trait, double yy, const arma::rowvec&,
                             const arma::rowvec&, const arma::rowvec& residual) const {
    return yy - transformed_score_norm_squared(static_cast<arma::uword>(trait)) +
      arma::dot(residual, residual);
  }

  inline double genetic_variance(int trait, const arma::rowvec& effects,
                                 const arma::rowvec& score,
                                 const arma::rowvec& residual,
                                 double n) const {
    return fitted_quadratic(trait, effects, score, residual) / n;
  }

  inline double rebuild_and_measure_drift(int trait, const arma::rowvec& score,
                                          const arma::rowvec& effects,
                                          arma::rowvec& residual) const {
    arma::rowvec rebuilt;
    rebuild(trait, score, effects, rebuilt);
    if (rebuilt.n_elem != residual.n_elem)
      throw std::invalid_argument("low-rank residual has an inconsistent size.");
    double max_abs_drift = 0.0;
    const double* rebuilt_ptr = rebuilt.memptr();
    const double* residual_ptr = residual.memptr();
    for (arma::uword k = 0; k < residual.n_elem; ++k) {
      max_abs_drift = std::max(max_abs_drift,
        std::abs(residual_ptr[k] - rebuilt_ptr[k]));
    }
    residual = std::move(rebuilt);
    return max_abs_drift;
  }

  inline void materialize_residual(int, const arma::rowvec& residual,
                                   arma::rowvec& marker_residual) const {
    marker_residual.zeros(static_cast<arma::uword>(marker_count));
    for (const BlockLowRankBlock& block : blocks) {
      for (int local = 0; local < block.size; ++local) {
        double value = 0.0;
        for (int k = 0; k < block.rank; ++k) {
          value += block.q(k, local) *
            residual(block.residual_offset + static_cast<arma::uword>(k));
        }
        marker_residual(static_cast<arma::uword>(block.start + local)) = value;
      }
    }
  }
};

inline void validate_block_low_rank_storage(const BlockLowRankStorage& storage) {
  if (storage.marker_count == 0 || storage.blocks.empty())
    throw std::invalid_argument("low-rank storage requires markers and blocks.");
  if (storage.block_of.size() != storage.marker_count ||
      storage.local_of.size() != storage.marker_count ||
      storage.diagonal.n_elem != storage.marker_count)
    throw std::invalid_argument("low-rank marker mappings have inconsistent sizes.");
  std::size_t expected_start = 0;
  arma::uword expected_offset = 0;
  for (std::size_t group = 0; group < storage.blocks.size(); ++group) {
    const BlockLowRankBlock& block = storage.blocks[group];
    if (block.start != static_cast<int>(expected_start) || block.size <= 0 ||
        block.rank <= 0 || block.rank > block.size ||
        block.residual_offset != expected_offset)
      throw std::invalid_argument("low-rank blocks must be positive, contiguous, and shape-consistent.");
    if (block.factor.size() != static_cast<std::size_t>(block.rank) *
                               static_cast<std::size_t>(block.size) ||
        block.transformed_score.n_rows != static_cast<arma::uword>(storage.trait_count) ||
        block.transformed_score.n_cols != static_cast<arma::uword>(block.rank))
      throw std::invalid_argument("low-rank factor or transformed-score shape is inconsistent.");
    for (int local = 0; local < block.size; ++local) {
      const std::size_t marker = expected_start + static_cast<std::size_t>(local);
      if (storage.block_of[marker] != static_cast<int>(group) ||
          storage.local_of[marker] != local ||
          !std::isfinite(storage.diagonal(static_cast<arma::uword>(marker))) ||
          storage.diagonal(static_cast<arma::uword>(marker)) <= 0.0)
        throw std::invalid_argument("low-rank marker mapping or runtime diagonal is invalid.");
    }
    expected_start += static_cast<std::size_t>(block.size);
    expected_offset += static_cast<arma::uword>(block.rank);
  }
  if (expected_start != storage.marker_count || expected_offset != storage.reduced_dimension)
    throw std::invalid_argument("low-rank blocks do not cover the operator domain.");
}

// Development-only synthetic operator used by the native hot-path benchmark.
// It bypasses factor construction deliberately so benchmark timings contain
// only sampling-state operations.
inline BlockLowRankStorage make_block_low_rank_hot_path_fixture(int markers,
                                                                int rank) {
  if (markers <= 0 || rank <= 0 || rank > markers)
    throw std::invalid_argument("benchmark markers/rank are inconsistent.");
  BlockLowRankStorage storage;
  storage.marker_count = static_cast<std::size_t>(markers);
  storage.trait_count = 1;
  storage.reduced_dimension = static_cast<arma::uword>(rank);
  storage.block_of.assign(static_cast<std::size_t>(markers), 0);
  storage.local_of.resize(static_cast<std::size_t>(markers));
  storage.diagonal.zeros(static_cast<arma::uword>(markers));
  storage.transformed_score_norm_squared.zeros(1);
  BlockLowRankBlock block;
  block.start = 0;
  block.size = markers;
  block.rank = rank;
  block.residual_offset = 0;
  block.factor.resize(static_cast<std::size_t>(markers) *
                      static_cast<std::size_t>(rank));
  block.transformed_score.set_size(1, static_cast<arma::uword>(rank));
  const double inv_sqrt_rank = 1.0 / std::sqrt(static_cast<double>(rank));
  for (int local = 0; local < markers; ++local) {
    double diagonal = 0.0;
    float* factor_ptr = block.factor.data() +
      static_cast<std::size_t>(local) * static_cast<std::size_t>(rank);
    for (int k = 0; k < rank; ++k) {
      const double value =
        (std::sin(0.013 * static_cast<double>((local + 1) * (k + 1))) +
         0.25 * std::cos(0.007 * static_cast<double>(local + k + 2))) *
        inv_sqrt_rank;
      factor_ptr[k] = static_cast<float>(value);
      diagonal += static_cast<double>(factor_ptr[k]) *
        static_cast<double>(factor_ptr[k]);
    }
    storage.local_of[static_cast<std::size_t>(local)] = local;
    storage.diagonal(static_cast<arma::uword>(local)) = diagonal;
  }
  for (int k = 0; k < rank; ++k) {
    const double value = std::cos(0.017 * static_cast<double>(k + 1));
    block.transformed_score(0, static_cast<arma::uword>(k)) = value;
    storage.transformed_score_norm_squared(0) += value * value;
  }
  storage.blocks.push_back(std::move(block));
  validate_block_low_rank_storage(storage);
  return storage;
}

}  // namespace core
}  // namespace sblr

#endif
