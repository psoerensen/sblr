#ifndef SBLR_BLR_BLOCK_EIGEN_H
#define SBLR_BLR_BLOCK_EIGEN_H

#include <armadillo>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace core {

enum class BlockEigenFilterMode {
  hard_truncate = 0,
  ridge_fixed = 1,
  ridge_lw = 2
};

struct BlockEigenFilterSpec {
  BlockEigenFilterMode mode = BlockEigenFilterMode::hard_truncate;
  double tau = 0.01;
  double eta = 0.0;
  double mu_floor = 0.01;
};

struct BlockEigenBlockDiagnostics {
  int start = 0;
  int size = 0;
  int n_kept = 0;
  double mu_min = 0.0;
  double shrink = 0.0;
};

struct BlockEigenBlockStorage {
  int start = 0;
  int size = 0;
  std::vector<float> upper_triangle;

  inline std::size_t packed_index(int i, int j) const {
    const std::size_t ii = static_cast<std::size_t>(i);
    const std::size_t jj = static_cast<std::size_t>(j);
    const std::size_t ss = static_cast<std::size_t>(size);
    return ii * (2 * ss - ii + 1) / 2 + (jj - ii);
  }

  inline double symmetric_at(int a, int b) const {
    const int i = (a <= b) ? a : b;
    const int j = (a <= b) ? b : a;
    return static_cast<double>(upper_triangle[packed_index(i, j)]);
  }

  // Compatibility spellings retained while scalar translation units migrate.
  inline std::size_t pidx(int i, int j) const { return packed_index(i, j); }
  inline double sym_at(int a, int b) const { return symmetric_at(a, b); }
  inline std::vector<float>& tri() = delete;
};

struct BlockEigenView {
  std::size_t marker_count = 0;
  const std::vector<BlockEigenBlockStorage>* blocks = nullptr;
  const int* block_of = nullptr;
  const int* local_of = nullptr;
  const arma::rowvec* diagonal = nullptr;

  inline const arma::rowvec& diag() const { return *diagonal; }

  inline void apply_offdiag(int marker, double difference, arma::rowvec& residual) const {
    const int group = block_of[static_cast<std::size_t>(marker)];
    const int local = local_of[static_cast<std::size_t>(marker)];
    const BlockEigenBlockStorage& block = (*blocks)[static_cast<std::size_t>(group)];
    for (int j = 0; j < block.size; ++j) {
      if (j == local) continue;
      const double value = block.symmetric_at(local, j);
      if (value != 0.0) residual(static_cast<arma::uword>(block.start + j)) -= value * difference;
    }
  }

  inline void apply_offdiag(int marker, double difference, std::vector<double>& residual) const {
    const int group = block_of[static_cast<std::size_t>(marker)];
    const int local = local_of[static_cast<std::size_t>(marker)];
    const BlockEigenBlockStorage& block = (*blocks)[static_cast<std::size_t>(group)];
    for (int j = 0; j < block.size; ++j) {
      if (j == local) continue;
      const double value = block.symmetric_at(local, j);
      if (value != 0.0) residual[static_cast<std::size_t>(block.start + j)] -= value * difference;
    }
  }

  inline void rebuild(const arma::rowvec& wy, const arma::rowvec& effects,
                      arma::rowvec& residual) const {
    residual = wy;
    for (const BlockEigenBlockStorage& block : *blocks) {
      for (int i = 0; i < block.size; ++i) {
        const double effect = effects(static_cast<arma::uword>(block.start + i));
        if (effect == 0.0) continue;
        for (int j = 0; j < block.size; ++j) {
          residual(static_cast<arma::uword>(block.start + j)) -=
            block.symmetric_at(i, j) * effect;
        }
      }
    }
  }

  inline void rebuild(const std::vector<double>& wy, const std::vector<double>& effects,
                      std::vector<double>& residual) const {
    residual = wy;
    for (const BlockEigenBlockStorage& block : *blocks) {
      for (int i = 0; i < block.size; ++i) {
        const double effect = effects[static_cast<std::size_t>(block.start + i)];
        if (effect == 0.0) continue;
        for (int j = 0; j < block.size; ++j) {
          residual[static_cast<std::size_t>(block.start + j)] -=
            block.symmetric_at(i, j) * effect;
        }
      }
    }
  }
};

struct BlockEigenStorage {
  std::size_t marker_count = 0;
  std::vector<BlockEigenBlockStorage> blocks;
  std::vector<int> block_of;
  std::vector<int> local_of;
  arma::rowvec diagonal;

  inline BlockEigenView view() const {
    BlockEigenView result;
    result.marker_count = marker_count;
    result.blocks = &blocks;
    result.block_of = block_of.empty() ? nullptr : block_of.data();
    result.local_of = local_of.empty() ? nullptr : local_of.data();
    result.diagonal = &diagonal;
    return result;
  }

  inline const arma::rowvec& diag() const { return view().diag(); }
  inline void apply_offdiag(int i, double d, arma::rowvec& r) const { view().apply_offdiag(i, d, r); }
  inline void apply_offdiag(int i, double d, std::vector<double>& r) const { view().apply_offdiag(i, d, r); }
  inline void rebuild(const arma::rowvec& w, const arma::rowvec& b, arma::rowvec& r) const { view().rebuild(w, b, r); }
  inline void rebuild(const std::vector<double>& w, const std::vector<double>& b,
                      std::vector<double>& r) const { view().rebuild(w, b, r); }
};

inline std::size_t block_eigen_packed_length(int size) {
  if (size <= 0) throw std::invalid_argument("block-eigen block size must be positive.");
  const std::size_t s = static_cast<std::size_t>(size);
  if (s > (std::numeric_limits<std::size_t>::max() - s) / s) {
    throw std::overflow_error("block-eigen packed length overflow.");
  }
  return s * (s + 1) / 2;
}

inline void validate_block_eigen_view(const BlockEigenView& view) {
  if (view.marker_count == 0 || view.blocks == nullptr || view.blocks->empty())
    throw std::invalid_argument("block-eigen view requires markers and blocks.");
  if (view.block_of == nullptr || view.local_of == nullptr || view.diagonal == nullptr)
    throw std::invalid_argument("block-eigen view contains null mapping or diagonal data.");
  if (view.diagonal->n_elem != view.marker_count)
    throw std::invalid_argument("block-eigen diagonal length does not match marker count.");

  std::size_t expected_start = 0;
  for (std::size_t group = 0; group < view.blocks->size(); ++group) {
    const BlockEigenBlockStorage& block = (*view.blocks)[group];
    if (block.size <= 0 || block.start != static_cast<int>(expected_start))
      throw std::invalid_argument("block-eigen blocks must be positive, contiguous, and start at zero.");
    if (block.upper_triangle.size() != block_eigen_packed_length(block.size))
      throw std::invalid_argument("block-eigen packed triangle has the wrong length.");
    for (float value : block.upper_triangle) {
      if (!std::isfinite(static_cast<double>(value)))
        throw std::invalid_argument("block-eigen packed values must be finite.");
    }
    for (int local = 0; local < block.size; ++local) {
      const std::size_t marker = expected_start + static_cast<std::size_t>(local);
      if (marker >= view.marker_count || view.block_of[marker] != static_cast<int>(group) ||
          view.local_of[marker] != local)
        throw std::invalid_argument("block-eigen mappings do not identify their marker positions.");
      const double packed_diagonal = block.symmetric_at(local, local);
      const double runtime_diagonal = (*view.diagonal)(static_cast<arma::uword>(marker));
      if (!std::isfinite(runtime_diagonal) || runtime_diagonal <= 0.0 ||
          runtime_diagonal != packed_diagonal)
        throw std::invalid_argument("block-eigen runtime diagonal must equal the positive packed diagonal.");
    }
    expected_start += static_cast<std::size_t>(block.size);
  }
  if (expected_start != view.marker_count)
    throw std::invalid_argument("block-eigen blocks must cover the complete marker domain.");
}

inline void validate_block_eigen_storage(const BlockEigenStorage& storage) {
  if (storage.block_of.size() != storage.marker_count || storage.local_of.size() != storage.marker_count)
    throw std::invalid_argument("block-eigen mapping lengths do not match marker count.");
  validate_block_eigen_view(storage.view());
}

}  // namespace core
}  // namespace sblr

#endif
