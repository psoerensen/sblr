#ifndef SBLR_CORE_BLR_SPARSE_LD_CSR_H
#define SBLR_CORE_BLR_SPARSE_LD_CSR_H

#include <armadillo>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace core {

// Canonical owning storage for one immutable sparse-LD operator. Legacy
// payload names preserve compatibility with existing scalar consumers.
struct SparseLdCsrStorage {
  std::vector<std::uint64_t> ptr;
  std::vector<int> idx;
  std::vector<float> xij;

  std::uint64_t input_nnz = 0;
  std::uint64_t symmetric_nnz = 0;
  double max_abs_r = 0.0;
  double max_abs_xij = 0.0;
};

// Immutable borrowed view of one marker-order domain. The off-diagonal CSR
// payload and separate diagonal may be shared independently by callers.
struct SparseLdCsrView {
  std::size_t marker_count = 0;
  const std::uint64_t* row_ptr = nullptr;
  std::size_t row_ptr_size = 0;
  const int* column_index = nullptr;
  const float* offdiag_xij = nullptr;
  std::size_t nonzero_count = 0;
  const arma::rowvec* diagonal = nullptr;

  inline const arma::rowvec& diag() const { return *diagonal; }

  inline void apply_offdiag(int marker, double difference,
                            arma::rowvec& residual) const {
    const std::uint64_t start = row_ptr[static_cast<std::size_t>(marker)];
    const std::uint64_t end = row_ptr[static_cast<std::size_t>(marker + 1)];
    for (std::uint64_t position = start; position < end; ++position) {
      const int neighbour = column_index[static_cast<std::size_t>(position)];
      residual(static_cast<arma::uword>(neighbour)) -=
        static_cast<double>(offdiag_xij[static_cast<std::size_t>(position)]) *
        difference;
    }
  }

  inline void rebuild(const arma::rowvec& trait_wy,
                      const arma::rowvec& effects,
                      arma::rowvec& residual) const {
    residual = trait_wy;
    for (std::size_t marker = 0; marker < marker_count; ++marker) {
      const arma::uword marker_u = static_cast<arma::uword>(marker);
      const double effect = effects(marker_u);
      if (effect == 0.0) continue;
      residual(marker_u) -= (*diagonal)(marker_u) * effect;
      apply_offdiag(static_cast<int>(marker), effect, residual);
    }
  }
};

inline void validate_sparse_ld_csr_view(const SparseLdCsrView& view) {
  if (view.marker_count == 0) {
    throw std::invalid_argument("sparse-LD CSR marker_count must be positive");
  }
  if (view.row_ptr == nullptr ||
      view.row_ptr_size != view.marker_count + 1) {
    throw std::invalid_argument("sparse-LD CSR row_ptr is incomplete");
  }
  if (view.row_ptr[0] != 0) {
    throw std::invalid_argument("sparse-LD CSR row_ptr must start at zero");
  }
  for (std::size_t row = 0; row < view.marker_count; ++row) {
    if (view.row_ptr[row] > view.row_ptr[row + 1]) {
      throw std::invalid_argument("sparse-LD CSR row_ptr must be nondecreasing");
    }
  }
  if (view.row_ptr[view.marker_count] != view.nonzero_count) {
    throw std::invalid_argument("sparse-LD CSR row_ptr does not match nnz");
  }
  if (view.nonzero_count > 0 &&
      (view.column_index == nullptr || view.offdiag_xij == nullptr)) {
    throw std::invalid_argument("sparse-LD CSR off-diagonal payload is incomplete");
  }
  for (std::size_t row = 0; row < view.marker_count; ++row) {
    for (std::uint64_t position = view.row_ptr[row];
         position < view.row_ptr[row + 1]; ++position) {
      const std::size_t p = static_cast<std::size_t>(position);
      const int column = view.column_index[p];
      if (column < 0 || static_cast<std::size_t>(column) >= view.marker_count) {
        throw std::invalid_argument("sparse-LD CSR column index is out of range");
      }
      if (static_cast<std::size_t>(column) == row) {
        throw std::invalid_argument("sparse-LD CSR diagonal entries are not allowed");
      }
      if (!std::isfinite(static_cast<double>(view.offdiag_xij[p]))) {
        throw std::invalid_argument("sparse-LD CSR values must be finite");
      }
    }
  }
  if (view.diagonal == nullptr ||
      view.diagonal->n_elem != view.marker_count) {
    throw std::invalid_argument("sparse-LD CSR diagonal is incomplete");
  }
  for (std::size_t marker = 0; marker < view.marker_count; ++marker) {
    const double value = (*view.diagonal)(static_cast<arma::uword>(marker));
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument("sparse-LD CSR diagonal must be finite and positive");
    }
  }
}

}  // namespace core
}  // namespace sblr

#endif
