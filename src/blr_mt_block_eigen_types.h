#ifndef SBLR_BLR_MT_BLOCK_EIGEN_TYPES_H
#define SBLR_BLR_MT_BLOCK_EIGEN_TYPES_H

#include "blr_block_eigen.h"
#include "blr_mt_default_types.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace mt {

struct MtBlockEigenBundleView {
  std::size_t marker_count = 0;
  std::vector<sblr::core::BlockEigenView> trait_operator;
};

struct MtBlockEigenDataView {
  const std::vector<std::vector<double>>& wy;
  const std::vector<double>& yy;
  const std::vector<int>& n;
  MtBlockEigenBundleView operators;
};

inline void validate_mt_block_eigen_data(
    const MtBlockEigenDataView& data,
    const MtDefaultModelSpec& model) {
  const std::size_t traits = data.wy.size();
  const std::size_t markers = data.operators.marker_count;
  if (traits == 0 || markers == 0 || data.yy.size() != traits ||
      data.n.size() != traits || data.operators.trait_operator.size() != traits) {
    throw std::invalid_argument("mt block-eigen data dimensions are incomplete");
  }
  for (std::size_t trait = 0; trait < traits; ++trait) {
    if (data.wy[trait].size() != markers || data.n[trait] <= 0) {
      throw std::invalid_argument("mt block-eigen trait summary dimensions are inconsistent");
    }
    const auto& view = data.operators.trait_operator[trait];
    if (view.marker_count != markers) {
      throw std::invalid_argument("mt block-eigen trait marker count is inconsistent");
    }
    sblr::core::validate_block_eigen_view(view);
  }
  for (const auto& pattern : model.models) {
    if (pattern.size() != traits) {
      throw std::invalid_argument("mt block-eigen model pattern trait dimension is inconsistent");
    }
  }
  for (const auto& set : model.sets) {
    for (int marker : set) {
      if (marker < 0 || static_cast<std::size_t>(marker) >= markers) {
        throw std::invalid_argument("mt block-eigen set marker is out of range");
      }
    }
  }
}

}  // namespace mt
}  // namespace sblr

#endif
