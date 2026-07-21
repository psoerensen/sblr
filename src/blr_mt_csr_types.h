#ifndef SBLR_BLR_MT_CSR_TYPES_H
#define SBLR_BLR_MT_CSR_TYPES_H

#include "blr_mt_default_types.h"
#include "blr_sparse_ld_csr.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace mt {

struct MtSparseLdBundleView {
 std::size_t marker_count=0;
 std::vector<sblr::core::SparseLdCsrView> trait_ld;
};

struct MtCsrDataView {
 const std::vector<std::vector<double>>& wy;
 const std::vector<double>& yy;
 const std::vector<int>& n;
 MtSparseLdBundleView ld;
};

inline void validate_mt_csr_data(
 const MtCsrDataView& data,
 const MtDefaultModelSpec& model
) {
 const std::size_t traits=data.wy.size();
 const std::size_t markers=data.ld.marker_count;
 if (traits==0 || markers==0 || data.yy.size()!=traits ||
     data.n.size()!=traits || data.ld.trait_ld.size()!=traits) {
  throw std::invalid_argument("mt CSR data dimensions are incomplete");
 }
 for (std::size_t trait=0; trait<traits; ++trait) {
  if (data.wy[trait].size()!=markers || data.n[trait]<=0) {
   throw std::invalid_argument("mt CSR trait summary dimensions are inconsistent");
  }
  if (data.ld.trait_ld[trait].marker_count!=markers) {
   throw std::invalid_argument("mt CSR trait marker count is inconsistent");
  }
  sblr::core::validate_sparse_ld_csr_view(data.ld.trait_ld[trait]);
 }
 for (const auto& pattern : model.models) {
  if (pattern.size()!=traits) {
   throw std::invalid_argument("mt CSR model pattern trait dimension is inconsistent");
  }
 }
 for (const auto& set : model.sets) {
  for (int marker : set) {
   if (marker<0 || static_cast<std::size_t>(marker)>=markers) {
    throw std::invalid_argument("mt CSR set marker is out of range");
   }
  }
 }
}

}  // namespace mt
}  // namespace sblr

#endif
