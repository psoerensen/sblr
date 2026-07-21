#ifndef SBLR_BLR_MT_LD_ACCESS_H
#define SBLR_BLR_MT_LD_ACCESS_H

#include "blr_mt_csr_types.h"
#include "blr_mt_block_eigen_types.h"

namespace sblr {
namespace mt {

inline std::size_t mt_trait_count(const MtDefaultDataView& data) {
 return data.wy.size();
}
inline std::size_t mt_marker_count(const MtDefaultDataView& data) {
 return data.wy.empty() ? 0 : data.wy[0].size();
}
inline double mt_diagonal(const MtDefaultDataView& data, int trait, int marker) {
 return data.ww[static_cast<std::size_t>(trait)][static_cast<std::size_t>(marker)];
}
inline void mt_rebuild_residuals(
 const MtDefaultDataView& data,
 const std::vector<std::vector<double>>& effects,
 std::vector<std::vector<double>>& residuals
) {
 const int nt=static_cast<int>(mt_trait_count(data));
 const int m=static_cast<int>(mt_marker_count(data));
 for (int i=0; i<m; ++i) {
  for (int t=0; t<nt; ++t) residuals[t][i]=data.wy[t][i];
 }
 for (int i=0; i<m; ++i) {
  for (int t=0; t<nt; ++t) {
   if (effects[t][i]!=0.0) {
    for (std::size_t j=0; j<data.XXindices[i].size(); ++j) {
     residuals[t][data.XXindices[i][j]] -=
      data.XXvalues[t][i][j]*effects[t][i];
    }
   }
  }
 }
}
inline void mt_apply_marker_difference(
 const MtDefaultDataView& data, int trait, int marker, double difference,
 std::vector<double>& residual
) {
 const std::size_t nnz=data.XXindices[marker].size();
 for (std::size_t j=0; j<nnz; ++j) {
  residual[data.XXindices[marker][j]] -= data.XXvalues[trait][marker][j]*difference;
 }
}

inline std::size_t mt_trait_count(const MtCsrDataView& data) {
 return data.wy.size();
}
inline std::size_t mt_marker_count(const MtCsrDataView& data) {
 return data.ld.marker_count;
}
inline double mt_diagonal(const MtCsrDataView& data, int trait, int marker) {
 return data.ld.trait_ld[static_cast<std::size_t>(trait)].diag()(
  static_cast<arma::uword>(marker));
}
inline void mt_rebuild_residuals(
 const MtCsrDataView& data,
 const std::vector<std::vector<double>>& effects,
 std::vector<std::vector<double>>& residuals
) {
 for (std::size_t trait=0; trait<data.wy.size(); ++trait) {
  data.ld.trait_ld[trait].rebuild(data.wy[trait], effects[trait], residuals[trait]);
 }
}
inline void mt_apply_marker_difference(
 const MtCsrDataView& data, int trait, int marker, double difference,
 std::vector<double>& residual
) {
 residual[static_cast<std::size_t>(marker)] -=
  mt_diagonal(data, trait, marker)*difference;
 data.ld.trait_ld[static_cast<std::size_t>(trait)].apply_offdiag(
  marker, difference, residual);
}

inline std::size_t mt_trait_count(const MtBlockEigenDataView& data) {
 return data.wy.size();
}
inline std::size_t mt_marker_count(const MtBlockEigenDataView& data) {
 return data.operators.marker_count;
}
inline double mt_diagonal(const MtBlockEigenDataView& data, int trait, int marker) {
 return data.operators.trait_operator[static_cast<std::size_t>(trait)].diag()(
  static_cast<arma::uword>(marker));
}
inline void mt_rebuild_residuals(
 const MtBlockEigenDataView& data,
 const std::vector<std::vector<double>>& effects,
 std::vector<std::vector<double>>& residuals
) {
 for (std::size_t trait=0; trait<data.wy.size(); ++trait) {
  data.operators.trait_operator[trait].rebuild(
   data.wy[trait], effects[trait], residuals[trait]);
 }
}
inline void mt_apply_marker_difference(
 const MtBlockEigenDataView& data, int trait, int marker, double difference,
 std::vector<double>& residual
) {
 residual[static_cast<std::size_t>(marker)] -=
  mt_diagonal(data, trait, marker)*difference;
 data.operators.trait_operator[static_cast<std::size_t>(trait)].apply_offdiag(
  marker, difference, residual);
}

}  // namespace mt
}  // namespace sblr

#endif
