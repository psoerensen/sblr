#ifndef SBLR_BLR_MT_BED_TYPES_H
#define SBLR_BLR_MT_BED_TYPES_H

#include "blr_bed_family_types.h"
#include "blr_mt_default_types.h"

#include <armadillo>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace sblr {
namespace mt {

struct MtBedMarkerMap {
 double value[4]={0.0, 0.0, 0.0, 0.0};
 double xx=0.0;
};

template <class PackedGenotype>
struct MtBedDataView {
 sblr::core::BedPackedGenotypeView<PackedGenotype> genotype;
 const std::vector<MtBedMarkerMap>& marker_maps;
 const arma::mat& phenotype;
 const arma::mat& marker_wy;
 const std::vector<int>& marker_order;
};

struct MtBedInitialState {
 std::vector<std::vector<double>> beta;
 std::vector<std::vector<double>> b;
 std::vector<std::vector<int>> state;
 std::vector<int> component;
 arma::mat B;
 arma::mat E;
 std::vector<double> pi;
};

struct MtBedExecutionSpec {
 bool updateB=false;
 bool updateE=false;
 bool updatePi=false;
 std::string residual_covariance;
 int nit=0;
 int nburn=0;
 int nthin=1;
 std::uint32_t seed=1;
 int method=4;
};

struct MtBedCoreDiagnostics {
 std::size_t marker_cholesky_jitter_attempts=0;
 double marker_cholesky_max_increment=0.0;
 std::size_t full_e_updates=0;
 std::size_t diagonal_e_updates=0;
};

struct MtBedMarkerKernelModel {
 arma::mat precision;
 arma::mat lower;
 arma::vec rhs;
 arma::vec mean;
 arma::mat covariance;
 double log_weight=0.0;
};

struct MtBedMarkerKernelResult {
 std::vector<MtBedMarkerKernelModel> models;
 std::vector<double> log_weight;
 std::vector<double> probability;
};

struct MtBedCoreResult {
 MtDefaultCoreResult posterior;
 MtBedCoreDiagnostics diagnostics;
};

}  // namespace mt
}  // namespace sblr

#endif
