#ifndef SBLR_BLR_MT_BED_ACCESS_H
#define SBLR_BLR_MT_BED_ACCESS_H

#include "blr_mt_bed_types.h"
#include "packed_bed.h"

#include <armadillo>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace mt {

template <class PackedGenotype>
inline void validate_mt_bed_genotype_view(
 const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype
) {
 if (genotype.marker_count==0)
  throw std::invalid_argument("mtblr_bed_internal: marker_count must be positive.");
 if (genotype.sample_count<=1)
  throw std::invalid_argument("mtblr_bed_internal: sample_count must exceed one.");
 if (genotype.packed_markers==nullptr)
  throw std::invalid_argument("mtblr_bed_internal: packed genotype pointer is null.");
 const std::size_t expected=(genotype.sample_count+3u)/4u;
 if (genotype.bytes_per_marker!=expected ||
     genotype.stride<genotype.bytes_per_marker ||
     genotype.packed_size<genotype.marker_count*genotype.stride)
  throw std::invalid_argument("mtblr_bed_internal: invalid packed genotype dimensions.");
}

template <class PackedGenotype>
inline void decode_mt_bed_marker(
 const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype,
 const MtBedMarkerMap& map,
 std::size_t marker,
 arma::vec& workspace
) {
 if (marker>=genotype.marker_count)
  throw std::out_of_range("decode_mt_bed_marker: marker is out of range.");
 if (workspace.n_elem!=genotype.sample_count)
  throw std::invalid_argument("decode_mt_bed_marker: workspace length is invalid.");
 const std::uint8_t* packed=
  genotype.packed_markers+marker*genotype.stride;
 for (std::size_t sample=0; sample<genotype.sample_count; ++sample)
  workspace(static_cast<arma::uword>(sample))=
   map.value[get_bed_code(packed, static_cast<int>(sample))];
}

template <class PackedGenotype>
inline std::vector<MtBedMarkerMap> build_mt_bed_marker_maps(
 const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype,
 const std::vector<double>& af
) {
 validate_mt_bed_genotype_view(genotype);
 if (af.size()!=genotype.marker_count)
  throw std::invalid_argument("mtblr_bed_internal: af length must equal marker count.");
 std::vector<MtBedMarkerMap> maps(genotype.marker_count);
 arma::vec workspace(genotype.sample_count);
 for (std::size_t marker=0; marker<genotype.marker_count; ++marker) {
  const double p=af[marker];
  if (!std::isfinite(p) || p<=0.0 || p>=1.0)
   throw std::invalid_argument(
    "mtblr_bed_internal: allele frequencies must be finite and in (0, 1).");
  const double denominator=std::sqrt(2.0*p*(1.0-p));
  maps[marker].value[0]=(2.0-2.0*p)/denominator;
  maps[marker].value[1]=0.0;
  maps[marker].value[2]=(1.0-2.0*p)/denominator;
  maps[marker].value[3]=(0.0-2.0*p)/denominator;
  decode_mt_bed_marker(genotype, maps[marker], marker, workspace);
  long double xx=0.0L;
  for (std::size_t sample=0; sample<genotype.sample_count; ++sample) {
   const long double x=workspace(static_cast<arma::uword>(sample));
   xx+=x*x;
  }
  maps[marker].xx=static_cast<double>(xx);
  if (!std::isfinite(maps[marker].xx) || maps[marker].xx<=0.0)
   throw std::invalid_argument(
    "mtblr_bed_internal: every selected marker must have finite positive x'x.");
 }
 return maps;
}

template <class PackedGenotype>
inline arma::mat compute_mt_bed_marker_wy(
 const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype,
 const std::vector<MtBedMarkerMap>& maps,
 const arma::mat& phenotype,
 arma::vec& workspace
) {
 arma::mat wy(genotype.marker_count, phenotype.n_cols, arma::fill::zeros);
 for (std::size_t marker=0; marker<genotype.marker_count; ++marker) {
  decode_mt_bed_marker(genotype, maps[marker], marker, workspace);
  for (arma::uword trait=0; trait<phenotype.n_cols; ++trait) {
   long double score=0.0L;
   for (std::size_t sample=0; sample<genotype.sample_count; ++sample)
    score+=static_cast<long double>(
     workspace(static_cast<arma::uword>(sample)))*
     static_cast<long double>(
      phenotype(static_cast<arma::uword>(sample),trait));
   wy(static_cast<arma::uword>(marker),trait)=static_cast<double>(score);
  }
 }
 return wy;
}

inline std::vector<int> compute_mt_bed_marker_order(
 const arma::mat& marker_wy,
 const std::vector<MtBedMarkerMap>& maps
) {
 const std::size_t m=maps.size();
 std::vector<double> score(m, 0.0);
 for (std::size_t marker=0; marker<m; ++marker)
  for (arma::uword trait=0; trait<marker_wy.n_cols; ++trait) {
   const double ratio=
    marker_wy(static_cast<arma::uword>(marker),trait)/maps[marker].xx;
   score[marker]=std::max(score[marker], ratio*ratio);
  }
 std::vector<int> order(m);
 std::iota(order.begin(), order.end(), 0);
 std::sort(order.begin(), order.end(), [&](int left, int right) {
  if (score[static_cast<std::size_t>(left)]!=
      score[static_cast<std::size_t>(right)])
   return score[static_cast<std::size_t>(left)]>
          score[static_cast<std::size_t>(right)];
  return left<right;
 });
 return order;
}

inline std::vector<std::vector<double>> mt_bed_trait_major(
 const arma::mat& marker_by_trait
) {
 std::vector<std::vector<double>> value(
  marker_by_trait.n_cols,
  std::vector<double>(marker_by_trait.n_rows, 0.0));
 for (arma::uword trait=0; trait<marker_by_trait.n_cols; ++trait)
  for (arma::uword marker=0; marker<marker_by_trait.n_rows; ++marker)
   value[trait][marker]=marker_by_trait(marker,trait);
 return value;
}

}  // namespace mt
}  // namespace sblr

#endif
