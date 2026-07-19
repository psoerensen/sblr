#ifndef SBLR_BLR_BED_FAMILY_TYPES_H
#define SBLR_BLR_BED_FAMILY_TYPES_H

#include <cstddef>
#include <cstdint>

namespace sblr { namespace core {

struct BedFamilyTaskIndex {
 int job;
 int trait;
 int chain;
};

inline BedFamilyTaskIndex make_bed_family_task_index(
 int job,
 int trait_count
) {
 return BedFamilyTaskIndex{job,job % trait_count,job / trait_count};
}

inline std::uint64_t resolve_bed_family_logical_chain_seed(
 int seed,
 int trait,
 int chain
) {
 return static_cast<unsigned int>(
  seed+1000003*(trait+1)+9176*(chain+1)
 );
}

// Borrowed immutable view of fit-owned, already decoded SNP-major storage.
template <class PackedGenotype>
struct BedPackedGenotypeView {
 const PackedGenotype& storage;
 const std::uint8_t* packed_markers;
 std::size_t packed_size;
 std::size_t marker_count;
 std::size_t sample_count;
 std::size_t bytes_per_marker;
 std::size_t stride;
};

}} // namespace sblr::core

#endif
