#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_RNG_H
#define SBLR_BLR_BED_SCHEDULED_BAYESC_RNG_H

#include <cstdint>
#include <random>

namespace sblr { namespace core {

// Fit-bounded RNG state owned by one logical scheduled packed-BED BayesC chain.
struct BedScheduledBayesCChainRng {
 std::mt19937 engine;
 std::normal_distribution<double> normal;
 std::uniform_real_distribution<double> uniform;

 explicit BedScheduledBayesCChainRng(std::uint64_t seed)
  : engine(static_cast<std::mt19937::result_type>(seed)),
    normal(0.0, 1.0),
    uniform(0.0, 1.0) {}
};

} }

#endif
