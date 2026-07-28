#ifndef SBLR_BLR_MT_BED_CONVERGENCE_TYPES_H
#define SBLR_BLR_MT_BED_CONVERGENCE_TYPES_H

#include <string>
#include <vector>

namespace sblr {
namespace mt {

struct MtBedConvergenceQuantity {
 std::string group;
 int trait=-1;
 int trait2=-1;
 int marker=-1;
 int component=-1;
 int pattern=-1;
 int annotation=-1;
 int stick=-1;
 int tier=1;
 bool updated=false;
 bool derived=false;
 bool structural=false;
};

struct MtBedConvergenceTraceBundle {
 int nchains=0;
 int postburn_draws=0;
 std::vector<MtBedConvergenceQuantity> quantities;
 // Quantity-major, then chain-major, with iteration varying fastest.
 std::vector<double> values;
};

}  // namespace mt
}  // namespace sblr

#endif
