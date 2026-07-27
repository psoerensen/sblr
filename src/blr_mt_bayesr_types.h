#ifndef SBLR_BLR_MT_BAYESR_TYPES_H
#define SBLR_BLR_MT_BAYESR_TYPES_H

#include <armadillo>
#include <cstddef>
#include <string>
#include <vector>

namespace sblr {
namespace mt {

struct MtJointStateSpec {
 std::vector<std::vector<int>> patterns;
 std::vector<int> component;
 std::vector<double> multiplier;
 std::vector<std::string> names;
 int component_count=0;
 bool scaled=false;
};

struct MtJointMarkerKernelState {
 arma::mat precision;
 arma::mat lower;
 arma::vec rhs;
 arma::vec mean;
 double log_weight=0.0;
};

struct MtJointMarkerKernelResult {
 std::vector<MtJointMarkerKernelState> states;
 std::vector<double> probability;
};

}  // namespace mt
}  // namespace sblr

#endif
