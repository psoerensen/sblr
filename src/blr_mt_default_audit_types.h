#ifndef SBLR_BLR_MT_DEFAULT_AUDIT_TYPES_H
#define SBLR_BLR_MT_DEFAULT_AUDIT_TYPES_H

#include <cstddef>
#include <string>
#include <vector>

namespace sblr {
namespace audit {

// Audit-only vocabulary: intentionally not included by production code.
struct MtDataSpec {
  std::size_t samples;
  std::size_t markers;
  std::size_t traits;
  std::string design_orientation;
  std::string phenotype_orientation;
  std::vector<std::string> marker_order;
  std::vector<std::string> trait_order;
};

struct MtModelSpec {
  std::vector<std::vector<int> > selection_patterns;
  std::string covariance_policy;
  std::string probability_policy;
  bool update_marker_covariance;
  bool update_residual_covariance;
  bool update_pattern_probabilities;
};

struct MtCovariancePriorSpec {
  double marker_degrees_of_freedom;
  double residual_degrees_of_freedom;
  std::vector<double> marker_scale_column_major;
  std::vector<double> residual_scale_column_major;
};

struct MtExecutionAuditSpec {
  int iterations;
  int burnin;
  int thinning;
  int resolved_seed;
  int cores;
  bool single_chain;
};

struct MtOwnershipAuditSpec {
  std::string design_owner;
  std::string phenotype_owner;
  std::string workspace_owner;
  std::string rng_owner;
  std::string result_owner;
};

struct MtChainResultVocabulary {
  std::size_t marker_trait_values;
  std::size_t trace_trait_values;
  std::size_t covariance_values;
};

struct MtExecutionResultVocabulary {
  std::vector<std::string> positional_fields;
  bool legacy_positional;
};

}  // namespace audit
}  // namespace sblr

#endif
