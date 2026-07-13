#ifndef SBLR_CORE_BLR_RESULT_H
#define SBLR_CORE_BLR_RESULT_H

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace core {

// Canonical flattened storage is row-major with the documented logical shape.
struct MarkerSummary {
  std::size_t marker_count = 0;
  std::size_t trait_count = 0;
  std::vector<double> effect_mean;  // markers x traits
  std::vector<double> pip;          // markers x traits
};

struct StateSummary {
  std::vector<double> final_effect;  // markers x traits
  std::vector<int> final_state;      // markers x traits
};

struct TraceSummary {
  std::size_t retained_samples = 0;
  std::size_t parameter_count = 0;
  std::vector<double> values;  // retained samples x parameter dimension
};

struct VarianceSummary {
  std::size_t trait_count = 0;
  std::vector<double> trait_covariance;     // traits x traits
  std::vector<double> residual_covariance;  // traits x traits
};

struct DiagnosticsSummary {
  std::size_t retained_samples = 0;
};

struct ChainSummary {
  int chain_index = 0;
  MarkerSummary marker;
};

struct OptionalResultDimensions {
  bool has_component_probability = false;
  std::size_t component_count = 0;  // markers x components x traits
  bool has_pattern_probability = false;
  std::size_t pattern_count = 0;    // markers x patterns
};

struct BlrResult {
  MarkerSummary marker;
  StateSummary state;
  TraceSummary trace;
  VarianceSummary variance;
  DiagnosticsSummary diagnostics;
  std::vector<ChainSummary> chains;
  OptionalResultDimensions optional;
};

inline std::size_t checked_product(std::size_t left, std::size_t right,
                                   const char* field) {
  if (left != 0 && right > static_cast<std::size_t>(-1) / left) {
    throw std::invalid_argument(field);
  }
  return left * right;
}

inline void validate_marker_summary(const MarkerSummary& marker) {
  if (marker.marker_count == 0) {
    throw std::invalid_argument("result$marker_count must be positive");
  }
  if (marker.trait_count == 0) {
    throw std::invalid_argument("result$trait_count must be positive");
  }
  const std::size_t expected = checked_product(
    marker.marker_count, marker.trait_count,
    "result marker dimensions overflow"
  );
  if (marker.effect_mean.size() != expected) {
    throw std::invalid_argument(
      "result$marker_effect_length must equal markers * traits"
    );
  }
  if (marker.pip.size() != expected) {
    throw std::invalid_argument(
      "result$marker_pip_length must equal markers * traits"
    );
  }
}

inline void validate_blr_result(const BlrResult& result) {
  validate_marker_summary(result.marker);
  const std::size_t marker_trait = checked_product(
    result.marker.marker_count, result.marker.trait_count,
    "result marker dimensions overflow"
  );
  if (result.state.final_effect.size() != marker_trait) {
    throw std::invalid_argument(
      "result$final_effect_length must equal markers * traits"
    );
  }
  if (result.state.final_state.size() != marker_trait) {
    throw std::invalid_argument(
      "result$final_state_length must equal markers * traits"
    );
  }
  const std::size_t trace_size = checked_product(
    result.trace.retained_samples, result.trace.parameter_count,
    "result trace dimensions overflow"
  );
  if (result.trace.values.size() != trace_size) {
    throw std::invalid_argument(
      "result$trace_length must equal retained samples * parameter dimension"
    );
  }
  if (result.variance.trait_count != result.marker.trait_count) {
    throw std::invalid_argument(
      "result variance trait count must equal marker trait count"
    );
  }
  const std::size_t covariance_size = checked_product(
    result.marker.trait_count, result.marker.trait_count,
    "result covariance dimensions overflow"
  );
  if (result.variance.trait_covariance.size() != covariance_size) {
    throw std::invalid_argument(
      "result$trait_covariance_length must equal traits * traits"
    );
  }
  if (result.variance.residual_covariance.size() != covariance_size) {
    throw std::invalid_argument(
      "result$residual_covariance_length must equal traits * traits"
    );
  }
  if (result.diagnostics.retained_samples != result.trace.retained_samples) {
    throw std::invalid_argument(
      "result diagnostic retained samples must equal trace retained samples"
    );
  }
  if (result.optional.has_component_probability &&
      result.optional.component_count == 0) {
    throw std::invalid_argument(
      "result component count must be positive when present"
    );
  }
  if (!result.optional.has_component_probability &&
      result.optional.component_count != 0) {
    throw std::invalid_argument(
      "result component count must be zero when absent"
    );
  }
  if (result.optional.has_pattern_probability &&
      result.optional.pattern_count == 0) {
    throw std::invalid_argument(
      "result pattern count must be positive when present"
    );
  }
  if (!result.optional.has_pattern_probability &&
      result.optional.pattern_count != 0) {
    throw std::invalid_argument(
      "result pattern count must be zero when absent"
    );
  }
  for (std::size_t index = 0; index < result.chains.size(); ++index) {
    const ChainSummary& chain = result.chains[index];
    if (chain.chain_index != static_cast<int>(index)) {
      throw std::invalid_argument(
        "result chain indices must be contiguous and zero-based"
      );
    }
    validate_marker_summary(chain.marker);
    if (chain.marker.marker_count != result.marker.marker_count ||
        chain.marker.trait_count != result.marker.trait_count) {
      throw std::invalid_argument(
        "result chain marker dimensions must match aggregate dimensions"
      );
    }
  }
}

}  // namespace core
}  // namespace sblr

#endif
