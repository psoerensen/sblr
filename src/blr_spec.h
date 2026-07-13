#ifndef SBLR_CORE_BLR_SPEC_H
#define SBLR_CORE_BLR_SPEC_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "blr_csr_contract.h"

namespace sblr {
namespace core {

enum class DataRepresentation { Csr };
enum class DesignContract { IndependentTraits };
enum class DataScaling { StandardizedGenotype };
enum class KernelFamily { Scalar };
enum class ModelFamily { BayesC };
enum class StateSpace { Binary };
enum class ProbabilityPolicyKind { GlobalBinary };
enum class ScalePolicyKind { Unit };
enum class TraitCovarianceKind { ScalarIndependent };
enum class ResidualCovarianceKind { ScalarIndependent };
enum class OperatorKind { Csr };
enum class BackendKind { UnscheduledCsrBayesC };

struct DataSpec {
  DataRepresentation representation = DataRepresentation::Csr;
  DesignContract design = DesignContract::IndependentTraits;
  std::size_t marker_count = 0;
  std::size_t trait_count = 0;
  std::vector<std::string> marker_ids;
  std::vector<std::string> trait_ids;
  std::vector<int> sample_size;
  DataScaling scaling = DataScaling::StandardizedGenotype;
  CsrResourceSpec csr;
};

struct ModelSpec {
  KernelFamily kernel = KernelFamily::Scalar;
  ModelFamily family = ModelFamily::BayesC;
  StateSpace state = StateSpace::Binary;
  ProbabilityPolicyKind probability = ProbabilityPolicyKind::GlobalBinary;
  ScalePolicyKind scale = ScalePolicyKind::Unit;
  TraitCovarianceKind trait_covariance =
    TraitCovarianceKind::ScalarIndependent;
  ResidualCovarianceKind residual_covariance =
    ResidualCovarianceKind::ScalarIndependent;
};

struct McmcControl {
  int nit = 0;
  int nburn = 0;
  int nthin = 1;
  int nchains = 1;
  int ncores = 1;
  int seed = 1;
  bool has_explicit_chain_seeds = false;
  std::vector<int> chain_seeds;
};

struct OutputSpec {
  bool marker_mean = true;
  bool marker_pip = true;
  bool parameter_traces = true;
  bool final_state = true;
  bool keep_chain_summaries = false;
};

struct ExecutionSpec {
  OperatorKind operator_kind = OperatorKind::Csr;
  BackendKind backend = BackendKind::UnscheduledCsrBayesC;
  bool scheduled = false;
};

struct ResolvedSpec {
  std::string schema_name = "blr_resolved_spec";
  int schema_version = 1;
  DataSpec data;
  ModelSpec model;
  McmcControl mcmc;
  OutputSpec output;
  ExecutionSpec execution;
};

inline bool has_empty_or_duplicate_ids(const std::vector<std::string>& ids) {
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i].empty()) return true;
    if (std::find(ids.begin(), ids.begin() + i, ids[i]) != ids.begin() + i) {
      return true;
    }
  }
  return false;
}

inline void validate_resolved_spec(const ResolvedSpec& spec) {
  if (spec.schema_name != "blr_resolved_spec") {
    throw std::invalid_argument("schema$name must be blr_resolved_spec");
  }
  if (spec.schema_version != 1) {
    throw std::invalid_argument("schema$version must be 1");
  }
  if (spec.data.marker_count == 0) {
    throw std::invalid_argument("data$n_markers must be positive");
  }
  if (spec.data.trait_count == 0) {
    throw std::invalid_argument("data$n_traits must be positive");
  }
  if (spec.data.marker_ids.size() != spec.data.marker_count ||
      has_empty_or_duplicate_ids(spec.data.marker_ids)) {
    throw std::invalid_argument(
      "data$marker_ids must be unique and match data$n_markers"
    );
  }
  if (spec.data.trait_ids.size() != spec.data.trait_count ||
      has_empty_or_duplicate_ids(spec.data.trait_ids)) {
    throw std::invalid_argument(
      "data$trait_ids must be unique and match data$n_traits"
    );
  }
  if (spec.data.sample_size.size() != spec.data.trait_count) {
    throw std::invalid_argument(
      "data$sample_size must contain one value per trait"
    );
  }
  for (int value : spec.data.sample_size) {
    if (value <= 0) {
      throw std::invalid_argument("data$sample_size values must be positive");
    }
  }
  validate_csr_resource(spec.data.csr);
  if (spec.data.csr.marker_count != spec.data.marker_count) {
    throw std::invalid_argument(
      "data$csr$marker_count must equal data$n_markers"
    );
  }
  if (spec.mcmc.nit <= 0) {
    throw std::invalid_argument("mcmc$nit must be positive");
  }
  if (spec.mcmc.nburn < 0) {
    throw std::invalid_argument("mcmc$nburn must be non-negative");
  }
  if (spec.mcmc.nthin <= 0) {
    throw std::invalid_argument("mcmc$nthin must be positive");
  }
  if (spec.mcmc.nchains <= 0) {
    throw std::invalid_argument("mcmc$nchains must be positive");
  }
  if (spec.mcmc.ncores <= 0) {
    throw std::invalid_argument("mcmc$ncores must be positive");
  }
  if (spec.mcmc.has_explicit_chain_seeds) {
    if (spec.mcmc.chain_seeds.size() !=
        static_cast<std::size_t>(spec.mcmc.nchains)) {
      throw std::invalid_argument(
        "mcmc$chain_seeds must contain one seed per chain"
      );
    }
  } else if (!spec.mcmc.chain_seeds.empty()) {
    throw std::invalid_argument(
      "mcmc$chain_seeds must be empty when explicit seeds are absent"
    );
  }
  if (spec.execution.scheduled) {
    throw std::invalid_argument(
      "execution$scheduled must be false in Phase 1"
    );
  }
}

}  // namespace core
}  // namespace sblr

#endif
