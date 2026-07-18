#ifndef SBLR_BLR_BED_BAYESR_AUDIT_TYPES_H
#define SBLR_BLR_BED_BAYESR_AUDIT_TYPES_H

// Binding-neutral Phase 13A vocabulary. These contracts describe the current
// packed-BED BayesR implementation; they do not invoke or replace the sampler.

#include "blr_scheduled_execution_types.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr {
namespace audit {

struct BedBayesRComponentSpec {
  std::size_t null_component = 0;
  std::vector<double> scales;
  std::vector<double> initial_probabilities;
  std::vector<double> dirichlet_alpha;

  void validate() const {
    const std::size_t k = scales.size();
    if (k < 2) throw std::invalid_argument("BayesR requires null plus at least one non-null component");
    if (null_component != 0) throw std::invalid_argument("BayesR null component index must be zero");
    if (initial_probabilities.size() != k || dirichlet_alpha.size() != k)
      throw std::invalid_argument("BayesR component vectors must have identical lengths");
    if (scales[0] != 0.0) throw std::invalid_argument("BayesR null component scale must be zero");
    double total = 0.0;
    for (std::size_t i = 0; i < k; ++i) {
      if (!std::isfinite(scales[i]) || (i > 0 && scales[i] <= 0.0))
        throw std::invalid_argument("BayesR non-null scales must be finite and positive");
      if (!std::isfinite(initial_probabilities[i]) || initial_probabilities[i] < 0.0)
        throw std::invalid_argument("BayesR probabilities must be finite and nonnegative");
      if (!std::isfinite(dirichlet_alpha[i]) || dirichlet_alpha[i] <= 0.0)
        throw std::invalid_argument("BayesR Dirichlet parameters must be finite and positive");
      total += initial_probabilities[i];
    }
    if (std::abs(total - 1.0) > 1e-12)
      throw std::invalid_argument("BayesR probabilities must sum to one");
  }
};

struct BedBayesRSchedulerSpec {
  sblr::core::ScheduledSweepControl sweep;
  sblr::core::NullSkipControl null_skip;
  sblr::core::CandidateControl candidate;
  bool skip_nulls_burnin_only = false;

  void validate() const {
    sweep.validate();
    null_skip.validate();
    candidate.validate();
  }
};

struct BedBayesRExecutionAuditSpec {
  int iterations = 0;
  int burnin = 0;
  int thinning = 1;
  int chains = 1;
  int cores = 1;
  std::uint64_t seed = 0;
  BedBayesRComponentSpec components;
  BedBayesRSchedulerSpec scheduler;

  void validate() const {
    if (iterations <= 0 || burnin < 0 || thinning <= 0)
      throw std::invalid_argument("BayesR iteration controls are invalid");
    if (chains <= 0 || cores <= 0)
      throw std::invalid_argument("BayesR chain and core counts must be positive");
    components.validate();
    scheduler.validate();
  }
};

struct BedBayesROwnershipVocabulary {
  const char* genotype_owner = "fit";
  const char* genotype_access = "borrowed immutable";
  const char* chain_state_owner = "logical chain";
  const char* rng_owner = "logical chain";
  const char* distribution_owner = "logical chain";
  const char* file_lifetime = "fit-local decode only";
  bool per_chain_genotype_copy = false;
  bool mcmc_disk_access = false;
};

struct BedBayesRChainResultVocabulary {
  bool marker_means = true;
  bool marker_pip = true;
  bool component_probabilities = true;
  bool final_effects_and_assignments = true;
  bool variance_traces = true;
  bool cpo = true;
  bool retained_count = true;
  bool timing_and_failure = true;
};

struct BedBayesRExecutionResultVocabulary {
  bool chain_mean_sd_min_max = true;
  bool component_probability_aggregation = true;
  bool trace_aggregation = true;
  bool final_state_aggregation = true;
  bool cpo_aggregation = true;
};

} // namespace audit
} // namespace sblr

#endif
