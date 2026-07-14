#ifndef SBLR_CORE_BLR_SCALAR_EXECUTION_H
#define SBLR_CORE_BLR_SCALAR_EXECUTION_H

#include "st_chain_utils.h"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr {
namespace core {

// Trait-major scalar execution identity. This is the ordering already used by
// ordinary CSR BayesC and CSR BayesR: task = trait * nchains + chain.
struct ScalarChainTask {
  std::size_t trait_index = 0;
  std::size_t chain_index = 0;
  std::size_t task_index = 0;
};

inline std::vector<ScalarChainTask> make_scalar_chain_tasks(
  std::size_t trait_count,
  std::size_t chain_count
) {
  if (trait_count == 0) {
    throw std::invalid_argument("scalar task trait_count must be positive");
  }
  if (chain_count == 0) {
    throw std::invalid_argument("scalar task chain_count must be positive");
  }
  if (trait_count > std::numeric_limits<std::size_t>::max() / chain_count) {
    throw std::overflow_error("scalar task count overflows size_t");
  }
  const std::size_t task_count = trait_count * chain_count;
  if (task_count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("scalar task count exceeds the execution index range");
  }

  std::vector<ScalarChainTask> tasks;
  tasks.reserve(task_count);
  for (std::size_t task = 0; task < task_count; ++task) {
    ScalarChainTask identity;
    identity.trait_index = task / chain_count;
    identity.chain_index = task % chain_count;
    identity.task_index = task;
    tasks.push_back(identity);
  }
  return tasks;
}

inline void validate_scalar_chain_seeds(
  std::size_t chain_count,
  const std::vector<int>& explicit_chain_seeds
) {
  if (!explicit_chain_seeds.empty() &&
      explicit_chain_seeds.size() != chain_count) {
    throw std::invalid_argument(
      "scalar explicit chain seeds must match chain_count"
    );
  }
}

// Delegates to the established seed helpers so integer conversion and stream
// mapping remain exactly the same as the current BayesC and BayesR backends.
inline unsigned int resolve_scalar_chain_seed(
  int seed,
  std::size_t chain_count,
  const std::vector<int>& explicit_chain_seeds,
  const ScalarChainTask& task
) {
  validate_scalar_chain_seeds(chain_count, explicit_chain_seeds);
  if (task.chain_index >= chain_count) {
    throw std::invalid_argument("scalar task chain index is out of range");
  }
  if (!explicit_chain_seeds.empty()) {
    return stblr_seed_with_chain_base(
      explicit_chain_seeds[task.chain_index],
      static_cast<int>(task.trait_index)
    );
  }
  if (chain_count == 1) {
    return stblr_trait_seed(seed, static_cast<int>(task.trait_index));
  }
  return stblr_chain_seed(
    seed,
    static_cast<int>(task.trait_index),
    static_cast<int>(task.chain_index)
  );
}

inline bool scalar_iteration_is_retained(
  int iteration,
  int burnin,
  int thinning
) {
  return iteration >= burnin && (iteration - burnin) % thinning == 0;
}

// Lightweight execution metadata shared by scalar chain implementations.
// It owns no operator or chain state. A physical chain continues to own its
// RNG engine and every stateful normal, uniform, chi-square, gamma, discrete,
// categorical, or proposal distribution used by its model-specific kernel.
struct ScalarChainExecutionStatus {
  ScalarChainTask task;
  unsigned int seed = 0u;
  bool failed = false;
  std::string failure_message;
  double elapsed_seconds = 0.0;
  double retained_samples = 0.0;
};

}  // namespace core
}  // namespace sblr

#endif
