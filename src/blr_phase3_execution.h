#ifndef SBLR_BLR_PHASE3_EXECUTION_H
#define SBLR_BLR_PHASE3_EXECUTION_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

struct BlrPhase3ExecutionContract {
  int seed_contract_version = 0;
  int retention_contract_version = 0;
  int scheduler_version = 0;
  std::vector<std::uint32_t> task_seeds;
  std::vector<std::string> task_ids;
  std::vector<int> retained_transition_indices;

  bool active() const {
    return seed_contract_version == 1 ||
      retention_contract_version == 1 || scheduler_version == 1;
  }

  bool retained(int post_burn_iteration) const {
    if (retention_contract_version != 1) return false;
    return std::binary_search(
      retained_transition_indices.begin(),
      retained_transition_indices.end(), post_burn_iteration);
  }
};

inline std::uint64_t blr_phase3_fnv1a64(const std::string& value) {
  std::uint64_t hash = UINT64_C(0xcbf29ce484222325);
  for (unsigned char byte : value) {
    hash ^= static_cast<std::uint64_t>(byte);
    hash *= UINT64_C(0x00000100000001b3);
  }
  return hash;
}

inline std::uint64_t blr_phase3_splitmix64(std::uint64_t value) {
  value += UINT64_C(0x9e3779b97f4a7c15);
  value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
  return value ^ (value >> 31);
}

inline std::uint32_t blr_phase3_seed_v1(
    std::uint32_t user_seed,
    int mode_code,
    const std::string& identity,
    std::uint32_t chain_index) {
  std::uint64_t mixed = static_cast<std::uint64_t>(user_seed);
  mixed = blr_phase3_splitmix64(mixed ^ UINT64_C(1));
  mixed = blr_phase3_splitmix64(
    mixed ^ static_cast<std::uint64_t>(mode_code));
  mixed = blr_phase3_splitmix64(mixed ^ blr_phase3_fnv1a64(identity));
  mixed = blr_phase3_splitmix64(
    mixed ^ static_cast<std::uint64_t>(chain_index));
  return static_cast<std::uint32_t>(mixed ^ (mixed >> 32));
}

inline BlrPhase3ExecutionContract parse_blr_phase3_execution_contract(
    Rcpp::Nullable<Rcpp::List> value,
    int expected_tasks,
    int sampling_iterations) {
  BlrPhase3ExecutionContract out;
  if (value.isNull()) return out;

  Rcpp::List input(value);
  const std::vector<std::string> required = {
    "seed_contract_version", "retention_contract_version",
    "scheduler_version", "task_seeds", "task_ids",
    "retained_transition_indices"
  };
  Rcpp::CharacterVector supplied = input.names();
  if (supplied.size() != static_cast<R_xlen_t>(required.size())) {
    throw std::runtime_error(
      "execution_contract must contain exactly the Phase 3 fields.");
  }
  for (const std::string& name : required) {
    if (!input.containsElementNamed(name.c_str())) {
      throw std::runtime_error(
        "execution_contract is missing required field '" + name + "'.");
    }
  }

  out.seed_contract_version = Rcpp::as<int>(input["seed_contract_version"]);
  out.retention_contract_version =
    Rcpp::as<int>(input["retention_contract_version"]);
  out.scheduler_version = Rcpp::as<int>(input["scheduler_version"]);
  if (out.seed_contract_version != 1 ||
      out.retention_contract_version != 1 ||
      out.scheduler_version != 1) {
    throw std::runtime_error(
      "execution_contract versions must all equal 1 for Phase 3 execution.");
  }

  Rcpp::NumericVector seeds = input["task_seeds"];
  Rcpp::CharacterVector task_ids = input["task_ids"];
  if (seeds.size() != expected_tasks || task_ids.size() != expected_tasks) {
    throw std::runtime_error(
      "execution_contract task seeds and IDs must match the logical task count.");
  }
  out.task_seeds.reserve(static_cast<std::size_t>(expected_tasks));
  out.task_ids.reserve(static_cast<std::size_t>(expected_tasks));
  for (int task = 0; task < expected_tasks; ++task) {
    const double seed = seeds[task];
    if (!std::isfinite(seed) || seed < 0.0 ||
        seed > 4294967295.0 || seed != std::floor(seed)) {
      throw std::runtime_error(
        "execution_contract task seeds must be exact uint32 values.");
    }
    if (task_ids[task] == NA_STRING) {
      throw std::runtime_error(
        "execution_contract task IDs must be nonmissing strings.");
    }
    const std::string task_id = Rcpp::as<std::string>(task_ids[task]);
    if (task_id.empty()) {
      throw std::runtime_error(
        "execution_contract task IDs must be nonempty strings.");
    }
    out.task_seeds.push_back(static_cast<std::uint32_t>(seed));
    out.task_ids.push_back(task_id);
  }
  std::vector<std::string> sorted_ids = out.task_ids;
  std::sort(sorted_ids.begin(), sorted_ids.end());
  if (std::adjacent_find(sorted_ids.begin(), sorted_ids.end()) !=
      sorted_ids.end()) {
    throw std::runtime_error("execution_contract task IDs must be unique.");
  }

  Rcpp::IntegerVector retained = input["retained_transition_indices"];
  int previous = 0;
  for (R_xlen_t index = 0; index < retained.size(); ++index) {
    const int iteration = retained[index];
    if (iteration == NA_INTEGER || iteration <= previous ||
        iteration < 1 || iteration > sampling_iterations) {
      throw std::runtime_error(
        "execution_contract retained indices must be strictly increasing "
        "post-burn iterations within sampling_iterations.");
    }
    out.retained_transition_indices.push_back(iteration);
    previous = iteration;
  }
  return out;
}

inline std::uint32_t blr_phase3_task_seed(
    const BlrPhase3ExecutionContract& contract,
    int task,
    std::uint32_t legacy_seed) {
  if (contract.seed_contract_version == 1) {
    return contract.task_seeds[static_cast<std::size_t>(task)];
  }
  return legacy_seed;
}

inline bool blr_phase3_iteration_is_retained(
    const BlrPhase3ExecutionContract& contract,
    int iteration,
    int burn_in,
    int thin_interval) {
  if (iteration < burn_in) return false;
  if (contract.retention_contract_version == 1) {
    return contract.retained(iteration - burn_in + 1);
  }
  return ((iteration - burn_in) % thin_interval) == 0;
}

inline int blr_phase3_configured_workers(int requested, int tasks) {
  return std::max(1, std::min(std::max(1, requested), std::max(1, tasks)));
}

inline Rcpp::List blr_phase3_worker_diagnostics(
    const BlrPhase3ExecutionContract& contract,
    int requested_cores,
    int configured_workers,
    const std::vector<int>& worker_ids,
    const std::vector<int>& team_sizes) {
  if (!contract.active()) return R_NilValue;
  int actual_team_size = 1;
  for (int value : team_sizes) actual_team_size = std::max(actual_team_size, value);
#ifdef _OPENMP
  const bool openmp_available = true;
  const int runtime_maximum = omp_get_max_threads();
#else
  const bool openmp_available = false;
  const int runtime_maximum = 1;
#endif
  return Rcpp::List::create(
    Rcpp::Named("requested_cores") = requested_cores,
    Rcpp::Named("configured_workers") = configured_workers,
    Rcpp::Named("actual_team_size") = actual_team_size,
    Rcpp::Named("task_worker_ids") = worker_ids,
    Rcpp::Named("scheduler_version") = contract.scheduler_version,
    Rcpp::Named("logical_task_order") = contract.task_ids,
    Rcpp::Named("openmp_available") = openmp_available,
    Rcpp::Named("runtime_maximum_workers") = runtime_maximum,
    Rcpp::Named("diagnostics_rng_draws") = 0
  );
}

#endif
