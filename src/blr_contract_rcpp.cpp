#include <Rcpp.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "blr_scalar_execution.h"
#include "blr_phase3_execution.h"

// [[Rcpp::export]]
double blr_phase3_seed_v1_internal(
    double user_seed,
    std::string analysis_mode,
    Rcpp::Nullable<Rcpp::CharacterVector> trait_id,
    double chain_index) {
  if (!std::isfinite(user_seed) || user_seed < 0.0 ||
      user_seed > 4294967295.0 || user_seed != std::floor(user_seed) ||
      !std::isfinite(chain_index) || chain_index < 0.0 ||
      chain_index > 4294967295.0 || chain_index != std::floor(chain_index)) {
    throw std::runtime_error("seed inputs must be exact uint32 values.");
  }
  int mode_code = 0;
  std::string identity;
  if (analysis_mode == "single_trait") {
    mode_code = 1;
    identity = "sblr:single_trait";
  } else if (analysis_mode == "independent_traits") {
    mode_code = 2;
    if (trait_id.isNull()) {
      throw std::runtime_error("independent_traits requires a trait ID.");
    }
    Rcpp::CharacterVector trait_value(trait_id);
    if (trait_value.size() != 1 || trait_value[0] == NA_STRING) {
      throw std::runtime_error("independent_traits requires one trait ID.");
    }
    identity = Rcpp::as<std::string>(trait_value[0]);
    if (identity.empty()) {
      throw std::runtime_error("independent_traits requires a nonempty trait ID.");
    }
  } else if (analysis_mode == "joint_multitrait") {
    mode_code = 3;
    identity = "sblr:joint_multitrait";
  } else {
    throw std::runtime_error("unknown analysis mode.");
  }
  return static_cast<double>(blr_phase3_seed_v1(
    static_cast<std::uint32_t>(user_seed), mode_code, identity,
    static_cast<std::uint32_t>(chain_index)));
}

// Internal binding for exact logical-task seed resolution.
// [[Rcpp::export]]
Rcpp::IntegerVector blr_scalar_seeds_cpp(
  int seed,
  int ntraits,
  int nchains,
  Rcpp::IntegerVector chain_seeds
) {
  using namespace sblr::core;
  if (ntraits < 0 || nchains < 0) {
    throw std::invalid_argument("seed dimensions must be non-negative");
  }
  const std::vector<int> explicit_seeds =
    Rcpp::as<std::vector<int>>(chain_seeds);
  const std::vector<ScalarChainTask> tasks = make_scalar_chain_tasks(
    static_cast<std::size_t>(ntraits), static_cast<std::size_t>(nchains)
  );
  validate_scalar_chain_seeds(static_cast<std::size_t>(nchains), explicit_seeds);
  Rcpp::IntegerVector resolved(tasks.size());
  for (std::size_t index = 0; index < tasks.size(); ++index) {
    resolved[index] = static_cast<int>(resolve_scalar_chain_seed(
      seed, static_cast<std::size_t>(nchains), explicit_seeds, tasks[index]
    ));
  }
  return resolved;
}
