#include <Rcpp.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "blr_scalar_execution.h"

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
