#ifndef SBLR_BLR_PHASE4A_CHENG_MT_H
#define SBLR_BLR_PHASE4A_CHENG_MT_H

#include "blr_mt_bed_access.h"
#include "blr_mt_covariance_rng.h"

#include <armadillo>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr {
namespace phase4a {

constexpr std::size_t maximum_trait_count = 12u;
using ActivityPattern = std::vector<int>;
using ActivityPatterns = std::vector<ActivityPattern>;
using ProbabilityVector = std::vector<double>;

inline void validate_activity_patterns(
    const ActivityPatterns& patterns, std::size_t trait_count) {
  if (trait_count < 2u || trait_count > maximum_trait_count) {
    throw std::invalid_argument(
      "Complete Cheng activity-pattern enumeration supports T in [2, 12].");
  }
  const std::size_t expected_count = std::size_t{1u} << trait_count;
  if (patterns.size() != expected_count) {
    throw std::invalid_argument(
      "The Cheng activity-pattern count must equal 2^T.");
  }
  for (std::size_t state = 0u; state < expected_count; ++state) {
    if (patterns[state].size() != trait_count) {
      throw std::invalid_argument(
        "Every Cheng activity pattern must contain one bit per trait.");
    }
    for (std::size_t trait = 0u; trait < trait_count; ++trait) {
      const int expected = static_cast<int>((state >> trait) & 1u);
      if (patterns[state][trait] != expected) {
        throw std::invalid_argument(
          "Cheng activity patterns must be complete and canonical with the first trait changing fastest.");
      }
    }
  }
}

inline std::size_t checked_size_product(
    std::initializer_list<std::size_t> factors, const char* component) {
  std::size_t value = 1u;
  for (const std::size_t factor : factors) {
    if (value != 0u && factor > std::numeric_limits<std::size_t>::max() / value) {
      throw std::overflow_error(std::string(
        "Cheng native allocation overflow before sampling in ") + component + ".");
    }
    value *= factor;
  }
  return value;
}

inline std::size_t checked_size_sum(
    std::initializer_list<std::size_t> values, const char* component) {
  std::size_t total = 0u;
  for (const std::size_t value : values) {
    if (value > std::numeric_limits<std::size_t>::max() - total) {
      throw std::overflow_error(std::string(
        "Cheng native allocation overflow before sampling in ") + component + ".");
    }
    total += value;
  }
  return total;
}

inline void validate_native_allocation_dimensions(
    std::size_t trait_count, std::size_t pattern_count,
    std::size_t marker_count, std::size_t observation_count,
    std::size_t chain_count, std::size_t retained_count,
    std::size_t convergence_count, std::size_t scale_count = 0u) {
  checked_size_product({pattern_count, trait_count}, "pattern metadata");
  checked_size_product({marker_count, trait_count}, "marker-trait state");
  checked_size_product({observation_count, trait_count},
    "observation-trait state");
  checked_size_product({chain_count, retained_count, marker_count, trait_count},
    "retained effect output");
  checked_size_product({chain_count, retained_count, observation_count,
    trait_count}, "retained prediction output");
  checked_size_product({chain_count, convergence_count, pattern_count},
    "convergence probability output");
  checked_size_product({chain_count, convergence_count, trait_count,
    trait_count}, "convergence covariance output");
  checked_size_product({chain_count, pattern_count},
    "compact pattern occupancy diagnostics");
  if (scale_count > 0u) {
    const std::size_t nonnull_states = checked_size_product(
      {pattern_count - 1u, scale_count}, "pattern-by-scale state workspace");
    if (nonnull_states == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error(
        "Cheng native allocation overflow before sampling in pattern-by-scale state workspace.");
    }
    checked_size_product({marker_count, scale_count},
      "marker-component posterior output");
    checked_size_product({chain_count, retained_count, marker_count},
      "retained component assignments");
    checked_size_product({chain_count, retained_count, scale_count},
      "retained scale probability output");
    checked_size_product({chain_count, convergence_count, scale_count},
      "convergence scale probability output");
    checked_size_product({chain_count, scale_count},
      "compact scale occupancy diagnostics");
  }
}

inline arma::mat require_spd(
    const arma::mat& value, std::size_t trait_count, const char* what) {
  if (value.n_rows != trait_count || value.n_cols != trait_count ||
      !value.is_finite()) {
    throw std::invalid_argument(std::string(what) +
      " must be a finite T x T matrix.");
  }
  const arma::mat symmetric = 0.5 * (value + value.t());
  arma::mat lower;
  if (arma::norm(value - value.t(), "inf") > 1e-12 ||
      !arma::chol(lower, symmetric, "lower")) {
    throw std::invalid_argument(std::string(what) +
      " must be symmetric positive definite.");
  }
  return symmetric;
}

inline arma::mat inverse_spd(const arma::mat& value) {
  arma::mat lower;
  if (!arma::chol(lower, value, "lower")) {
    throw std::runtime_error("An SPD solve received a non-SPD matrix.");
  }
  arma::mat intermediate;
  arma::mat out;
  const arma::mat identity = arma::eye(value.n_rows, value.n_cols);
  if (!arma::solve(intermediate, arma::trimatl(lower), identity) ||
      !arma::solve(out, arma::trimatu(lower.t()), intermediate)) {
    throw std::runtime_error("A Cholesky SPD solve failed.");
  }
  return 0.5 * (out + out.t());
}

inline arma::vec solve_spd(const arma::mat& value, const arma::vec& rhs) {
  arma::mat lower;
  if (!arma::chol(lower, value, "lower")) {
    throw std::runtime_error("An SPD solve received a non-SPD matrix.");
  }
  arma::vec intermediate;
  arma::vec out;
  if (!arma::solve(intermediate, arma::trimatl(lower), rhs) ||
      !arma::solve(out, arma::trimatu(lower.t()), intermediate)) {
    throw std::runtime_error("A Cholesky SPD solve failed.");
  }
  return out;
}

inline arma::mat solve_spd(const arma::mat& value, const arma::mat& rhs) {
  arma::mat lower;
  if (!arma::chol(lower, value, "lower")) {
    throw std::runtime_error("An SPD solve received a non-SPD matrix.");
  }
  arma::mat intermediate;
  arma::mat out;
  if (!arma::solve(intermediate, arma::trimatl(lower), rhs) ||
      !arma::solve(out, arma::trimatu(lower.t()), intermediate)) {
    throw std::runtime_error("A Cholesky SPD solve failed.");
  }
  return out;
}

inline double logdet_spd(const arma::mat& value) {
  arma::mat lower;
  if (!arma::chol(lower, value, "lower")) {
    throw std::runtime_error("A log determinant received a non-SPD matrix.");
  }
  return 2.0 * arma::accu(arma::log(lower.diag()));
}

inline arma::uvec active_indices(const ActivityPattern& pattern) {
  arma::uvec out(static_cast<arma::uword>(std::count(
    pattern.begin(), pattern.end(), 1)));
  arma::uword position = 0u;
  for (std::size_t trait = 0u; trait < pattern.size(); ++trait) {
    if (pattern[trait] == 1) out(position++) = trait;
  }
  return out;
}

inline arma::uvec inactive_indices(const ActivityPattern& pattern) {
  arma::uvec out(static_cast<arma::uword>(std::count(
    pattern.begin(), pattern.end(), 0)));
  arma::uword position = 0u;
  for (std::size_t trait = 0u; trait < pattern.size(); ++trait) {
    if (pattern[trait] == 0) out(position++) = trait;
  }
  return out;
}

struct PatternKernel {
  ProbabilityVector probability;
  ProbabilityVector log_weight;
  std::vector<arma::vec> active_mean;
  std::vector<arma::mat> active_covariance;

  explicit PatternKernel(std::size_t pattern_count)
      : probability(pattern_count, 0.0),
        log_weight(pattern_count, 0.0),
        active_mean(pattern_count),
        active_covariance(pattern_count) {}
};

// Phase 7 factorized activity-pattern by positive-scale mixture. The native
// Markov state stores base effects beta ~ N(0, Vb). The likelihood-facing
// completed effect is theta = sqrt(gamma_k * q_j) * beta.
struct ScaleMixture {
  bool enabled = false;
  ProbabilityVector scales;
  ProbabilityVector probability;
  ProbabilityVector dirichlet_prior;
  ProbabilityVector marker_multiplier;
};

inline void validate_scale_mixture(
    const ScaleMixture& mixture, std::size_t marker_count) {
  if (!mixture.enabled) return;
  const std::size_t scale_count = mixture.scales.size();
  if (scale_count == 0u || mixture.probability.size() != scale_count ||
      mixture.dirichlet_prior.size() != scale_count ||
      mixture.marker_multiplier.size() != marker_count) {
    throw std::invalid_argument(
      "Pattern-by-scale mixture dimensions are inconsistent.");
  }
  double total = 0.0;
  for (std::size_t scale = 0u; scale < scale_count; ++scale) {
    if (!std::isfinite(mixture.scales[scale]) || mixture.scales[scale] <= 0.0 ||
        (scale > 0u && mixture.scales[scale] <= mixture.scales[scale - 1u]) ||
        !std::isfinite(mixture.probability[scale]) ||
        mixture.probability[scale] <= 0.0 ||
        !std::isfinite(mixture.dirichlet_prior[scale]) ||
        mixture.dirichlet_prior[scale] <= 0.0) {
      throw std::invalid_argument(
        "Positive scales must be finite and strictly ordered, with positive scale probabilities and Dirichlet shapes.");
    }
    total += mixture.probability[scale];
  }
  for (double value : mixture.marker_multiplier) {
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument(
        "Marker variance multipliers must be finite and positive.");
    }
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::invalid_argument("Scale probabilities do not normalize.");
  }
}

struct PatternScaleKernel {
  ProbabilityVector probability;
  ProbabilityVector log_weight;
  std::vector<std::size_t> pattern;
  std::vector<int> scale;
  std::vector<arma::vec> active_mean;
  std::vector<arma::mat> active_covariance;
};

constexpr std::size_t pattern_scale_armadillo_container_slot_bytes = 256u;
constexpr std::size_t pattern_scale_kernel_container_bytes = 256u;
static_assert(sizeof(arma::vec) <= pattern_scale_armadillo_container_slot_bytes,
  "The Phase 7 memory contract must cover an Armadillo vector container.");
static_assert(sizeof(arma::mat) <= pattern_scale_armadillo_container_slot_bytes,
  "The Phase 7 memory contract must cover an Armadillo matrix container.");
static_assert(sizeof(PatternScaleKernel) <= pattern_scale_kernel_container_bytes,
  "The Phase 7 memory contract must cover the pattern-scale kernel object.");

struct PatternScaleWorkspaceAllocation {
  std::size_t nonnull_states = 0u;
  std::size_t state_count = 0u;
  std::size_t state_table_bytes = 0u;
  std::size_t active_container_bytes = 0u;
  std::size_t active_numeric_bytes = 0u;
  std::size_t kernel_container_bytes = 0u;
  std::size_t total_bytes = 0u;
};

inline bool uses_pattern_scale_workspace(const ScaleMixture& mixture) {
  return mixture.enabled && !(mixture.scales.size() == 1u &&
    mixture.scales[0u] == 1.0 &&
    std::all_of(mixture.marker_multiplier.begin(),
                mixture.marker_multiplier.end(),
                [](double value) { return value == 1.0; }));
}

inline PatternScaleWorkspaceAllocation pattern_scale_workspace_allocation(
    const ActivityPatterns& patterns, std::size_t scale_count,
    std::size_t concurrent_chains) {
  if (patterns.empty() || scale_count == 0u || concurrent_chains == 0u) {
    throw std::invalid_argument(
      "Pattern-scale workspace dimensions must be positive before sampling.");
  }
  validate_activity_patterns(patterns, patterns.front().size());
  PatternScaleWorkspaceAllocation out;
  out.nonnull_states = checked_size_product(
    {patterns.size() - 1u, scale_count},
    "pattern-scale non-null state count");
  out.state_count = checked_size_sum(
    {out.nonnull_states, 1u}, "pattern-scale state count");
  std::size_t active_moment_elements = 0u;
  for (const ActivityPattern& pattern : patterns) {
    const std::size_t active = static_cast<std::size_t>(
      std::accumulate(pattern.begin(), pattern.end(), 0));
    if (active == 0u) continue;
    active_moment_elements = checked_size_sum(
      {active_moment_elements, active,
       checked_size_product({active, active},
         "pattern-scale active covariance elements")},
      "pattern-scale active moment elements");
  }
  const std::size_t concurrent_states = checked_size_product(
    {concurrent_chains, out.state_count},
    "pattern-scale concurrent state count");
  out.state_table_bytes = checked_size_product(
    {concurrent_states, 28u}, "pattern-scale state tables");
  out.active_container_bytes = checked_size_product(
    {concurrent_states, 2u, pattern_scale_armadillo_container_slot_bytes},
    "pattern-scale active containers");
  out.active_numeric_bytes = checked_size_product(
    {concurrent_chains, scale_count, active_moment_elements, sizeof(double)},
    "pattern-scale active numeric payloads");
  out.kernel_container_bytes = checked_size_product(
    {concurrent_chains, pattern_scale_kernel_container_bytes},
    "pattern-scale kernel containers");
  out.total_bytes = checked_size_sum(
    {out.state_table_bytes, out.active_container_bytes,
     out.active_numeric_bytes, out.kernel_container_bytes},
    "pattern-scale conditional workspace total");
  return out;
}

inline void guard_pattern_scale_workspace_allocation(
    const ActivityPatterns& patterns, std::size_t scale_count,
    std::size_t concurrent_chains, double native_memory_limit_bytes) {
  if (std::isnan(native_memory_limit_bytes) || native_memory_limit_bytes < 0.0 ||
      (std::isinf(native_memory_limit_bytes) && native_memory_limit_bytes < 0.0)) {
    throw std::invalid_argument(
      "Pattern-scale native memory limit must be nonnegative or positive Inf.");
  }
  const PatternScaleWorkspaceAllocation allocation =
    pattern_scale_workspace_allocation(patterns, scale_count, concurrent_chains);
  if (std::isfinite(native_memory_limit_bytes) &&
      static_cast<long double>(allocation.total_bytes) >
        static_cast<long double>(native_memory_limit_bytes)) {
    throw std::runtime_error(
      "Pattern-scale conditional workspace allocation guard failed before sampling: "
      "T=" + std::to_string(patterns.front().size()) +
      ", K=" + std::to_string(scale_count) +
      ", nonnull_states=" + std::to_string(allocation.nonnull_states) +
      ", state_count=" + std::to_string(allocation.state_count) +
      ", concurrent_chains=" + std::to_string(concurrent_chains) +
      ", pattern_scale_conditional_workspace=" +
        std::to_string(allocation.total_bytes) +
      ", native_memory_limit_bytes=" +
        std::to_string(native_memory_limit_bytes) + ".");
  }
}

inline PatternScaleKernel pattern_scale_kernel_information(
    const arma::vec& score,
    const arma::mat& likelihood_precision,
    const arma::mat& marker_covariance,
    const ProbabilityVector& pattern_probability,
    const ProbabilityVector& scale_probability,
    const ProbabilityVector& scales,
    double marker_multiplier,
    const ActivityPatterns& patterns) {
  const std::size_t trait_count = score.n_elem;
  validate_activity_patterns(patterns, trait_count);
  if (!score.is_finite() || !likelihood_precision.is_finite() ||
      likelihood_precision.n_rows != trait_count ||
      likelihood_precision.n_cols != trait_count ||
      marker_covariance.n_rows != trait_count ||
      marker_covariance.n_cols != trait_count ||
      pattern_probability.size() != patterns.size() || scales.empty() ||
      scale_probability.size() != scales.size() ||
      !std::isfinite(marker_multiplier) || marker_multiplier <= 0.0) {
    throw std::invalid_argument(
      "Pattern-by-scale marker information has inconsistent dimensions.");
  }
  const std::size_t nonnull_states = checked_size_product(
    {patterns.size() - 1u, scales.size()}, "pattern-by-scale state workspace");
  if (nonnull_states == std::numeric_limits<std::size_t>::max()) {
    throw std::overflow_error(
      "Cheng native allocation overflow before sampling in pattern-by-scale state workspace.");
  }
  const std::size_t state_count = nonnull_states + 1u;
  PatternScaleKernel out;
  out.probability.assign(state_count, 0.0);
  out.log_weight.assign(state_count, 0.0);
  out.pattern.assign(state_count, 0u);
  out.scale.assign(state_count, -1);
  out.active_mean.resize(state_count);
  out.active_covariance.resize(state_count);
  if (!std::isfinite(pattern_probability[0u]) ||
      pattern_probability[0u] <= 0.0) {
    throw std::invalid_argument("The null pattern probability must be positive.");
  }
  out.log_weight[0u] = std::log(pattern_probability[0u]);
  std::size_t state = 1u;
  for (std::size_t pattern = 1u; pattern < patterns.size(); ++pattern) {
    if (!std::isfinite(pattern_probability[pattern]) ||
        pattern_probability[pattern] <= 0.0) {
      throw std::invalid_argument("Activity-pattern probabilities must be positive.");
    }
    const arma::uvec active = active_indices(patterns[pattern]);
    const arma::mat prior = marker_covariance.submat(active, active);
    const arma::mat prior_precision = inverse_spd(prior);
    const arma::mat active_likelihood =
      likelihood_precision.submat(active, active);
    const arma::vec active_score = score.elem(active);
    for (std::size_t scale = 0u; scale < scales.size(); ++scale, ++state) {
      if (!std::isfinite(scales[scale]) || scales[scale] <= 0.0 ||
          !std::isfinite(scale_probability[scale]) ||
          scale_probability[scale] <= 0.0) {
        throw std::invalid_argument(
          "Scale values and probabilities must be finite and positive.");
      }
      const double multiplier = scales[scale] * marker_multiplier;
      const double root = std::sqrt(multiplier);
      const arma::mat precision = prior_precision +
        multiplier * active_likelihood;
      const arma::vec scaled_score = root * active_score;
      const arma::vec mean = solve_spd(precision, scaled_score);
      out.pattern[state] = pattern;
      out.scale[state] = static_cast<int>(scale);
      out.active_mean[state] = mean;
      out.active_covariance[state] = inverse_spd(precision);
      out.log_weight[state] = std::log(pattern_probability[pattern]) +
        std::log(scale_probability[scale]) - 0.5 * logdet_spd(prior) -
        0.5 * logdet_spd(precision) +
        0.5 * arma::dot(scaled_score, mean);
    }
  }
  const double maximum = *std::max_element(
    out.log_weight.begin(), out.log_weight.end());
  double total = 0.0;
  for (std::size_t index = 0u; index < state_count; ++index) {
    out.probability[index] = std::exp(out.log_weight[index] - maximum);
    total += out.probability[index];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error("Pattern-by-scale weights did not normalize.");
  }
  for (double& value : out.probability) value /= total;
  return out;
}

inline PatternKernel pattern_kernel(
    const arma::vec& score,
    double marker_sum_squares,
    const arma::mat& marker_covariance,
    const arma::mat& residual_precision,
    const ProbabilityVector& pattern_probability,
    const ActivityPatterns& patterns) {
  const std::size_t trait_count = score.n_elem;
  validate_activity_patterns(patterns, trait_count);
  if (!score.is_finite() || !std::isfinite(marker_sum_squares) ||
      marker_sum_squares <= 0.0 ||
      marker_covariance.n_rows != trait_count ||
      marker_covariance.n_cols != trait_count ||
      residual_precision.n_rows != trait_count ||
      residual_precision.n_cols != trait_count ||
      pattern_probability.size() != patterns.size()) {
    throw std::invalid_argument(
      "The Cheng marker score, covariance, or pattern dimensions are invalid.");
  }
  PatternKernel out(patterns.size());
  const arma::vec full_score = residual_precision * score;
  const arma::mat marker_precision = inverse_spd(marker_covariance);

  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    const double prior_mass = pattern_probability[state];
    if (!std::isfinite(prior_mass) || prior_mass <= 0.0) {
      throw std::invalid_argument(
        "Cheng activity-pattern probabilities must be finite and positive.");
    }
    out.log_weight[state] = std::log(prior_mass);
    const arma::uvec active = active_indices(patterns[state]);
    if (active.n_elem == 0u) continue;

    if (active.n_elem == 1u) {
      const arma::uword coordinate = active(0u);
      const double prior_variance = marker_covariance(coordinate, coordinate);
      const double precision = 1.0 / prior_variance +
        marker_sum_squares * residual_precision(coordinate, coordinate);
      if (!std::isfinite(precision) || precision <= 0.0) {
        throw std::runtime_error("A Cheng scalar precision is not positive.");
      }
      const double covariance = 1.0 / precision;
      const double mean = covariance * full_score(coordinate);
      out.active_mean[state] = arma::vec(1u);
      out.active_mean[state](0u) = mean;
      out.active_covariance[state] = arma::mat(1u, 1u);
      out.active_covariance[state](0u, 0u) = covariance;
      out.log_weight[state] += -0.5 * std::log(prior_variance) -
        0.5 * std::log(precision) + 0.5 * full_score(coordinate) * mean;
      continue;
    }

    if (active.n_elem == trait_count) {
      const arma::mat precision = marker_precision +
        marker_sum_squares * residual_precision;
      const arma::vec mean = solve_spd(precision, full_score);
      out.active_mean[state] = mean;
      out.active_covariance[state] = inverse_spd(precision);
      out.log_weight[state] += -0.5 * logdet_spd(marker_covariance) -
        0.5 * logdet_spd(precision) + 0.5 * arma::dot(full_score, mean);
      continue;
    }

    const arma::mat prior = marker_covariance.submat(active, active);
    const arma::mat precision = inverse_spd(prior) +
      marker_sum_squares * residual_precision.submat(active, active);
    const arma::vec active_score = full_score.elem(active);
    const arma::vec mean = solve_spd(precision, active_score);
    out.active_mean[state] = mean;
    out.active_covariance[state] = inverse_spd(precision);
    out.log_weight[state] += -0.5 * logdet_spd(prior) -
      0.5 * logdet_spd(precision) + 0.5 * arma::dot(active_score, mean);
  }

  const double maximum = *std::max_element(
    out.log_weight.begin(), out.log_weight.end());
  double total = 0.0;
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    out.probability[state] = std::exp(out.log_weight[state] - maximum);
    total += out.probability[state];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error("Cheng pattern weights did not normalize.");
  }
  for (double& value : out.probability) value /= total;
  return out;
}

inline std::size_t draw_pattern(
    const ProbabilityVector& probability, std::mt19937& rng) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double draw = uniform(rng);
  double cumulative = 0.0;
  for (std::size_t state = 0u; state < probability.size(); ++state) {
    cumulative += probability[state];
    if (draw < cumulative || state + 1u == probability.size()) return state;
  }
  return probability.size() - 1u;
}

inline arma::vec draw_active(
    const arma::vec& mean,
    const arma::mat& covariance,
    std::mt19937& rng) {
  arma::mat lower;
  if (!arma::chol(lower, covariance, "lower")) {
    throw std::runtime_error(
      "A Cheng active-effect covariance is not positive definite.");
  }
  std::normal_distribution<double> normal(0.0, 1.0);
  arma::vec z(mean.n_elem);
  for (arma::uword index = 0u; index < mean.n_elem; ++index) {
    z(index) = normal(rng);
  }
  return mean + lower * z;
}

inline arma::rowvec complete_latent(
    const ActivityPattern& pattern,
    const arma::vec& active_value,
    const arma::mat& marker_covariance,
    std::mt19937& rng) {
  const std::size_t trait_count = pattern.size();
  const arma::uvec active = active_indices(pattern);
  const arma::uvec inactive = inactive_indices(pattern);
  arma::rowvec completed(trait_count);
  completed.fill(std::numeric_limits<double>::quiet_NaN());
  if (active.n_elem == 0u) return completed;
  if (active_value.n_elem != active.n_elem) {
    throw std::invalid_argument(
      "A Cheng active draw does not match its activity pattern.");
  }
  for (arma::uword index = 0u; index < active.n_elem; ++index) {
    completed(active(index)) = active_value(index);
  }
  if (inactive.n_elem == 0u) return completed;

  std::normal_distribution<double> normal(0.0, 1.0);
  if (active.n_elem == 1u && inactive.n_elem == 1u) {
    const arma::uword a = active(0u);
    const arma::uword i = inactive(0u);
    const double coefficient = marker_covariance(i, a) /
      marker_covariance(a, a);
    const double conditional_variance = marker_covariance(i, i) -
      marker_covariance(i, a) * coefficient;
    if (!std::isfinite(conditional_variance) ||
        conditional_variance <= 0.0) {
      throw std::runtime_error(
        "A Cheng conditional-completion variance is not positive.");
    }
    completed(i) = coefficient * active_value(0u) +
      std::sqrt(conditional_variance) * normal(rng);
    return completed;
  }

  const arma::mat covariance_aa = marker_covariance.submat(active, active);
  const arma::mat covariance_ia = marker_covariance.submat(inactive, active);
  const arma::mat covariance_ai = marker_covariance.submat(active, inactive);
  const arma::mat covariance_ii = marker_covariance.submat(inactive, inactive);
  const arma::vec conditional_mean = covariance_ia *
    solve_spd(covariance_aa, active_value);
  arma::mat conditional_covariance = covariance_ii - covariance_ia *
    solve_spd(covariance_aa, covariance_ai);
  conditional_covariance = 0.5 *
    (conditional_covariance + conditional_covariance.t());
  arma::mat lower;
  if (!conditional_covariance.is_finite() ||
      !arma::chol(lower, conditional_covariance, "lower")) {
    throw std::runtime_error(
      "A Cheng conditional-completion covariance is not positive definite.");
  }
  arma::vec z(inactive.n_elem);
  for (arma::uword index = 0u; index < inactive.n_elem; ++index) {
    z(index) = normal(rng);
  }
  const arma::vec inactive_value = conditional_mean + lower * z;
  for (arma::uword index = 0u; index < inactive.n_elem; ++index) {
    completed(inactive(index)) = inactive_value(index);
  }
  return completed;
}

inline ProbabilityVector draw_dirichlet(
    const ProbabilityVector& shape, std::mt19937& rng) {
  ProbabilityVector out(shape.size(), 0.0);
  double total = 0.0;
  for (std::size_t state = 0u; state < shape.size(); ++state) {
    if (!std::isfinite(shape[state]) || shape[state] <= 0.0) {
      throw std::invalid_argument("Dirichlet shapes must be finite and positive.");
    }
    std::gamma_distribution<double> gamma(shape[state], 1.0);
    out[state] = gamma(rng);
    total += out[state];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error("A Cheng Dirichlet draw was invalid.");
  }
  for (double& value : out) value /= total;
  return out;
}

struct ChainResult {
  std::vector<arma::mat> realised_draws;
  std::vector<arma::mat> latent_draws;
  std::vector<std::vector<int>> state_draws;
  std::vector<arma::mat> covariance_draws;
  std::vector<arma::mat> residual_covariance_draws;
  std::vector<ProbabilityVector> probability_draws;
  std::vector<std::vector<int>> component_draws;
  std::vector<ProbabilityVector> scale_probability_draws;
  std::vector<arma::mat> prediction_draws;
  std::vector<arma::mat> convergence_covariance;
  std::vector<arma::mat> convergence_residual_covariance;
  std::vector<ProbabilityVector> convergence_probability;
  std::vector<ProbabilityVector> convergence_scale_probability;
  std::vector<int> convergence_active_count;
  arma::mat final_realised;
  arma::mat final_latent;
  std::vector<int> final_state;
  arma::mat final_covariance;
  arma::mat final_residual_covariance;
  ProbabilityVector final_probability;
  std::vector<int> final_component;
  ProbabilityVector final_scale_probability;
  arma::mat final_prediction;
  arma::mat last_covariance_statistic;
  arma::mat last_covariance_scale;
  double last_covariance_df = 0.0;
  int last_active_count = 0;
  arma::mat last_residual_statistic;
  arma::mat last_residual_scale;
  double last_residual_df = 0.0;
  int residual_covariance_update_count = 0;
  std::vector<std::uint64_t> pattern_occupancy_count;
  std::uint64_t pattern_change_count = 0u;
  std::vector<std::uint64_t> scale_occupancy_count;
  std::uint64_t scale_change_count = 0u;
};

template <class PackedGenotype>
ChainResult run_chain(
    const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype,
    const std::vector<sblr::mt::MtBedMarkerMap>& marker_maps,
    const arma::mat& phenotype,
    const ActivityPatterns& patterns,
    const arma::mat& initial_residual_covariance,
    const arma::mat& initial_marker_covariance,
    const ProbabilityVector& initial_probability,
    const ProbabilityVector& dirichlet_prior,
    double prior_df,
    const arma::mat& prior_scale,
    bool update_covariance,
    bool update_probability,
    bool update_residual_covariance,
    double residual_prior_df,
    const arma::mat& residual_prior_scale,
    int burn_in_iterations,
    int sampling_iterations,
    const std::vector<int>& retained_transition_indices,
    std::uint32_t seed,
    const ScaleMixture& scale_mixture = ScaleMixture()) {
  const std::size_t trait_count = phenotype.n_cols;
  validate_activity_patterns(patterns, trait_count);
  const std::size_t pattern_count = patterns.size();
  if (genotype.marker_count != marker_maps.size() ||
      phenotype.n_rows != genotype.sample_count || !phenotype.is_finite() ||
      initial_probability.size() != pattern_count ||
      dirichlet_prior.size() != pattern_count) {
    throw std::invalid_argument("Cheng MT BED data dimensions are inconsistent.");
  }
  if (burn_in_iterations < 0 || sampling_iterations < 1 ||
      burn_in_iterations >
        std::numeric_limits<int>::max() - sampling_iterations) {
    throw std::invalid_argument(
      "Cheng iteration counts overflow before sampling.");
  }
  if (!std::isfinite(prior_df) || prior_df <= trait_count - 1u) {
    throw std::invalid_argument(
      "Marker inverse-Wishart degrees of freedom must exceed T - 1.");
  }
  arma::mat residual_covariance = require_spd(
    initial_residual_covariance, trait_count, "initial residual covariance");
  arma::mat marker_covariance = require_spd(
    initial_marker_covariance, trait_count, "initial marker covariance");
  const arma::mat covariance_prior_scale = require_spd(
    prior_scale, trait_count, "marker-covariance prior scale");
  arma::mat residual_precision = inverse_spd(residual_covariance);
  arma::mat residual_covariance_prior_scale;
  if (update_residual_covariance) {
    if (!std::isfinite(residual_prior_df) ||
        residual_prior_df <= trait_count - 1u) {
      throw std::invalid_argument(
        "Residual inverse-Wishart degrees of freedom must exceed T - 1.");
    }
    residual_covariance_prior_scale = require_spd(
      residual_prior_scale, trait_count, "residual-covariance prior scale");
  }

  const std::size_t marker_count = genotype.marker_count;
  validate_scale_mixture(scale_mixture, marker_count);
  const bool scaled = scale_mixture.enabled;
  const bool cheng_reduction =
      scaled && scale_mixture.scales.size() == 1u &&
      scale_mixture.scales[0] == 1.0 &&
      std::all_of(scale_mixture.marker_multiplier.begin(),
                  scale_mixture.marker_multiplier.end(),
                  [](double value) { return value == 1.0; });
  arma::mat realised(marker_count, trait_count, arma::fill::zeros);
  arma::mat latent(marker_count, trait_count);
  latent.fill(std::numeric_limits<double>::quiet_NaN());
  arma::mat base;
  if (scaled) {
    base.set_size(marker_count, trait_count);
    base.fill(std::numeric_limits<double>::quiet_NaN());
  }
  std::vector<int> state(marker_count, 0);
  std::vector<int> component(marker_count, -1);
  arma::mat residual = phenotype;
  arma::vec marker_workspace(genotype.sample_count, arma::fill::zeros);
  ProbabilityVector probability = initial_probability;
  ProbabilityVector scale_probability = scale_mixture.probability;
  if (scaled) {
    const double total = std::accumulate(
      scale_probability.begin(), scale_probability.end(), 0.0);
    for (double& value : scale_probability) value /= total;
  }
  std::mt19937 rng(seed);

  std::vector<int> retained_position(
    static_cast<std::size_t>(sampling_iterations) + 1u, -1);
  for (std::size_t index = 0u; index < retained_transition_indices.size(); ++index) {
    const int iteration = retained_transition_indices[index];
    if (iteration < 1 || iteration > sampling_iterations) {
      throw std::invalid_argument("A retained Cheng iteration is out of range.");
    }
    retained_position[static_cast<std::size_t>(iteration)] =
      static_cast<int>(index);
  }

  ChainResult out;
  out.realised_draws.resize(retained_transition_indices.size());
  out.latent_draws.resize(retained_transition_indices.size());
  out.state_draws.resize(retained_transition_indices.size());
  out.covariance_draws.resize(retained_transition_indices.size());
  if (update_residual_covariance) {
    out.residual_covariance_draws.resize(retained_transition_indices.size());
  }
  out.probability_draws.resize(retained_transition_indices.size());
  if (scaled) {
    out.component_draws.resize(retained_transition_indices.size());
    out.scale_probability_draws.resize(retained_transition_indices.size());
  }
  out.prediction_draws.resize(retained_transition_indices.size());
  out.convergence_covariance.reserve(static_cast<std::size_t>(sampling_iterations));
  if (update_residual_covariance) {
    out.convergence_residual_covariance.reserve(
      static_cast<std::size_t>(sampling_iterations));
  }
  out.convergence_probability.reserve(static_cast<std::size_t>(sampling_iterations));
  if (scaled) out.convergence_scale_probability.reserve(
    static_cast<std::size_t>(sampling_iterations));
  out.convergence_active_count.reserve(static_cast<std::size_t>(sampling_iterations));
  out.pattern_occupancy_count.assign(pattern_count, 0u);
  if (scaled) out.scale_occupancy_count.assign(
    scale_mixture.scales.size(), 0u);

  const int total_iterations = burn_in_iterations + sampling_iterations;
  for (int iteration = 0; iteration < total_iterations; ++iteration) {
    for (std::size_t marker = 0u; marker < marker_count; ++marker) {
      sblr::mt::decode_mt_bed_marker(
        genotype, marker_maps[marker], marker, marker_workspace);
      const int previous_state = state[marker];
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        for (std::size_t trait = 0u; trait < trait_count; ++trait) {
          residual(static_cast<arma::uword>(sample), trait) +=
            x * realised(static_cast<arma::uword>(marker), trait);
        }
      }
      arma::vec score(trait_count, arma::fill::zeros);
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        for (std::size_t trait = 0u; trait < trait_count; ++trait) {
          score(trait) += x * residual(static_cast<arma::uword>(sample), trait);
        }
      }
      std::size_t new_state = 0u;
      int new_component = -1;
      arma::vec active;
      if (!scaled || cheng_reduction) {
        const PatternKernel conditional = pattern_kernel(
          score, marker_maps[marker].xx, marker_covariance,
          residual_precision, probability, patterns);
        new_state = draw_pattern(conditional.probability, rng);
        if (new_state != 0u) {
          active = draw_active(conditional.active_mean[new_state],
                               conditional.active_covariance[new_state], rng);
        }
      } else {
        const arma::vec information_score = residual_precision * score;
        const arma::mat information_precision =
          marker_maps[marker].xx * residual_precision;
        const PatternScaleKernel conditional =
          pattern_scale_kernel_information(
            information_score, information_precision, marker_covariance,
            probability, scale_probability, scale_mixture.scales,
            scale_mixture.marker_multiplier[marker], patterns);
        const std::size_t joint = draw_pattern(conditional.probability, rng);
        new_state = conditional.pattern[joint];
        new_component = conditional.scale[joint];
        if (new_state != 0u) {
          active = draw_active(conditional.active_mean[joint],
                               conditional.active_covariance[joint], rng);
        }
      }
      state[marker] = static_cast<int>(new_state);
      const int previous_component = component[marker];
      if (cheng_reduction && new_state != 0u) new_component = 0;
      component[marker] = new_component;
      ++out.pattern_occupancy_count[new_state];
      if (new_component >= 0) {
        ++out.scale_occupancy_count[static_cast<std::size_t>(new_component)];
      }
      if (static_cast<std::size_t>(previous_state) != new_state) {
        ++out.pattern_change_count;
      }
      if (previous_component != new_component) ++out.scale_change_count;
      if (new_state == 0u) {
        realised.row(static_cast<arma::uword>(marker)).zeros();
        latent.row(static_cast<arma::uword>(marker)).fill(
          std::numeric_limits<double>::quiet_NaN());
        if (scaled) {
          base.row(static_cast<arma::uword>(marker)).fill(
            std::numeric_limits<double>::quiet_NaN());
        }
      } else {
        const arma::rowvec completed = complete_latent(
          patterns[new_state], active, marker_covariance, rng);
        if (scaled) base.row(static_cast<arma::uword>(marker)) = completed;
        const double root = scaled ? std::sqrt(
          scale_mixture.scales[static_cast<std::size_t>(new_component)] *
          scale_mixture.marker_multiplier[marker]) : 1.0;
        latent.row(static_cast<arma::uword>(marker)) = root * completed;
        for (std::size_t trait = 0u; trait < trait_count; ++trait) {
          realised(static_cast<arma::uword>(marker), trait) =
            patterns[new_state][trait] * root * completed(trait);
        }
      }
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        for (std::size_t trait = 0u; trait < trait_count; ++trait) {
          residual(static_cast<arma::uword>(sample), trait) -=
            x * realised(static_cast<arma::uword>(marker), trait);
        }
      }
    }

    std::vector<int> counts(pattern_count, 0);
    for (int value : state) ++counts[static_cast<std::size_t>(value)];
    if (update_probability) {
      ProbabilityVector shape(pattern_count, 0.0);
      for (std::size_t index = 0u; index < pattern_count; ++index) {
        shape[index] = dirichlet_prior[index] + counts[index];
      }
      probability = draw_dirichlet(shape, rng);
    }
    if (scaled && scale_mixture.scales.size() > 1u) {
      ProbabilityVector shape(scale_mixture.scales.size(), 0.0);
      for (std::size_t index = 0u; index < shape.size(); ++index) {
        shape[index] = scale_mixture.dirichlet_prior[index];
      }
      for (int value : component) {
        if (value >= 0) ++shape[static_cast<std::size_t>(value)];
      }
      scale_probability = draw_dirichlet(shape, rng);
    }

    arma::mat statistic(trait_count, trait_count, arma::fill::zeros);
    int active_count = 0;
    for (std::size_t marker = 0u; marker < marker_count; ++marker) {
      if (state[marker] == 0) continue;
      const arma::rowvec value = scaled ?
        base.row(static_cast<arma::uword>(marker)) :
        latent.row(static_cast<arma::uword>(marker));
      statistic += value.t() * value;
      ++active_count;
    }
    const arma::mat posterior_scale = covariance_prior_scale + statistic;
    const double posterior_df = prior_df + active_count;
    if (update_covariance) {
      marker_covariance = sblr::mt::draw_inverse_wishart(
        posterior_df, posterior_scale, rng);
      marker_covariance = require_spd(
        marker_covariance, trait_count, "sampled marker covariance");
    }
    out.last_covariance_statistic = statistic;
    out.last_covariance_scale = posterior_scale;
    out.last_covariance_df = posterior_df;
    out.last_active_count = active_count;

    if (update_residual_covariance) {
      arma::mat residual_statistic = residual.t() * residual;
      residual_statistic = 0.5 * (residual_statistic + residual_statistic.t());
      const arma::mat residual_posterior_scale =
        residual_covariance_prior_scale + residual_statistic;
      const double residual_posterior_df = residual_prior_df +
        static_cast<double>(genotype.sample_count);
      residual_covariance = sblr::mt::draw_inverse_wishart(
        residual_posterior_df, residual_posterior_scale, rng);
      residual_covariance = require_spd(
        residual_covariance, trait_count, "sampled residual covariance");
      residual_precision = inverse_spd(residual_covariance);
      out.last_residual_statistic = residual_statistic;
      out.last_residual_scale = residual_posterior_scale;
      out.last_residual_df = residual_posterior_df;
      ++out.residual_covariance_update_count;
    }

    if (iteration >= burn_in_iterations) {
      const int post_burn = iteration - burn_in_iterations + 1;
      out.convergence_covariance.push_back(marker_covariance);
      if (update_residual_covariance) {
        out.convergence_residual_covariance.push_back(residual_covariance);
      }
      out.convergence_probability.push_back(probability);
      if (scaled) out.convergence_scale_probability.push_back(
        scale_probability);
      out.convergence_active_count.push_back(active_count);
      const int retained = retained_position[static_cast<std::size_t>(post_burn)];
      if (retained >= 0) {
        const std::size_t index = static_cast<std::size_t>(retained);
        out.realised_draws[index] = realised;
        out.latent_draws[index] = latent;
        out.state_draws[index] = state;
        out.covariance_draws[index] = marker_covariance;
        if (update_residual_covariance) {
          out.residual_covariance_draws[index] = residual_covariance;
        }
        out.probability_draws[index] = probability;
        if (scaled) {
          out.component_draws[index] = component;
          out.scale_probability_draws[index] = scale_probability;
        }
        out.prediction_draws[index] = phenotype - residual;
      }
    }
  }

  out.final_realised = realised;
  out.final_latent = latent;
  out.final_state = state;
  out.final_covariance = marker_covariance;
  out.final_residual_covariance = residual_covariance;
  out.final_probability = probability;
  out.final_component = component;
  out.final_scale_probability = scale_probability;
  out.final_prediction = phenotype - residual;
  return out;
}

}  // namespace phase4a
}  // namespace sblr

#endif
