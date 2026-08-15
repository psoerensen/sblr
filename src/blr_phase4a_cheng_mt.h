#ifndef SBLR_BLR_PHASE4A_CHENG_MT_H
#define SBLR_BLR_PHASE4A_CHENG_MT_H

#include "blr_mt_bed_access.h"
#include "blr_mt_covariance_rng.h"

#include <armadillo>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr {
namespace phase4a {

constexpr std::size_t trait_count = 2u;
constexpr std::size_t pattern_count = 4u;

inline const std::array<std::array<int, trait_count>, pattern_count>&
activity_patterns() {
  static const std::array<std::array<int, trait_count>, pattern_count> value{{
    {{0, 0}}, {{1, 0}}, {{0, 1}}, {{1, 1}}
  }};
  return value;
}

inline arma::mat require_spd(const arma::mat& value, const char* what) {
  if (value.n_rows != trait_count || value.n_cols != trait_count ||
      !value.is_finite()) {
    throw std::invalid_argument(std::string(what) +
      " must be a finite 2 x 2 matrix.");
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

inline double logdet_spd(const arma::mat& value) {
  arma::mat lower;
  if (!arma::chol(lower, value, "lower")) {
    throw std::runtime_error("A log determinant received a non-SPD matrix.");
  }
  return 2.0 * arma::accu(arma::log(lower.diag()));
}

struct PatternKernel {
  std::array<double, pattern_count> probability{{0.0, 0.0, 0.0, 0.0}};
  std::array<double, pattern_count> log_weight{{0.0, 0.0, 0.0, 0.0}};
  std::array<arma::vec, pattern_count> active_mean;
  std::array<arma::mat, pattern_count> active_covariance;
};

inline PatternKernel pattern_kernel(
    const arma::vec& score,
    double marker_sum_squares,
    const arma::mat& marker_covariance,
    const arma::mat& residual_precision,
    const std::array<double, pattern_count>& pattern_probability) {
  if (score.n_elem != trait_count || !score.is_finite() ||
      !std::isfinite(marker_sum_squares) || marker_sum_squares <= 0.0) {
    throw std::invalid_argument(
      "The Phase 4a marker score and sum of squares are invalid.");
  }
  PatternKernel out;
  const arma::vec full_score = residual_precision * score;
  const arma::mat marker_precision = inverse_spd(marker_covariance);

  for (std::size_t state = 0; state < pattern_count; ++state) {
    const double prior_mass = pattern_probability[state];
    if (!std::isfinite(prior_mass) || prior_mass <= 0.0) {
      throw std::invalid_argument(
        "Phase 4a activity-pattern probabilities must be finite and positive.");
    }
    out.log_weight[state] = std::log(prior_mass);
    if (state == 0u) continue;

    if (state == 1u || state == 2u) {
      const arma::uword active = state == 1u ? 0u : 1u;
      const double prior_variance = marker_covariance(active, active);
      const double precision = 1.0 / prior_variance +
        marker_sum_squares * residual_precision(active, active);
      if (!std::isfinite(precision) || precision <= 0.0) {
        throw std::runtime_error("A Phase 4a scalar precision is not positive.");
      }
      const double covariance = 1.0 / precision;
      const double mean = covariance * full_score(active);
      out.active_mean[state] = arma::vec(1u);
      out.active_mean[state](0u) = mean;
      out.active_covariance[state] = arma::mat(1u, 1u);
      out.active_covariance[state](0u, 0u) = covariance;
      out.log_weight[state] += -0.5 * std::log(prior_variance) -
        0.5 * std::log(precision) + 0.5 * full_score(active) * mean;
      continue;
    }

    const arma::mat precision = marker_precision +
      marker_sum_squares * residual_precision;
    const arma::vec mean = solve_spd(precision, full_score);
    out.active_mean[state] = mean;
    out.active_covariance[state] = inverse_spd(precision);
    out.log_weight[state] += -0.5 * logdet_spd(marker_covariance) -
      0.5 * logdet_spd(precision) + 0.5 * arma::dot(full_score, mean);
  }

  double maximum = out.log_weight[0u];
  for (std::size_t state = 1u; state < pattern_count; ++state) {
    maximum = std::max(maximum, out.log_weight[state]);
  }
  double total = 0.0;
  for (std::size_t state = 0u; state < pattern_count; ++state) {
    out.probability[state] = std::exp(out.log_weight[state] - maximum);
    total += out.probability[state];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error("Phase 4a pattern weights did not normalize.");
  }
  for (double& value : out.probability) value /= total;
  return out;
}

inline std::size_t draw_pattern(
    const std::array<double, pattern_count>& probability,
    std::mt19937& rng) {
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  const double draw = uniform(rng);
  double cumulative = 0.0;
  for (std::size_t state = 0u; state < pattern_count; ++state) {
    cumulative += probability[state];
    if (draw < cumulative || state + 1u == pattern_count) return state;
  }
  return pattern_count - 1u;
}

inline arma::vec draw_active(
    const arma::vec& mean,
    const arma::mat& covariance,
    std::mt19937& rng) {
  arma::mat lower;
  if (!arma::chol(lower, covariance, "lower")) {
    throw std::runtime_error(
      "A Phase 4a active-effect covariance is not positive definite.");
  }
  std::normal_distribution<double> normal(0.0, 1.0);
  arma::vec z(mean.n_elem);
  for (arma::uword index = 0u; index < mean.n_elem; ++index) {
    z(index) = normal(rng);
  }
  return mean + lower * z;
}

inline arma::rowvec complete_latent(
    std::size_t state,
    const arma::vec& active_value,
    const arma::mat& marker_covariance,
    std::mt19937& rng) {
  arma::rowvec completed(trait_count);
  completed.fill(std::numeric_limits<double>::quiet_NaN());
  if (state == 0u) return completed;
  if (state == 3u) {
    completed(0u) = active_value(0u);
    completed(1u) = active_value(1u);
    return completed;
  }
  const arma::uword active = state == 1u ? 0u : 1u;
  const arma::uword inactive = 1u - active;
  completed(active) = active_value(0u);
  const double coefficient = marker_covariance(inactive, active) /
    marker_covariance(active, active);
  const double conditional_variance = marker_covariance(inactive, inactive) -
    marker_covariance(inactive, active) * coefficient;
  if (!std::isfinite(conditional_variance) || conditional_variance <= 0.0) {
    throw std::runtime_error(
      "A Phase 4a conditional-completion variance is not positive.");
  }
  std::normal_distribution<double> normal(0.0, 1.0);
  completed(inactive) = coefficient * active_value(0u) +
    std::sqrt(conditional_variance) * normal(rng);
  return completed;
}

inline std::array<double, pattern_count> draw_dirichlet(
    const std::array<double, pattern_count>& shape,
    std::mt19937& rng) {
  std::array<double, pattern_count> out{{0.0, 0.0, 0.0, 0.0}};
  double total = 0.0;
  for (std::size_t state = 0u; state < pattern_count; ++state) {
    if (!std::isfinite(shape[state]) || shape[state] <= 0.0) {
      throw std::invalid_argument("Dirichlet shapes must be finite and positive.");
    }
    std::gamma_distribution<double> gamma(shape[state], 1.0);
    out[state] = gamma(rng);
    total += out[state];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error("A Phase 4a Dirichlet draw was invalid.");
  }
  for (double& value : out) value /= total;
  return out;
}

struct ChainResult {
  std::vector<arma::mat> realised_draws;
  std::vector<arma::mat> latent_draws;
  std::vector<std::vector<int>> state_draws;
  std::vector<arma::mat> covariance_draws;
  std::vector<std::array<double, pattern_count>> probability_draws;
  std::vector<arma::mat> prediction_draws;
  std::vector<arma::mat> convergence_covariance;
  std::vector<std::array<double, pattern_count>> convergence_probability;
  std::vector<int> convergence_active_count;
  arma::mat final_realised;
  arma::mat final_latent;
  std::vector<int> final_state;
  arma::mat final_covariance;
  std::array<double, pattern_count> final_probability{{0.0, 0.0, 0.0, 0.0}};
  arma::mat final_prediction;
  arma::mat last_covariance_statistic;
  arma::mat last_covariance_scale;
  double last_covariance_df = 0.0;
  int last_active_count = 0;
  std::array<std::array<std::uint64_t, pattern_count>, pattern_count>
    transition_count{};
};

template <class PackedGenotype>
ChainResult run_chain(
    const sblr::core::BedPackedGenotypeView<PackedGenotype>& genotype,
    const std::vector<sblr::mt::MtBedMarkerMap>& marker_maps,
    const arma::mat& phenotype,
    const arma::mat& fixed_residual_covariance,
    const arma::mat& initial_marker_covariance,
    const std::array<double, pattern_count>& initial_probability,
    const std::array<double, pattern_count>& dirichlet_prior,
    double prior_df,
    const arma::mat& prior_scale,
    bool update_covariance,
    bool update_probability,
    int burn_in_iterations,
    int sampling_iterations,
    const std::vector<int>& retained_transition_indices,
    std::uint32_t seed) {
  if (genotype.marker_count != marker_maps.size() ||
      phenotype.n_rows != genotype.sample_count ||
      phenotype.n_cols != trait_count || !phenotype.is_finite()) {
    throw std::invalid_argument("Phase 4a BED data dimensions are inconsistent.");
  }
  if (!std::isfinite(prior_df) || prior_df <= 1.0) {
    throw std::invalid_argument(
      "Phase 4a inverse-Wishart degrees of freedom must exceed T - 1.");
  }
  const arma::mat residual_covariance =
    require_spd(fixed_residual_covariance, "fixed residual covariance");
  arma::mat marker_covariance =
    require_spd(initial_marker_covariance, "initial marker covariance");
  const arma::mat covariance_prior_scale =
    require_spd(prior_scale, "marker-covariance prior scale");
  const arma::mat residual_precision = inverse_spd(residual_covariance);

  const std::size_t marker_count = genotype.marker_count;
  arma::mat realised(marker_count, trait_count, arma::fill::zeros);
  arma::mat latent(marker_count, trait_count);
  latent.fill(std::numeric_limits<double>::quiet_NaN());
  std::vector<int> state(marker_count, 0);
  arma::mat residual = phenotype;
  arma::vec marker_workspace(genotype.sample_count, arma::fill::zeros);
  std::array<double, pattern_count> probability = initial_probability;
  std::mt19937 rng(seed);

  std::vector<int> retained_position(
    static_cast<std::size_t>(sampling_iterations) + 1u, -1);
  for (std::size_t index = 0u; index < retained_transition_indices.size(); ++index) {
    const int iteration = retained_transition_indices[index];
    if (iteration < 1 || iteration > sampling_iterations) {
      throw std::invalid_argument("A retained Phase 4a iteration is out of range.");
    }
    retained_position[static_cast<std::size_t>(iteration)] =
      static_cast<int>(index);
  }

  ChainResult out;
  out.realised_draws.resize(retained_transition_indices.size());
  out.latent_draws.resize(retained_transition_indices.size());
  out.state_draws.resize(retained_transition_indices.size());
  out.covariance_draws.resize(retained_transition_indices.size());
  out.probability_draws.resize(retained_transition_indices.size());
  out.prediction_draws.resize(retained_transition_indices.size());
  out.convergence_covariance.reserve(static_cast<std::size_t>(sampling_iterations));
  out.convergence_probability.reserve(static_cast<std::size_t>(sampling_iterations));
  out.convergence_active_count.reserve(static_cast<std::size_t>(sampling_iterations));

  const int total_iterations = burn_in_iterations + sampling_iterations;
  for (int iteration = 0; iteration < total_iterations; ++iteration) {
    for (std::size_t marker = 0u; marker < marker_count; ++marker) {
      sblr::mt::decode_mt_bed_marker(
        genotype, marker_maps[marker], marker, marker_workspace);
      const int previous_state = state[marker];
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        residual(static_cast<arma::uword>(sample), 0u) +=
          x * realised(static_cast<arma::uword>(marker), 0u);
        residual(static_cast<arma::uword>(sample), 1u) +=
          x * realised(static_cast<arma::uword>(marker), 1u);
      }
      arma::vec score(trait_count, arma::fill::zeros);
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        score(0u) += x * residual(static_cast<arma::uword>(sample), 0u);
        score(1u) += x * residual(static_cast<arma::uword>(sample), 1u);
      }
      const PatternKernel conditional = pattern_kernel(
        score, marker_maps[marker].xx, marker_covariance,
        residual_precision, probability);
      const std::size_t new_state = draw_pattern(conditional.probability, rng);
      state[marker] = static_cast<int>(new_state);
      ++out.transition_count[static_cast<std::size_t>(previous_state)][new_state];
      if (new_state == 0u) {
        realised.row(static_cast<arma::uword>(marker)).zeros();
        latent.row(static_cast<arma::uword>(marker)).fill(
          std::numeric_limits<double>::quiet_NaN());
      } else {
        const arma::vec active = draw_active(
          conditional.active_mean[new_state],
          conditional.active_covariance[new_state], rng);
        const arma::rowvec completed = complete_latent(
          new_state, active, marker_covariance, rng);
        latent.row(static_cast<arma::uword>(marker)) = completed;
        const auto& pattern = activity_patterns()[new_state];
        realised(static_cast<arma::uword>(marker), 0u) =
          pattern[0u] * completed(0u);
        realised(static_cast<arma::uword>(marker), 1u) =
          pattern[1u] * completed(1u);
      }
      for (std::size_t sample = 0u; sample < genotype.sample_count; ++sample) {
        const double x = marker_workspace(static_cast<arma::uword>(sample));
        residual(static_cast<arma::uword>(sample), 0u) -=
          x * realised(static_cast<arma::uword>(marker), 0u);
        residual(static_cast<arma::uword>(sample), 1u) -=
          x * realised(static_cast<arma::uword>(marker), 1u);
      }
    }

    std::array<int, pattern_count> counts{{0, 0, 0, 0}};
    for (int value : state) ++counts[static_cast<std::size_t>(value)];
    if (update_probability) {
      std::array<double, pattern_count> shape;
      for (std::size_t index = 0u; index < pattern_count; ++index) {
        shape[index] = dirichlet_prior[index] + counts[index];
      }
      probability = draw_dirichlet(shape, rng);
    }

    arma::mat statistic(trait_count, trait_count, arma::fill::zeros);
    int active_count = 0;
    for (std::size_t marker = 0u; marker < marker_count; ++marker) {
      if (state[marker] == 0) continue;
      const arma::rowvec value = latent.row(static_cast<arma::uword>(marker));
      statistic += value.t() * value;
      ++active_count;
    }
    const arma::mat posterior_scale = covariance_prior_scale + statistic;
    const double posterior_df = prior_df + active_count;
    if (update_covariance) {
      marker_covariance = sblr::mt::draw_inverse_wishart(
        posterior_df, posterior_scale, rng);
      marker_covariance = require_spd(marker_covariance,
                                      "sampled marker covariance");
    }
    out.last_covariance_statistic = statistic;
    out.last_covariance_scale = posterior_scale;
    out.last_covariance_df = posterior_df;
    out.last_active_count = active_count;

    if (iteration >= burn_in_iterations) {
      const int post_burn = iteration - burn_in_iterations + 1;
      out.convergence_covariance.push_back(marker_covariance);
      out.convergence_probability.push_back(probability);
      out.convergence_active_count.push_back(active_count);
      const int retained = retained_position[static_cast<std::size_t>(post_burn)];
      if (retained >= 0) {
        const std::size_t index = static_cast<std::size_t>(retained);
        out.realised_draws[index] = realised;
        out.latent_draws[index] = latent;
        out.state_draws[index] = state;
        out.covariance_draws[index] = marker_covariance;
        out.probability_draws[index] = probability;
        out.prediction_draws[index] = phenotype - residual;
      }
    }
  }

  out.final_realised = realised;
  out.final_latent = latent;
  out.final_state = state;
  out.final_covariance = marker_covariance;
  out.final_probability = probability;
  out.final_prediction = phenotype - residual;
  return out;
}

}  // namespace phase4a
}  // namespace sblr

#endif
