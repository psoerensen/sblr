#ifndef SBLR_BLR_SUMMARY_MT_CHENG_H
#define SBLR_BLR_SUMMARY_MT_CHENG_H

#include "blr_phase4a_cheng_mt.h"
#include "blr_sparse_ld_csr.h"

#include <armadillo>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sblr {
namespace summary_mt {

enum class OperatorKind { csr, block_eigen };

struct BlockFactor {
  std::vector<int> marker;
  arma::mat factor;
};

// One immutable provider-local cross-product operator. CSR uses the shared
// Phase 2 view contract. Block eigen retains U sqrt(Lambda) factors and never
// materializes the reconstructed block matrix.
struct OperatorResource {
  OperatorKind kind = OperatorKind::csr;
  std::size_t marker_count = 0u;
  sblr::core::SparseLdCsrStorage csr;
  arma::rowvec diagonal;
  std::vector<BlockFactor> blocks;
  std::vector<int> block_of;
  std::vector<int> local_in_block;

  inline sblr::core::SparseLdCsrView csr_view() const {
    sblr::core::SparseLdCsrView view;
    view.marker_count = marker_count;
    view.row_ptr = csr.ptr.data();
    view.row_ptr_size = csr.ptr.size();
    view.column_index = csr.idx.empty() ? nullptr : csr.idx.data();
    view.offdiag_xij = csr.xij.empty() ? nullptr : csr.xij.data();
    view.nonzero_count = csr.idx.size();
    view.diagonal = &diagonal;
    return view;
  }

  inline double diag(std::size_t marker) const {
    return diagonal(static_cast<arma::uword>(marker));
  }

  inline void apply_difference(
      std::size_t marker, double difference,
      std::vector<double>& residual) const {
    if (kind == OperatorKind::csr) {
      residual[marker] -= diag(marker) * difference;
      csr_view().apply_offdiag(static_cast<int>(marker), difference, residual);
      return;
    }
    const int group = block_of.at(marker);
    const int local = local_in_block.at(marker);
    const BlockFactor& block = blocks.at(static_cast<std::size_t>(group));
    const arma::rowvec selected = block.factor.row(
      static_cast<arma::uword>(local));
    for (std::size_t row = 0u; row < block.marker.size(); ++row) {
      const double value = arma::dot(
        block.factor.row(static_cast<arma::uword>(row)), selected);
      residual[static_cast<std::size_t>(block.marker[row])] -=
        value * difference;
    }
  }
};

inline void validate_operator_resource(const OperatorResource& resource) {
  if (resource.marker_count == 0u ||
      resource.diagonal.n_elem != resource.marker_count ||
      !resource.diagonal.is_finite() || arma::any(resource.diagonal <= 0.0)) {
    throw std::invalid_argument(
      "A summary operator requires a finite positive diagonal.");
  }
  if (resource.kind == OperatorKind::csr) {
    sblr::core::validate_sparse_ld_csr_view(resource.csr_view());
    return;
  }
  if (resource.blocks.empty() ||
      resource.block_of.size() != resource.marker_count ||
      resource.local_in_block.size() != resource.marker_count) {
    throw std::invalid_argument(
      "A block-eigen summary operator has incomplete block metadata.");
  }
  std::vector<int> seen(resource.marker_count, 0);
  for (std::size_t group = 0u; group < resource.blocks.size(); ++group) {
    const BlockFactor& block = resource.blocks[group];
    if (block.marker.empty() || block.factor.n_rows != block.marker.size() ||
        block.factor.n_cols == 0u || !block.factor.is_finite()) {
      throw std::invalid_argument(
        "A block-eigen summary factor has invalid dimensions.");
    }
    for (std::size_t local = 0u; local < block.marker.size(); ++local) {
      const int marker = block.marker[local];
      if (marker < 0 || static_cast<std::size_t>(marker) >= resource.marker_count ||
          seen[static_cast<std::size_t>(marker)]++ != 0 ||
          resource.block_of[static_cast<std::size_t>(marker)] !=
            static_cast<int>(group) ||
          resource.local_in_block[static_cast<std::size_t>(marker)] !=
            static_cast<int>(local)) {
        throw std::invalid_argument(
          "Block-eigen summary blocks must partition the local marker domain.");
      }
    }
  }
}

struct Provider {
  std::string id;
  std::size_t trait = 0u;
  std::size_t resource = 0u;
  std::vector<std::size_t> local_to_global;
  std::vector<double> score;
  double residual_scale = 1.0;
  double sample_size = 0.0;
};

inline void validate_provider(
    const Provider& provider,
    const std::vector<OperatorResource>& resources,
    std::size_t trait_count,
    std::size_t global_marker_count) {
  if (provider.id.empty() || provider.trait >= trait_count ||
      provider.resource >= resources.size() ||
      provider.local_to_global.size() != provider.score.size() ||
      provider.score.size() != resources[provider.resource].marker_count ||
      !std::isfinite(provider.residual_scale) ||
      provider.residual_scale <= 0.0 ||
      !std::isfinite(provider.sample_size) || provider.sample_size <= 0.0) {
    throw std::invalid_argument(
      "A summary likelihood provider has invalid dimensions or scale.");
  }
  std::vector<int> seen(global_marker_count, 0);
  for (std::size_t local = 0u; local < provider.score.size(); ++local) {
    const std::size_t global = provider.local_to_global[local];
    if (global >= global_marker_count || seen[global]++ != 0 ||
        !std::isfinite(provider.score[local])) {
      throw std::invalid_argument(
        "A summary provider map or score is invalid.");
    }
  }
}

struct ProviderMarkerReference {
  std::size_t provider = 0u;
  std::size_t local_marker = 0u;
};

inline sblr::phase4a::PatternKernel pattern_kernel(
    const arma::vec& score,
    const arma::vec& likelihood_diagonal,
    const arma::mat& marker_covariance,
    const sblr::phase4a::ProbabilityVector& pattern_probability,
    const sblr::phase4a::ActivityPatterns& patterns) {
  const std::size_t trait_count = score.n_elem;
  if (likelihood_diagonal.n_elem != trait_count ||
      marker_covariance.n_rows != trait_count ||
      marker_covariance.n_cols != trait_count ||
      pattern_probability.size() != patterns.size() ||
      !score.is_finite() || !likelihood_diagonal.is_finite() ||
      arma::any(likelihood_diagonal < 0.0)) {
    throw std::invalid_argument(
      "Summary Cheng marker conditional dimensions are invalid.");
  }
  sblr::phase4a::PatternKernel out(patterns.size());
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    const double prior_mass = pattern_probability[state];
    if (!std::isfinite(prior_mass) || prior_mass <= 0.0) {
      throw std::invalid_argument(
        "Summary Cheng pattern probabilities must be finite and positive.");
    }
    out.log_weight[state] = std::log(prior_mass);
    const arma::uvec active = sblr::phase4a::active_indices(patterns[state]);
    if (active.n_elem == 0u) continue;
    const arma::mat prior = marker_covariance.submat(active, active);
    arma::mat precision = sblr::phase4a::inverse_spd(prior);
    precision.diag() += likelihood_diagonal.elem(active);
    const arma::vec active_score = score.elem(active);
    const arma::vec mean = sblr::phase4a::solve_spd(precision, active_score);
    out.active_mean[state] = mean;
    out.active_covariance[state] = sblr::phase4a::inverse_spd(precision);
    out.log_weight[state] +=
      -0.5 * sblr::phase4a::logdet_spd(prior) -
      0.5 * sblr::phase4a::logdet_spd(precision) +
      0.5 * arma::dot(active_score, mean);
  }
  const double maximum = *std::max_element(
    out.log_weight.begin(), out.log_weight.end());
  double total = 0.0;
  for (std::size_t state = 0u; state < patterns.size(); ++state) {
    out.probability[state] = std::exp(out.log_weight[state] - maximum);
    total += out.probability[state];
  }
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::runtime_error(
      "Summary Cheng pattern weights did not normalize.");
  }
  for (double& value : out.probability) value /= total;
  return out;
}

struct ChainResult {
  std::vector<arma::mat> realised_draws;
  std::vector<arma::mat> latent_draws;
  std::vector<std::vector<int>> state_draws;
  std::vector<arma::mat> covariance_draws;
  std::vector<sblr::phase4a::ProbabilityVector> probability_draws;
  std::vector<arma::mat> convergence_covariance;
  std::vector<sblr::phase4a::ProbabilityVector> convergence_probability;
  std::vector<int> convergence_active_count;
  arma::mat final_realised;
  arma::mat final_latent;
  std::vector<int> final_state;
  arma::mat final_covariance;
  sblr::phase4a::ProbabilityVector final_probability;
  std::vector<std::vector<double>> final_provider_residual_score;
  arma::mat last_covariance_statistic;
  arma::mat last_covariance_scale;
  double last_covariance_df = 0.0;
  int last_active_count = 0;
  std::vector<std::uint64_t> pattern_occupancy_count;
  std::uint64_t pattern_change_count = 0u;
};

inline ChainResult run_chain(
    const std::vector<OperatorResource>& resources,
    const std::vector<Provider>& providers,
    const std::vector<std::vector<ProviderMarkerReference>>& marker_providers,
    std::size_t marker_count,
    std::size_t trait_count,
    const sblr::phase4a::ActivityPatterns& patterns,
    const arma::mat& initial_marker_covariance,
    const sblr::phase4a::ProbabilityVector& initial_probability,
    const sblr::phase4a::ProbabilityVector& dirichlet_prior,
    double prior_df,
    const arma::mat& prior_scale,
    bool update_covariance,
    bool update_probability,
    int burn_in_iterations,
    int sampling_iterations,
    const std::vector<int>& retained_transition_indices,
    std::uint32_t seed) {
  sblr::phase4a::validate_activity_patterns(patterns, trait_count);
  if (marker_count == 0u || marker_providers.size() != marker_count ||
      providers.empty() || initial_probability.size() != patterns.size() ||
      dirichlet_prior.size() != patterns.size() ||
      burn_in_iterations < 0 || sampling_iterations < 1 ||
      burn_in_iterations >
        std::numeric_limits<int>::max() - sampling_iterations ||
      !std::isfinite(prior_df) || prior_df <= trait_count - 1u) {
    throw std::invalid_argument(
      "Summary Cheng chain inputs are inconsistent.");
  }
  arma::mat marker_covariance = sblr::phase4a::require_spd(
    initial_marker_covariance, trait_count, "initial marker covariance");
  const arma::mat covariance_prior_scale = sblr::phase4a::require_spd(
    prior_scale, trait_count, "marker-covariance prior scale");
  arma::mat realised(marker_count, trait_count, arma::fill::zeros);
  arma::mat latent(marker_count, trait_count);
  latent.fill(std::numeric_limits<double>::quiet_NaN());
  std::vector<int> state(marker_count, 0);
  sblr::phase4a::ProbabilityVector probability = initial_probability;
  std::vector<std::vector<double>> residual(providers.size());
  for (std::size_t provider = 0u; provider < providers.size(); ++provider) {
    residual[provider] = providers[provider].score;
  }
  std::mt19937 rng(seed);

  std::vector<int> retained_position(
    static_cast<std::size_t>(sampling_iterations) + 1u, -1);
  for (std::size_t index = 0u; index < retained_transition_indices.size(); ++index) {
    const int iteration = retained_transition_indices[index];
    if (iteration < 1 || iteration > sampling_iterations) {
      throw std::invalid_argument(
        "A retained summary Cheng iteration is out of range.");
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
  out.convergence_covariance.reserve(static_cast<std::size_t>(sampling_iterations));
  out.convergence_probability.reserve(static_cast<std::size_t>(sampling_iterations));
  out.convergence_active_count.reserve(static_cast<std::size_t>(sampling_iterations));
  out.pattern_occupancy_count.assign(patterns.size(), 0u);

  const int total_iterations = burn_in_iterations + sampling_iterations;
  for (int iteration = 0; iteration < total_iterations; ++iteration) {
    for (std::size_t marker = 0u; marker < marker_count; ++marker) {
      const int previous_state = state[marker];
      arma::vec score(trait_count, arma::fill::zeros);
      arma::vec diagonal(trait_count, arma::fill::zeros);
      for (const ProviderMarkerReference& reference : marker_providers[marker]) {
        const Provider& provider = providers[reference.provider];
        const OperatorResource& resource = resources[provider.resource];
        const double old_effect = realised(
          static_cast<arma::uword>(marker),
          static_cast<arma::uword>(provider.trait));
        resource.apply_difference(
          reference.local_marker, -old_effect, residual[reference.provider]);
        score(static_cast<arma::uword>(provider.trait)) +=
          residual[reference.provider][reference.local_marker] /
          provider.residual_scale;
        diagonal(static_cast<arma::uword>(provider.trait)) +=
          resource.diag(reference.local_marker) / provider.residual_scale;
      }
      const sblr::phase4a::PatternKernel conditional = pattern_kernel(
        score, diagonal, marker_covariance, probability, patterns);
      const std::size_t new_state = sblr::phase4a::draw_pattern(
        conditional.probability, rng);
      state[marker] = static_cast<int>(new_state);
      ++out.pattern_occupancy_count[new_state];
      if (static_cast<std::size_t>(previous_state) != new_state) {
        ++out.pattern_change_count;
      }
      if (new_state == 0u) {
        realised.row(static_cast<arma::uword>(marker)).zeros();
        latent.row(static_cast<arma::uword>(marker)).fill(
          std::numeric_limits<double>::quiet_NaN());
      } else {
        const arma::vec active = sblr::phase4a::draw_active(
          conditional.active_mean[new_state],
          conditional.active_covariance[new_state], rng);
        const arma::rowvec completed = sblr::phase4a::complete_latent(
          patterns[new_state], active, marker_covariance, rng);
        latent.row(static_cast<arma::uword>(marker)) = completed;
        for (std::size_t trait = 0u; trait < trait_count; ++trait) {
          realised(static_cast<arma::uword>(marker),
                   static_cast<arma::uword>(trait)) =
            patterns[new_state][trait] * completed(trait);
        }
      }
      for (const ProviderMarkerReference& reference : marker_providers[marker]) {
        const Provider& provider = providers[reference.provider];
        const double new_effect = realised(
          static_cast<arma::uword>(marker),
          static_cast<arma::uword>(provider.trait));
        resources[provider.resource].apply_difference(
          reference.local_marker, new_effect, residual[reference.provider]);
      }
    }

    std::vector<int> counts(patterns.size(), 0);
    for (int value : state) ++counts[static_cast<std::size_t>(value)];
    if (update_probability) {
      sblr::phase4a::ProbabilityVector shape(patterns.size(), 0.0);
      for (std::size_t index = 0u; index < patterns.size(); ++index) {
        shape[index] = dirichlet_prior[index] + counts[index];
      }
      probability = sblr::phase4a::draw_dirichlet(shape, rng);
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
      marker_covariance = sblr::phase4a::require_spd(
        sblr::mt::draw_inverse_wishart(posterior_df, posterior_scale, rng),
        trait_count, "sampled marker covariance");
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
      }
    }
  }
  out.final_realised = realised;
  out.final_latent = latent;
  out.final_state = state;
  out.final_covariance = marker_covariance;
  out.final_probability = probability;
  out.final_provider_residual_score = std::move(residual);
  return out;
}

}  // namespace summary_mt
}  // namespace sblr

#endif
