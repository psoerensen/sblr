#ifndef SBLR_BLR_BLOCK_RESIDUAL_POLICY_H
#define SBLR_BLR_BLOCK_RESIDUAL_POLICY_H

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <RcppArmadillo.h>

namespace sblr { namespace core {

enum class BlockResidualPolicy {
  GlobalProjectedLegacy = 0,
  GctbBlock = 1,
  FixedBlock = 2
};

enum class BlockVeMode {
  FixVe = 0,
  SamVe = 1,
  MixVe = 2,
  AllMixVe = 3
};

struct BlockResidualControl {
  BlockResidualPolicy policy = BlockResidualPolicy::GlobalProjectedLegacy;
  BlockVeMode mode = BlockVeMode::FixVe;
  double resam_threshold = 1.1;
  double minimum_ve_ratio = 0.7;
  bool keep_history = false;

  inline bool uses_block_variance() const {
    return policy == BlockResidualPolicy::GctbBlock ||
      policy == BlockResidualPolicy::FixedBlock;
  }
};

inline BlockResidualPolicy parse_block_residual_policy(const std::string& value) {
  if (value == "global_projected_legacy")
    return BlockResidualPolicy::GlobalProjectedLegacy;
  if (value == "gctb_block") return BlockResidualPolicy::GctbBlock;
  if (value == "fixed_block") return BlockResidualPolicy::FixedBlock;
  throw std::invalid_argument("unknown block residual policy: " + value);
}

inline BlockVeMode parse_block_ve_mode(const std::string& value) {
  if (value == "fixVe") return BlockVeMode::FixVe;
  if (value == "samVe") return BlockVeMode::SamVe;
  if (value == "mixVe") return BlockVeMode::MixVe;
  if (value == "allMixVe") return BlockVeMode::AllMixVe;
  throw std::invalid_argument("unknown block Ve mode: " + value);
}

inline const char* block_residual_policy_name(BlockResidualPolicy value) {
  switch (value) {
    case BlockResidualPolicy::GctbBlock: return "gctb_block";
    case BlockResidualPolicy::FixedBlock: return "fixed_block";
    default: return "global_projected_legacy";
  }
}

inline const char* block_ve_mode_name(BlockVeMode value) {
  switch (value) {
    case BlockVeMode::SamVe: return "samVe";
    case BlockVeMode::MixVe: return "mixVe";
    case BlockVeMode::AllMixVe: return "allMixVe";
    default: return "fixVe";
  }
}

inline void validate_block_residual_control(const BlockResidualControl& control) {
  if (!std::isfinite(control.resam_threshold) || control.resam_threshold <= 0.0)
    throw std::invalid_argument("resam_thresh must be finite and positive.");
  if (!std::isfinite(control.minimum_ve_ratio) ||
      control.minimum_ve_ratio <= 0.0)
    throw std::invalid_argument("minimum_ve_ratio must be finite and positive.");
  if (control.policy == BlockResidualPolicy::FixedBlock &&
      control.mode != BlockVeMode::FixVe)
    throw std::invalid_argument("fixed_block requires block_ve_mode = 'fixVe'.");
}

inline bool block_ve_eligible(double genetic_variance, double effect_ss,
                              double threshold) {
  if (genetic_variance < 1e-8 || effect_ss < 1e-8) return false;
  return effect_ss / genetic_variance > threshold;
}

inline bool block_ve_ratio_accepted(double ratio, double minimum) {
  return std::isfinite(ratio) && ratio > minimum;
}

inline bool block_ve_selected(BlockVeMode mode, int iteration,
                              double genetic_variance, double effect_ss,
                              double threshold, unsigned char& permanent) {
  switch (mode) {
    case BlockVeMode::FixVe:
      return false;
    case BlockVeMode::SamVe:
      return true;
    case BlockVeMode::AllMixVe:
      return block_ve_eligible(genetic_variance, effect_ss, threshold);
    case BlockVeMode::MixVe:
      if (iteration < 50) return false;
      if (iteration == 50) {
        permanent = block_ve_eligible(
          genetic_variance, effect_ss, threshold) ? 1u : 0u;
      }
      return permanent != 0u;
  }
  return false;
}

struct BlockResidualChainState {
  arma::vec value;
  arma::vec posterior_sum;
  arma::uvec resampled;
  arma::uvec reset_to_phenotype;
  std::vector<unsigned char> mix_selected;
  arma::mat history;
  double posterior_draws = 0.0;
};

inline BlockResidualChainState make_block_residual_chain_state(
    int blocks, int trace_length, double phenotype_variance,
    bool keep_history) {
  if (blocks <= 0 || trace_length <= 0 ||
      !std::isfinite(phenotype_variance) || phenotype_variance <= 0.0)
    throw std::invalid_argument("invalid block residual state dimensions or phenotype variance.");
  BlockResidualChainState state;
  state.value.set_size(static_cast<arma::uword>(blocks));
  state.value.fill(phenotype_variance);
  state.posterior_sum.zeros(static_cast<arma::uword>(blocks));
  state.resampled.zeros(static_cast<arma::uword>(blocks));
  state.reset_to_phenotype.zeros(static_cast<arma::uword>(blocks));
  state.mix_selected.assign(static_cast<std::size_t>(blocks), 0u);
  if (keep_history) state.history.zeros(trace_length, blocks);
  return state;
}

template <class Operator>
inline double block_residual_variance_for_marker(
    const Operator& op, int marker, double legacy_variance,
    const BlockResidualControl& control,
    const BlockResidualChainState& state) {
  if (!control.uses_block_variance()) return legacy_variance;
  const int block = op.marker_block(marker);
  if (block < 0 || block >= static_cast<int>(state.value.n_elem))
    throw std::runtime_error("marker has no valid block residual variance.");
  return state.value(static_cast<arma::uword>(block));
}

template <class Operator>
inline void update_block_residual_variance(
    const Operator& op, int trait, int iteration, bool retain,
    double sample_size, double phenotype_variance, double residual_df,
    const arma::rowvec& effects, const arma::rowvec& residual,
    const BlockResidualControl& control, std::mt19937& generator,
    BlockResidualChainState& state) {
  if (!control.uses_block_variance()) return;
  const int blocks = op.block_count();
  if (blocks <= 0 || static_cast<int>(state.value.n_elem) != blocks)
    throw std::runtime_error("block residual state does not match the operator.");
  const double prior_scale = (residual_df - 2.0) / residual_df *
    phenotype_variance;
  if (!std::isfinite(prior_scale) || prior_scale <= 0.0)
    throw std::runtime_error("block residual prior scale is invalid.");

  for (int block = 0; block < blocks; ++block) {
    const arma::uword bu = static_cast<arma::uword>(block);
    const double genetic_variance = op.block_genetic_variance(
      trait, block, effects, residual, sample_size);
    const double effect_ss = op.block_effect_ss(block, effects);
    unsigned char& permanent = state.mix_selected[static_cast<std::size_t>(block)];
    const bool selected = control.policy == BlockResidualPolicy::GctbBlock &&
      block_ve_selected(control.mode, iteration, genetic_variance, effect_ss,
                        control.resam_threshold, permanent);
    if (selected) {
      const double residual_ss = op.block_residual_norm_squared(block, residual);
      const double scale = residual_ss + residual_df * prior_scale;
      const double degrees_freedom = residual_df + op.block_rank(block);
      if (!std::isfinite(scale) || scale <= 0.0 ||
          !std::isfinite(degrees_freedom) || degrees_freedom <= 0.0)
        throw std::runtime_error("official-compatible block residual draw is invalid.");
      std::chi_squared_distribution<double> rchisq(degrees_freedom);
      const double draw = scale / rchisq(generator);
      ++state.resampled(bu);
      if (draw > 0.0 && block_ve_ratio_accepted(
          draw / phenotype_variance, control.minimum_ve_ratio)) {
        state.value(bu) = draw;
      } else {
        state.value(bu) = phenotype_variance;
        ++state.reset_to_phenotype(bu);
      }
    } else {
      state.value(bu) = phenotype_variance;
    }
  }

  if (control.keep_history)
    state.history.row(static_cast<arma::uword>(iteration)) = state.value.t();
  if (retain) {
    state.posterior_sum += state.value;
    state.posterior_draws += 1.0;
  }
}

inline double mean_block_residual_variance(const BlockResidualChainState& state) {
  return arma::mean(state.value);
}

inline arma::vec posterior_mean_block_residual_variance(
    const BlockResidualChainState& state) {
  if (state.posterior_draws <= 0.0) return state.value;
  return state.posterior_sum / state.posterior_draws;
}

} }

#endif
