#ifndef ST_BAYESRC_PAIRWISE_ALLOCATION_H
#define ST_BAYESRC_PAIRWISE_ALLOCATION_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

struct StBayesRCPairState {
 int component_i = 0;
 int component_j = 0;
 double log_weight = -std::numeric_limits<double>::infinity();
 double probability = 0.0;
 double mean_i = 0.0;
 double mean_j = 0.0;
 double covariance_ii = 0.0;
 double covariance_ij = 0.0;
 double covariance_jj = 0.0;
};

inline std::vector<StBayesRCPairState> st_bayesrc_pairwise_conditional(
 const arma::rowvec& prior_i,
 const arma::rowvec& prior_j,
 const arma::vec& gamma,
 double marker_variance,
 double prior_scale_i,
 double prior_scale_j,
 double residual_variance,
 double diagonal_i,
 double diagonal_j,
 double cross_product,
 double score_i,
 double score_j
) {
 const int component_count = static_cast<int>(gamma.n_elem);
 if (component_count < 2 || prior_i.n_elem != gamma.n_elem ||
     prior_j.n_elem != gamma.n_elem || gamma(0) != 0.0 ||
     !prior_i.is_finite() || !prior_j.is_finite() ||
     arma::any(prior_i <= 0.0) || arma::any(prior_j <= 0.0) ||
     !std::isfinite(marker_variance) || marker_variance <= 0.0 ||
     !std::isfinite(prior_scale_i) || prior_scale_i <= 0.0 ||
     !std::isfinite(prior_scale_j) || prior_scale_j <= 0.0 ||
     !std::isfinite(residual_variance) || residual_variance <= 0.0 ||
     !std::isfinite(diagonal_i) || diagonal_i <= 0.0 ||
     !std::isfinite(diagonal_j) || diagonal_j <= 0.0 ||
     !std::isfinite(cross_product) || !std::isfinite(score_i) ||
     !std::isfinite(score_j)) {
  throw std::invalid_argument("BayesRC pairwise conditional inputs are invalid.");
 }
 const double determinant_g = diagonal_i * diagonal_j -
  cross_product * cross_product;
 const double determinant_tolerance = 1e-12 *
  std::max(1.0, diagonal_i * diagonal_j);
 if (determinant_g < -determinant_tolerance) {
  throw std::invalid_argument("BayesRC pair operator submatrix is not positive semidefinite.");
 }
 for (int component = 1; component < component_count; ++component) {
  if (!std::isfinite(gamma(static_cast<arma::uword>(component))) ||
      gamma(static_cast<arma::uword>(component)) <= 0.0) {
   throw std::invalid_argument("BayesRC active component scales must be positive finite.");
  }
 }

 const double inv_ve = 1.0 / residual_variance;
 const double h_i = score_i * inv_ve;
 const double h_j = score_j * inv_ve;
 std::vector<StBayesRCPairState> states;
 states.reserve(static_cast<std::size_t>(component_count * component_count));
 double maximum = -std::numeric_limits<double>::infinity();

 for (int ci = 0; ci < component_count; ++ci) {
  for (int cj = 0; cj < component_count; ++cj) {
   StBayesRCPairState state;
   state.component_i = ci;
   state.component_j = cj;
   state.log_weight = std::log(prior_i(static_cast<arma::uword>(ci))) +
    std::log(prior_j(static_cast<arma::uword>(cj)));
   const bool active_i = ci > 0;
   const bool active_j = cj > 0;
   if (active_i && !active_j) {
    const double variance_i = marker_variance * prior_scale_i *
     gamma(static_cast<arma::uword>(ci));
    const double precision = diagonal_i * inv_ve + 1.0 / variance_i;
    state.covariance_ii = 1.0 / precision;
    state.mean_i = state.covariance_ii * h_i;
    state.log_weight += -0.5 * (std::log(variance_i) +
     std::log(precision)) + 0.5 * h_i * state.mean_i;
   } else if (!active_i && active_j) {
    const double variance_j = marker_variance * prior_scale_j *
     gamma(static_cast<arma::uword>(cj));
    const double precision = diagonal_j * inv_ve + 1.0 / variance_j;
    state.covariance_jj = 1.0 / precision;
    state.mean_j = state.covariance_jj * h_j;
    state.log_weight += -0.5 * (std::log(variance_j) +
     std::log(precision)) + 0.5 * h_j * state.mean_j;
   } else if (active_i && active_j) {
    const double variance_i = marker_variance * prior_scale_i *
     gamma(static_cast<arma::uword>(ci));
    const double variance_j = marker_variance * prior_scale_j *
     gamma(static_cast<arma::uword>(cj));
    const double p11 = diagonal_i * inv_ve + 1.0 / variance_i;
    const double p22 = diagonal_j * inv_ve + 1.0 / variance_j;
    const double p12 = cross_product * inv_ve;
    const double determinant = p11 * p22 - p12 * p12;
    if (!std::isfinite(determinant) || determinant <= 0.0) {
     throw std::runtime_error("BayesRC pair posterior precision is not positive definite.");
    }
    state.covariance_ii = p22 / determinant;
    state.covariance_ij = -p12 / determinant;
    state.covariance_jj = p11 / determinant;
    state.mean_i = state.covariance_ii * h_i + state.covariance_ij * h_j;
    state.mean_j = state.covariance_ij * h_i + state.covariance_jj * h_j;
    state.log_weight += -0.5 * (std::log(variance_i) +
     std::log(variance_j) + std::log(determinant)) +
     0.5 * (h_i * state.mean_i + h_j * state.mean_j);
   }
   maximum = std::max(maximum, state.log_weight);
   states.push_back(state);
  }
 }
 double total = 0.0;
 for (StBayesRCPairState& state : states) {
  state.probability = std::exp(state.log_weight - maximum);
  total += state.probability;
 }
 if (!std::isfinite(total) || total <= 0.0)
  throw std::runtime_error("BayesRC pairwise probabilities cannot be normalized.");
 for (StBayesRCPairState& state : states) state.probability /= total;
 return states;
}

struct StBayesRCPairDraw {
 int component_i = 0;
 int component_j = 0;
 double effect_i = 0.0;
 double effect_j = 0.0;
};

inline StBayesRCPairDraw st_bayesrc_draw_pairwise_conditional(
 const std::vector<StBayesRCPairState>& states,
 std::mt19937& generator
) {
 if (states.empty())
  throw std::invalid_argument("BayesRC pairwise state list is empty.");
 std::uniform_real_distribution<double> uniform(0.0, 1.0);
 const double u = uniform(generator);
 double cumulative = 0.0;
 const StBayesRCPairState* selected = &states.back();
 for (const auto& state : states) {
  cumulative += state.probability;
  if (u <= cumulative) { selected = &state; break; }
 }
 std::normal_distribution<double> normal(0.0, 1.0);
 StBayesRCPairDraw draw;
 draw.component_i = selected->component_i;
 draw.component_j = selected->component_j;
 if (draw.component_i > 0 && draw.component_j > 0) {
  const double l11 = std::sqrt(selected->covariance_ii);
  const double l21 = selected->covariance_ij / l11;
  const double conditional_variance = selected->covariance_jj - l21 * l21;
  if (!std::isfinite(conditional_variance) || conditional_variance <= 0.0)
   throw std::runtime_error("BayesRC pair covariance is not positive definite.");
  const double z1 = normal(generator);
  const double z2 = normal(generator);
  draw.effect_i = selected->mean_i + l11 * z1;
  draw.effect_j = selected->mean_j + l21 * z1 +
   std::sqrt(conditional_variance) * z2;
 } else if (draw.component_i > 0) {
  draw.effect_i = selected->mean_i +
   std::sqrt(selected->covariance_ii) * normal(generator);
 } else if (draw.component_j > 0) {
  draw.effect_j = selected->mean_j +
   std::sqrt(selected->covariance_jj) * normal(generator);
 }
 return draw;
}

#endif
