#ifndef ST_BAYESRC_INFORMATION_DIAGNOSTICS_H
#define ST_BAYESRC_INFORMATION_DIAGNOSTICS_H

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

// Development-only diagnostic state for Phase 4D.  These quantities are
// deterministic functions of probabilities and states already present in the
// Gibbs sampler.  They must never be used to alter a transition.
struct StBayesRCInformationDiagnosticState {
 bool enabled = false;
 arma::mat current_prior_probability;
 arma::mat current_rb_probability;
 arma::mat prior_probability_accum;
 arma::mat rb_probability_accum;
 arma::mat prior_component_count_trace;
 arma::mat rb_component_count_trace;
 arma::mat hard_component_count_trace;
 arma::mat hard_stick_trace;
 arma::mat soft_stick_trace;
 arma::mat hard_annotation_trace;
 arma::mat soft_annotation_trace;
 arma::mat information_gain_trace;
 double retained_count = 0.0;
};

inline void st_bayesrc_information_initialize(
 StBayesRCInformationDiagnosticState& state,
 bool enabled,
 int marker_count,
 int component_count,
 int annotation_count,
 int retained_iterations
) {
 state.enabled = enabled;
 if (!enabled) return;
 const int stick_count = component_count - 1;
 state.current_prior_probability = arma::mat(
  marker_count, component_count, arma::fill::zeros);
 state.current_rb_probability = arma::mat(
  marker_count, component_count, arma::fill::zeros);
 state.prior_probability_accum = arma::mat(
  marker_count, component_count, arma::fill::zeros);
 state.rb_probability_accum = arma::mat(
  marker_count, component_count, arma::fill::zeros);
 state.prior_component_count_trace = arma::mat(
  retained_iterations, component_count, arma::fill::zeros);
 state.rb_component_count_trace = arma::mat(
  retained_iterations, component_count, arma::fill::zeros);
 state.hard_component_count_trace = arma::mat(
  retained_iterations, component_count, arma::fill::zeros);
 state.hard_stick_trace = arma::mat(
  retained_iterations, 3 * stick_count, arma::fill::zeros);
 state.soft_stick_trace = arma::mat(
  retained_iterations, 3 * stick_count, arma::fill::zeros);
 state.hard_annotation_trace = arma::mat(
  retained_iterations, 3 * stick_count * annotation_count, arma::fill::zeros);
 state.soft_annotation_trace = arma::mat(
  retained_iterations, 3 * stick_count * annotation_count, arma::fill::zeros);
 // entropy mean/median/p90, KL mean/median/p90, TV mean/median/p90.
 state.information_gain_trace = arma::mat(
  retained_iterations, 9, arma::fill::zeros);
}

inline void st_bayesrc_information_capture_probability(
 const std::vector<double>& log_probability,
 const arma::rowvec& prior_probability,
 arma::uword marker,
 StBayesRCInformationDiagnosticState& state
) {
 if (!state.enabled) return;
 if (log_probability.size() != state.current_rb_probability.n_cols ||
     prior_probability.n_elem != state.current_prior_probability.n_cols ||
     marker >= state.current_rb_probability.n_rows) {
  throw std::invalid_argument("invalid SBayesRC information diagnostic dimensions");
 }
 double maximum = -std::numeric_limits<double>::infinity();
 for (double value : log_probability) maximum = std::max(maximum, value);
 if (!std::isfinite(maximum))
  throw std::runtime_error("SBayesRC diagnostic log probabilities are invalid");
 double total = 0.0;
 for (std::size_t component = 0; component < log_probability.size(); ++component) {
  const double weight = std::exp(log_probability[component] - maximum);
  state.current_rb_probability(marker, static_cast<arma::uword>(component)) = weight;
  total += weight;
 }
 if (!std::isfinite(total) || total <= 0.0)
  throw std::runtime_error("SBayesRC diagnostic probabilities cannot be normalized");
 state.current_rb_probability.row(marker) /= total;
 state.current_prior_probability.row(marker) = prior_probability;
}

inline double st_bayesrc_information_quantile(
 std::vector<double> value,
 double probability
) {
 if (value.empty()) return std::numeric_limits<double>::quiet_NaN();
 std::sort(value.begin(), value.end());
 const double position = probability * static_cast<double>(value.size() - 1u);
 const std::size_t lower = static_cast<std::size_t>(std::floor(position));
 const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
 const double fraction = position - static_cast<double>(lower);
 return value[lower] + fraction * (value[upper] - value[lower]);
}

inline void st_bayesrc_information_capture_iteration(
 const arma::Row<int>& component,
 const arma::mat& annotation,
 int retained_index,
 bool posterior_retain,
 StBayesRCInformationDiagnosticState& state
) {
 if (!state.enabled) return;
 const arma::uword row = static_cast<arma::uword>(retained_index);
 const arma::uword marker_count = component.n_elem;
 const arma::uword component_count = state.current_rb_probability.n_cols;
 const arma::uword stick_count = component_count - 1u;
 const arma::uword annotation_count = annotation.n_cols;
 if (row >= state.prior_component_count_trace.n_rows ||
     annotation.n_rows != marker_count ||
     state.current_rb_probability.n_rows != marker_count) {
  throw std::invalid_argument("invalid SBayesRC information iteration state");
 }

 state.prior_component_count_trace.row(row) =
  arma::sum(state.current_prior_probability, 0);
 state.rb_component_count_trace.row(row) =
  arma::sum(state.current_rb_probability, 0);
 for (arma::uword marker = 0; marker < marker_count; ++marker) {
  const int hard = component(marker);
  if (hard < 0 || hard >= static_cast<int>(component_count))
   throw std::runtime_error("invalid hard component in information diagnostic");
  state.hard_component_count_trace(row, static_cast<arma::uword>(hard)) += 1.0;
 }

 for (arma::uword stick = 0; stick < stick_count; ++stick) {
  arma::vec hard_eligible(marker_count, arma::fill::zeros);
  arma::vec hard_success(marker_count, arma::fill::zeros);
  arma::vec soft_eligible(marker_count, arma::fill::zeros);
  arma::vec soft_success(marker_count, arma::fill::zeros);
  for (arma::uword marker = 0; marker < marker_count; ++marker) {
   const int hard = component(marker);
   hard_eligible(marker) = hard >= static_cast<int>(stick) ? 1.0 : 0.0;
   hard_success(marker) = hard > static_cast<int>(stick) ? 1.0 : 0.0;
   soft_eligible(marker) = arma::accu(state.current_rb_probability.row(marker).cols(
    stick, component_count - 1u));
   soft_success(marker) = arma::accu(state.current_rb_probability.row(marker).cols(
    stick + 1u, component_count - 1u));
  }
  const double hard_e = arma::accu(hard_eligible);
  const double hard_s = arma::accu(hard_success);
  const double soft_e = arma::accu(soft_eligible);
  const double soft_s = arma::accu(soft_success);
  const arma::uword offset = 3u * stick;
  state.hard_stick_trace(row, offset) = hard_e;
  state.hard_stick_trace(row, offset + 1u) = hard_s;
  state.hard_stick_trace(row, offset + 2u) = hard_e > 0.0 ? hard_s / hard_e : 0.0;
  state.soft_stick_trace(row, offset) = soft_e;
  state.soft_stick_trace(row, offset + 1u) = soft_s;
  state.soft_stick_trace(row, offset + 2u) = soft_e > 0.0 ? soft_s / soft_e : 0.0;

  const arma::vec hard_q_weight = hard_eligible;
  const arma::vec soft_q_weight = soft_eligible;
  for (arma::uword annotation_index = 0;
       annotation_index < annotation_count; ++annotation_index) {
   const arma::vec values = annotation.col(annotation_index);
   const arma::uword base = (3u * stick) * annotation_count + annotation_index;
   state.hard_annotation_trace(row, base) = arma::dot(values, hard_eligible);
   state.hard_annotation_trace(row, base + annotation_count) =
    arma::dot(values, hard_success);
   state.hard_annotation_trace(row, base + 2u * annotation_count) =
    arma::dot(values % values, hard_q_weight);
   state.soft_annotation_trace(row, base) = arma::dot(values, soft_eligible);
   state.soft_annotation_trace(row, base + annotation_count) =
    arma::dot(values, soft_success);
   state.soft_annotation_trace(row, base + 2u * annotation_count) =
    arma::dot(values % values, soft_q_weight);
  }
 }

 std::vector<double> entropy(marker_count), divergence(marker_count), variation(marker_count);
 for (arma::uword marker = 0; marker < marker_count; ++marker) {
  double h = 0.0, kl = 0.0, tv = 0.0;
  for (arma::uword k = 0; k < component_count; ++k) {
   const double posterior = state.current_rb_probability(marker, k);
   const double prior = std::max(state.current_prior_probability(marker, k), 1e-300);
   if (posterior > 0.0) {
    h -= posterior * std::log(posterior);
    kl += posterior * std::log(posterior / prior);
   }
   tv += std::abs(posterior - prior);
  }
  entropy[marker] = h;
  divergence[marker] = std::max(0.0, kl);
  variation[marker] = 0.5 * tv;
 }
 const auto mean = [](const std::vector<double>& x) {
  double total = 0.0;
  for (double value : x) total += value;
  return total / static_cast<double>(x.size());
 };
 state.information_gain_trace(row, 0u) = mean(entropy);
 state.information_gain_trace(row, 1u) = st_bayesrc_information_quantile(entropy, 0.5);
 state.information_gain_trace(row, 2u) = st_bayesrc_information_quantile(entropy, 0.9);
 state.information_gain_trace(row, 3u) = mean(divergence);
 state.information_gain_trace(row, 4u) = st_bayesrc_information_quantile(divergence, 0.5);
 state.information_gain_trace(row, 5u) = st_bayesrc_information_quantile(divergence, 0.9);
 state.information_gain_trace(row, 6u) = mean(variation);
 state.information_gain_trace(row, 7u) = st_bayesrc_information_quantile(variation, 0.5);
 state.information_gain_trace(row, 8u) = st_bayesrc_information_quantile(variation, 0.9);

 if (posterior_retain) {
  state.prior_probability_accum += state.current_prior_probability;
  state.rb_probability_accum += state.current_rb_probability;
  state.retained_count += 1.0;
 }
}

#endif
