#ifndef SBLR_BLR_CSR_LOGVAR_BAYESC_CORE_IMPL_H
#define SBLR_BLR_CSR_LOGVAR_BAYESC_CORE_IMPL_H

#include "blr_csr_logvar_bayesc_types.h"
#include "blr_csr_bayesc_core_impl.h"

namespace sblr {
namespace logvar {

class CsrLogvarBayesCPolicy {
 public:
  CsrLogvarBayesCPolicy(
    const CsrLogvarBayesCPolicyInput& input,
    int task,
    int trait,
    int marker_count
  ) : input_(input), output_(input.outputs->at(static_cast<std::size_t>(task))) {
    output_.theta = input.theta_initial->col(static_cast<arma::uword>(trait));
    output_.prior_scale = calculate_prior_scale(
      *input.annotation, output_.theta, &output_.diagnostics).t();
    if (output_.prior_scale.n_elem != static_cast<arma::uword>(marker_count)) {
      throw std::invalid_argument("BayesC-LV prior scale marker count mismatch.");
    }
    output_.theta_trace.zeros(
      static_cast<arma::uword>(input.trace_count), output_.theta.n_elem);
    output_.theta_sum.zeros(output_.theta.n_elem);
    output_.prior_scale_sum.zeros(output_.prior_scale.n_elem);
  }

  bool provides_prior_scale() const noexcept { return true; }
  const arma::rowvec& prior_scale() const noexcept {
    return output_.prior_scale;
  }

  void after_vb_update(
    const arma::rowvec& effects,
    const arma::Row<int>& state,
    double marker_variance,
    std::mt19937& generator,
    int
  ) {
    if (!input_.update_theta) return;
    const arma::vec effect = effects.t();
    const arma::ivec active_state = arma::conv_to<arma::ivec>::from(state.t());
    const bool empty = arma::accu(active_state > 0) == 0;
    const auto likelihood = [&](const arma::vec& value) {
      return theta_log_likelihood_bayesc(
        value, *input_.annotation, effect, active_state, marker_variance);
    };
    output_.theta = elliptical_slice_update(
      output_.theta, input_.theta_prior_sd, empty, likelihood, generator,
      output_.diagnostics);
    output_.prior_scale = calculate_prior_scale(
      *input_.annotation, output_.theta, &output_.diagnostics).t();
  }

  void capture(int iteration) {
    output_.theta_trace.row(static_cast<arma::uword>(iteration)) =
      output_.theta.t();
  }

  void retain(int) {
    output_.theta_sum += output_.theta;
    output_.prior_scale_sum += output_.prior_scale;
    output_.retained_samples += 1.0;
  }

  void finish() noexcept {}

 private:
  const CsrLogvarBayesCPolicyInput& input_;
  CsrLogvarBayesCChainOutput& output_;
};

class CsrLogvarBayesCPolicyFactory {
 public:
  explicit CsrLogvarBayesCPolicyFactory(
    const CsrLogvarBayesCPolicyInput& input
  ) : input_(input) {}

  CsrLogvarBayesCPolicy make(int task, int trait, int, int marker_count) {
    return CsrLogvarBayesCPolicy(input_, task, trait, marker_count);
  }

 private:
  const CsrLogvarBayesCPolicyInput& input_;
};

}  // namespace logvar
}  // namespace sblr

#endif
