#ifndef SBLR_BLR_CSR_LOGVAR_BAYESR_CORE_IMPL_H
#define SBLR_BLR_CSR_LOGVAR_BAYESR_CORE_IMPL_H

#include "blr_csr_bayesr_policy.h"
#include "blr_csr_logvar_bayesr_types.h"

namespace sblr {
namespace logvar {

class CsrLogvarBayesRPolicy final : public CsrBayesRPolicy {
public:
 CsrLogvarBayesRPolicy(
  const CsrLogvarBayesRPolicyInput& input, int task, int trait,
  int marker_count)
  : input_(input), output_(input.outputs->at(static_cast<std::size_t>(task))) {
  output_.theta = input.theta_initial->col(static_cast<arma::uword>(trait));
  output_.prior_scale = calculate_prior_scale(
   *input.annotation, output_.theta, &output_.diagnostics).t();
  if (output_.prior_scale.n_elem != static_cast<arma::uword>(marker_count))
   throw std::invalid_argument("BayesR-LV prior scale marker count mismatch.");
  output_.theta_trace.zeros(
   static_cast<arma::uword>(input.trace_count), output_.theta.n_elem);
  output_.theta_sum.zeros(output_.theta.n_elem);
  output_.prior_scale_sum.zeros(output_.prior_scale.n_elem);
 }

 bool provides_prior_scale() const noexcept override { return true; }
 const arma::rowvec& prior_scale() const noexcept override {
  return output_.prior_scale;
 }

 void after_vb_update(
  const arma::rowvec& effects, const arma::Row<int>& component,
  double marker_variance, const arma::vec& mixture_var,
  std::mt19937& generator, int) override {
  if (!input_.update_theta) return;
  const arma::vec effect = effects.t();
  const arma::ivec component_vec =
   arma::conv_to<arma::ivec>::from(component.t());
  const bool empty = arma::accu(component_vec > 0) == 0;
  const auto likelihood = [&](const arma::vec& value) {
   return theta_log_likelihood_bayesr(
    value, *input_.annotation, effect, component_vec, marker_variance,
    mixture_var);
  };
  output_.theta = elliptical_slice_update(
   output_.theta, input_.theta_prior_sd, empty, likelihood, generator,
   output_.diagnostics);
  output_.prior_scale = calculate_prior_scale(
   *input_.annotation, output_.theta, &output_.diagnostics).t();
 }

 void capture(int iteration) override {
  output_.theta_trace.row(static_cast<arma::uword>(iteration)) =
   output_.theta.t();
 }
 void retain(int) override {
  output_.theta_sum += output_.theta;
  output_.prior_scale_sum += output_.prior_scale;
  output_.retained_samples += 1.0;
 }
 void finish() override {}

private:
 const CsrLogvarBayesRPolicyInput& input_;
 CsrLogvarBayesRChainOutput& output_;
};

class CsrLogvarBayesRPolicyFactory final : public CsrBayesRPolicyFactory {
public:
 explicit CsrLogvarBayesRPolicyFactory(
  const CsrLogvarBayesRPolicyInput& input) : input_(input) {}
 CsrBayesRPolicyHandle make(int task, int trait, int, int marker_count) override {
  return CsrBayesRPolicyHandle(std::unique_ptr<CsrBayesRPolicy>(
   new CsrLogvarBayesRPolicy(input_, task, trait, marker_count)));
 }
private:
 const CsrLogvarBayesRPolicyInput& input_;
};

}  // namespace logvar
}  // namespace sblr

#endif
