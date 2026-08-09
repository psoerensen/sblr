#ifndef SBLR_BLR_CSR_BAYESR_POLICY_H
#define SBLR_BLR_CSR_BAYESR_POLICY_H

#include <RcppArmadillo.h>

#include <memory>
#include <random>
#include <stdexcept>
#include <utility>

class CsrBayesRPolicy {
public:
 virtual ~CsrBayesRPolicy() = default;
 virtual bool provides_prior_scale() const noexcept = 0;
 virtual const arma::rowvec& prior_scale() const = 0;
 virtual void after_vb_update(
  const arma::rowvec&, const arma::Row<int>&, double,
  const arma::vec&, std::mt19937&, int) = 0;
 virtual void capture(int) = 0;
 virtual void retain(int) = 0;
 virtual void finish() = 0;
};

class CsrBayesRPolicyHandle {
 std::unique_ptr<CsrBayesRPolicy> policy_;
public:
 explicit CsrBayesRPolicyHandle(std::unique_ptr<CsrBayesRPolicy> policy)
  : policy_(std::move(policy)) {
  if (!policy_) throw std::invalid_argument("BayesR policy factory returned null");
 }
 bool provides_prior_scale() const noexcept {
  return policy_->provides_prior_scale();
 }
 const arma::rowvec& prior_scale() const { return policy_->prior_scale(); }
 void after_vb_update(
  const arma::rowvec& b, const arma::Row<int>& component, double vb,
  const arma::vec& mixture_var, std::mt19937& gen, int iteration) {
  policy_->after_vb_update(b, component, vb, mixture_var, gen, iteration);
 }
 void capture(int iteration) { policy_->capture(iteration); }
 void retain(int iteration) { policy_->retain(iteration); }
 void finish() { policy_->finish(); }
};

class CsrBayesRPolicyFactory {
public:
 virtual ~CsrBayesRPolicyFactory() = default;
 virtual CsrBayesRPolicyHandle make(int task, int trait, int chain, int m) = 0;
};

#endif
