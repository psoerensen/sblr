#ifndef SBLR_BLR_CSR_BAYESC_POLICY_H
#define SBLR_BLR_CSR_BAYESC_POLICY_H

#include <RcppArmadillo.h>

#include <memory>
#include <random>
#include <stdexcept>
#include <utility>

class CsrBayesCPolicy {
 public:
  virtual ~CsrBayesCPolicy() = default;
  virtual bool provides_prior_scale() const noexcept = 0;
  virtual const arma::rowvec& prior_scale() const = 0;
  virtual void after_vb_update(
    const arma::rowvec&, const arma::Row<int>&, double,
    std::mt19937&, int) = 0;
  virtual void capture(int) = 0;
  virtual void retain(int) = 0;
  virtual void finish() = 0;
};

class CsrBayesCPolicyHandle {
 public:
  explicit CsrBayesCPolicyHandle(std::unique_ptr<CsrBayesCPolicy> policy)
    : policy_(std::move(policy)) {
    if (!policy_) {
      throw std::invalid_argument("BayesC policy factory returned null");
    }
  }
  bool provides_prior_scale() const noexcept {
    return policy_->provides_prior_scale();
  }
  const arma::rowvec& prior_scale() const { return policy_->prior_scale(); }
  void after_vb_update(
    const arma::rowvec& b, const arma::Row<int>& d, double vb,
    std::mt19937& gen, int iteration
  ) {
    policy_->after_vb_update(b, d, vb, gen, iteration);
  }
  void capture(int iteration) { policy_->capture(iteration); }
  void retain(int iteration) { policy_->retain(iteration); }
  void finish() { policy_->finish(); }

 private:
  std::unique_ptr<CsrBayesCPolicy> policy_;
};

class CsrBayesCPolicyFactory {
 public:
  virtual ~CsrBayesCPolicyFactory() = default;
  virtual CsrBayesCPolicyHandle make(int task, int trait, int chain, int m) = 0;
};

class CsrBayesCNoOpErasedPolicy : public CsrBayesCPolicy {
 public:
  bool provides_prior_scale() const noexcept override { return false; }
  const arma::rowvec& prior_scale() const override {
    throw std::logic_error("ordinary BayesC has no policy-owned prior scale");
  }
  void after_vb_update(
    const arma::rowvec&, const arma::Row<int>&, double,
    std::mt19937&, int) noexcept override {}
  void capture(int) noexcept override {}
  void retain(int) noexcept override {}
  void finish() noexcept override {}
};

inline CsrBayesCPolicyHandle make_bayesc_noop_policy() {
  return CsrBayesCPolicyHandle(
    std::unique_ptr<CsrBayesCPolicy>(new CsrBayesCNoOpErasedPolicy())
  );
}

#endif
