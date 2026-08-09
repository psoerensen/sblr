// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "blr_block_eigen_rcpp_adapter.h"
#include "blr_csr_logvar_bayesc_core_impl.h"
#include "st_logvar_annotation_rcpp.h"

namespace {

class BlockLogvarBayesCPolicy final : public CsrBayesCPolicy {
 public:
  BlockLogvarBayesCPolicy(
    const sblr::logvar::CsrLogvarBayesCPolicyInput& input,
    int task, int trait, int marker_count
  ) : policy_(input, task, trait, marker_count) {}

  bool provides_prior_scale() const noexcept override {
    return policy_.provides_prior_scale();
  }
  const arma::rowvec& prior_scale() const override {
    return policy_.prior_scale();
  }
  void after_vb_update(
    const arma::rowvec& b, const arma::Row<int>& d, double vb,
    std::mt19937& gen, int iteration
  ) override {
    policy_.after_vb_update(b, d, vb, gen, iteration);
  }
  void capture(int iteration) override { policy_.capture(iteration); }
  void retain(int iteration) override { policy_.retain(iteration); }
  void finish() override { policy_.finish(); }

 private:
  sblr::logvar::CsrLogvarBayesCPolicy policy_;
};

class BlockLogvarBayesCPolicyFactory final : public CsrBayesCPolicyFactory {
 public:
  explicit BlockLogvarBayesCPolicyFactory(
    const sblr::logvar::CsrLogvarBayesCPolicyInput& input
  ) : input_(input) {}

  CsrBayesCPolicyHandle make(
    int task, int trait, int, int marker_count
  ) override {
    return CsrBayesCPolicyHandle(std::unique_ptr<CsrBayesCPolicy>(
      new BlockLogvarBayesCPolicy(input_, task, trait, marker_count)
    ));
  }

 private:
  const sblr::logvar::CsrLogvarBayesCPolicyInput& input_;
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_block_eigen_logvar_bayesc(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  double adjE,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  double pi_prior_a,
  double pi_prior_b,
  int ncores,
  int seed,
  int nchains,
  bool keep_chains,
  std::vector<int> chain_seeds,
  Rcpp::IntegerVector convergence_markers,
  bool convergence_b,
  bool convergence_d,
  Rcpp::CharacterVector bed_files,
  int n_bed,
  Rcpp::List cls,
  Rcpp::Nullable<Rcpp::IntegerVector> rows,
  Rcpp::NumericVector af,
  Rcpp::IntegerVector block_start,
  std::string eigen_filter,
  double eigen_tau,
  double eigen_eta,
  std::string representation,
  double eigen_prop,
  int low_rank_residual_rebuild_every,
  arma::mat annotation,
  arma::mat theta_init,
  double theta_prior_sd,
  bool updateTheta
) {
  const int nt = static_cast<int>(wy.size());
  if (nt <= 0) {
    throw std::invalid_argument("BayesC-LV block fit requires at least one trait.");
  }
  const int m = static_cast<int>(wy[0].size());
  if (m <= 0 || annotation.n_rows != static_cast<arma::uword>(m) ||
      annotation.n_cols == 0 || !annotation.is_finite()) {
    throw std::invalid_argument(
      "BayesC-LV block annotation must be a finite m by p matrix.");
  }
  const int p = static_cast<int>(annotation.n_cols);
  if (theta_init.n_rows != static_cast<arma::uword>(p) ||
      theta_init.n_cols != static_cast<arma::uword>(nt) ||
      !theta_init.is_finite()) {
    throw std::invalid_argument(
      "BayesC-LV block theta_init must be finite p by nt.");
  }
  if (!std::isfinite(theta_prior_sd) || theta_prior_sd <= 0.0) {
    throw std::invalid_argument(
      "BayesC-LV block theta_prior_sd must be positive finite.");
  }

  const int task_count = nt * nchains;
  std::vector<sblr::logvar::CsrLogvarBayesCChainOutput> outputs(
    static_cast<std::size_t>(task_count));
  sblr::logvar::CsrLogvarBayesCPolicyInput policy_input;
  policy_input.annotation = &annotation;
  policy_input.theta_initial = &theta_init;
  policy_input.theta_prior_sd = theta_prior_sd;
  policy_input.update_theta = updateTheta;
  policy_input.trace_count = nit + nburn;
  policy_input.burnin = nburn;
  policy_input.thinning = nthin;
  policy_input.outputs = &outputs;
  BlockLogvarBayesCPolicyFactory policy_factory(policy_input);

  Rcpp::List raw = stblr_cpg_omp_csr_block_eigen_with_policy(
    std::move(wy), std::move(ww), std::move(yy), std::move(b_init),
    std::move(d_init), use_d_init, std::move(r_init), use_r_init,
    rebuild_r_before_updateE, "", std::move(B), std::move(E),
    std::move(ssb_prior), std::move(sse_prior), std::move(pi), nub, nue,
    updateB, updateE, updatePi, adjE, std::move(n), nit, nburn, nthin,
    pi_prior_a, pi_prior_b, ncores, seed, nchains, keep_chains,
    std::move(chain_seeds), false, 0.05, 0.8, 50, 1, R_NilValue,
    false, 0.0, Rcpp::NumericVector::create(-3.0, 2.0), 0.35,
    R_NilValue, convergence_markers, convergence_b, convergence_d,
    bed_files, n_bed, cls, rows, af, block_start, std::move(eigen_filter),
    eigen_tau, eigen_eta, std::move(representation), eigen_prop,
    low_rank_residual_rebuild_every, &policy_factory
  );
  return decorate_logvar_raw(
    raw, outputs, m, nt, p, nit + nburn, nchains, theta_prior_sd,
    updateTheta, "sbayesc_logvar", "block_eigen_logvar_bayesc"
  );
}
