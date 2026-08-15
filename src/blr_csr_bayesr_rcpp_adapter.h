#ifndef SBLR_BLR_CSR_BAYESR_RCPP_ADAPTER_H
#define SBLR_BLR_CSR_BAYESR_RCPP_ADAPTER_H

#include <RcppArmadillo.h>

#include "blr_csr_bayesr_policy.h"
#include "blr_phase3_execution.h"

Rcpp::List stblr_cpg_omp_csr_bayesr_with_policy(
 std::vector<std::vector<double>> wy,
 std::vector<std::vector<double>> ww,
 std::vector<double> yy,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<double>> comp_init,
 bool use_comp_init,
 std::vector<std::vector<double>> r_init,
 bool use_r_init,
 bool rebuild_r_before_updateE,
 std::string ld_prefix,
 arma::mat B,
 arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<double> pi,
 std::vector<double> mixture_var,
 std::vector<double> alpha,
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
 int ncores,
 int seed,
 int nchains,
 bool keep_chains,
 std::vector<int> chain_seeds,
 int updateE_start,
 int updateE_every,
 bool updateLDswap,
 double ld_swap_prob,
 double ld_swap_r2,
 int ld_swap_max_friends,
 int ld_swap_moves,
 Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_prior_scale,
 bool estimate_maf_effect_s,
 double maf_effect_s_init,
 Rcpp::NumericVector maf_effect_s_prior,
 double maf_effect_s_proposal_sd,
 Rcpp::Nullable<Rcpp::NumericVector> maf_effect_s_log_h,
 Rcpp::IntegerVector convergence_markers,
 bool convergence_probability,
 bool convergence_b,
 bool convergence_d,
 bool convergence_component,
 CsrBayesRPolicyFactory* policy_factory,
 const BlrPhase3ExecutionContract& execution_contract
);

#endif
