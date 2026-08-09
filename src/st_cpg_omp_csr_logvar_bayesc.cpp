// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "st_csr_common.h"
#include "st_ld_operator.h"
#include "st_chain_utils.h"
#include "blr_csr_bayesc_operator_adapter.h"
#include "blr_csr_bayesc_rcpp_adapter.h"
#include "blr_csr_logvar_bayesc_core_impl.h"

// [[Rcpp::export]]
Rcpp::List stblr_cpg_omp_csr_logvar_bayesc(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  std::vector<std::vector<double>> d_init,
  bool use_d_init,
  std::vector<std::vector<double>> r_init,
  bool use_r_init,
  bool rebuild_r_before_updateE,
  std::string ld_prefix,
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
  bool updateLDswap,
  double ld_swap_prob,
  double ld_swap_r2,
  int ld_swap_max_friends,
  int ld_swap_moves,
  arma::mat annotation,
  arma::mat theta_init,
  double theta_prior_sd,
  bool updateTheta,
  Rcpp::IntegerVector convergence_markers,
  bool convergence_b,
  bool convergence_d
) {
  const int nt = static_cast<int>(wy.size());
  if (nt <= 0) throw std::invalid_argument("BayesC-LV requires at least one trait.");
  const int m = static_cast<int>(wy[0].size());
  if (m <= 0 || nit <= 0 || nburn < 0 || nthin <= 0 || nchains <= 0 ||
      ncores <= 0) {
    throw std::invalid_argument("BayesC-LV MCMC dimensions are invalid.");
  }
  if (!chain_seeds.empty() && static_cast<int>(chain_seeds.size()) != nchains) {
    throw std::invalid_argument("BayesC-LV chain_seeds must have length nchains.");
  }
  if (annotation.n_rows != static_cast<arma::uword>(m) ||
      annotation.n_cols == 0 || !annotation.is_finite()) {
    throw std::invalid_argument("BayesC-LV annotation must be a finite m by p matrix.");
  }
  const int p = static_cast<int>(annotation.n_cols);
  if (theta_init.n_rows != static_cast<arma::uword>(p) ||
      theta_init.n_cols != static_cast<arma::uword>(nt) ||
      !theta_init.is_finite()) {
    throw std::invalid_argument("BayesC-LV theta_init must be finite p by nt.");
  }
  if (!std::isfinite(theta_prior_sd) || theta_prior_sd <= 0.0) {
    throw std::invalid_argument("BayesC-LV theta_prior_sd must be positive finite.");
  }
  if (!std::isfinite(ld_swap_prob) || ld_swap_prob < 0.0 ||
      ld_swap_prob > 1.0 || !std::isfinite(ld_swap_r2) || ld_swap_r2 < 0.0 ||
      ld_swap_r2 > 1.0 || ld_swap_max_friends <= 0 || ld_swap_moves < 0) {
    throw std::invalid_argument("BayesC-LV LD-swap controls are invalid.");
  }
  if (static_cast<int>(ww.size()) != nt || static_cast<int>(yy.size()) != nt ||
      static_cast<int>(b_init.size()) != nt || static_cast<int>(n.size()) != nt ||
      static_cast<int>(ssb_prior.size()) != nt ||
      static_cast<int>(sse_prior.size()) != nt || pi.size() != 2) {
    throw std::invalid_argument("BayesC-LV trait/prior dimensions are inconsistent.");
  }
  if (B.n_rows != static_cast<arma::uword>(nt) ||
      B.n_cols != static_cast<arma::uword>(nt) ||
      E.n_rows != static_cast<arma::uword>(nt) ||
      E.n_cols != static_cast<arma::uword>(nt)) {
    throw std::invalid_argument("BayesC-LV B and E must be nt by nt.");
  }

  arma::mat wy_mat(nt, m, arma::fill::zeros);
  arma::mat ww_mat(nt, m, arma::fill::zeros);
  arma::mat b_mat(nt, m, arma::fill::zeros);
  arma::vec yy_vec(nt, arma::fill::zeros);
  arma::mat ssb_mat(nt, nt, arma::fill::zeros);
  arma::mat sse_mat(nt, nt, arma::fill::zeros);
  for (int trait = 0; trait < nt; ++trait) {
    if (static_cast<int>(wy[trait].size()) != m ||
        static_cast<int>(ww[trait].size()) != m ||
        static_cast<int>(b_init[trait].size()) != m ||
        static_cast<int>(ssb_prior[trait].size()) != nt ||
        static_cast<int>(sse_prior[trait].size()) != nt) {
      throw std::invalid_argument("BayesC-LV marker/prior dimensions are inconsistent.");
    }
    yy_vec(static_cast<arma::uword>(trait)) = yy[trait];
    for (int marker = 0; marker < m; ++marker) {
      wy_mat(trait, marker) = wy[trait][marker];
      ww_mat(trait, marker) = ww[trait][marker];
      b_mat(trait, marker) = b_init[trait][marker];
      if (!std::isfinite(ww_mat(trait, marker)) || ww_mat(trait, marker) <= 0.0) {
        throw std::invalid_argument("BayesC-LV ww must be positive finite.");
      }
    }
    for (int other = 0; other < nt; ++other) {
      ssb_mat(trait, other) = ssb_prior[trait][other];
      sse_mat(trait, other) = sse_prior[trait][other];
    }
  }
  for (int trait = 1; trait < nt; ++trait) {
    if (n[trait] != n[0]) {
      throw std::invalid_argument("BayesC-LV shared CSR requires equal sample sizes.");
    }
    for (int marker = 0; marker < m; ++marker) {
      const double tolerance = 1e-8 * std::max(1.0, std::abs(ww_mat(0, marker)));
      if (std::abs(ww_mat(trait, marker) - ww_mat(0, marker)) > tolerance) {
        throw std::invalid_argument("BayesC-LV shared CSR requires equal ww across traits.");
      }
    }
  }
  if (use_d_init && static_cast<int>(d_init.size()) != nt) {
    throw std::invalid_argument("BayesC-LV d_init trait count mismatch.");
  }
  if (use_r_init && static_cast<int>(r_init.size()) != nt) {
    throw std::invalid_argument("BayesC-LV r_init trait count mismatch.");
  }

  std::vector<double> xx(static_cast<std::size_t>(m));
  arma::rowvec xx_row(m);
  for (int marker = 0; marker < m; ++marker) {
    xx[static_cast<std::size_t>(marker)] = ww_mat(0, marker);
    xx_row(static_cast<arma::uword>(marker)) = ww_mat(0, marker);
  }
  STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);
  CsrOperator op(ld, xx_row);
  LDLDFriends friends;
  if (updateLDswap) {
    friends = build_ld_swap_friends_st_csr(
      m, ld, xx, ld_swap_r2, ld_swap_max_friends);
  } else {
    friends.ptr.assign(static_cast<std::size_t>(m) + 1, 0);
  }

  std::vector<int> order(static_cast<std::size_t>(m));
  std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
  for (int marker = 0; marker < m; ++marker) {
    order[static_cast<std::size_t>(marker)] = marker;
    for (int trait = 0; trait < nt; ++trait) {
      const double bhat = wy_mat(trait, marker) / ww_mat(trait, marker);
      x2[static_cast<std::size_t>(marker)] = std::max(
        x2[static_cast<std::size_t>(marker)], bhat * bhat);
    }
  }
  std::sort(order.begin(), order.end(), [&](int left, int right) {
    return x2[static_cast<std::size_t>(left)] > x2[static_cast<std::size_t>(right)];
  });
  const std::vector<int> convergence = Rcpp::as<std::vector<int>>(
    convergence_markers);

  sblr::core::ResolvedSpec specification;
  specification.data.marker_count = static_cast<std::size_t>(m);
  specification.data.trait_count = static_cast<std::size_t>(nt);
  for (int marker = 0; marker < m; ++marker) {
    specification.data.marker_ids.push_back("marker_" + std::to_string(marker + 1));
  }
  for (int trait = 0; trait < nt; ++trait) {
    specification.data.trait_ids.push_back("trait_" + std::to_string(trait + 1));
  }
  specification.data.sample_size = n;
  specification.data.csr.resource_id = ld_prefix;
  specification.data.csr.marker_count = static_cast<std::size_t>(m);
  specification.mcmc.nit = nit;
  specification.mcmc.nburn = nburn;
  specification.mcmc.nthin = nthin;
  specification.mcmc.nchains = nchains;
  specification.mcmc.ncores = ncores;
  specification.mcmc.seed = seed;
  specification.mcmc.has_explicit_chain_seeds = !chain_seeds.empty();
  specification.mcmc.chain_seeds = chain_seeds;
  specification.output.keep_chain_summaries = keep_chains;

  sblr::core::CsrBayesCExecutionInput input;
  input.specification = std::move(specification);
  input.data.marker_count = static_cast<std::size_t>(m);
  input.data.trait_count = static_cast<std::size_t>(nt);
  input.data.ld = op.view();
  input.data.wy = &wy_mat;
  input.data.yy = &yy_vec;
  input.data.sample_size = &n;
  input.priors.marker_variance = &B;
  input.priors.residual_variance = &E;
  input.priors.marker_scale_prior = &ssb_mat;
  input.priors.residual_scale_prior = &sse_mat;
  input.priors.inclusion_probability = pi;
  input.priors.marker_degrees_freedom = nub;
  input.priors.residual_degrees_freedom = nue;
  input.priors.inclusion_prior_active = pi_prior_a;
  input.priors.inclusion_prior_null = pi_prior_b;
  input.initial.effects = &b_mat;
  input.initial.inclusion = &d_init;
  input.initial.residual = &r_init;
  input.initial.use_inclusion = use_d_init;
  input.initial.use_residual = use_r_init;
  input.controls.nit = nit;
  input.controls.nburn = nburn;
  input.controls.nthin = nthin;
  input.controls.nchains = nchains;
  input.controls.ncores = ncores;
  input.controls.seed = seed;
  input.controls.chain_seeds = chain_seeds;
  input.controls.keep_chains = keep_chains;
  input.controls.update_marker_variance = updateB;
  input.controls.update_residual_variance = updateE;
  input.controls.update_inclusion_probability = updatePi;
  input.controls.rebuild_residual_before_update = rebuild_r_before_updateE;
  input.controls.residual_adjustment = adjE;
  input.controls.update_ld_swap = updateLDswap;
  input.controls.ld_swap_probability = ld_swap_prob;
  input.controls.ld_swap_moves = ld_swap_moves;
  input.controls.convergence_markers = convergence;
  input.controls.convergence_b = convergence_b;
  input.controls.convergence_d = convergence_d;
  input.output.keep_chains = keep_chains;
  input.ld_friends.row_ptr = friends.ptr.data();
  input.ld_friends.row_ptr_size = friends.ptr.size();
  input.ld_friends.index = friends.idx.empty() ? nullptr : friends.idx.data();
  input.ld_friends.friend_count = friends.idx.size();
  input.marker_order = &order;

  const int ntasks = nt * nchains;
  std::vector<sblr::logvar::CsrLogvarBayesCChainOutput> policy_outputs(
    static_cast<std::size_t>(ntasks));
  sblr::logvar::CsrLogvarBayesCPolicyInput policy_input;
  policy_input.annotation = &annotation;
  policy_input.theta_initial = &theta_init;
  policy_input.theta_prior_sd = theta_prior_sd;
  policy_input.update_theta = updateTheta;
  policy_input.trace_count = nit + nburn;
  policy_input.burnin = nburn;
  policy_input.thinning = nthin;
  policy_input.outputs = &policy_outputs;
  sblr::logvar::CsrLogvarBayesCPolicyFactory policy_factory(policy_input);
  const sblr::core::CsrBayesCResult result =
    sblr::core::run_csr_bayesc_engine(input, policy_factory);

  const CsrBayesCRawConversionContext conversion = {
    m, nt, nit, nburn, nthin, ncores, nchains, keep_chains,
    pi_prior_a, pi_prior_b, updateLDswap, false, false, &n, &convergence
  };
  Rcpp::List raw = stblr_csr_bayesc_result_to_raw(result, conversion);
  Rcpp::List meta = raw["meta"];
  meta["model"] = "sbayesc_logvar";
  meta["backend"] = "csr_logvar_bayesc";
  meta["prior_type"] = "annotation_log_variance";
  meta["n_annotations"] = p;

  Rcpp::NumericMatrix theta_mean(p, nt);
  Rcpp::NumericMatrix q_mean(m, nt);
  Rcpp::NumericMatrix q_chain_mean(m, ntasks);
  Rcpp::NumericVector theta_trace(
    static_cast<R_xlen_t>(nit + nburn) * p * ntasks);
  theta_trace.attr("dim") = Rcpp::IntegerVector::create(nit + nburn, p, ntasks);
  double updates = 0.0, evaluations = 0.0, contractions = 0.0;
  std::size_t max_evaluations = 0, max_contractions = 0;
  double min_log_q = std::numeric_limits<double>::infinity();
  double max_log_q = -std::numeric_limits<double>::infinity();
  for (int trait = 0; trait < nt; ++trait) {
    for (int chain = 0; chain < nchains; ++chain) {
      const int task = trait * nchains + chain;
      const auto& current = policy_outputs[static_cast<std::size_t>(task)];
      const double retained = std::max(current.retained_samples, 1.0);
      for (int column = 0; column < p; ++column) {
        theta_mean(column, trait) += current.theta_sum(column) /
          retained / static_cast<double>(nchains);
        for (int iteration = 0; iteration < nit + nburn; ++iteration) {
          theta_trace[iteration + (nit + nburn) *
            (column + p * task)] = current.theta_trace(iteration, column);
        }
      }
      for (int marker = 0; marker < m; ++marker) {
        const double chain_q = current.prior_scale_sum(marker) / retained;
        q_chain_mean(marker, task) = chain_q;
        q_mean(marker, trait) += chain_q / static_cast<double>(nchains);
      }
      updates += current.diagnostics.theta_updates;
      evaluations += current.diagnostics.likelihood_evaluations;
      contractions += current.diagnostics.bracket_contractions;
      max_evaluations = std::max(
        max_evaluations, current.diagnostics.max_likelihood_evaluations);
      max_contractions = std::max(
        max_contractions, current.diagnostics.max_bracket_contractions);
      min_log_q = std::min(min_log_q, current.diagnostics.min_log_q);
      max_log_q = std::max(max_log_q, current.diagnostics.max_log_q);
    }
  }
  Rcpp::NumericMatrix variance_ratio(p, nt);
  for (int trait = 0; trait < nt; ++trait) {
    for (int column = 0; column < p; ++column) {
      variance_ratio(column, trait) = std::exp(theta_mean(column, trait));
    }
  }
  raw["annotation"] = Rcpp::List::create(
    Rcpp::Named("theta") = theta_mean,
    Rcpp::Named("theta_trace") = theta_trace,
    Rcpp::Named("variance_ratio") = variance_ratio,
    Rcpp::Named("marker_prior_scale") = q_mean,
    Rcpp::Named("marker_prior_scale_chain") = q_chain_mean,
    Rcpp::Named("theta_prior_sd") = theta_prior_sd,
    Rcpp::Named("update_theta") = updateTheta
  );
  Rcpp::List diagnostics = raw["diagnostics"];
  diagnostics["logvar"] = Rcpp::List::create(
    Rcpp::Named("theta_updates") = updates,
    Rcpp::Named("mean_likelihood_evaluations_per_update") =
      updates > 0.0 ? evaluations / updates : 0.0,
    Rcpp::Named("max_likelihood_evaluations") = max_evaluations,
    Rcpp::Named("mean_bracket_contractions") =
      updates > 0.0 ? contractions / updates : 0.0,
    Rcpp::Named("max_bracket_contractions") = max_contractions,
    Rcpp::Named("min_log_q") = min_log_q,
    Rcpp::Named("max_log_q") = max_log_q
  );
  raw["diagnostics"] = diagnostics;
  return raw;
}
