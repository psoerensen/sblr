#ifndef ST_BAYESRC_ANNOTATION_PRIOR_H
#define ST_BAYESRC_ANNOTATION_PRIOR_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include <armadillo>
#include "blr_normal_probability.h"

inline double st_bayesrc_safe_pnorm(double x) {
 double p = sblr::core::StandardNormalProbability::cdf(x);
 if (!std::isfinite(p)) p = (x > 0.0) ? 1.0 : 0.0;
 return std::min(std::max(p, 1e-12), 1.0 - 1e-12);
}

inline double st_bayesrc_safe_qnorm(double p) {
 p = std::min(std::max(p, 1e-12), 1.0 - 1e-12);
 return sblr::core::StandardNormalProbability::quantile(p);
}

inline double st_bayesrc_sample_truncated_normal_std(
  double mu,
  bool positive,
  std::mt19937& gen
) {
 if (!std::isfinite(mu)) {
  throw std::runtime_error(
   "BayesRC truncated-normal location must be finite.");
 }

 // Work with Z ~ N(0, 1), Z > -mu. For a non-positive lower bound,
 // direct rejection accepts at least half of proposals. For a positive bound,
 // use the exact exponential-envelope rejection kernel; unlike clipped
 // inverse-CDF sampling, it remains valid in the far tail.
 const auto sample_lower = [&gen](double location) {
  const double lower = -location;
  std::normal_distribution<double> normal(0.0, 1.0);
  if (lower <= 0.0) {
   double z = 0.0;
   do z = normal(gen); while (z <= lower);
   const double draw = location + z;
   return draw > 0.0 ? draw : std::nextafter(0.0, 1.0);
  }

  const double rate = 0.5 * (lower + std::hypot(lower, 2.0));
  std::exponential_distribution<double> exponential(rate);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  for (;;) {
   const double z = lower + exponential(gen);
   const double log_acceptance = -0.5 * (z - rate) * (z - rate);
   const double u = std::max(uniform(gen),
    std::numeric_limits<double>::min());
   if (std::log(u) <= log_acceptance) {
    const double draw = location + z;
    return draw > 0.0 ? draw : std::nextafter(0.0, 1.0);
   }
  }
 };

 return positive ? sample_lower(mu) : -sample_lower(-mu);
}

inline arma::mat st_bayesrc_compute_tempered_snp_pi(
  const arma::mat& annotation_matrix,
  const arma::mat& annot_alpha,
  const arma::vec& baseline_intercept,
  double coupling,
  double pi_floor
) {
 const int m = static_cast<int>(annotation_matrix.n_rows);
 const int nstep = static_cast<int>(annot_alpha.n_cols);
 const int ncomponent = nstep + 1;
 if (baseline_intercept.n_elem != static_cast<arma::uword>(nstep) ||
     !baseline_intercept.is_finite() || !std::isfinite(coupling) ||
     coupling < 0.0 || coupling > 1.0) {
  throw std::runtime_error("BayesRC tempered coupling inputs are invalid.");
 }

 arma::mat p(m, nstep, arma::fill::zeros);
 arma::mat snp_pi(m, ncomponent, arma::fill::zeros);

 for (int j = 0; j < nstep; ++j) {
  arma::vec eta = (1.0 - coupling) * baseline_intercept(
   static_cast<arma::uword>(j)) +
   coupling * annotation_matrix * annot_alpha.col(static_cast<arma::uword>(j));
  for (int i = 0; i < m; ++i) {
   p(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) =
    st_bayesrc_safe_pnorm(eta(static_cast<arma::uword>(i)));
  }
 }

 for (int i = 0; i < m; ++i) {
  double prod_prev = 1.0;

  for (int k = 0; k < ncomponent; ++k) {
   double val;

   if (k < nstep) {
    const double pk = p(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
    val = prod_prev * (1.0 - pk);
    prod_prev *= pk;
   } else {
    val = prod_prev;
   }

   snp_pi(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) =
    std::max(val, pi_floor);
  }

  const double s = arma::accu(snp_pi.row(static_cast<arma::uword>(i)));
  if (!std::isfinite(s) || s <= 0.0) {
   throw std::runtime_error("compute_snp_pi_from_alpha: invalid row probability sum.");
  }
  snp_pi.row(static_cast<arma::uword>(i)) /= s;
 }

 return snp_pi;
}

inline arma::mat st_bayesrc_compute_snp_pi(
  const arma::mat& annotation_matrix,
  const arma::mat& annot_alpha,
  double pi_floor
) {
 return st_bayesrc_compute_tempered_snp_pi(
  annotation_matrix, annot_alpha,
  arma::vec(annot_alpha.n_cols, arma::fill::zeros), 1.0, pi_floor);
}

inline double st_bayesrc_tempered_log_allocation_prior(
  const arma::mat& annotation_matrix,
  const arma::mat& annot_alpha,
  const arma::vec& baseline_intercept,
  double coupling,
  double pi_floor,
  const arma::Row<int>& component
) {
 const arma::mat probability = st_bayesrc_compute_tempered_snp_pi(
  annotation_matrix, annot_alpha, baseline_intercept, coupling, pi_floor);
 if (component.n_elem != probability.n_rows) {
  throw std::runtime_error(
   "BayesRC tempered allocation state has the wrong marker count.");
 }
 double result = 0.0;
 for (arma::uword marker = 0; marker < component.n_elem; ++marker) {
  const int state = component(marker);
  if (state < 0 || state >= static_cast<int>(probability.n_cols)) {
   throw std::runtime_error(
    "BayesRC tempered allocation state contains an invalid component.");
  }
  result += std::log(probability(marker, static_cast<arma::uword>(state)));
 }
 return result;
}

inline void st_bayesrc_build_step_indicators(
  const arma::Row<int>& component,
  arma::Mat<int>& step_indicators
) {
 const int m = static_cast<int>(component.n_elem);
 const int nstep = static_cast<int>(step_indicators.n_cols);

 for (int i = 0; i < m; ++i) {
  const int ci = component(static_cast<arma::uword>(i));
  for (int j = 0; j < nstep; ++j) {
   step_indicators(static_cast<arma::uword>(i), static_cast<arma::uword>(j)) =
    (ci > j) ? 1 : 0;
  }
 }
}

// Resolved at the R boundary.  The intercept prior is deliberately separate
// from alpha_init and from the hierarchical prior on non-intercept effects.
struct StBayesRCInterceptPrior {
 bool legacy_flat;
 bool update_variance;
 int allocation_updates_per_cycle;
 int annotation_updates_per_cycle;
 bool coupling_tempering;
 int coupling_swap_every;
 bool px_sandwich;
 double px_log_scale_sd;
 arma::vec mean;
 arma::vec precision;
};

inline StBayesRCInterceptPrior st_bayesrc_parse_intercept_prior(
  const arma::mat& resolved,
  int nstep
) {
 if ((resolved.n_rows != 3 && resolved.n_rows != 4 &&
      resolved.n_rows != 6 && resolved.n_rows != 8 &&
      resolved.n_rows != 10) ||
     resolved.n_cols != static_cast<arma::uword>(nstep) ||
     !resolved.is_finite()) {
  throw std::runtime_error(
   "BayesRC resolved annotation control must be a finite 3, 4, 6, 8, or 10 by nstep matrix.");
 }
 const bool legacy_flat = resolved(0, 0) == 1.0;
 const bool update_variance = resolved.n_rows == 3 || resolved(3, 0) == 1.0;
 const int allocation_updates = resolved.n_rows >= 6
  ? static_cast<int>(resolved(4, 0)) : 1;
 const int annotation_updates = resolved.n_rows >= 6
  ? static_cast<int>(resolved(5, 0)) : 1;
 const bool coupling_tempering = resolved.n_rows >= 8 && resolved(6, 0) == 1.0;
 const int coupling_swap_every = resolved.n_rows >= 8
  ? static_cast<int>(resolved(7, 0)) : 0;
 const bool px_sandwich = resolved.n_rows == 10 && resolved(8, 0) == 1.0;
 const double px_log_scale_sd = resolved.n_rows == 10 ? resolved(9, 0) : 0.0;
 if (allocation_updates <= 0 || annotation_updates <= 0 ||
     (resolved.n_rows >= 6 &&
      (resolved(4, 0) != static_cast<double>(allocation_updates) ||
       resolved(5, 0) != static_cast<double>(annotation_updates)))) {
  throw std::runtime_error(
   "BayesRC diagnostic kernel-update counts must be positive integers.");
 }
 if (resolved.n_rows >= 8 &&
     ((!coupling_tempering && resolved(6, 0) != 0.0) ||
      (coupling_tempering && coupling_swap_every <= 0) ||
      resolved(7, 0) != static_cast<double>(coupling_swap_every))) {
  throw std::runtime_error(
   "BayesRC coupling-tempering controls are invalid.");
 }
 if (resolved.n_rows == 10 &&
     ((!px_sandwich && resolved(8, 0) != 0.0) ||
      (px_sandwich && (!std::isfinite(px_log_scale_sd) ||
                       px_log_scale_sd <= 0.0)))) {
  throw std::runtime_error("BayesRC PX sandwich controls are invalid.");
 }
 for (int j = 0; j < nstep; ++j) {
  if (resolved(0, static_cast<arma::uword>(j)) != (legacy_flat ? 1.0 : 0.0)) {
   throw std::runtime_error("BayesRC intercept prior type must be consistent across sticks.");
  }
  if (!legacy_flat && resolved(2, static_cast<arma::uword>(j)) <= 0.0) {
   throw std::runtime_error("BayesRC proper intercept prior precisions must be positive.");
  }
  if (resolved.n_rows == 4 &&
      resolved(3, static_cast<arma::uword>(j)) !=
       (update_variance ? 1.0 : 0.0)) {
   throw std::runtime_error(
    "BayesRC annotation-variance update policy must be consistent across sticks.");
  }
  if (resolved.n_rows >= 6 &&
      (resolved(3, static_cast<arma::uword>(j)) !=
        (update_variance ? 1.0 : 0.0) ||
       resolved(4, static_cast<arma::uword>(j)) != allocation_updates ||
       resolved(5, static_cast<arma::uword>(j)) != annotation_updates)) {
   throw std::runtime_error(
    "BayesRC diagnostic annotation controls must be consistent across sticks.");
  }
  if (resolved.n_rows >= 8 &&
      (resolved(6, static_cast<arma::uword>(j)) !=
        (coupling_tempering ? 1.0 : 0.0) ||
       resolved(7, static_cast<arma::uword>(j)) != coupling_swap_every)) {
   throw std::runtime_error(
    "BayesRC coupling-tempering controls must be consistent across sticks.");
  }
  if (resolved.n_rows == 10 &&
      (resolved(8, static_cast<arma::uword>(j)) !=
        (px_sandwich ? 1.0 : 0.0) ||
       resolved(9, static_cast<arma::uword>(j)) != px_log_scale_sd)) {
   throw std::runtime_error(
    "BayesRC PX sandwich controls must be consistent across sticks.");
  }
 }
 return StBayesRCInterceptPrior{
  legacy_flat, update_variance, allocation_updates, annotation_updates,
  coupling_tempering, coupling_swap_every, px_sandwich, px_log_scale_sd,
  resolved.row(1).t(), resolved.row(2).t()
 };
}

struct StBayesRCAnnotationUpdateDiagnostics {
 std::vector<int> eligible;
 std::vector<int> continuation;
 std::vector<int> prior_only;
 std::vector<int> px_attempted;
 std::vector<int> px_accepted;
 std::vector<double> px_abs_log_scale;
 std::vector<double> px_alpha_jump;
};

struct StBayesRCScalarConditional {
 double mean;
 double variance;
};

inline StBayesRCScalarConditional st_bayesrc_scalar_conditional(
 double likelihood_diagonal,
 double likelihood_rhs,
 double prior_mean,
 double prior_precision
) {
 if (!std::isfinite(likelihood_diagonal) || likelihood_diagonal < 0.0 ||
     !std::isfinite(likelihood_rhs) || !std::isfinite(prior_mean) ||
     !std::isfinite(prior_precision) || prior_precision < 0.0 ||
     likelihood_diagonal + prior_precision <= 0.0) {
  throw std::runtime_error("BayesRC alpha conditional parameters are invalid.");
 }
 const double variance = 1.0 / (likelihood_diagonal + prior_precision);
 return StBayesRCScalarConditional{
  variance * (likelihood_rhs + prior_precision * prior_mean), variance
 };
}

inline double st_bayesrc_sample_sigma_sq_alpha(
 double sum_squares,
 int coefficient_count,
 double prior_a,
 double prior_b,
 std::mt19937& gen
) {
 const double df = static_cast<double>(coefficient_count) + prior_a;
 const double scale = sum_squares + prior_b;
 if (!std::isfinite(scale) || scale <= 0.0 || !std::isfinite(df) || df <= 0.0) {
  throw std::runtime_error("BayesRC sigmaSqAlpha conditional parameters are invalid.");
 }
 std::chi_squared_distribution<double> rchisq(df);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 const double value = scale / chi2;
 if (!std::isfinite(value) || value <= 0.0) {
  throw std::runtime_error("BayesRC sigmaSqAlpha update produced an invalid value.");
 }
 return value;
}

inline StBayesRCAnnotationUpdateDiagnostics
st_bayesrc_update_tempered_annotation_prior(
  const arma::mat& annotation_matrix,
  const arma::Row<int>& component,
  arma::mat& annot_alpha,
  arma::vec& annot_sigma_sq_alpha,
  const StBayesRCInterceptPrior& intercept_prior,
  double annot_sigma_sq_alpha_a,
  double annot_sigma_sq_alpha_b,
  const arma::vec& baseline_intercept,
  double coupling,
  std::mt19937& gen
) {
 const int m = static_cast<int>(annotation_matrix.n_rows);
 const int n_annot = static_cast<int>(annotation_matrix.n_cols);
 const int nstep = static_cast<int>(annot_alpha.n_cols);
 if (baseline_intercept.n_elem != static_cast<arma::uword>(nstep) ||
     !baseline_intercept.is_finite() || !std::isfinite(coupling) ||
     coupling < 0.0 || coupling > 1.0) {
  throw std::runtime_error("BayesRC tempered annotation inputs are invalid.");
 }

 arma::Mat<int> step_indicators(m, nstep, arma::fill::zeros);
 st_bayesrc_build_step_indicators(component, step_indicators);
 StBayesRCAnnotationUpdateDiagnostics diagnostics{
 std::vector<int>(static_cast<std::size_t>(nstep), 0),
 std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<double>(static_cast<std::size_t>(nstep), 0.0),
  std::vector<double>(static_cast<std::size_t>(nstep), 0.0)
 };

 for (int j = 0; j < nstep; ++j) {
  std::vector<int> idx;
  idx.reserve(static_cast<std::size_t>(m));

  for (int i = 0; i < m; ++i) {
   if (j == 0 || step_indicators(static_cast<arma::uword>(i),
                                  static_cast<arma::uword>(j - 1)) > 0) {
    idx.push_back(i);
   }
  }

  const int nj = static_cast<int>(idx.size());
  diagnostics.eligible[static_cast<std::size_t>(j)] = nj;
  int n_continuation = 0;
  for (int i : idx) {
   n_continuation += step_indicators(static_cast<arma::uword>(i),
                                     static_cast<arma::uword>(j));
  }
  diagnostics.continuation[static_cast<std::size_t>(j)] = n_continuation;

  if (idx.empty()) {
   if (intercept_prior.legacy_flat) {
    throw std::runtime_error(
     "BayesRC legacy flat intercept prior is improper for an empty eligible stick.");
   }
   diagnostics.prior_only[static_cast<std::size_t>(j)] = 1;
   std::normal_distribution<double> intercept_norm(
    intercept_prior.mean(static_cast<arma::uword>(j)),
    1.0 / std::sqrt(intercept_prior.precision(static_cast<arma::uword>(j))));
   annot_alpha(0, static_cast<arma::uword>(j)) = intercept_norm(gen);
   const double sig = annot_sigma_sq_alpha(static_cast<arma::uword>(j));
   if (!std::isfinite(sig) || sig <= 0.0) {
    throw std::runtime_error("BayesRC sigmaSqAlpha must remain positive and finite.");
   }
   std::normal_distribution<double> coefficient_norm(0.0, std::sqrt(sig));
   for (int k = 1; k < n_annot; ++k) {
    annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) =
     coefficient_norm(gen);
   }
  } else if (intercept_prior.legacy_flat &&
             (n_continuation == 0 || n_continuation == nj)) {
   throw std::runtime_error(
    "BayesRC legacy flat intercept prior is unsafe under complete separation; use the proper normal prior.");
  }

  arma::vec mu(nj, arma::fill::zeros);

  for (int ii = 0; ii < nj; ++ii) {
   const int i = idx[static_cast<std::size_t>(ii)];
   double s = (1.0 - coupling) * baseline_intercept(
    static_cast<arma::uword>(j));
   for (int k = 0; k < n_annot; ++k) {
    s += coupling * annotation_matrix(
     static_cast<arma::uword>(i), static_cast<arma::uword>(k)) *
     annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j));
   }
   mu(static_cast<arma::uword>(ii)) = s;
  }

  arma::vec latent(nj, arma::fill::zeros);

  for (int ii = 0; ii < nj; ++ii) {
   const int i = idx[static_cast<std::size_t>(ii)];
   const bool positive = step_indicators(static_cast<arma::uword>(i),
                                          static_cast<arma::uword>(j)) > 0;
   latent(static_cast<arma::uword>(ii)) = st_bayesrc_sample_truncated_normal_std(
    mu(static_cast<arma::uword>(ii)), positive, gen
   );
  }

  if (intercept_prior.px_sandwich && idx.empty()) {
   // With no likelihood contribution the branch above has already drawn the
   // complete coefficient vector from its conditional prior.  There is no
   // latent scale on which to apply a sandwich move; finish with the unchanged
   // sigmaSqAlpha full conditional.
   double ss = 0.0;
   int ncoef = 0;
   for (int k = 1; k < n_annot; ++k) {
    const double ak = annot_alpha(static_cast<arma::uword>(k),
                                  static_cast<arma::uword>(j));
    ss += ak * ak;
    ++ncoef;
   }
   if (intercept_prior.update_variance) {
    annot_sigma_sq_alpha(static_cast<arma::uword>(j)) =
     st_bayesrc_sample_sigma_sq_alpha(
      ss, ncoef, annot_sigma_sq_alpha_a, annot_sigma_sq_alpha_b, gen);
   }
   continue;
  }

  if (intercept_prior.px_sandwich) {
   if (coupling != 1.0) {
    throw std::runtime_error(
     "BayesRC PX sandwich is only defined for the ordinary coupling endpoint.");
   }
   arma::mat design(static_cast<arma::uword>(nj),
                    static_cast<arma::uword>(n_annot), arma::fill::zeros);
   for (int ii = 0; ii < nj; ++ii) {
    design.row(static_cast<arma::uword>(ii)) = annotation_matrix.row(
     static_cast<arma::uword>(idx[static_cast<std::size_t>(ii)]));
   }
   arma::vec prior_mean(n_annot, arma::fill::zeros);
   arma::vec prior_precision(n_annot, arma::fill::zeros);
   if (intercept_prior.legacy_flat) {
    throw std::runtime_error(
     "BayesRC PX sandwich requires a proper intercept prior.");
   }
   prior_mean(0) = intercept_prior.mean(static_cast<arma::uword>(j));
   prior_precision(0) = intercept_prior.precision(static_cast<arma::uword>(j));
   const double sig = annot_sigma_sq_alpha(static_cast<arma::uword>(j));
   if (!std::isfinite(sig) || sig <= 0.0) {
    throw std::runtime_error("BayesRC sigmaSqAlpha must remain positive and finite.");
   }
   for (int k = 1; k < n_annot; ++k) prior_precision(k) = 1.0 / sig;

   const arma::mat precision = design.t() * design + arma::diagmat(prior_precision);
   arma::mat chol_precision;
   if (!arma::chol(chol_precision, precision)) {
    throw std::runtime_error("BayesRC PX alpha precision is not positive definite.");
   }
   const arma::vec xtz = design.t() * latent;
   const arma::vec prior_rhs = prior_precision % prior_mean;
   const arma::vec solved_xtz = arma::solve(
    arma::trimatu(chol_precision),
    arma::solve(arma::trimatl(chol_precision.t()), xtz));
   const arma::vec solved_prior_rhs = arma::solve(
    arma::trimatu(chol_precision),
    arma::solve(arma::trimatl(chol_precision.t()), prior_rhs));
   // Evaluate z'z - z'X P^{-1} X'z through its positive residual-plus-prior
   // identity.  The direct subtraction loses precision for small, nearly
   // saturated eligible sets even though P is strictly positive definite.
   const arma::vec projected_residual = latent - design * solved_xtz;
   const double quadratic_a =
    arma::dot(projected_residual, projected_residual) +
    arma::dot(prior_precision % solved_xtz, solved_xtz);
   const double linear_b = arma::dot(xtz, solved_prior_rhs);
   if (!std::isfinite(quadratic_a) || quadratic_a <= 0.0 ||
       !std::isfinite(linear_b)) {
    throw std::runtime_error("BayesRC PX latent marginal is invalid.");
   }

   std::normal_distribution<double> log_scale_proposal(
    0.0, intercept_prior.px_log_scale_sd);
   const double log_scale = log_scale_proposal(gen);
   double log_ratio = -std::numeric_limits<double>::infinity();
   if (std::isfinite(log_scale) && log_scale < 350.0) {
    const double scale = std::exp(log_scale);
    const double scale_sq = scale * scale;
    if (std::isfinite(scale_sq)) {
     log_ratio = static_cast<double>(nj) * log_scale - 0.5 *
      (quadratic_a * (scale_sq - 1.0) -
       2.0 * linear_b * (scale - 1.0));
    }
   }
   diagnostics.px_attempted[static_cast<std::size_t>(j)] = 1;
   std::uniform_real_distribution<double> uniform(0.0, 1.0);
   const double log_uniform = std::log(std::max(
    uniform(gen), std::numeric_limits<double>::min()));
   if (log_uniform < std::min(0.0, log_ratio)) {
    const double scale = std::exp(log_scale);
    latent *= scale;
    diagnostics.px_accepted[static_cast<std::size_t>(j)] = 1;
    diagnostics.px_abs_log_scale[static_cast<std::size_t>(j)] =
     std::abs(log_scale);
   }

   const arma::vec alpha_old = annot_alpha.col(static_cast<arma::uword>(j));
   const arma::vec conditional_rhs = design.t() * latent + prior_rhs;
   const arma::vec conditional_mean = arma::solve(
    arma::trimatu(chol_precision),
    arma::solve(arma::trimatl(chol_precision.t()), conditional_rhs));
   arma::vec standard_normal(n_annot, arma::fill::zeros);
   std::normal_distribution<double> normal(0.0, 1.0);
   for (int k = 0; k < n_annot; ++k) standard_normal(k) = normal(gen);
   const arma::vec alpha_new = conditional_mean + arma::solve(
    arma::trimatu(chol_precision), standard_normal);
   if (!alpha_new.is_finite()) {
    throw std::runtime_error("BayesRC PX alpha update produced a non-finite state.");
   }
   annot_alpha.col(static_cast<arma::uword>(j)) = alpha_new;
   diagnostics.px_alpha_jump[static_cast<std::size_t>(j)] =
    arma::norm(alpha_new - alpha_old);

   double ss = 0.0;
   int ncoef = 0;
   for (int k = 1; k < n_annot; ++k) {
    ss += alpha_new(static_cast<arma::uword>(k)) *
      alpha_new(static_cast<arma::uword>(k));
    ++ncoef;
   }
   if (intercept_prior.update_variance) {
    annot_sigma_sq_alpha(static_cast<arma::uword>(j)) =
     st_bayesrc_sample_sigma_sq_alpha(
      ss, ncoef, annot_sigma_sq_alpha_a, annot_sigma_sq_alpha_b, gen);
   }
   continue;
  }

  arma::vec resid = latent - mu;

  for (int k = 0; k < n_annot; ++k) {
   const bool intercept = k == 0;
   double diag_k = 0.0;
   double rhs = 0.0;
   const double old = annot_alpha(static_cast<arma::uword>(k),
                                  static_cast<arma::uword>(j));

   for (int ii = 0; ii < nj; ++ii) {
    const int i = idx[static_cast<std::size_t>(ii)];
    const double x = coupling * annotation_matrix(
     static_cast<arma::uword>(i), static_cast<arma::uword>(k));
    diag_k += x * x;
    rhs += x * resid(static_cast<arma::uword>(ii));
   }

   rhs += diag_k * old;
   if (diag_k <= 0.0) {
    // A non-intercept column can be identically zero on the current
    // stick-eligible subset even when the stick itself is non-empty.  Its
    // likelihood is then absent, so its full conditional is its hierarchical
    // prior.  Leaving the old value in place creates an artificial absorbing
    // direction and contaminates the following sigmaSqAlpha update.
    if (intercept && nj == 0) {
     // The prior-only branch above has already drawn this intercept.
     continue;
    }
    if (intercept) {
     if (intercept_prior.legacy_flat) {
      throw std::runtime_error(
       "BayesRC flat intercept prior is improper at zero coupling.");
     }
     std::normal_distribution<double> intercept_norm(
      intercept_prior.mean(static_cast<arma::uword>(j)),
      1.0 / std::sqrt(intercept_prior.precision(static_cast<arma::uword>(j))));
     annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) =
      intercept_norm(gen);
    } else {
     const double sig = annot_sigma_sq_alpha(static_cast<arma::uword>(j));
     if (!std::isfinite(sig) || sig <= 0.0) {
      throw std::runtime_error("BayesRC sigmaSqAlpha must remain positive and finite.");
     }
     std::normal_distribution<double> coefficient_norm(0.0, std::sqrt(sig));
     annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) =
      coefficient_norm(gen);
    }
    continue;
   }

   double prior_mean = 0.0;
   double prior_precision = 0.0;
   if (intercept && intercept_prior.legacy_flat) {
    prior_precision = 0.0;
   } else if (intercept) {
    prior_precision = intercept_prior.precision(static_cast<arma::uword>(j));
    prior_mean = intercept_prior.mean(static_cast<arma::uword>(j));
   } else {
    const double sig = annot_sigma_sq_alpha(static_cast<arma::uword>(j));
    if (!std::isfinite(sig) || sig <= 0.0) {
     throw std::runtime_error("BayesRC sigmaSqAlpha must remain positive and finite.");
    }
    prior_precision = 1.0 / sig;
   }

   const auto conditional = st_bayesrc_scalar_conditional(
    diag_k, rhs, prior_mean, prior_precision);
   std::normal_distribution<double> norm(
    conditional.mean, std::sqrt(conditional.variance));
   const double annot_alpha_new = norm(gen);

   annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) =
    annot_alpha_new;

   const double diff_old_new = old - annot_alpha_new;
   if (diff_old_new != 0.0) {
    for (int ii = 0; ii < nj; ++ii) {
     const int i = idx[static_cast<std::size_t>(ii)];
     const double x = coupling * annotation_matrix(
      static_cast<arma::uword>(i), static_cast<arma::uword>(k));
     resid(static_cast<arma::uword>(ii)) += x * diff_old_new;
    }
   }
  }

  double ss = 0.0;
  int ncoef = 0;
  for (int k = 0; k < n_annot; ++k) {
   if (k == 0) continue;
   const double ak = annot_alpha(static_cast<arma::uword>(k),
                                 static_cast<arma::uword>(j));
   ss += ak * ak;
   ++ncoef;
  }

  if (intercept_prior.update_variance) {
   annot_sigma_sq_alpha(static_cast<arma::uword>(j)) =
    st_bayesrc_sample_sigma_sq_alpha(
     ss, ncoef, annot_sigma_sq_alpha_a, annot_sigma_sq_alpha_b, gen);
  }
 }
 return diagnostics;
}

inline StBayesRCAnnotationUpdateDiagnostics st_bayesrc_update_annotation_prior(
  const arma::mat& annotation_matrix,
  const arma::Row<int>& component,
  arma::mat& annot_alpha,
  arma::vec& annot_sigma_sq_alpha,
  const StBayesRCInterceptPrior& intercept_prior,
  double annot_sigma_sq_alpha_a,
  double annot_sigma_sq_alpha_b,
  std::mt19937& gen
) {
 return st_bayesrc_update_tempered_annotation_prior(
  annotation_matrix, component, annot_alpha, annot_sigma_sq_alpha,
  intercept_prior, annot_sigma_sq_alpha_a, annot_sigma_sq_alpha_b,
  intercept_prior.mean, 1.0, gen);
}

#endif
