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

inline arma::mat st_bayesrc_compute_snp_pi(
  const arma::mat& annotation_matrix,
  const arma::mat& annot_alpha,
  double pi_floor
) {
 const int m = static_cast<int>(annotation_matrix.n_rows);
 const int nstep = static_cast<int>(annot_alpha.n_cols);
 const int ncomponent = nstep + 1;

 arma::mat p(m, nstep, arma::fill::zeros);
 arma::mat snp_pi(m, ncomponent, arma::fill::zeros);

 for (int j = 0; j < nstep; ++j) {
  arma::vec eta = annotation_matrix * annot_alpha.col(static_cast<arma::uword>(j));
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
 arma::vec mean;
 arma::vec precision;
};

inline StBayesRCInterceptPrior st_bayesrc_parse_intercept_prior(
  const arma::mat& resolved,
  int nstep
) {
 if ((resolved.n_rows != 3 && resolved.n_rows != 4) ||
     resolved.n_cols != static_cast<arma::uword>(nstep) ||
     !resolved.is_finite()) {
  throw std::runtime_error(
   "BayesRC resolved intercept prior must be a finite 3 or 4 by nstep matrix.");
 }
 const bool legacy_flat = resolved(0, 0) == 1.0;
 const bool update_variance = resolved.n_rows == 3 || resolved(3, 0) == 1.0;
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
 }
 return StBayesRCInterceptPrior{
  legacy_flat, update_variance, resolved.row(1).t(), resolved.row(2).t()
 };
}

struct StBayesRCAnnotationUpdateDiagnostics {
 std::vector<int> eligible;
 std::vector<int> continuation;
 std::vector<int> prior_only;
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
 const int m = static_cast<int>(annotation_matrix.n_rows);
 const int n_annot = static_cast<int>(annotation_matrix.n_cols);
 const int nstep = static_cast<int>(annot_alpha.n_cols);

 arma::Mat<int> step_indicators(m, nstep, arma::fill::zeros);
 st_bayesrc_build_step_indicators(component, step_indicators);
 StBayesRCAnnotationUpdateDiagnostics diagnostics{
  std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<int>(static_cast<std::size_t>(nstep), 0),
  std::vector<int>(static_cast<std::size_t>(nstep), 0)
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
   double s = 0.0;
   for (int k = 0; k < n_annot; ++k) {
    s += annotation_matrix(static_cast<arma::uword>(i), static_cast<arma::uword>(k)) *
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

  arma::vec resid = latent - mu;

  for (int k = 0; k < n_annot; ++k) {
   const bool intercept = k == 0;
   double diag_k = 0.0;
   double rhs = 0.0;
   const double old = annot_alpha(static_cast<arma::uword>(k),
                                  static_cast<arma::uword>(j));

   for (int ii = 0; ii < nj; ++ii) {
    const int i = idx[static_cast<std::size_t>(ii)];
    const double x = annotation_matrix(static_cast<arma::uword>(i),
                                       static_cast<arma::uword>(k));
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
     throw std::runtime_error(
      "BayesRC non-empty stick has an unavailable intercept column.");
    }
    const double sig = annot_sigma_sq_alpha(static_cast<arma::uword>(j));
    if (!std::isfinite(sig) || sig <= 0.0) {
     throw std::runtime_error("BayesRC sigmaSqAlpha must remain positive and finite.");
    }
    std::normal_distribution<double> coefficient_norm(0.0, std::sqrt(sig));
    annot_alpha(static_cast<arma::uword>(k), static_cast<arma::uword>(j)) =
     coefficient_norm(gen);
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
     const double x = annotation_matrix(static_cast<arma::uword>(i),
                                        static_cast<arma::uword>(k));
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

#endif
