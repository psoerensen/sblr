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

inline void st_bayesrc_update_annotation_prior(
  const arma::mat& annotation_matrix,
  const arma::Row<int>& component,
  arma::mat& annot_alpha,
  arma::vec& annot_sigma_sq_alpha,
  bool intercept_flat,
  double annot_sigma_sq_alpha_a,
  double annot_sigma_sq_alpha_b,
  std::mt19937& gen
) {
 const int m = static_cast<int>(annotation_matrix.n_rows);
 const int n_annot = static_cast<int>(annotation_matrix.n_cols);
 const int nstep = static_cast<int>(annot_alpha.n_cols);

 arma::Mat<int> step_indicators(m, nstep, arma::fill::zeros);
 st_bayesrc_build_step_indicators(component, step_indicators);

 for (int j = 0; j < nstep; ++j) {
  std::vector<int> idx;
  idx.reserve(static_cast<std::size_t>(m));

  for (int i = 0; i < m; ++i) {
   if (j == 0 || step_indicators(static_cast<arma::uword>(i),
                                  static_cast<arma::uword>(j - 1)) > 0) {
    idx.push_back(i);
   }
  }

  if (idx.empty()) continue;

  const int nj = static_cast<int>(idx.size());
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
   const bool flat_prior = intercept_flat && (k == 0);
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
   if (diag_k <= 0.0) continue;

   double inv_lhs;
   if (flat_prior) {
    inv_lhs = 1.0 / diag_k;
   } else {
    const double sig = std::max(
     annot_sigma_sq_alpha(static_cast<arma::uword>(j)), 1e-12
    );
    inv_lhs = 1.0 / (diag_k + 1.0 / sig);
   }

   const double mean = inv_lhs * rhs;
   const double sd = std::sqrt(inv_lhs);
   std::normal_distribution<double> norm(mean, sd);
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
   if (intercept_flat && k == 0) continue;
   const double ak = annot_alpha(static_cast<arma::uword>(k),
                                 static_cast<arma::uword>(j));
   ss += ak * ak;
   ++ncoef;
  }

  if (ncoef > 0) {
   const double df = static_cast<double>(ncoef) + annot_sigma_sq_alpha_a;
   const double scale = ss + annot_sigma_sq_alpha_b;
   std::chi_squared_distribution<double> rchisq(df);
   const double chi2 = std::max(rchisq(gen), 1e-300);
   annot_sigma_sq_alpha(static_cast<arma::uword>(j)) =
    std::max(scale / chi2, 1e-12);
  }
 }
}

#endif
