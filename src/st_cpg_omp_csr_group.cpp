// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"
#include "st_csr_common.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// =============================================================================
// STBLR summary-stat CSR with group-level annotation priors
// =============================================================================
//
// Exported function:
//
//   stblr_cpg_omp_csr_group_annot(...)
//
// Purpose:
//   Group-level annotation prior model for summary-statistic CSR BayesC.
//
// Model:
//   Each marker i belongs to exactly one group g_i in {0, ..., G-1}.
//
//   pi_g,t ~ Beta(pi_group_prior_a[g], pi_group_prior_b[g])
//   d_i,t | g_i ~ Bernoulli(pi_g_i,t)
//
//   b_i,t | d_i,t = 1, g_i ~ N(0, vb_t * lambda_g_i,t)
//
// where lambda_g is a group-specific variance multiplier.
//
// Options:
//   updatePi = TRUE:
//     updates group-specific pi_g using marker inclusion counts by group.
//
//   updateGroupVb = TRUE:
//     updates group-specific variance multipliers lambda_g by inverse-chi-square
//     style conditional updates, normalized to have weighted mean 1 across
//     groups. This keeps vb_t interpretable as the global scale.
//
//   updateB = TRUE:
//     updates global vb_t conditional on current lambda_g.
//
// This is intentionally a full-sweep CSR sampler. No sparse scheduling.
//
// Return layout, aligned with current CSR BLR conventions plus group outputs:
//
//   0  bm
//   1  dm
//   2  wy
//   3  r
//   4  b
//   5  d
//   6  marker index
//   7  vbs
//   8  vgs
//   9  ves
//   10 covb
//   11 covg
//   12 cove
//   13 vb
//   14 vg
//   15 ve
//   16 final global pi, c(pi0, pi1), where pi1 is marker-weighted mean pi_g
//   17 posterior mean global pi, c(pi0, pi1)
//   18 diagnostics: reserved/log-CPO placeholders + runtime, length 4
//   19 diagnostics: nsamples, n, length 2
//   20 vle trace
//   21 vld trace
//   22 posterior mean group pi, flattened ngroup x nt by trait slot
//   23 posterior mean group variance multiplier, flattened ngroup x nt by trait slot
//   24 posterior mean group inclusion counts, length ngroup by trait slot
//   25 group sizes, length ngroup by trait slot
//
// R indexing:
//   raw_fit[[23]] = group_pi
//   raw_fit[[24]] = group_vb_multiplier
//   raw_fit[[25]] = group_nincluded
//   raw_fit[[26]] = group_size
//
// Neutral model check:
//   If all markers have group = 1 in R, ngroup = 1, and group_vb_multiplier = 1,
//   this reduces to the ordinary CSR BayesC sampler, apart from harmless extra
//   bookkeeping.
// =============================================================================

inline double clamp_prob_group(double x) {
 if (!std::isfinite(x)) {
  throw std::runtime_error("clamp_prob_group: probability is NaN/Inf.");
 }
 return std::min(std::max(x, 1e-300), 1.0 - 1e-12);
}

inline double computeLE_ST_csr_group(
  int m,
  const arma::rowvec& b,
  const arma::rowvec& ww,
  int n
) {
 double vle = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double bi = b(iu);

  if (bi != 0.0) {
   vle += ww(iu) * bi * bi;
  }
 }

 return vle / static_cast<double>(n);
}

inline void sampleBetaC_ST_csr_group(
  int i,
  double pi1_i,
  double vb_t,
  double vb_mult_i,
  double vei_i,
  const arma::rowvec& ww,
  arma::rowvec& r,
  arma::rowvec& b,
  arma::Row<int>& d,
  const STLDCSR& ld,
  std::mt19937& gen
) {
 const arma::uword iu = static_cast<arma::uword>(i);
 const double wi = ww(iu);

 if (!std::isfinite(wi) || wi <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid ww value.");
 }
 if (!std::isfinite(vb_t) || vb_t <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid global vb.");
 }
 if (!std::isfinite(vb_mult_i) || vb_mult_i <= 0.0) {
  throw std::runtime_error("sampleBetaC_ST_csr_group: invalid group vb multiplier.");
 }

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm01(0.0, 1.0);

 pi1_i = clamp_prob_group(pi1_i);
 const double pi0_i = std::max(1.0 - pi1_i, 1e-300);
 const double vbi = std::max(vb_t * vb_mult_i, 1e-12);
 const double vei_safe = std::max(vei_i, 1e-300);

 const double score = r(iu) + wi * b(iu);
 const double denom = std::max(vei_safe + wi * vbi, 1e-300);

 const double logBF =
  0.5 * std::log(vei_safe / denom) +
  0.5 * score * score * vbi / (vei_safe * denom);

 const double logp1 = std::log(pi1_i) + logBF;
 const double logp0 = std::log(pi0_i);
 const double delta_log = logp0 - logp1;

 double p1 = 0.0;
 if (delta_log > 35.0) p1 = 0.0;
 else if (delta_log < -35.0) p1 = 1.0;
 else p1 = 1.0 / (1.0 + std::exp(delta_log));

 const int di = (runif(gen) < p1) ? 1 : 0;
 double b_new = 0.0;

 if (di == 1) {
  const double lhs = wi + vei_safe / vbi;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b(iu);

 if (diff != 0.0) {
  r(iu) -= wi * diff;

  const uint64_t start = ld.ptr[static_cast<std::size_t>(i)];
  const uint64_t end   = ld.ptr[static_cast<std::size_t>(i + 1)];

  for (uint64_t p = start; p < end; ++p) {
   const int j = ld.idx[static_cast<std::size_t>(p)];
   r(static_cast<arma::uword>(j)) -=
    static_cast<double>(ld.xij[static_cast<std::size_t>(p)]) * diff;
  }
 }

 b(iu) = b_new;
 d(iu) = di;
}

inline void sampleB_ST_csr_group(
  int m,
  double nub,
  double& vb,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  const arma::rowvec& group_vb_multiplier,
  double ssb_prior,
  std::mt19937& gen
) {
 double ssb_scaled = 0.0;
 double dfb = 0.0;

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   const int g = group(iu);
   const double mult = group_vb_multiplier(static_cast<arma::uword>(g));

   if (!std::isfinite(mult) || mult <= 0.0) {
    throw std::runtime_error("sampleB_ST_csr_group: invalid group vb multiplier.");
   }

   ssb_scaled += b(iu) * b(iu) / mult;
   dfb += 1.0;
  }
 }

 const double scale = ssb_scaled + nub * ssb_prior;
 if (!std::isfinite(scale) || scale <= 0.0) {
  throw std::runtime_error("sampleB_ST_csr_group: invalid scale.");
 }

 std::chi_squared_distribution<double> rchisq(dfb + nub);
 const double chi2 = std::max(rchisq(gen), 1e-300);
 vb = std::max(scale / chi2, 1e-12);
}

inline void samplePiGroups_ST_csr_group(
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  arma::rowvec& group_pi,
  const arma::rowvec& prior_a,
  const arma::rowvec& prior_b,
  int ngroup,
  std::mt19937& gen
) {
 arma::rowvec c1 = prior_a;
 arma::rowvec c0 = prior_b;

 if (static_cast<int>(prior_a.n_elem) != ngroup ||
     static_cast<int>(prior_b.n_elem) != ngroup ||
     static_cast<int>(group_pi.n_elem) != ngroup) {
  throw std::runtime_error("samplePiGroups_ST_csr_group: invalid group prior dimensions.");
 }

 for (arma::uword i = 0; i < d.n_elem; ++i) {
  const int g = group(i);
  if (g < 0 || g >= ngroup) {
   throw std::runtime_error("samplePiGroups_ST_csr_group: group index out of range.");
  }

  if (d(i) > 0) c1(static_cast<arma::uword>(g)) += 1.0;
  else          c0(static_cast<arma::uword>(g)) += 1.0;
 }

 for (int g = 0; g < ngroup; ++g) {
  const arma::uword gu = static_cast<arma::uword>(g);
  if (!std::isfinite(c1(gu)) || c1(gu) <= 0.0 ||
      !std::isfinite(c0(gu)) || c0(gu) <= 0.0) {
      throw std::runtime_error("samplePiGroups_ST_csr_group: invalid Beta posterior shape.");
  }

  std::gamma_distribution<double> rg0(c0(gu), 1.0);
  std::gamma_distribution<double> rg1(c1(gu), 1.0);

  const double g0 = std::max(rg0(gen), 1e-300);
  const double g1 = std::max(rg1(gen), 1e-300);
  group_pi(gu) = g1 / (g0 + g1);
 }
}

inline void sampleGroupVbMultipliers_ST_csr_group(
  int m,
  int ngroup,
  double nub_group,
  double ssb_group_prior,
  double vb_global,
  const arma::rowvec& b,
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  const arma::rowvec& group_size,
  arma::rowvec& group_vb_multiplier,
  bool normalize_group_vb,
  std::mt19937& gen
) {
 if (!std::isfinite(vb_global) || vb_global <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: invalid global vb.");
 }
 if (!std::isfinite(nub_group) || nub_group <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: nub_group must be positive.");
 }
 if (!std::isfinite(ssb_group_prior) || ssb_group_prior <= 0.0) {
  throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: ssb_group_prior must be positive.");
 }

 arma::rowvec ssb(ngroup, arma::fill::zeros);
 arma::rowvec df(ngroup, arma::fill::zeros);

 for (int i = 0; i < m; ++i) {
  const arma::uword iu = static_cast<arma::uword>(i);
  if (d(iu) > 0) {
   const int g = group(iu);
   if (g < 0 || g >= ngroup) {
    throw std::runtime_error("sampleGroupVbMultipliers_ST_csr_group: group index out of range.");
   }
   const arma::uword gu = static_cast<arma::uword>(g);
   ssb(gu) += b(iu) * b(iu) / vb_global;
   df(gu) += 1.0;
  }
 }

 for (int g = 0; g < ngroup; ++g) {
  const arma::uword gu = static_cast<arma::uword>(g);
  const double scale = ssb(gu) + nub_group * ssb_group_prior;
  std::chi_squared_distribution<double> rchisq(df(gu) + nub_group);
  const double chi2 = std::max(rchisq(gen), 1e-300);
  group_vb_multiplier(gu) = std::max(scale / chi2, 1e-12);
 }

 if (normalize_group_vb) {
  double weighted_mean = 0.0;
  double total_size = 0.0;

  for (int g = 0; g < ngroup; ++g) {
   const arma::uword gu = static_cast<arma::uword>(g);
   weighted_mean += group_size(gu) * group_vb_multiplier(gu);
   total_size += group_size(gu);
  }

  if (total_size > 0.0) weighted_mean /= total_size;

  if (std::isfinite(weighted_mean) && weighted_mean > 0.0) {
   group_vb_multiplier /= weighted_mean;
  }
 }
}

inline double marker_weighted_mean_group_value(
  const arma::rowvec& x,
  const arma::rowvec& group_size
) {
 if (x.n_elem != group_size.n_elem) {
  throw std::runtime_error("marker_weighted_mean_group_value: dimension mismatch.");
 }

 double s = 0.0;
 double n = 0.0;

 for (arma::uword g = 0; g < x.n_elem; ++g) {
  s += x(g) * group_size(g);
  n += group_size(g);
 }

 if (n <= 0.0) return NA_REAL;
 return s / n;
}

inline arma::rowvec count_group_inclusions(
  const arma::Row<int>& d,
  const arma::Row<int>& group,
  int ngroup
) {
 arma::rowvec out(ngroup, arma::fill::zeros);

 for (arma::uword i = 0; i < d.n_elem; ++i) {
  if (d(i) > 0) {
   const int g = group(i);
   if (g < 0 || g >= ngroup) {
    throw std::runtime_error("count_group_inclusions: group index out of range.");
   }
   out(static_cast<arma::uword>(g)) += 1.0;
  }
 }

 return out;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> stblr_cpg_omp_csr_group_annot(
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
  std::vector<int> group_index,
  int ngroup,
  std::vector<std::vector<double>> group_pi_init,
  std::vector<double> pi_group_prior_a,
  std::vector<double> pi_group_prior_b,
  std::vector<std::vector<double>> group_vb_multiplier_init,
  bool updateGroupVb,
  double nub_group,
  double ssb_group_prior,
  bool normalize_group_vb,
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
  int seed
) {
 const int nt = static_cast<int>(wy.size());

 if (nt <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: nt must be positive.");
 const int m = static_cast<int>(wy[0].size());
 if (m <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: m must be positive.");

 if (nit <= 0 || nburn < 0 || nthin <= 0) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: invalid nit/nburn/nthin.");
 }

 if (ngroup <= 0) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: ngroup must be positive.");
 if (static_cast<int>(group_index.size()) != m) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: group_index must have length m.");
 }

 if (pi.size() != 2) throw std::runtime_error("stblr_cpg_omp_csr_group_annot: pi must have length 2.");
 if ((int)ww.size() != nt || (int)b_init.size() != nt || (int)yy.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: inconsistent trait dimensions.");
 }
 if ((int)ssb_prior.size() != nt || (int)sse_prior.size() != nt) {
  throw std::runtime_error("stblr_cpg_omp_csr_group_annot: priors must be nt x nt.");
 }
 if ((int)B.n_rows != nt || (int)B.n_cols != nt) throw std::runtime_error("B must be nt x nt.");
 if ((int)E.n_rows != nt || (int)E.n_cols != nt) throw std::runtime_error("E must be nt x nt.");

 if (static_cast<int>(group_pi_init.size()) != nt ||
     static_cast<int>(group_vb_multiplier_init.size()) != nt) {
  throw std::runtime_error("group_pi_init and group_vb_multiplier_init must have length nt.");
 }
 for (int t = 0; t < nt; ++t) {
  if (static_cast<int>(group_pi_init[static_cast<std::size_t>(t)].size()) != ngroup ||
      static_cast<int>(group_vb_multiplier_init[static_cast<std::size_t>(t)].size()) != ngroup) {
   throw std::runtime_error("group initial vectors must have length ngroup for each trait.");
  }
 }
 if (static_cast<int>(pi_group_prior_a.size()) != ngroup ||
     static_cast<int>(pi_group_prior_b.size()) != ngroup) {
  throw std::runtime_error("pi_group_prior_a and pi_group_prior_b must have length ngroup.");
 }

 arma::Row<int> group(m, arma::fill::zeros);
 arma::rowvec group_size(ngroup, arma::fill::zeros);

 for (int i = 0; i < m; ++i) {
  const int g0 = group_index[static_cast<std::size_t>(i)];
  if (g0 < 0 || g0 >= ngroup) {
   throw std::runtime_error("group_index must be 0-based and in [0, ngroup-1].");
  }
  group(static_cast<arma::uword>(i)) = g0;
  group_size(static_cast<arma::uword>(g0)) += 1.0;
 }
 for (int g = 0; g < ngroup; ++g) {
  if (group_size(static_cast<arma::uword>(g)) <= 0.0) {
   throw std::runtime_error("Every group must contain at least one marker.");
  }
 }

 arma::rowvec prior_a(ngroup, arma::fill::zeros);
 arma::rowvec prior_b(ngroup, arma::fill::zeros);
 for (int g = 0; g < ngroup; ++g) {
  const double a = pi_group_prior_a[static_cast<std::size_t>(g)];
  const double b = pi_group_prior_b[static_cast<std::size_t>(g)];
  if (!std::isfinite(a) || a <= 0.0 || !std::isfinite(b) || b <= 0.0) {
   throw std::runtime_error("group Beta prior shapes must be finite and positive.");
  }
  prior_a(static_cast<arma::uword>(g)) = a;
  prior_b(static_cast<arma::uword>(g)) = b;
 }

 for (int t = 0; t < nt; ++t) {
  if ((int)wy[t].size() != m || (int)ww[t].size() != m || (int)b_init[t].size() != m) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot: inconsistent marker dimensions.");
  }
 }
 if (use_r_init) {
  if (static_cast<int>(r_init.size()) != nt) throw std::runtime_error("r_init must have length nt.");
  for (int t = 0; t < nt; ++t) if (static_cast<int>(r_init[t].size()) != m) throw std::runtime_error("r_init[t] must have length m.");
 }
 if (use_d_init) {
  if (static_cast<int>(d_init.size()) != nt) throw std::runtime_error("d_init must have length nt.");
  for (int t = 0; t < nt; ++t) if (static_cast<int>(d_init[t].size()) != m) throw std::runtime_error("d_init[t] must have length m.");
 }

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);
 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);
 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(static_cast<arma::uword>(t)) = yy[t];
  for (int i = 0; i < m; ++i) {
   wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = wy[t][i];
   ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = ww[t][i];
   b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) = b_init[t][i];
  }
  if ((int)ssb_prior[t].size() != nt || (int)sse_prior[t].size() != nt) {
   throw std::runtime_error("priors must be nt x nt.");
  }
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = ssb_prior[t][t2];
   sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t2)) = sse_prior[t][t2];
  }
 }

 for (int t = 1; t < nt; ++t) {
  if (n[t] != n[0]) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot: shared LD scaling assumes equal n across traits.");
  }
 }
 for (int t = 1; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   const double w0 = ww_mat(0, static_cast<arma::uword>(i));
   const double wt = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   const double tol = 1e-8 * std::max(1.0, std::abs(w0));
   if (!std::isfinite(w0) || !std::isfinite(wt) || w0 <= 0.0 || wt <= 0.0) {
    throw std::runtime_error("ww contains invalid value before LD pre-scaling.");
   }
   if (std::abs(w0 - wt) > tol) {
    throw std::runtime_error("ww differs across traits; pre-scaled shared ST LD is invalid.");
   }
  }
 }

 std::vector<double> xx(static_cast<std::size_t>(m), 0.0);
 for (int i = 0; i < m; ++i) {
  const double wi = ww_mat(0, static_cast<arma::uword>(i));
  if (!std::isfinite(wi) || wi <= 0.0) throw std::runtime_error("ww contains invalid value in trait 0.");
  xx[static_cast<std::size_t>(i)] = wi;
 }
 STLDCSR ld = read_and_build_st_ld_csr(ld_prefix, m, xx);

 std::vector<double> x2(static_cast<std::size_t>(m), 0.0);
 std::vector<int> order(static_cast<std::size_t>(m));
 for (int i = 0; i < m; ++i) {
  double best = 0.0;
  for (int t = 0; t < nt; ++t) {
   const double wi = ww_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   if (wi > 0.0) {
    const double bhat = wy_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i)) / wi;
    best = std::max(best, bhat * bhat);
   }
  }
  x2[static_cast<std::size_t>(i)] = best;
  order[static_cast<std::size_t>(i)] = i;
 }
 std::sort(order.begin(), order.end(), [&](int a, int b) {
  return x2[static_cast<std::size_t>(a)] > x2[static_cast<std::size_t>(b)];
 });

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);

 arma::mat group_pi_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_vb_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_nincluded_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_size_mat(nt, ngroup, arma::fill::zeros);

 std::vector<int> failed(static_cast<std::size_t>(nt), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(nt));
 std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
 std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);

 int nthreads = 1;
#ifdef _OPENMP
 omp_set_dynamic(0);
 nthreads = std::max(1, std::min(ncores, nt));
 omp_set_num_threads(nthreads);
 Rcpp::Rcout << "STBLR group annotation CSR OpenMP requested threads = "
             << nthreads << ", omp_get_max_threads = " << omp_get_max_threads()
             << ", num procs = " << omp_get_num_procs() << "\n";
#endif

 Rcpp::Rcout << "STBLR group annotation CSR: ngroup=" << ngroup
             << ", updatePi=" << updatePi
             << ", updateGroupVb=" << updateGroupVb
             << ", updateB=" << updateB
             << ", normalize_group_vb=" << normalize_group_vb
             << "\n";

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int t = 0; t < nt; ++t) {
#ifdef _OPENMP
  const double wall_start = omp_get_wtime();
  thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
#else
  const double wall_start = 0.0;
  thread_used[static_cast<std::size_t>(t)] = 0;
#endif

  try {
   std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));

   arma::rowvec r_t(m, arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   if (use_d_init) {
    for (int i = 0; i < m; ++i) d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
   } else {
    for (int i = 0; i < m; ++i) d_t(static_cast<arma::uword>(i)) = b_t(static_cast<arma::uword>(i)) != 0.0 ? 1 : 0;
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    if (!r_t.is_finite()) throw std::runtime_error("r_init contains NaN/Inf.");
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   arma::rowvec group_pi_t(ngroup, arma::fill::zeros);
   arma::rowvec group_vb_multiplier_t(ngroup, arma::fill::ones);
   for (int g = 0; g < ngroup; ++g) {
    const double pg = group_pi_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(g)];
    const double mg = group_vb_multiplier_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(g)];
    if (!std::isfinite(pg) || pg <= 0.0 || pg >= 1.0) throw std::runtime_error("group_pi_init contains invalid value.");
    if (!std::isfinite(mg) || mg <= 0.0) throw std::runtime_error("group_vb_multiplier_init contains invalid value.");
    group_pi_t(static_cast<arma::uword>(g)) = pg;
    group_vb_multiplier_t(static_cast<arma::uword>(g)) = mg;
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_group(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);

   arma::rowvec group_pi_accum(ngroup, arma::fill::zeros);
   arma::rowvec group_vb_accum(ngroup, arma::fill::zeros);
   arma::rowvec group_nincluded_accum(ngroup, arma::fill::zeros);

   double nsamples_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const arma::uword iu = static_cast<arma::uword>(i);
     const int g = group(iu);
     const arma::uword gu = static_cast<arma::uword>(g);

     sampleBetaC_ST_csr_group(
      i,
      group_pi_t(gu),
      vb_t,
      group_vb_multiplier_t(gu),
      vei_t,
      ww_t,
      r_t,
      b_t,
      d_t,
      ld,
      gen_t
     );
    }

    if (updateB) {
     sampleB_ST_csr_group(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      group,
      group_vb_multiplier_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );
    }

    if (updateGroupVb) {
     sampleGroupVbMultipliers_ST_csr_group(
      m,
      ngroup,
      nub_group,
      ssb_group_prior,
      vb_t,
      b_t,
      d_t,
      group,
      group_size,
      group_vb_multiplier_t,
      normalize_group_vb,
      gen_t
     );
    }

    if (updateE) {
     if (rebuild_r_before_updateE) rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
     sampleE_ST_csr(
      m,
      nue,
      ve_t,
      b_t,
      wy_t,
      r_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t)),
      n[t],
       gen_t
     );
    }

    if (updatePi) {
     samplePiGroups_ST_csr_group(d_t, group, group_pi_t, prior_a, prior_b, ngroup, gen_t);
    }

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_group(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;
    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vg_t)) throw std::runtime_error("vg became NaN/Inf.");
    if (!std::isfinite(vle_t)) throw std::runtime_error("vle became NaN/Inf.");
    if (!std::isfinite(vld_t)) throw std::runtime_error("vld became NaN/Inf.");
    if (!std::isfinite(vei_t) || vei_t <= 0.0) throw std::runtime_error("vei became invalid.");

    const double global_pi = marker_weighted_mean_group_value(group_pi_t, group_size);

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    pis_t(static_cast<arma::uword>(it)) = global_pi;
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;
     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }
     group_pi_accum += group_pi_t;
     group_vb_accum += group_vb_multiplier_t;
     group_nincluded_accum += count_group_inclusions(d_t, group, ngroup);
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;
   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   group_pi_accum /= nsamples_t;
   group_vb_accum /= nsamples_t;
   group_nincluded_accum /= nsamples_t;

   bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
   dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
   b_mat.row(static_cast<arma::uword>(t)) = b_t;
   r_mat.row(static_cast<arma::uword>(t)) = r_t;
   d_mat.row(static_cast<arma::uword>(t)) = d_t;
   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
   vles_mat.row(static_cast<arma::uword>(t)) = vles_t;
   vlds_mat.row(static_cast<arma::uword>(t)) = vlds_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_pi(static_cast<arma::uword>(t)) = marker_weighted_mean_group_value(group_pi_t, group_size);
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;

   group_pi_mean.row(static_cast<arma::uword>(t)) = group_pi_accum;
   group_vb_mean.row(static_cast<arma::uword>(t)) = group_vb_accum;
   group_nincluded_mean.row(static_cast<arma::uword>(t)) = group_nincluded_accum;
   group_size_mat.row(static_cast<arma::uword>(t)) = group_size;

#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif

  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = e.what();
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  } catch (...) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = "unknown error";
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  }
 }

#ifdef _OPENMP
 for (int t = 0; t < nt; ++t) {
  Rcpp::Rcout << "trait " << t << " used thread " << thread_used[static_cast<std::size_t>(t)]
              << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)] << "\n";
 }
#endif

 for (int t = 0; t < nt; ++t) {
  if (failed[static_cast<std::size_t>(t)]) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot failed for trait " + std::to_string(t) + ": " + errors[static_cast<std::size_t>(t)]);
  }
 }

 std::vector<std::vector<std::vector<double>>> result(26);
 for (int k = 0; k < 26; ++k) result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int k = 0; k <= 6; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(m));
  for (int k = 7; k <= 9; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nit + nburn));
  for (int k = 10; k <= 15; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nt));
  result[16][ts].resize(2);
  result[17][ts].resize(2);
  result[18][ts].resize(4);
  result[19][ts].resize(2);
  result[20][ts].resize(static_cast<std::size_t>(nit + nburn));
  result[21][ts].resize(static_cast<std::size_t>(nit + nburn));
  result[22][ts].resize(static_cast<std::size_t>(ngroup));
  result[23][ts].resize(static_cast<std::size_t>(ngroup));
  result[24][ts].resize(static_cast<std::size_t>(ngroup));
  result[25][ts].resize(static_cast<std::size_t>(ngroup));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int i = 0; i < m; ++i) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);
   const std::size_t is = static_cast<std::size_t>(i);
   result[0][ts][is] = bm_mat(tu, iu);
   result[1][ts][is] = dm_mat(tu, iu);
   result[2][ts][is] = wy_mat(tu, iu);
   result[3][ts][is] = r_mat(tu, iu);
   result[4][ts][is] = b_mat(tu, iu);
   result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
   result[6][ts][is] = static_cast<double>(i);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int it = 0; it < nit + nburn; ++it) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword itu = static_cast<arma::uword>(it);
   const std::size_t its = static_cast<std::size_t>(it);
   result[7][ts][its] = vbs_mat(tu, itu);
   result[8][ts][its] = vgs_mat(tu, itu);
   result[9][ts][its] = ves_mat(tu, itu);
   result[20][ts][its] = vles_mat(tu, itu);
   result[21][ts][its] = vlds_mat(tu, itu);
  }
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  const std::size_t t1s = static_cast<std::size_t>(t1);
  for (int t2 = 0; t2 < nt; ++t2) {
   const std::size_t t2s = static_cast<std::size_t>(t2);
   result[10][t1s][t2s] = 0.0;
   result[11][t1s][t2s] = 0.0;
   result[12][t1s][t2s] = 0.0;
   result[13][t1s][t2s] = 0.0;
   result[14][t1s][t2s] = 0.0;
   result[15][t1s][t2s] = 0.0;
  }
  result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
  result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
  result[16][ts][1] = final_pi(static_cast<arma::uword>(t));

  double mean_pi = 0.0;
  int npi = 0;
  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   ++npi;
  }
  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi(static_cast<arma::uword>(t));
  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  result[18][ts][0] = 0.0;
  result[18][ts][1] = 0.0;
  result[18][ts][2] = trait_seconds[static_cast<std::size_t>(t)];
  result[18][ts][3] = trait_seconds[static_cast<std::size_t>(t)];
  result[19][ts][0] = nsamples_vec(static_cast<arma::uword>(t));
  result[19][ts][1] = static_cast<double>(n[t]);

  for (int g = 0; g < ngroup; ++g) {
   const std::size_t gs = static_cast<std::size_t>(g);
   result[22][ts][gs] = group_pi_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[23][ts][gs] = group_vb_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[24][ts][gs] = group_nincluded_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[25][ts][gs] = group_size_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
  }
 }

 return result;
}
