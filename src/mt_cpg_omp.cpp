// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

#include "cpg_samplers.h"
#include "distributions.h"

#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;


inline void sampleB_cpg_arma_omp(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const arma::mat& beta,
  const arma::mat& ssb_prior,
  std::mt19937& gen)
{
 if (nub <= nt - 1) {
  throw std::runtime_error("sampleB_cpg_arma_omp: nub must be > nt - 1.");
 }

 if (!beta.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp: beta contains NaN/Inf.");
 }

 arma::mat S_post = ssb_prior + beta * beta.t();

 // Numerical cleanup
 S_post = 0.5 * (S_post + S_post.t());

 if (!S_post.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp: S_post contains NaN/Inf.");
 }

 // Ensure SPD before inverse-Wishart draw
 arma::mat L;
 bool ok = arma::chol(L, S_post, "lower");

 double jitter = 1e-10;
 int tries = 0;

 while (!ok && tries < 8) {
  S_post.diag() += jitter;
  S_post = 0.5 * (S_post + S_post.t());

  ok = arma::chol(L, S_post, "lower");

  jitter *= 10.0;
  ++tries;
 }

 if (!ok) {
  throw std::runtime_error("sampleB_cpg_arma_omp: S_post is not SPD after jitter.");
 }

 const unsigned int df_post = static_cast<unsigned int>(nub + m);

 B = rinvwishart_cpp(df_post, S_post, gen);

 B = 0.5 * (B + B.t());
 B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);

 if (!B.is_finite()) {
  throw std::runtime_error("sampleB_cpg_arma_omp: sampled B contains NaN/Inf.");
 }
}


inline void sampleE_cpg_arma_omp(
  int nt,
  int nue,
  arma::mat& E,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const arma::mat& sse_prior,
  const arma::vec& yy,
  const std::vector<int>& n,
  std::mt19937& gen)
{
 E.zeros();

 for (int t = 0; t < nt; ++t) {
  double sse = arma::dot(b.row(t), r.row(t) + wy.row(t));
  sse = yy(t) - sse;

  const double scale = sse + nue * sse_prior(t, t);

  if (!std::isfinite(scale) || scale <= 0.0) {
   throw std::runtime_error("sampleE_cpg_arma_omp: invalid residual scale.");
  }


  std::chi_squared_distribution<double> rchisq(n[t] + nue);
  const double chi2 = rchisq(gen);

  E(t, t) = scale / chi2;
 }
}


inline void computeG_cpg_arma_omp(
  int nt,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const std::vector<int>& n,
  arma::mat& G)
{
 const arma::mat Q = wy - r;

 G.set_size(nt, nt);

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1; t2 < nt; ++t2) {
   const double denom = std::sqrt((double)n[t1] * (double)n[t2]);

   double s12 = arma::dot(b.row(t1), Q.row(t2));

   if (t1 != t2) {
    const double s21 = arma::dot(b.row(t2), Q.row(t1));
    s12 = 0.5 * (s12 + s21);
   }

   const double gij = s12 / denom;
   G(t1, t2) = gij;
   G(t2, t1) = gij;
  }
 }
}


inline void sampleBetaCPG_Mt_latent_arma(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel_local,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const arma::mat& ww,
  arma::mat& r,        // ✅ LIVE
  arma::mat& beta,
  arma::mat& b,        // ✅ LIVE
  arma::Mat<int>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 arma::vec rhs(nt), rhs_base(nt), z(nt), beta_new(nt), b_new(nt);
 arma::mat C(nt, nt), L(nt, nt);

 std::vector<double> loglik(nmodels), w(nmodels);

 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm(0.0, 1.0);

 // -----------------------------------
 // Base RHS from CURRENT state
 // -----------------------------------
 for (int t = 0; t < nt; ++t) {
  rhs_base(t) = Ei(t, t) * (r(t, i) + ww(t, i) * b(t, i));
 }

 // -----------------------------------
 // Model loop
 // -----------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -INFINITY;
   continue;
  }

  rhs.zeros();
  C = Bi;

  const auto& mk = models[k];

  for (int t = 0; t < nt; ++t) {
   if (mk[t]) {
    rhs(t) = rhs_base(t);
    C(t, t) += ww(t, i) * Ei(t, t);
   }
  }

  if (!arma::chol(L, C, "lower")) {
   loglik[k] = -INFINITY;
   continue;
  }

  double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  arma::vec y    = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean = arma::solve(arma::trimatu(L.t()), y);

  double quad = arma::dot(rhs, mean);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // -----------------------------------
 // Sample model
 // -----------------------------------
 double max_log = *std::max_element(loglik.begin(), loglik.end());

 double sumw = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_log) : 0.0;
  sumw += w[k];
 }

 double u = runif(gen) * sumw;

 int mselect = nmodels - 1;
 double cum = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  cum += w[k];
  if (u <= cum) {
   mselect = k;
   break;
  }
 }

 cmodel_local[mselect] += 1.0;

 const auto& msel = models[mselect];

 // -----------------------------------
 // Draw beta
 // -----------------------------------
 rhs.zeros();
 C = Bi;

 for (int t = 0; t < nt; ++t) {
  d(t, i) = msel[t];
  if (msel[t]) {
   rhs(t) = rhs_base(t);
   C(t, t) += ww(t, i) * Ei(t, t);
  }
 }

 arma::chol(L, C, "lower");

 arma::vec y    = arma::solve(arma::trimatl(L), rhs);
 arma::vec mean = arma::solve(arma::trimatu(L.t()), y);

 for (int t = 0; t < nt; ++t) z(t) = norm(gen);

 beta_new = mean + arma::solve(arma::trimatu(L.t()), z);

 b_new.zeros();
 for (int t = 0; t < nt; ++t)
  if (msel[t]) b_new(t) = beta_new(t);

  // -----------------------------------
  // Residual update (CRITICAL)
  // -----------------------------------
  const auto& idx = XXindices[i];
  const size_t nnz = idx.size();

  for (int t = 0; t < nt; ++t) {

   double diff = b_new(t) - b(t, i);

   if (diff != 0.0) {
    const auto& val = XXvalues[t][i];
    for (size_t j = 0; j < nnz; ++j) {
     r(t, idx[j]) -= val[j] * diff;
    }
   }
  }

  // store
  beta.col(i) = beta_new;
  b.col(i)    = b_new;
}

inline void rebuild_residual_from_b(
  int nt,
  int m,
  const arma::mat& wy,
  const arma::mat& b,
  arma::mat& r,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues)
{
 r = wy;

 for (int i = 0; i < m; ++i) {
  const auto& idx = XXindices[i];
  const size_t nnz = idx.size();

  for (int t = 0; t < nt; ++t) {
   const double bi = b(t, i);
   if (bi == 0.0) continue;

   const auto& val = XXvalues[t][i];
   for (size_t j = 0; j < nnz; ++j) {
    r(t, idx[j]) -= val[j] * bi;
   }
  }
 }
}

std::vector<std::vector<int>> build_independent_sets(
  int m,
  const std::vector<std::vector<int>>& XXindices)
{
 // Step 1: build inverted index (row → markers)
 std::unordered_map<int, std::vector<int>> row_to_markers;

 for (int i = 0; i < m; ++i) {
  for (int idx : XXindices[i]) {
   row_to_markers[idx].push_back(i);
  }
 }

 // Step 2: build adjacency list
 std::vector<std::unordered_set<int>> adj(m);

 for (const auto& kv : row_to_markers) {
  const std::vector<int>& markers = kv.second;

  for (size_t a = 0; a < markers.size(); ++a) {
   for (size_t b = a + 1; b < markers.size(); ++b) {
    int i = markers[a];
    int j = markers[b];
    adj[i].insert(j);
    adj[j].insert(i);
   }
  }
 }

 // Step 3: greedy coloring
 std::vector<int> color(m, -1);

 for (int i = 0; i < m; ++i) {

  std::unordered_set<int> used;

  for (int j : adj[i]) {
   if (color[j] != -1)
    used.insert(color[j]);
  }

  int c = 0;
  while (used.count(c)) ++c;

  color[i] = c;
 }

 // Step 4: collect sets
 int max_color = *std::max_element(color.begin(), color.end());

 std::vector<std::vector<int>> sets(max_color + 1);

 for (int i = 0; i < m; ++i) {
  sets[color[i]].push_back(i);
 }

 return sets;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtblr_cpg_omp(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b_init,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<int>>& sets,
  arma::mat B,
  arma::mat E,
  std::vector<std::vector<double>> ssb_prior,
  std::vector<std::vector<double>> sse_prior,
  std::vector<std::vector<int>> models,
  std::vector<double> pi,
  double nub,
  double nue,
  bool updateB,
  bool updateE,
  bool updatePi,
  std::vector<int> n,
  int nit,
  int nburn,
  int nthin,
  int seed,
  int method)
{
 const int nt = wy.size();
 const int m  = wy[0].size();
 const int nmodels = models.size();

 double nsamples = 0.0;

 arma::mat wy_mat(nt, m, arma::fill::zeros);
 arma::mat ww_mat(nt, m, arma::fill::zeros);
 arma::mat b_mat(nt, m, arma::fill::zeros);
 arma::mat beta_mat(nt, m, arma::fill::zeros);
 arma::mat r_mat(nt, m, arma::fill::zeros);

 arma::vec yy_vec(nt, arma::fill::zeros);
 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 arma::mat sse_prior_mat(nt, nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  yy_vec(t) = yy[t];

  for (int i = 0; i < m; ++i) {
   wy_mat(t, i) = wy[t][i];
   ww_mat(t, i) = ww[t][i];
   b_mat(t, i)  = b_init[t][i];
  }

  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }
 }

 // Initial residual: r = wy - XX b
 rebuild_residual_from_b(nt, m, wy_mat, b_mat, r_mat, XXindices, XXvalues);

 arma::Mat<int> d_mat(nt, m, arma::fill::zeros);

 arma::mat Bi, Ei, G(nt, nt, arma::fill::zeros);

 if (!arma::inv_sympd(Bi, B)) {
  throw std::runtime_error("Initial Bi inversion failed.");
 }

 if (!arma::inv_sympd(Ei, E)) {
  throw std::runtime_error("Initial Ei inversion failed.");
 }

 std::vector<double> cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));

 std::vector<std::vector<double>> ves(nt, std::vector<double>(nit + nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit + nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit + nburn, 0.0));

 std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));

 std::mt19937 gen(seed);

 std::vector<std::vector<int>> auto_sets;

 // if (sets.empty()) {
 //  auto_sets = build_independent_sets(m, XXindices);
 // } else {
 //  auto_sets = sets;
 // }
 auto_sets = build_independent_sets(m, XXindices);
 // Dense/conflicted LD: too many colors -> sequential is faster
 const bool use_serial =
  auto_sets.size() > static_cast<size_t>(m / 4);

#ifdef _OPENMP
 omp_set_num_threads(std::min(omp_get_max_threads(), 8));
#endif

 for (int it = 0; it < nit + nburn; ++it) {

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

  // ----------------------------
  // Sample B / E
  // ----------------------------
  if (updateB && method == 4) {
   sampleB_cpg_arma_omp(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);

   if (!arma::inv_sympd(Bi, B)) {
    throw std::runtime_error("Bi inversion failed.");
   }
  }

  if (updateE) {
   sampleE_cpg_arma_omp(
    nt, nue, E,
    b_mat, wy_mat, r_mat,
    sse_prior_mat, yy_vec, n,
    gen
   );

   if (!arma::inv_sympd(Ei, E)) {
    throw std::runtime_error("Ei inversion failed.");
   }
  }

  // ----------------------------
  // Marker updates
  // Sequential over sets, parallel within each set.
  // Requires markers inside each set to be residual-write independent.
  // ----------------------------
  if (method == 4) {

   if (use_serial) {

    for (int i = 0; i < m; i++) {

     // Only update if marker i is in the current set
      sampleBetaCPG_Mt_latent_arma(
       i,
       nt,
       nmodels,
       models,
       cmodel,
       pi,
       Ei,
       Bi,
       ww_mat,
       r_mat,
       beta_mat,
       b_mat,
       d_mat,
       XXindices,
       XXvalues,
       gen
      );

    }


   } else {
#ifdef _OPENMP
#pragma omp parallel
#endif
{
#ifdef _OPENMP
 std::mt19937 gen_local(
   seed + 100000 * it + omp_get_thread_num()
 );
#else
 std::mt19937 gen_local(
   seed + 100000 * it
 );
#endif

 std::vector<double> cmodel_local(nmodels, 0.0);

 // Loop over independent sets (sequential)
 for (size_t s = 0; s < auto_sets.size(); ++s) {

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
  for (int kk = 0; kk < (int)auto_sets[s].size(); ++kk) {

   const int i = auto_sets[s][kk];

   sampleBetaCPG_Mt_latent_arma(
    i,
    nt,
    nmodels,
    models,
    cmodel_local,
    pi,
    Ei,
    Bi,
    ww_mat,
    r_mat,
    beta_mat,
    b_mat,
    d_mat,
    XXindices,
    XXvalues,
    gen_local
   );
  }
 }

#ifdef _OPENMP
#pragma omp critical
#endif
{
 for (int k = 0; k < nmodels; ++k) {
  cmodel[k] += cmodel_local[k];
 }
}
}
   }


  }


  // ----------------------------
  // Sample pi
  // ----------------------------
  if (updatePi && method == 4) {
   samplePi_cpg(cmodel, pi, gen);

   if (it > nburn) {
    for (int k = 0; k < nmodels; ++k) {
     pis[k] += pi[k];
    }
   }
  }

  // ----------------------------
  // Update G
  // ----------------------------
  computeG_cpg_arma_omp(nt, b_mat, wy_mat, r_mat, n, G);

  for (int t = 0; t < nt; ++t) {
   vbs[t][it] = B(t, t);
   ves[t][it] = E(t, t);
   vgs[t][it] = G(t, t);
  }

  if (it > nburn) {
   for (int t1 = 0; t1 < nt; ++t1) {
    for (int t2 = 0; t2 < nt; ++t2) {
     cvbm[t1][t2] += B(t1, t2);
     cvem[t1][t2] += E(t1, t2);
     cvgm[t1][t2] += G(t1, t2);
    }
   }
  }

  if ((it > nburn) && (it % nthin == 0)) {
   nsamples += 1.0;

   for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
     if (d_mat(t, i) > 0) {
      dm[t][i] += 1.0;
      bm[t][i] += b_mat(t, i);
     }
    }
   }
  }
 }

 // ----------------------------
 // Build result
 // ----------------------------
 std::vector<std::vector<std::vector<double>>> result(20);

 for (int k = 0; k < 20; ++k) {
  result[k].resize(nt);
 }

 for (int t = 0; t < nt; ++t) {
  result[0][t].resize(m);
  result[1][t].resize(m);
  result[2][t].resize(m);
  result[3][t].resize(m);
  result[4][t].resize(m);
  result[5][t].resize(m);
  result[6][t].resize(m);

  result[7][t].resize(nit + nburn);
  result[8][t].resize(nit + nburn);
  result[9][t].resize(nit + nburn);

  result[10][t].resize(nt);
  result[11][t].resize(nt);
  result[12][t].resize(nt);
  result[13][t].resize(nt);
  result[14][t].resize(nt);
  result[15][t].resize(nt);

  result[16][t].resize(nmodels);
  result[17][t].resize(nmodels);

  result[18][t].resize(4);
  result[19][t].resize(2);
 }

 const double denom = std::max(nsamples, 1.0);

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   result[0][t][i] = bm[t][i] / denom;
   result[1][t][i] = dm[t][i] / denom;
   result[2][t][i] = wy_mat(t, i);
   result[3][t][i] = r_mat(t, i);
   result[4][t][i] = b_mat(t, i);
   result[5][t][i] = d_mat(t, i);
   result[6][t][i] = i;
  }
 }

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < nit + nburn; ++i) {
   result[7][t][i] = vbs[t][i];
   result[8][t][i] = vgs[t][i];
   result[9][t][i] = ves[t][i];
  }
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = 0; t2 < nt; ++t2) {
   result[10][t1][t2] = cvbm[t1][t2] / denom;
   result[11][t1][t2] = cvgm[t1][t2] / denom;
   result[12][t1][t2] = cvem[t1][t2] / denom;
   result[13][t1][t2] = B(t1, t2);
   result[14][t1][t2] = G(t1, t2);
   result[15][t1][t2] = E(t1, t2);
  }
 }

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < nmodels; ++i) {
   result[16][t][i] = pi[i];
   result[17][t][i] = pis[i] / denom;
  }
 }

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < 4; ++i) {
   result[18][t][i] = 0.0;
  }
  for (int i = 0; i < 2; ++i) {
   result[19][t][i] = 0.0;
  }
 }

 return result;
}

// // [[Rcpp::export]]
// std::vector<std::vector<std::vector<double>>> mtblr_cpg_omp(
//   std::vector<std::vector<double>> wy,
//   std::vector<std::vector<double>> ww,
//   std::vector<double> yy,
//   std::vector<std::vector<double>> b_init,
//   const std::vector<std::vector<std::vector<double>>>& XXvalues,
//   const std::vector<std::vector<int>>& XXindices,
//   const std::vector<std::vector<int>>& sets,
//   arma::mat B,
//   arma::mat E,
//   std::vector<std::vector<double>> ssb_prior,
//   std::vector<std::vector<double>> sse_prior,
//   std::vector<std::vector<int>> models,
//   std::vector<double> pi,
//   double nub,
//   double nue,
//   bool updateB,
//   bool updateE,
//   bool updatePi,
//   std::vector<int> n,
//   int nit,
//   int nburn,
//   int nthin,
//   int seed,
//   int method)
// {
//  const int nt = wy.size();
//  const int m  = wy[0].size();
//  const int nmodels = models.size();
//
//  double nsamples = 0.0;
//
//  // ----------------------------
//  // Convert to arma
//  // ----------------------------
//  arma::mat wy_mat(nt, m, fill::zeros);
//  arma::mat ww_mat(nt, m, fill::zeros);
//  arma::mat b_mat(nt, m, fill::zeros);
//  arma::mat beta_mat(nt, m, fill::zeros);
//  arma::mat r_mat(nt, m, fill::zeros);
//
//  arma::vec yy_vec(nt, fill::zeros);
//  arma::mat ssb_prior_mat(nt, nt, fill::zeros);
//  arma::mat sse_prior_mat(nt, nt, fill::zeros);
//
//  for (int t = 0; t < nt; ++t) {
//   yy_vec(t) = yy[t];
//   for (int i = 0; i < m; ++i) {
//    wy_mat(t, i) = wy[t][i];
//    ww_mat(t, i) = ww[t][i];
//    b_mat(t, i)  = b_init[t][i];
//   }
//   for (int t2 = 0; t2 < nt; ++t2) {
//    ssb_prior_mat(t, t2) = ssb_prior[t][t2];
//    sse_prior_mat(t, t2) = sse_prior[t][t2];
//   }
//  }
//
//  // initial residual
//  rebuild_residual_from_b(nt, m, wy_mat, b_mat, r_mat, XXindices, XXvalues);
//
//  arma::Mat<int> d_mat(nt, m, fill::zeros);
//
//  arma::mat Bi, Ei, G(nt, nt, fill::zeros);
//  if (!arma::inv_sympd(Bi, B)) throw std::runtime_error("Initial Bi inversion failed.");
//  if (!arma::inv_sympd(Ei, E)) throw std::runtime_error("Initial Ei inversion failed.");
//
//  std::vector<double> cmodel(nmodels), pis(nmodels, 0.0);
//
//  // summaries
//  std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
//  std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
//  std::vector<std::vector<double>> ves(nt, std::vector<double>(nit + nburn, 0.0));
//  std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit + nburn, 0.0));
//  std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit + nburn, 0.0));
//  std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
//  std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
//  std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
//
//  std::mt19937 gen(seed);
//
//  for (int it = 0; it < nit + nburn; ++it) {
//
//   std::fill(cmodel.begin(), cmodel.end(), 1.0);
//
//   // ----------------------------
//   // Sample B / E from current state
//   // ----------------------------
//   if (updateB && method == 4) {
//    sampleB_cpg_arma_omp(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);
//    if (!arma::inv_sympd(Bi, B)) throw std::runtime_error("Bi inversion failed.");
//   }
//
//   if (updateE) {
//    sampleE_cpg_arma_omp(nt, nue, E, b_mat, wy_mat, r_mat, sse_prior_mat, yy_vec, n, gen);
//    if (!arma::inv_sympd(Ei, E)) throw std::runtime_error("Ei inversion failed.");
//   }
//
//   // ----------------------------
//   // Freeze state for parallel block updates
//   // ----------------------------
//   const arma::mat r_old    = r_mat;
//   const arma::mat b_old    = b_mat;
//   const arma::mat beta_old = beta_mat; // optional if needed later
//
//   // ----------------------------
//   // Parallel block updates
//   // ----------------------------
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
// {
// #ifdef _OPENMP
//  std::mt19937 gen_local(seed + 1000 * it + omp_get_thread_num());
// #else
//  std::mt19937 gen_local(seed + 1000 * it);
// #endif
//  std::vector<double> cmodel_local(nmodels, 0.0);
//
// #ifdef _OPENMP
// #pragma omp for schedule(dynamic)
// #endif
//  for (size_t s = 0; s < sets.size(); ++s) {
//   for (int i : sets[s]) {
//    sampleBetaCPG_Mt_latent_blocked(
//     i, nt, nmodels,
//     models,
//     cmodel_local, pi,
//     Ei, Bi,
//     ww_mat,
//     r_old, b_old,
//     beta_mat, b_mat, d_mat,
//     gen_local
//    );
//   }
//  }
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//  for (int k = 0; k < nmodels; ++k) {
//   cmodel[k] += cmodel_local[k];
//  }
// }
// }
//
// // ----------------------------
// // Rebuild residual AFTER all block updates
// // ----------------------------
// rebuild_residual_from_b(nt, m, wy_mat, b_mat, r_mat, XXindices, XXvalues);
//
// // ----------------------------
// // Sample pi
// // ----------------------------
// if (updatePi && method == 4) {
//  samplePi_cpg(cmodel, pi, gen);
//  if (it > nburn) {
//   for (int k = 0; k < nmodels; ++k) pis[k] += pi[k];
//  }
// }
//
// // ----------------------------
// // Update G
// // ----------------------------
// computeG_cpg_arma_omp(nt, b_mat, wy_mat, r_mat, n, G);
//
// for (int t = 0; t < nt; ++t) {
//  vbs[t][it] = B(t, t);
//  ves[t][it] = E(t, t);
//  vgs[t][it] = G(t, t);
// }
//
// if (it > nburn) {
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    cvbm[t1][t2] += B(t1, t2);
//    cvem[t1][t2] += E(t1, t2);
//    cvgm[t1][t2] += G(t1, t2);
//   }
//  }
// }
//
// if ((it > nburn) && (it % nthin == 0)) {
//  nsamples += 1.0;
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    if (d_mat(t, i) > 0) {
//     dm[t][i] += 1.0;
//     bm[t][i] += b_mat(t, i);
//    }
//   }
//  }
// }
//  }
//
//  // ----------------------------
//  // Build result in your old format
//  // ----------------------------
//  std::vector<std::vector<std::vector<double>>> result(20);
//  for (int k = 0; k < 20; ++k) result[k].resize(nt);
//
//  for (int t = 0; t < nt; ++t) {
//   result[0][t].resize(m);
//   result[1][t].resize(m);
//   result[2][t].resize(m);
//   result[3][t].resize(m);
//   result[4][t].resize(m);
//   result[5][t].resize(m);
//   result[6][t].resize(m);
//   result[7][t].resize(nit + nburn);
//   result[8][t].resize(nit + nburn);
//   result[9][t].resize(nit + nburn);
//   result[10][t].resize(nt);
//   result[11][t].resize(nt);
//   result[12][t].resize(nt);
//   result[13][t].resize(nt);
//   result[14][t].resize(nt);
//   result[15][t].resize(nt);
//   result[16][t].resize(nmodels);
//   result[17][t].resize(nmodels);
//   result[18][t].resize(4);
//   result[19][t].resize(2);
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < m; ++i) {
//    result[0][t][i] = bm[t][i] / std::max(nsamples, 1.0);
//    result[1][t][i] = dm[t][i] / std::max(nsamples, 1.0);
//    result[2][t][i] = wy_mat(t, i);
//    result[3][t][i] = r_mat(t, i);
//    result[4][t][i] = b_mat(t, i);
//    result[5][t][i] = d_mat(t, i);
//    result[6][t][i] = i;
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < nit + nburn; ++i) {
//    result[7][t][i] = vbs[t][i];
//    result[8][t][i] = vgs[t][i];
//    result[9][t][i] = ves[t][i];
//   }
//  }
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = 0; t2 < nt; ++t2) {
//    result[10][t1][t2] = cvbm[t1][t2] / std::max(nsamples, 1.0);
//    result[11][t1][t2] = cvgm[t1][t2] / std::max(nsamples, 1.0);
//    result[12][t1][t2] = cvem[t1][t2] / std::max(nsamples, 1.0);
//    result[13][t1][t2] = B(t1, t2);
//    result[14][t1][t2] = G(t1, t2);
//    result[15][t1][t2] = E(t1, t2);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < nmodels; ++i) {
//    result[16][t][i] = pi[i];
//    result[17][t][i] = pis[i] / std::max(nsamples, 1.0);
//   }
//  }
//
//  for (int t = 0; t < nt; ++t) {
//   for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
//   for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
//  }
//
//  return result;
// }

// inline void sampleBetaCPG_Mt_latent_blocked(
//   int i,
//   int nt,
//   int nmodels,
//   const std::vector<std::vector<int>>& models,
//   std::vector<double>& cmodel_local,
//   const std::vector<double>& pi,
//   const arma::mat& Ei,
//   const arma::mat& Bi,
//   const arma::mat& ww,
//   const arma::mat& r_old,
//   const arma::mat& b_old,
//   arma::mat& beta,
//   arma::mat& b,
//   arma::Mat<int>& d,
//   std::mt19937& gen)
// {
//  arma::vec rhs(nt), rhs_base(nt), z(nt), beta_new(nt), b_new(nt);
//  arma::mat C(nt, nt), L(nt, nt);
//
//  std::vector<double> loglik(nmodels), w(nmodels);
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  std::normal_distribution<double> norm(0.0, 1.0);
//
//  // -----------------------------------
//  // Base RHS from frozen state
//  // -----------------------------------
//  for (int t = 0; t < nt; ++t) {
//   rhs_base(t) = Ei(t, t) * (r_old(t, i) + ww(t, i) * b_old(t, i));
//  }
//
//  // -----------------------------------
//  // Model loop
//  // -----------------------------------
//  for (int k = 0; k < nmodels; ++k) {
//
//   if (pi[k] <= 0.0) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   rhs.zeros();
//   C = Bi;
//
//   const auto& mk = models[k];
//
//   for (int t = 0; t < nt; ++t) {
//    if (mk[t]) {
//     rhs(t) = rhs_base(t);
//     C(t, t) += ww(t, i) * Ei(t, t);
//    }
//   }
//
//   if (!arma::chol(L, C, "lower")) {
//    loglik[k] = -INFINITY;
//    continue;
//   }
//
//   const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
//
//   arma::vec y    = arma::solve(arma::trimatl(L), rhs);
//   arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//   const double quad = arma::dot(rhs, mean);
//
//   loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
//  }
//
//  // -----------------------------------
//  // Sample model
//  // -----------------------------------
//  const double max_log = *std::max_element(loglik.begin(), loglik.end());
//
//  double sumw = 0.0;
//  for (int k = 0; k < nmodels; ++k) {
//   w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_log) : 0.0;
//   sumw += w[k];
//  }
//
//  const double u = runif(gen) * sumw;
//
//  int mselect = nmodels - 1;
//  double cum = 0.0;
//  for (int k = 0; k < nmodels; ++k) {
//   cum += w[k];
//   if (u <= cum) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel_local[mselect] += 1.0;
//
//  const auto& msel = models[mselect];
//
//  // -----------------------------------
//  // Selected model draw
//  // -----------------------------------
//  rhs.zeros();
//  C = Bi;
//
//  for (int t = 0; t < nt; ++t) {
//   d(t, i) = msel[t];
//   if (msel[t]) {
//    rhs(t) = rhs_base(t);
//    C(t, t) += ww(t, i) * Ei(t, t);
//   }
//  }
//
//  arma::chol(L, C, "lower");
//
//  arma::vec y    = arma::solve(arma::trimatl(L), rhs);
//  arma::vec mean = arma::solve(arma::trimatu(L.t()), y);
//
//  for (int t = 0; t < nt; ++t) {
//   z(t) = norm(gen);
//  }
//
//  beta_new = mean + arma::solve(arma::trimatu(L.t()), z);
//
//  b_new.zeros();
//  for (int t = 0; t < nt; ++t) {
//   if (msel[t]) b_new(t) = beta_new(t);
//  }
//
//  // -----------------------------------
//  // Write only marker state
//  // -----------------------------------
//  beta.col(i) = beta_new;
//  b.col(i)    = b_new;
// }
