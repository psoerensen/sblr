// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

#include "cpg_samplers.h"
#include "distributions.h"

#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>

using namespace arma;


inline void sampleBetaCPG_Mt_latent_fast_arma(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const std::vector<std::vector<double>>& ww,
  arma::mat& r,
  arma::mat& beta,
  arma::mat& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 arma::vec rhs(nt), rhs_base(nt), z(nt), beta_new(nt), b_new(nt);
 arma::mat C(nt, nt), L(nt, nt);

 std::vector<double> loglik(nmodels), w(nmodels);

 std::uniform_real_distribution<double> runif(0.0,1.0);
 std::normal_distribution<double> norm(0.0,1.0);

 // rhs_base
 for (int t=0; t<nt; ++t)
  rhs_base(t) = Ei(t,t) * (r(t,i) + ww[t][i] * b(t,i));

 // ----- MODEL LOOP -----
 for (int k=0; k<nmodels; ++k) {

  if (pi[k] <= 0.0) { loglik[k] = -INFINITY; continue; }

  rhs.zeros();
  C = Bi;

  const auto& mk = models[k];

  for (int t=0; t<nt; ++t) {
   if (mk[t]) {
    rhs(t) = rhs_base(t);
    C(t,t) += ww[t][i] * Ei(t,t);
   }
  }

  if (!arma::chol(L,C,"lower")) { loglik[k] = -INFINITY; continue; }

  double logdet = 2.0 * arma::sum(log(L.diag()));

  arma::vec y = solve(trimatl(L), rhs);
  arma::vec mean = solve(trimatu(L.t()), y);

  double quad = dot(rhs, mean);

  loglik[k] = std::log(pi[k]) - 0.5*logdet + 0.5*quad;
 }

 // ----- SAMPLE MODEL -----
 double max_log = *std::max_element(loglik.begin(), loglik.end());

 double sum=0.0;
 for (int k=0;k<nmodels;k++){
  w[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k]-max_log) : 0.0;
  sum += w[k];
 }

 double u = runif(gen)*sum;

 int mselect = nmodels-1;
 double cum=0.0;
 for (int k=0;k<nmodels;k++){
  cum += w[k];
  if (u<=cum){ mselect=k; break; }
 }

 cmodel[mselect] += 1.0;

 const auto& msel = models[mselect];

 for (int t=0;t<nt;t++) d[t][i]=msel[t];

 // ----- BUILD SELECTED MODEL -----
 rhs.zeros();
 C = Bi;

 for (int t=0;t<nt;t++){
  if (msel[t]){
   rhs(t)=rhs_base(t);
   C(t,t)+=ww[t][i]*Ei(t,t);
  }
 }

 arma::chol(L,C,"lower");

 arma::vec y = solve(trimatl(L), rhs);
 arma::vec mean = solve(trimatu(L.t()), y);

 for (int t=0;t<nt;t++) z(t)=norm(gen);

 beta_new = mean + solve(trimatu(L.t()), z);

 // effective effects
 b_new.zeros();
 for (int t=0;t<nt;t++)
  if (msel[t]) b_new(t)=beta_new(t);

  // ----- RESIDUAL UPDATE -----
  const auto& idx = XXindices[i];
  const size_t nnz = idx.size();

  for (int t=0;t<nt;t++){
   double diff = b_new(t) - b(t,i);

   if (diff!=0.0){
    const auto& val = XXvalues[t][i];
    for (size_t j=0;j<nnz;j++){
     r(t, idx[j]) -= val[j]*diff;
    }
   }
  }

  // store
  beta.col(i) = beta_new;
  b.col(i)    = b_new;
}



// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtblr_cpg_arma(
  std::vector<std::vector<double>> wy,
  std::vector<std::vector<double>> ww,
  std::vector<double> yy,
  std::vector<std::vector<double>> b,
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
 // --------------------------------------------------
 // Dimensions
 // --------------------------------------------------
 const int nt = static_cast<int>(wy.size());
 const int m  = static_cast<int>(wy[0].size());
 const int nmodels = static_cast<int>(models.size());
 double nsamples = 0.0;

 // --------------------------------------------------
 // Convert dense inputs to arma objects
 // --------------------------------------------------
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
   b_mat(t, i)  = b[t][i];
   r_mat(t, i)  = wy[t][i];
  }
  for (int t2 = 0; t2 < nt; ++t2) {
   ssb_prior_mat(t, t2) = ssb_prior[t][t2];
   sse_prior_mat(t, t2) = sse_prior[t][t2];
  }
 }

 // --------------------------------------------------
 // Allocate remaining state / summaries
 // --------------------------------------------------
 std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));

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

 std::vector<double> x2t(m, 0.0);
 std::vector<int> order(m);

 // --------------------------------------------------
 // Initialize residuals: r = wy - W'W b
 // --------------------------------------------------
 for (int i = 0; i < m; ++i) {
  for (int t = 0; t < nt; ++t) {
   if (b_mat(t, i) != 0.0) {
    const size_t nnz = XXindices[i].size();
    const std::vector<double>& val = XXvalues[t][i];
    for (size_t j = 0; j < nnz; ++j) {
     r_mat(t, XXindices[i][j]) -= val[j] * b_mat(t, i);
    }
   }
  }
 }

 // --------------------------------------------------
 // Ranking statistic
 // --------------------------------------------------
 for (int i = 0; i < m; ++i) {
  double best = 0.0;
  for (int t = 0; t < nt; ++t) {
   const double x2 = (wy_mat(t, i) / ww_mat(t, i)) * (wy_mat(t, i) / ww_mat(t, i));
   if (x2 > best) best = x2;
  }
  x2t[i] = best;
 }

 std::iota(order.begin(), order.end(), 0);
 std::sort(order.begin(), order.end(),
           [&](int i1, int i2) { return x2t[i1] > x2t[i2]; });

 // --------------------------------------------------
 // Marker -> set lookup
 // --------------------------------------------------
 std::vector<int> marker_to_set(m, -1);
 for (size_t s = 0; s < sets.size(); ++s) {
  for (int idx : sets[s]) {
   marker_to_set[idx] = static_cast<int>(s);
  }
 }

 // --------------------------------------------------
 // Precision / covariance objects
 // --------------------------------------------------
 arma::mat Bi(nt, nt, arma::fill::zeros);
 arma::mat Ei(nt, nt, arma::fill::zeros);
 arma::mat G(nt, nt, arma::fill::zeros);

 if (!arma::inv_sympd(Bi, B)) {
  throw std::runtime_error("Initial Bi inversion failed.");
 }
 if (!arma::inv_sympd(Ei, E)) {
  throw std::runtime_error("Initial Ei inversion failed.");
 }

 std::mt19937 gen(seed);

 // --------------------------------------------------
 // Gibbs sampler
 // --------------------------------------------------
 for (int it = 0; it < nit + nburn; ++it) {

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

  // -----------------------------------------------
  // Sample marker effects (BayesCP-G latent)
  // -----------------------------------------------
  if (method == 4) {
   for (size_t s = 0; s < sets.size(); ++s) {

    // NOTE:
    // This preserves your current logic:
    // B is sampled once per set.
    // If this was not intended, move it outside this loop.
    sampleB_cpg_arma(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);

    if (!arma::inv_sympd(Bi, B)) {
     throw std::runtime_error("Bi inversion failed.");
    }

    for (int isort = 0; isort < m; ++isort) {
     const int i = order[isort];

     if (marker_to_set[i] == static_cast<int>(s)) {
      sampleBetaCPG_Mt_latent_fast_arma(
       i, nt, nmodels,
       models,
       cmodel, pi,
       Ei, Bi,
       ww,
       r_mat, beta_mat, b_mat, d,
       XXindices, XXvalues,
       gen
      );
     }
    }
   }
  }

  // -----------------------------------------------
  // Store marker summaries
  // -----------------------------------------------
  for (int t = 0; t < nt; ++t) {
   for (int i = 0; i < m; ++i) {
    if (d[t][i] > 0) {
     if ((it > nburn) && (it % nthin == 0)) {
      dm[t][i] += 1.0;
      bm[t][i] += b_mat(t, i);
     }
    }
   }
  }

  // -----------------------------------------------
  // Sample pi
  // -----------------------------------------------
  if (updatePi && method == 4) {
   samplePi_cpg(cmodel, pi, gen);
   for (int k = 0; k < nmodels; ++k) {
    if (it > nburn) pis[k] += pi[k];
   }
  }

  // -----------------------------------------------
  // Sample marker covariance B
  // -----------------------------------------------
  if (updateB && method == 4) {
   sampleB_cpg_arma(nt, m, nub, B, beta_mat, ssb_prior_mat, gen);
   for (int t = 0; t < nt; ++t) {
    vbs[t][it] = B(t, t);
   }
   for (int t1 = 0; t1 < nt; ++t1) {
    for (int t2 = 0; t2 < nt; ++t2) {
     if (it > nburn) cvbm[t1][t2] += B(t1, t2);
    }
   }
  }

  // -----------------------------------------------
  // Update genetic covariance G
  // -----------------------------------------------
  computeG_cpg_arma(nt, m, b_mat, wy_mat, r_mat, n, G);
  for (int t = 0; t < nt; ++t) {
   vgs[t][it] = G(t, t);
  }
  for (int t1 = 0; t1 < nt; ++t1) {
   for (int t2 = 0; t2 < nt; ++t2) {
    if (it > nburn) cvgm[t1][t2] += G(t1, t2);
   }
  }

  // -----------------------------------------------
  // Sample residual covariance E
  // -----------------------------------------------
  if (updateE) {
   sampleE_cpg_arma(nt, m, nue, E, b_mat, wy_mat, r_mat, sse_prior_mat, yy_vec, n, gen);

   if (!arma::inv_sympd(Ei, E)) {
    throw std::runtime_error("Ei inversion failed.");
   }

   for (int t = 0; t < nt; ++t) {
    ves[t][it] = E(t, t);
   }
   for (int t1 = 0; t1 < nt; ++t1) {
    for (int t2 = 0; t2 < nt; ++t2) {
     if (it > nburn) cvem[t1][t2] += E(t1, t2);
    }
   }
  }

  if ((it > nburn) && (it % nthin == 0)) {
   nsamples += 1.0;
  }
 }

 // --------------------------------------------------
 // Build output
 // --------------------------------------------------
 std::vector<std::vector<std::vector<double>>> result(20);

 for (int k = 0; k < 20; ++k) result[k].resize(nt);

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

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   result[0][t][i] = bm[t][i] / nsamples;
   result[1][t][i] = dm[t][i] / nsamples;
   result[2][t][i] = wy_mat(t, i);
   result[3][t][i] = r_mat(t, i);
   result[4][t][i] = b_mat(t, i);
   result[5][t][i] = d[t][i];
   result[6][t][i] = order[i];
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
   result[10][t1][t2] = cvbm[t1][t2] / nsamples;
   result[11][t1][t2] = cvgm[t1][t2] / nsamples;
   result[12][t1][t2] = cvem[t1][t2] / nsamples;
   result[13][t1][t2] = B(t1, t2);
   result[14][t1][t2] = G(t1, t2);
   result[15][t1][t2] = E(t1, t2);
  }
 }

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < nmodels; ++i) {
   result[16][t][i] = pi[i];
   result[17][t][i] = pis[i] / nsamples;
  }
 }

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < 4; ++i) result[18][t][i] = 0.0;
  for (int i = 0; i < 2; ++i) result[19][t][i] = 0.0;
 }

 return result;
}

