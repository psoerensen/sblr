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


void sampleBetaCPG_Mt_latent_fast(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const std::vector<std::vector<double>>& ww,
  std::vector<std::vector<double>>& r,
  std::vector<std::vector<double>>& beta,
  std::vector<std::vector<double>>& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // -----------------------------------
 // PREALLOCATE
 // -----------------------------------
 arma::vec rhs(nt);
 arma::vec rhs_base(nt);
 arma::mat C(nt, nt);
 arma::mat L(nt, nt);
 arma::vec z(nt);
 arma::vec beta_new(nt);
 arma::vec b_new(nt);

 std::vector<double> loglik(nmodels);
 std::vector<double> w(nmodels);  // weights

 // RNG
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 std::normal_distribution<double> norm(0.0, 1.0);

 // -----------------------------------
 // PRECOMPUTE rhs_base
 // -----------------------------------
 for (int t = 0; t < nt; ++t) {
  rhs_base(t) = Ei(t,t) * (r[t][i] + ww[t][i] * b[t][i]);
 }

 // -----------------------------------
 // MODEL LOOP
 // -----------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -INFINITY;
   continue;
  }

  rhs.zeros();
  C = Bi;

  const std::vector<int>& mk = models[k];

  for (int t = 0; t < nt; ++t) {
   if (mk[t]) {
    rhs(t) = rhs_base(t);
    C(t,t) += ww[t][i] * Ei(t,t);
   }
  }

  if (!arma::chol(L, C, "lower")) {
   loglik[k] = -INFINITY;
   continue;
  }

  // log|C|
  double logdet = 0.0;
  for (int t = 0; t < nt; ++t)
   logdet += std::log(L(t,t));
  logdet *= 2.0;

  // solve
  arma::vec y    = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean = arma::solve(arma::trimatu(L.t()), y);

  double quad = arma::dot(rhs, mean);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // -----------------------------------
 // SAMPLE MODEL (log-sum-exp)
 // -----------------------------------
 double max_log = *std::max_element(loglik.begin(), loglik.end());

 double sum = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   w[k] = std::exp(loglik[k] - max_log);
   sum += w[k];
  } else {
   w[k] = 0.0;
  }
 }

 double u = runif(gen) * sum;

 int mselect = nmodels - 1;
 double cum = 0.0;

 for (int k = 0; k < nmodels; ++k) {
  cum += w[k];
  if (u <= cum) {
   mselect = k;
   break;
  }
 }

 // -----------------------------------
 // STORE MODEL
 // -----------------------------------
 cmodel[mselect] += 1.0;

 const std::vector<int>& msel = models[mselect];

 for (int t = 0; t < nt; ++t)
  d[t][i] = msel[t];

 // -----------------------------------
 // BUILD SELECTED MODEL
 // -----------------------------------
 rhs.zeros();
 C = Bi;

 for (int t = 0; t < nt; ++t) {
  if (msel[t]) {
   rhs(t) = rhs_base(t);
   C(t,t) += ww[t][i] * Ei(t,t);
  }
 }

 arma::chol(L, C, "lower");

 arma::vec y    = arma::solve(arma::trimatl(L), rhs);
 arma::vec mean = arma::solve(arma::trimatu(L.t()), y);

 for (int t = 0; t < nt; ++t)
  z(t) = norm(gen);

 beta_new = mean + arma::solve(arma::trimatu(L.t()), z);

 // -----------------------------------
 // EFFECTIVE EFFECT
 // -----------------------------------
 b_new.zeros();
 for (int t = 0; t < nt; ++t)
  if (msel[t])
   b_new(t) = beta_new(t);

  // -----------------------------------
  // RESIDUAL UPDATE (hotspot)
  // -----------------------------------
  const size_t nnz = XXindices[i].size();

  for (int t = 0; t < nt; ++t) {
   double diff = b_new(t) - b[t][i];

   if (diff != 0.0) {
    const std::vector<int>& idx = XXindices[i];
    const std::vector<double>& val = XXvalues[t][i];

    for (size_t j = 0; j < nnz; ++j) {
     r[t][idx[j]] -= val[j] * diff;
    }
   }
  }

  // -----------------------------------
  // STORE RESULTS
  // -----------------------------------
  for (int t = 0; t < nt; ++t) {
   beta[t][i] = beta_new(t);
   b[t][i]    = b_new(t);
  }
}



// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtblr_cpg(   std::vector<std::vector<double>> wy,
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
                                                        int method) {

 // Define local variables
 int nt = wy.size();
 int m = wy[0].size();
 int nmodels = models.size();
 double nsamples=0.0;

 std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
 std::vector<std::vector<double>> beta(nt, std::vector<double>(m, 0));

 std::vector<double> cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> ves(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));

 std::vector<double> x2t(m);
 std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
 std::vector<int> order(m);

 // Initialize variables
 for (int i = 0; i < m; i++) {
  for (int t = 0; t < nt; t++) {
   r[t][i] = wy[t][i];
   x2[t][i] = (wy[t][i]/ww[t][i])*(wy[t][i]/ww[t][i]);
  }
 }

 // Wy - W'Wb
 for ( int i = 0; i < m; i++) {
  for (int t = 0; t < nt; t++) {
   if (b[t][i]!= 0.0) {
    for (size_t j = 0; j < XXindices[i].size(); j++) {
     r[t][XXindices[i][j]]=r[t][XXindices[i][j]] - XXvalues[t][i][j]*b[t][i];
    }
   }
  }
 }

 for ( int i = 0; i < m; i++) {
  x2t[i] = 0.0;
  for ( int t = 0; t < nt; t++) {
   if(x2[t][i]>x2t[i]) {
    x2t[i] = x2[t][i];
   }
  }
 }


 // Establish order of markers as they are entered into the model
 std::iota(order.begin(), order.end(), 0);
 std::sort(  std::begin(order),
             std::end(order),
             [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );

 std::vector<int> marker_to_set(m, -1);

 for (size_t s = 0; s < sets.size(); ++s) {
  for (int idx : sets[s]) {
   marker_to_set[idx] = s;
  }
 }


 // Initialize (co)variance matrices
 arma::mat Bi = arma::inv(B);
 arma::mat Ei = arma::inv(E);
 arma::mat G(nt,nt, fill::zeros);

 // Start Gibbs sampler
 //std::random_device rd;
 std::mt19937 gen(seed);

 for ( int it = 0; it < nit+nburn; it++) {

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

  //Sample marker effects (BayesC)
  if (method==4) {
   for (size_t s = 0; s < sets.size(); ++s) {

    // Step 1: Sample B using only markers in the current set
    sampleB_cpg(nt, m, nub, B, beta, ssb_prior, gen);

    // Step 2: Invert B to get Bi
    if (!arma::inv_sympd(Bi, B)) throw std::runtime_error("Bi inversion failed");

    // Step 3: Sample marker effects for markers in this set
    // Convert current set to a std::unordered_set for fast lookup
    //std::unordered_set<int> current_set(set.begin(), set.end());

    for (int isort = 0; isort < m; isort++) {
     int i = order[isort];

     // Only update if marker i is in the current set
     if (marker_to_set[i] == s) {
      sampleBetaCPG_Mt_latent_fast(i, nt,
                                   nmodels, models, cmodel, pi,
                                   Ei, Bi,
                                   ww, r, beta, b, d,
                                   XXindices, XXvalues,
                                   gen);

     }
    }

   }
  }
  // Store values
  for (int t = 0; t < nt; t++) {
   for ( int i = 0; i < m; i++) {
    //if (d[t][i] == 1) {
    if (d[t][i] > 0) {
     if ((it > nburn) && (it % nthin == 0)) {
      dm[t][i] = dm[t][i] + 1.0;
      bm[t][i] = bm[t][i] + b[t][i];
     }
    }
   }
  }

  // Sample pi for Bayes C
  if(updatePi && method==4) {
   samplePi_cpg(cmodel, pi, gen);
   for (int k = 0; k<nmodels ; k++) {
    if(it>nburn) pis[k] = pis[k] + pi[k];
   }
  }

  // Sample marker variance
  if(updateB && method==4) {
   sampleB_cpg(nt, m, nub, B, beta, ssb_prior, gen);
   for (int t = 0; t < nt; t++) {
    vbs[t][it] = B(t,t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if(it>nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
    }
   }
  }

  //Update genetic variance
  computeG_cpg(nt, m, b, wy, r, n, G);
  for (int t = 0; t < nt; t++) {
   vgs[t][it] = G(t,t);
  }
  for (int t1 = 0; t1 < nt; t1++  ) {
   for (int t2 = 0; t2 < nt; t2++) {
    if(it>nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1,t2);
   }
  }

  // Sample residual variance
  if(updateE) {
   sampleE_cpg(nt, m, nue, E, b, wy, r, sse_prior, yy, n, gen);
   if (!arma::inv_sympd(Ei, E)) throw std::runtime_error("Ei inversion failed");
   for (int t = 0; t < nt; t++) {
    ves[t][it] = E(t,t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if(it>nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
    }
   }
  }

  if ( (it > nburn) && (it % nthin == 0) ) {
   nsamples = nsamples + 1.0;
  }

 }
 // Summarize results
 std::vector<std::vector<std::vector<double>>> result;
 result.resize(20);

 result[0].resize(nt);
 result[1].resize(nt);
 result[2].resize(nt);
 result[3].resize(nt);
 result[4].resize(nt);
 result[5].resize(nt);
 result[6].resize(nt);
 result[7].resize(nt);
 result[8].resize(nt);
 result[9].resize(nt);
 result[10].resize(nt);
 result[11].resize(nt);
 result[12].resize(nt);
 result[13].resize(nt);
 result[14].resize(nt);
 result[15].resize(nt);
 result[16].resize(nt);
 result[17].resize(nt);
 result[18].resize(nt);
 result[19].resize(nt);

 for (int t=0; t < nt; t++) {
  result[0][t].resize(m);
  result[1][t].resize(m);
  result[2][t].resize(m);
  result[3][t].resize(m);
  result[4][t].resize(m);
  result[5][t].resize(m);
  result[6][t].resize(m);
  result[7][t].resize(nit+nburn);
  result[8][t].resize(nit+nburn);
  result[9][t].resize(nit+nburn);
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

 for (int t=0; t < nt; t++) {
  for (int i=0; i < m; i++) {
   result[0][t][i] = bm[t][i]/nsamples;
   result[1][t][i] = dm[t][i]/nsamples;
   result[2][t][i] = wy[t][i];
   result[3][t][i] = r[t][i];
   result[4][t][i] = b[t][i];
   result[5][t][i] = d[t][i];
   result[6][t][i] = order[i];
  }
 }

 for (int t=0; t < nt; t++) {
  for (int i=0; i < nit+nburn; i++) {
   result[7][t][i] = vbs[t][i];
   result[8][t][i] = vgs[t][i];
   result[9][t][i] = ves[t][i];
  }
 }
 for (int t1=0; t1 < nt; t1++) {
  for (int t2=0; t2 < nt; t2++) {
   result[10][t1][t2] = cvbm[t1][t2]/nsamples;
   result[11][t1][t2] = cvgm[t1][t2]/nsamples;
   result[12][t1][t2] = cvem[t1][t2]/nsamples;
   result[13][t1][t2] = B(t1,t2);
   result[14][t1][t2] = G(t1,t2);
   result[15][t1][t2] = E(t1,t2);
  }
 }
 for (int t=0; t < nt; t++) {
  for (int i=0; i < nmodels; i++) {
   result[16][t][i] = pi[i];
   result[17][t][i] = pis[i]/nsamples;
  }
 }
 for (int t=0; t < nt; t++) {
  for (int i=0; i < 4; i++) {
   result[18][t][i] = 0.0;
  }
  for (int i=0; i < 2; i++) {
   result[19][t][i] = 0.0;
  }
 }
 return result;
}
