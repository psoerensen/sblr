// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "cpg_samplers.h"
#include "distributions.h"

#include <cmath>
#include <stdexcept>
#include <numeric>
#include <algorithm>


void sampleB_cpg_arma(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const arma::mat& beta,
  const arma::mat& ssb_prior,
  std::mt19937& gen)
{
 if (nub <= nt-1)
  throw std::runtime_error("nub must be > nt-1");

 arma::mat S_post = ssb_prior;

 // BLAS-3 (FAST)
 S_post += beta * beta.t();

 unsigned int df_post = nub + m;

 B = rinvwishart_cpp(df_post, S_post, gen);

 B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);
}


void sampleE_cpg_arma(
  int nt,
  int m,
  int nue,
  arma::mat& E,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const arma::mat& sse_prior,
  const arma::vec& yy,
  const std::vector<int>& n,
  std::mt19937& gen
)
{
 arma::vec Se_diag(nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {

  // SSE = b * (r + wy)
  double sse = arma::dot(b.row(t), r.row(t) + wy.row(t));

  sse = yy(t) - sse;

  Se_diag(t) = sse + nue * sse_prior(t,t);
 }

 // Sample diagonal of E
 for (int t = 0; t < nt; ++t) {
  std::chi_squared_distribution<double> rchisq(n[t] + nue);
  double chi2 = rchisq(gen);
  E(t,t) = Se_diag(t) / chi2;
 }
}


void computeG_cpg_arma(
  int nt,
  int m,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const std::vector<int>& n,
  arma::mat& G
)
{
 // Precompute Q = wy - r
 arma::mat Q = wy - r;

 G.set_size(nt, nt);

 for (int t1 = 0; t1 < nt; ++t1) {

  for (int t2 = t1; t2 < nt; ++t2) {

   double denom = std::sqrt((double)n[t1] * (double)n[t2]);

   // BLAS dot product
   double s12 = arma::dot(b.row(t1), Q.row(t2));

   if (t1 != t2) {
    double s21 = arma::dot(b.row(t2), Q.row(t1));
    s12 = 0.5 * (s12 + s21);
   }

   double gij = s12 / denom;

   G(t1,t2) = gij;
   G(t2,t1) = gij;
  }
 }
}


// Sample pi
void samplePi_cpg(std::vector<double>& cmodel,
                  std::vector<double>& pi,
                  std::mt19937& gen) {
 // Iterate over elements of cmodel
 for (size_t k = 0; k < cmodel.size(); k++) {
  // Create a gamma distribution with shape parameter cmodel[k] and scale parameter 1.0
  std::gamma_distribution<double> rgamma(cmodel[k], 1.0);
  // Generate a random gamma-distributed value using the provided random number generator gen
  double rg = rgamma(gen);
  // Store the generated value in the pi vector
  pi[k] = rg;
 }
 // Calculate the sum of all values in the pi vector
 double psum = std::accumulate(pi.begin(), pi.end(), 0.0);
 // Normalize the pi vector by dividing each element by the sum
 for (size_t k = 0; k < cmodel.size(); k++) {
  pi[k] = pi[k] / psum;
 }
 // Reset all elements in the cmodel vector to 1.0
 std::fill(cmodel.begin(), cmodel.end(), 1.0);
}




// Proper Gibbs update for B under:
//   beta_i ~ N(0, B)
//   B ~ InvWishart(nub, ssb_prior)
//
// Posterior:
//   B | beta ~ InvWishart(nub + m, ssb_prior + sum_i beta_i beta_i')

void sampleB_cpg(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const std::vector<std::vector<double>>& beta,
  const std::vector<std::vector<double>>& ssb_prior,
  std::mt19937& gen)
{
 if (nub <= nt - 1) {
  throw std::runtime_error("sampleB_latent: nub must be > nt - 1.");
 }

 // -----------------------------------
 // Build prior scale matrix S0
 // -----------------------------------
 arma::mat S_post(nt, nt);

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = 0; t2 < nt; ++t2) {
   S_post(t1, t2) = ssb_prior[t1][t2];
  }
 }

 // -----------------------------------
 // Accumulate beta beta'
 // -----------------------------------
 // IMPORTANT: avoid allocating arma::vec inside loop
 arma::vec beta_i(nt);

 for (int i = 0; i < m; ++i) {

  // fill vector (cheap, contiguous write)
  for (int t = 0; t < nt; ++t) {
   beta_i[t] = beta[t][i];
  }

  // rank-1 update
  S_post += beta_i * beta_i.t();
 }

 // -----------------------------------
 // Sample from Inv-Wishart
 // -----------------------------------
 const unsigned int df_post = static_cast<unsigned int>(nub + m);

 B = rinvwishart_cpp(df_post, S_post, gen);

 // -----------------------------------
 // Minimal cleanup (cheap + safe)
 // -----------------------------------
 B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);
}



// Sample residual covariance matrix E
void sampleE_cpg(
  int nt,
  int m,
  int nue,
  arma::mat& E,
  const std::vector<std::vector<double>>& b,
  const std::vector<std::vector<double>>& wy,
  const std::vector<std::vector<double>>& r,
  const std::vector<std::vector<double>>& sse_prior,
  const std::vector<double>& yy,
  const std::vector<int>& n,
  std::mt19937& gen
) {
 // Initialize the Se matrix to store the sum of squared residuals
 arma::mat Se(nt, nt, arma::fill::zeros);

 // Calculate the sum of squared residuals for each time point
 for (int t1 = 0; t1 < nt; t1++) {
  double sse = 0.0; // Initialize the sum of squared residuals for t1

  // Calculate the sum of squared residuals using b, wy, and r vectors
  for (int i = 0; i < m; i++) {
   sse += b[t1][i] * (r[t1][i] + wy[t1][i]);
  }

  // Calculate the adjusted sum of squared residuals
  sse = yy[t1] - sse;

  // Store the adjusted sum of squared residuals in the Se matrix
  Se(t1, t1) = sse + nue * sse_prior[t1][t1];
 }

 // Generate chi-squared random variables and update the E matrix
 for (int t = 0; t < nt; t++) {
  std::chi_squared_distribution<double> rchisq(n[t] + nue);
  double chi2 = rchisq(gen);
  E(t, t) = Se(t, t) / chi2;
 }

}


// Compute genetic covariance matrix G
// Interprets: wy[t] - r[t]  ≈  (X'X) b[t]
void computeG_cpg(
  int nt,
  int m,
  const std::vector<std::vector<double>>& b,
  const std::vector<std::vector<double>>& wy,
  const std::vector<std::vector<double>>& r,
  const std::vector<int>& n,
  arma::mat& G
) {
 // Basic dimension checks (defensive)
 if ((int)b.size() != nt || (int)wy.size() != nt || (int)r.size() != nt || (int)n.size() != nt) {
  throw std::runtime_error("computeG: nt does not match sizes of b/wy/r/n.");
 }
 for (int t = 0; t < nt; ++t) {
  if ((int)b[t].size() != m || (int)wy[t].size() != m || (int)r[t].size() != m) {
   throw std::runtime_error("computeG: m does not match marker dimension for some trait.");
  }
  if (n[t] <= 0) {
   throw std::runtime_error("computeG: n[t] must be positive.");
  }
 }

 // Ensure correct size and reset
 G.set_size(nt, nt);
 G.zeros();

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1; t2 < nt; ++t2) {

   // denom = sqrt(n1 * n2) matches your existing normalization
   const double denom = std::sqrt((double)n[t1] * (double)n[t2]);

   long double s12 = 0.0L;
   for (int i = 0; i < m; ++i) {
    const double q2 = wy[t2][i] - r[t2][i]; // ≈ (X'X b_t2)_i
    s12 += (long double)b[t1][i] * (long double)q2;
   }

   // Symmetrize for off-diagonals to reduce floating point drift
   if (t1 != t2) {
    long double s21 = 0.0L;
    for (int i = 0; i < m; ++i) {
     const double q1 = wy[t1][i] - r[t1][i]; // ≈ (X'X b_t1)_i
     s21 += (long double)b[t2][i] * (long double)q1;
    }
    s12 = 0.5L * (s12 + s21);
   }

   const double gij = (double)(s12 / (long double)denom);
   G(t1, t2) = gij;
   G(t2, t1) = gij;
  }
 }
}
