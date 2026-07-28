// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <memory>
#include <sstream>

#include "blr_mt_ld_access.h"
#include "blr_mt_covariance_rng.h"
#include "st_csr_common.h"
#include "st_block_eigen.h"
#include "st_block_eigen_rcpp.h"

using namespace Rcpp;
using namespace arma;

// Matrix multiplication helper
arma::mat mmult(const arma::mat& A, const arma::mat& B) {
 return A * B;
}


arma::mat mvrnormARMA(const arma::mat& sigma) {
 int ncols = sigma.n_cols;
 arma::mat Y = arma::randn(1, ncols);
 return Y * arma::chol(sigma);
}

arma::mat rwishart(unsigned int df, const arma::mat& S) {
 unsigned int m = S.n_rows;
 arma::mat Z(m, m, fill::zeros);

 // Fill the diagonal with sqrt of chi-squared draws
 for (unsigned int i = 0; i < m; i++) {
  Z(i, i) = std::sqrt(R::rchisq(df - i));
 }

 // Fill the lower triangle with standard normal values
 for (unsigned int j = 0; j < m; j++) {
  for (unsigned int i = j + 1; i < m; i++) {
   Z(i, j) = R::rnorm(0.0, 1.0);
  }
 }

 arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
 return C.t() * C;
}

arma::mat riwishart(unsigned int df, const arma::mat& S) {
 return arma::inv(rwishart(df, arma::inv(S)));
}

arma::mat shrink_B_from_Sb(const arma::mat& Sb,
                           const arma::mat& B_prior,
                           double diag_weight = 1.0,
                           double corr_shrink = 0.5,
                           double eps = 1e-8) {
 int nt = Sb.n_rows;

 arma::vec s(nt);
 for (int t = 0; t < nt; ++t) {
  double v = Sb(t, t) + diag_weight * B_prior(t, t);
  s(t) = std::sqrt(std::max(v, eps));
 }

 arma::mat R(nt, nt, arma::fill::eye);

 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1 + 1; t2 < nt; ++t2) {
   double denom = s(t1) * s(t2);
   double r = 0.0;
   if (denom > eps) {
    r = Sb(t1, t2) / denom;
   }

   // shrink toward 0
   r *= corr_shrink;

   // clamp for stability
   r = std::max(-0.95, std::min(0.95, r));

   R(t1, t2) = r;
   R(t2, t1) = r;
  }
 }

 arma::mat D = arma::diagmat(s);
 arma::mat B = D * R * D;

 // small ridge for SPD stability
 B.diag() += 1e-8;

 return 0.5 * (B + B.t());
}

// Draw W ~ Wishart(df, V)
// where V is the scale matrix
arma::mat rwishart(unsigned int df,
                   const arma::mat& V,
                   std::mt19937& gen) {
 return sblr::mt::draw_wishart(static_cast<double>(df), V, gen);
}

// Draw Sigma ~ InvWishart(df, S)
// using: if W ~ Wishart(df, S^{-1}), then Sigma = W^{-1}
arma::mat rinvwishart(unsigned int df,
                      const arma::mat& S,
                      std::mt19937& gen) {
 return sblr::mt::draw_inverse_wishart(
  static_cast<double>(df), S, gen);
}

// Proper Gibbs update for B under:
//   beta_i ~ N(0, B)
//   B ~ InvWishart(nub, ssb_prior)
//
// Posterior:
//   B | beta ~ InvWishart(nub + m, ssb_prior + sum_i beta_i beta_i')

void sampleB_latent(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const std::vector<std::vector<double>>& beta,   // latent effects [nt][m]
  const std::vector<std::vector<double>>& ssb_prior,
  std::mt19937& gen) {

 if (nub <= nt - 1) {
  throw std::runtime_error("sampleB_latent: nub must be > nt - 1.");
 }

 // Prior scale matrix S0
 arma::mat S0(nt, nt, arma::fill::zeros);
 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = 0; t2 < nt; ++t2) {
   S0(t1, t2) = ssb_prior[t1][t2];
  }
 }

 // Posterior scale:
 // S_post = S0 + sum_i beta_i beta_i'
 arma::mat S_post = S0;

 for (int i = 0; i < m; ++i) {
  arma::vec beta_i(nt);
  for (int t = 0; t < nt; ++t) {
   beta_i(t) = beta[t][i];
  }
  S_post += beta_i * beta_i.t();
 }

 S_post = 0.5 * (S_post + S_post.t());

 const unsigned int df_post = static_cast<unsigned int>(nub + m);

 B = rinvwishart(df_post, S_post, gen);

 // small numeric cleanup
 B = 0.5 * (B + B.t());
 B.diag() = arma::clamp(B.diag(), 1e-12, arma::datum::inf);
}

struct ModelCache {
 arma::mat U;     // Cholesky factor (upper) of C_k
 double logdet;   // log |C_k|
};

void precomputeModelCache(int nt,
                          int nmodels,
                          const std::vector<std::vector<int>>& models,
                          const arma::mat& Ei,
                          const arma::mat& Bi,
                          double ww_const,
                          const std::vector<int>& allowed_traits,
                          std::vector<ModelCache>& cache)
{
 cache.resize(nmodels);

 for(int k = 0; k < nmodels; k++)
 {
  bool valid = true;

  for(int t = 0; t < nt; t++)
  {
   if(models[k][t] == 1 && !allowed_traits[t])
   {
    valid = false;
    break;
   }
  }

  if(!valid)
  {
   cache[k].logdet = std::numeric_limits<double>::infinity();
   cache[k].U.reset();
   continue;
  }

  arma::mat C = Bi;

  for(int t = 0; t < nt; t++)
  {
   if(models[k][t] == 1)
    C(t,t) += ww_const * Ei(t,t);
  }

  arma::mat U;

  bool ok = arma::chol(U, C);

  if(!ok)
  {
   C.diag() += 1e-8 * arma::mean(C.diag());
   arma::chol(U, C);
  }

  double logdet = 2.0 * arma::sum(log(U.diag()));

  cache[k].U = U;
  cache[k].logdet = logdet;
 }
}

double quadratic_form(const arma::mat& U,
                      const std::vector<double>& rhs)
{
 arma::vec r(rhs);

 arma::vec y = arma::solve(arma::trimatu(U), r);

 return arma::dot(y, y);
}

arma::vec posterior_mean(const arma::mat& U,
                         const std::vector<double>& rhs)
{
 arma::vec r(rhs);

 arma::vec y = arma::solve(arma::trimatl(U.t()), r);
 arma::vec x = arma::solve(arma::trimatu(U), y);

 return x;
}

arma::vec sample_posterior(const arma::mat& U,
                           std::mt19937& gen)
{
 int nt = U.n_rows;

 arma::vec z(nt);

 std::normal_distribution<double> rnorm(0.0,1.0);

 for(int t = 0; t < nt; t++)
  z(t) = rnorm(gen);

 arma::vec sample = arma::solve(arma::trimatu(U), z);

 return sample;
}

arma::mat getBi_local(int i,
                      const arma::mat& B,
                      const std::vector<std::vector<double>>& b,
                      const std::vector<std::vector<int>>& d,
                      std::string method = "trace",
                      double lambda = -1.0) {
  int nt = b.size();
  arma::mat Sb_local(nt, nt, fill::zeros);
  bool any_active = false;

  // Compute local covariance matrix (outer product for included traits)
  for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = t1; t2 < nt; t2++) {
      if (d[t1][i] == 1 && d[t2][i] == 1) {
        double val = b[t1][i] * b[t2][i];
        Sb_local(t1, t2) = val;
        Sb_local(t2, t1) = val;
        any_active = true;
      }
    }
  }

  // Ensure positive definiteness
  for (int t = 0; t < nt; t++) {
    Sb_local(t, t) += 1e-6;
  }

  // Compute alpha based on selected shrinkage method
  double alpha = 0.5;
  if (!any_active) {
    return inv(B); // fall back to global
  }

  if (method == "trace") {
    double trSb = trace(Sb_local);
    double trB = trace(B);
    if (lambda < 0) lambda = trB / 10.0;
    alpha = trSb / (trSb + lambda);
  } else if (method == "norm") {
    double normSb = norm(Sb_local, "fro");
    double normB = norm(B, "fro");
    if (lambda < 0) lambda = normB / 10.0;
    alpha = normSb / (normSb + lambda);
  }

  alpha = std::min(std::max(alpha, 0.05), 0.95); // clip to [0.05, 0.95]

  // Shrink toward global covariance B
  arma::mat B_local = alpha * B + (1.0 - alpha) * Sb_local;
  return inv(B_local);
}

// Sample pi
void samplePi(std::vector<double>& cmodel,
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

void sampleBset_old(int nt,
             int m,
             int nub,
             arma::mat& B,
             const std::vector<std::vector<int>>& d,
             const std::vector<std::vector<double>>& b,
             const std::vector<std::vector<double>>& ssb_prior,
             const std::vector<int>& subset, // new argument
             std::mt19937& gen) {

 arma::mat Sb(nt, nt, arma::fill::zeros);
 arma::mat dfB(nt, nt, arma::fill::zeros);

 for (int t1 = 0; t1 < nt; t1++) {
  for (int t2 = t1; t2 < nt; t2++) {
   double ssb = 0.0;
   double dfb = 0.0;
   for (int idx = 0; idx < subset.size(); idx++) {
    int i = subset[idx];
    if (d[t1][i] == 1 && d[t2][i] == 1) {
     ssb += b[t1][i] * b[t2][i];
     dfb += 1.0;
    }
   }
   Sb(t1, t2) = Sb(t2, t1) = ssb;
   dfB(t1, t2) = dfB(t2, t1) = dfb;
  }
 }

 // Convert ssb_prior to arma::mat
 arma::mat ssb_prior_mat(nt, nt);
 for (int i = 0; i < nt; ++i) {
  for (int j = 0; j < nt; ++j) {
   ssb_prior_mat(i, j) = ssb_prior[i][j];
  }
 }

 // Posterior scale matrix
 arma::mat S = Sb + nub * ssb_prior_mat;

 // Degrees of freedom
 int df = static_cast<int>(arma::accu(dfB.diag()) / nt + nub);

 // Sample B from inverse Wishart
 B = riwishart(df, S);

 // Ensure symmetry
 B = 0.5 * (B + B.t());
}

void sampleB_old(int nt,
             int m,
             int nub,
             arma::mat& B,
             const std::vector<std::vector<int>>& d,
             const std::vector<std::vector<double>>& b,
             const std::vector<std::vector<double>>& ssb_prior,
             std::mt19937& gen) {

 arma::mat Sb(nt, nt, arma::fill::zeros);
 arma::mat dfB(nt, nt, arma::fill::zeros);

 for (int t1 = 0; t1 < nt; t1++) {
  for (int t2 = t1; t2 < nt; t2++) {
   double ssb = 0.0;
   double dfb = 0.0;
   for (int i = 0; i < m; i++) {
    if (d[t1][i] == 1 && d[t2][i] == 1) {
     ssb += b[t1][i] * b[t2][i];
     dfb += 1.0;
    }
   }
   Sb(t1, t2) = Sb(t2, t1) = ssb;
   dfB(t1, t2) = dfB(t2, t1) = dfb;
  }
 }

 // Convert ssb_prior to arma::mat
 arma::mat ssb_prior_mat(nt, nt);
 for (int i = 0; i < nt; ++i) {
  for (int j = 0; j < nt; ++j) {
   ssb_prior_mat(i, j) = ssb_prior[i][j];
  }
 }

 // Posterior scale matrix
 arma::mat S = Sb + nub * ssb_prior_mat;

 // Degrees of freedom
 int df = static_cast<int>(arma::accu(dfB.diag()) / nt + nub);

 // Sample B from inverse Wishart
 B = riwishart(df, S);

 // ensure sparsity
 double rho_shrink = 0.5;
 for(int t1=0;t1<nt;t1++){
  for(int t2=0;t2<nt;t2++){
   if(t1!=t2){
    B(t1,t2) *= rho_shrink;
   }
  }
 }
 //

 // Ensure symmetry
 B = 0.5 * (B + B.t());
}

void sampleBset(int nt,
                int m,
                int nub,
                arma::mat& B,
                const std::vector<std::vector<int>>& d,
                const std::vector<std::vector<double>>& b,
                const std::vector<std::vector<double>>& ssb_prior,
                const std::vector<int>& subset,
                std::mt19937& gen) {

 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 for (int t1 = 0; t1 < nt; ++t1)
  for (int t2 = 0; t2 < nt; ++t2)
   ssb_prior_mat(t1, t2) = ssb_prior[t1][t2];

 arma::vec var_draw(nt, arma::fill::zeros);
 arma::mat R(nt, nt, arma::fill::eye);

 // --------------------------------
 // 1. Trait-specific diagonal variances
 // --------------------------------
 for (int t = 0; t < nt; ++t) {

  double ss = 0.0;
  int n_active = 0;

  for (size_t idx = 0; idx < subset.size(); ++idx) {
   int i = subset[idx];

   if (d[t][i] == 1) {
    ss += b[t][i] * b[t][i];
    n_active++;
   }
  }

  double scale = ss + nub * ssb_prior_mat(t, t);
  double df = static_cast<double>(nub + n_active);

  std::chi_squared_distribution<double> rchisq(std::max(df, 1.0));
  double x = std::max(rchisq(gen), 1e-12);

  var_draw(t) = std::max(scale / x, 1e-8);
 }

 // --------------------------------
 // 2. Pairwise correlations
 // --------------------------------
 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1 + 1; t2 < nt; ++t2) {

   double s12 = 0.0;
   int n_shared = 0;

   for (size_t idx = 0; idx < subset.size(); ++idx) {

    int i = subset[idx];

    if (d[t1][i] == 1 && d[t2][i] == 1) {
     s12 += b[t1][i] * b[t2][i];
     n_shared++;
    }
   }

   double denom = std::sqrt(var_draw(t1) * var_draw(t2));
   double r_hat = 0.0;

   if (n_shared > 0 && denom > 1e-12) {

    double cov_hat = s12 / static_cast<double>(n_shared);

    double prior_cov = ssb_prior_mat(t1, t2);

    double w = static_cast<double>(n_shared) /
     static_cast<double>(n_shared + nub);

    cov_hat = w * cov_hat + (1.0 - w) * prior_cov;

    r_hat = cov_hat / denom;
   }

   double corr_shrink =
    static_cast<double>(n_shared) /
     static_cast<double>(n_shared + 20.0);

   double r = corr_shrink * r_hat;

   r = std::max(-0.95, std::min(0.95, r));

   R(t1, t2) = r;
   R(t2, t1) = r;
  }
 }

 // --------------------------------
 // 3. Reconstruct covariance matrix
 // --------------------------------
 arma::vec s = arma::sqrt(var_draw);
 arma::mat D = arma::diagmat(s);

 B = D * R * D;

 B.diag() += 1e-8;

 B = 0.5 * (B + B.t());

 // SPD safeguard
 arma::vec eigval;
 arma::mat eigvec;

 if (!arma::eig_sym(eigval, eigvec, B)) {
  B = arma::diagmat(arma::clamp(B.diag(), 1e-8, arma::datum::inf));
  return;
 }

 eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);

 B = eigvec * arma::diagmat(eigval) * eigvec.t();

 B = 0.5 * (B + B.t());
}

void sampleB(int nt,
             int m,
             int nub,
             arma::mat& B,
             const std::vector<std::vector<int>>& d,
             const std::vector<std::vector<double>>& b,
             const std::vector<std::vector<double>>& ssb_prior,
             std::mt19937& gen) {

 arma::mat ssb_prior_mat(nt, nt, arma::fill::zeros);
 for (int t1 = 0; t1 < nt; ++t1)
  for (int t2 = 0; t2 < nt; ++t2)
   ssb_prior_mat(t1, t2) = ssb_prior[t1][t2];

 arma::vec var_draw(nt, arma::fill::zeros);
 arma::mat R(nt, nt, arma::fill::eye);

 // ----------------------------
 // 1. Update diagonal variances
 // ----------------------------
 for (int t = 0; t < nt; ++t) {
  double ss = 0.0;
  int n_active = 0;

  for (int i = 0; i < m; ++i) {
   if (d[t][i] == 1) {
    ss += b[t][i] * b[t][i];
    n_active++;
   }
  }

  // posterior scale for variance component
  double scale = ss + nub * ssb_prior_mat(t, t);
  double df = static_cast<double>(nub + n_active);

  // inverse-chi-square style draw:
  // if X ~ ChiSq(df), then scale / X is scaled-inv-chi-square
  std::chi_squared_distribution<double> rchisq(std::max(df, 1.0));
  double x = std::max(rchisq(gen), 1e-12);

  var_draw(t) = std::max(scale / x, 1e-8);
 }

 // --------------------------------------------
 // 2. Update off-diagonal correlations pairwise
 // --------------------------------------------
 for (int t1 = 0; t1 < nt; ++t1) {
  for (int t2 = t1 + 1; t2 < nt; ++t2) {

   double s12 = 0.0;
   int n_shared = 0;

   for (int i = 0; i < m; ++i) {
    if (d[t1][i] == 1 && d[t2][i] == 1) {
     s12 += b[t1][i] * b[t2][i];
     n_shared++;
    }
   }

   double denom = std::sqrt(var_draw(t1) * var_draw(t2));
   double r_hat = 0.0;

   if (n_shared > 0 && denom > 1e-12) {
    // crude empirical covariance among shared markers
    double cov_hat = s12 / static_cast<double>(n_shared);

    // optional prior contribution
    double prior_cov = ssb_prior_mat(t1, t2);
    double w = static_cast<double>(n_shared) /
     static_cast<double>(n_shared + nub);

    cov_hat = w * cov_hat + (1.0 - w) * prior_cov;

    r_hat = cov_hat / denom;
   }

   // shrink correlation according to shared support
   double corr_shrink = static_cast<double>(n_shared) /
    static_cast<double>(n_shared + 20.0);

   double r = corr_shrink * r_hat;

   // clamp
   r = std::max(-0.95, std::min(0.95, r));

   R(t1, t2) = r;
   R(t2, t1) = r;
  }
 }

 // ----------------------------
 // 3. Reconstruct covariance B
 // ----------------------------
 arma::vec s = arma::sqrt(var_draw);
 arma::mat D = arma::diagmat(s);
 B = D * R * D;

 // small ridge for stability
 B.diag() += 1e-8;

 // enforce symmetry
 B = 0.5 * (B + B.t());

 // optional SPD fix via eigenvalue floor
 arma::vec eigval;
 arma::mat eigvec;
 if (!arma::eig_sym(eigval, eigvec, B)) {
  B = arma::diagmat(arma::clamp(B.diag(), 1e-8, arma::datum::inf));
  return;
 }

 eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
 B = eigvec * arma::diagmat(eigval) * eigvec.t();
 B = 0.5 * (B + B.t());
}

// void sampleB(int nt,
//              int m,
//              int nub,
//              arma::mat& B,
//              const std::vector<std::vector<int>>& d,
//              const std::vector<std::vector<double>>& b,
//              const std::vector<std::vector<double>>& ssb_prior,
//              std::mt19937& gen) {
//
//  arma::mat Sb(nt, nt, arma::fill::zeros);
//  arma::mat dfB(nt, nt, arma::fill::zeros);
//
//  double n_shared = 0.0;
//
//  // Build covariance from current effects
//  for (int t1 = 0; t1 < nt; t1++) {
//   for (int t2 = t1; t2 < nt; t2++) {
//
//    double ssb = 0.0;
//    double dfb = 0.0;
//
//    for (int i = 0; i < m; i++) {
//
//     int model_size = 0;
//     for (int tt = 0; tt < nt; ++tt)
//      model_size += d[tt][i];
//
//     if (d[t1][i] == 1 && d[t2][i] == 1) {
//
//      if (t1 == t2 || model_size >= 2) {
//       ssb += b[t1][i] * b[t2][i];
//       dfb += 1.0;
//
//       if (t1 != t2 && model_size >= 2)
//        n_shared += 1.0;
//      }
//     }
//    }
//
//    Sb(t1, t2) = Sb(t2, t1) = ssb;
//    dfB(t1, t2) = dfB(t2, t1) = dfb;
//   }
//  }
//
//  // Convert prior matrix
//  arma::mat ssb_prior_mat(nt, nt);
//  for (int i = 0; i < nt; ++i)
//   for (int j = 0; j < nt; ++j)
//    ssb_prior_mat(i, j) = ssb_prior[i][j];
//
//  arma::mat S = Sb + nub * ssb_prior_mat;
//
//  int df = static_cast<int>(arma::accu(dfB.diag()) / nt + nub);
//
//  // Inverse Wishart draw
//  arma::mat B_draw = riwishart(df, S);
//  B_draw = 0.5 * (B_draw + B_draw.t());
//
//  // Extract variances
//  arma::vec s(nt);
//  for (int t = 0; t < nt; ++t)
//   s(t) = std::sqrt(std::max(B_draw(t, t), 1e-8));
//
//  // Adaptive correlation shrinkage
//  double corr_shrink = n_shared / (n_shared + 20.0);
//
//  arma::mat R(nt, nt, arma::fill::eye);
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1 + 1; t2 < nt; ++t2) {
//
//    double denom = s(t1) * s(t2);
//    double r = 0.0;
//
//    if (denom > 1e-12)
//     r = B_draw(t1, t2) / denom;
//
//    // Shrink correlation
//    r *= corr_shrink;
//
//    // Clamp
//    r = std::max(-0.95, std::min(0.95, r));
//
//    R(t1, t2) = r;
//    R(t2, t1) = r;
//   }
//  }
//
//  arma::mat D = arma::diagmat(s);
//  B = D * R * D;
//
//  // Small ridge
//  B.diag() += 1e-8;
//
//  B = 0.5 * (B + B.t());
// }

// void sampleB(int nt,
//              int m,
//              int nub,
//              arma::mat& B,
//              const std::vector<std::vector<int>>& d,
//              const std::vector<std::vector<double>>& b,
//              const std::vector<std::vector<double>>& ssb_prior,
//              std::mt19937& gen) {
//
//  arma::mat Sb(nt, nt, arma::fill::zeros);
//  arma::mat dfB(nt, nt, arma::fill::zeros);
//
//  for (int t1 = 0; t1 < nt; t1++) {
//   for (int t2 = t1; t2 < nt; t2++) {
//    double ssb = 0.0;
//    double dfb = 0.0;
//
//    for (int i = 0; i < m; i++) {
//     int model_size = 0;
//     for (int tt = 0; tt < nt; ++tt) model_size += d[tt][i];
//
//     if (d[t1][i] == 1 && d[t2][i] == 1) {
//      if (t1 == t2 || model_size >= 2) {
//       ssb += b[t1][i] * b[t2][i];
//       dfb += 1.0;
//      }
//     }
//    }
//
//    Sb(t1, t2) = Sb(t2, t1) = ssb;
//    dfB(t1, t2) = dfB(t2, t1) = dfb;
//   }
//  }
//
//  arma::mat ssb_prior_mat(nt, nt);
//  for (int i = 0; i < nt; ++i) {
//   for (int j = 0; j < nt; ++j) {
//    ssb_prior_mat(i, j) = ssb_prior[i][j];
//   }
//  }
//
//  arma::mat S = Sb + nub * ssb_prior_mat;
//  int df = static_cast<int>(arma::accu(dfB.diag()) / nt + nub);
//
//  arma::mat B_draw = riwishart(df, S);
//  B_draw = 0.5 * (B_draw + B_draw.t());
//
//  // convert to D R D and shrink off-diagonal correlations
//  arma::vec s(nt);
//  for (int t = 0; t < nt; ++t) {
//   s(t) = std::sqrt(std::max(B_draw(t, t), 1e-8));
//  }
//
//  arma::mat R(nt, nt, arma::fill::eye);
//  double corr_shrink = 0.4;
//
//  for (int t1 = 0; t1 < nt; ++t1) {
//   for (int t2 = t1 + 1; t2 < nt; ++t2) {
//    double r = B_draw(t1, t2) / (s(t1) * s(t2));
//    r *= corr_shrink;
//    r = std::max(-0.95, std::min(0.95, r));
//    R(t1, t2) = r;
//    R(t2, t1) = r;
//   }
//  }
//
//  arma::mat D = arma::diagmat(s);
//  B = D * R * D;
//  B.diag() += 1e-8;
//  B = 0.5 * (B + B.t());
// }

// Sample residual covariance matrix E
void sampleE(
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

 // Calculate the inverse of the E matrix
 //arma::mat Ei = arma::inv(E);
}

// Sample residual covariance matrix E (full, with correlations)
// Hybrid approach:
//   - exact diagonal (from SSE)
//   - approximate off-diagonal (from residual projections)
//   - shrinkage + SPD enforcement

void sampleE_full(
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

 // ------------------------------------------------------------
 // 1) Diagonal update (exact, same as before)
 // ------------------------------------------------------------
 arma::vec var(nt, arma::fill::zeros);

 for (int t = 0; t < nt; t++) {

  double sse = 0.0;

  for (int i = 0; i < m; i++) {
   sse += b[t][i] * (r[t][i] + wy[t][i]);
  }

  sse = yy[t] - sse;

  double scale = sse + nue * sse_prior[t][t];
  double df = static_cast<double>(n[t] + nue);

  std::chi_squared_distribution<double> rchisq(df);
  double chi2 = std::max(rchisq(gen), 1e-12);

  var(t) = std::max(scale / chi2, 1e-8);
 }

 // ------------------------------------------------------------
 // 2) Build correlation matrix from residual projections
 // ------------------------------------------------------------
 arma::mat R(nt, nt, arma::fill::eye);

 for (int t1 = 0; t1 < nt; t1++) {
  for (int t2 = t1 + 1; t2 < nt; t2++) {

   double s12 = 0.0;
   int count = 0;

   // approximate cross-residual covariance using marker projections
   for (int i = 0; i < m; i++) {
    s12 += r[t1][i] * r[t2][i];
    count++;
   }

   double cov = 0.0;
   if (count > 0) {
    cov = s12 / static_cast<double>(count);
   }

   double denom = std::sqrt(var(t1) * var(t2));
   double corr = 0.0;

   if (denom > 1e-12) {
    corr = cov / denom;
   }

   // shrink correlation (important for stability)
   double shrink = static_cast<double>(count) /
    static_cast<double>(count + 50.0);

   corr *= shrink;

   // clamp to valid range
   corr = std::max(-0.95, std::min(0.95, corr));

   R(t1, t2) = corr;
   R(t2, t1) = corr;
  }
 }

 // ------------------------------------------------------------
 // 3) Reconstruct covariance matrix E
 // ------------------------------------------------------------
 arma::vec s = arma::sqrt(var);
 arma::mat D = arma::diagmat(s);

 E = D * R * D;

 // enforce symmetry
 E = 0.5 * (E + E.t());

 // ------------------------------------------------------------
 // 4) Ensure SPD (very important)
 // ------------------------------------------------------------
 arma::vec eigval;
 arma::mat eigvec;

 if (!arma::eig_sym(eigval, eigvec, E)) {
  // fallback to diagonal if decomposition fails
  E = arma::diagmat(arma::clamp(var, 1e-8, arma::datum::inf));
  return;
 }

 // floor eigenvalues
 eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);

 E = eigvec * arma::diagmat(eigval) * eigvec.t();
 E = 0.5 * (E + E.t());
}

// Sample residual covariance matrix E
// Principled version for current summaries:
//   - exact diagonal updates
//   - fixed residual correlation matrix Re
void sampleE(
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
  const arma::mat& Re,   // fixed residual correlation matrix
  std::mt19937& gen
) {
 arma::vec var(nt, arma::fill::zeros);

 // Exact traitwise variance updates
 for (int t = 0; t < nt; t++) {
  double sse = 0.0;

  for (int i = 0; i < m; i++) {
   sse += b[t][i] * (r[t][i] + wy[t][i]);
  }

  sse = yy[t] - sse;

  double scale = sse + nue * sse_prior[t][t];
  double df = static_cast<double>(n[t] + nue);

  std::chi_squared_distribution<double> rchisq(df);
  double chi2 = std::max(rchisq(gen), 1e-12);

  var(t) = std::max(scale / chi2, 1e-8);
 }

 // Reconstruct E = D Re D
 arma::vec sd = arma::sqrt(var);
 arma::mat D = arma::diagmat(sd);

 E = D * Re * D;
 E = 0.5 * (E + E.t());

 // small cleanup
 arma::vec eigval;
 arma::mat eigvec;
 if (!arma::eig_sym(eigval, eigvec, E)) {
  E = arma::diagmat(var);
  return;
 }

 eigval = arma::clamp(eigval, 1e-8, arma::datum::inf);
 E = eigvec * arma::diagmat(eigval) * eigvec.t();
 E = 0.5 * (E + E.t());
}

// Sample residual covariance matrix E
void sampleE_eigen(
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

 std::vector<int> nq(nt, 0);

 // Calculate the sum of squared residuals for each time point
 for (int t1 = 0; t1 < nt; t1++) {
  double sse_rot = 0.0; // Initialize the sum of squared residuals for t1

  // Calculate the sum of squared residuals using b, wy, and r vectors
  for (int i = 0; i < m; i++) {
   if (r[t1][i] != 0.0) {
    sse_rot += r[t1][i]*r[t1][i];
    nq[t1] += 1;
   }
  }

  // scale back to original likelihood scale: sse = n[t1] * sum(r^2)
  const double sse = static_cast<double>(n[t1]) * sse_rot;

  // Store the adjusted sum of squared residuals in the Se matrix
  Se(t1, t1) = sse + nue * sse_prior[t1][t1];
 }

 // Generate chi-squared random variables and update the E matrix
 for (int t = 0; t < nt; t++) {
  std::chi_squared_distribution<double> rchisq(nq[t] + nue);
  double chi2 = rchisq(gen);
  E(t, t) = Se(t, t) / chi2;
 }

}

// Compute genetic covariance matrix G
// Interprets: wy[t] - r[t]  ≈  (X'X) b[t]
void computeG(
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
    Rcpp::stop("computeG: nt does not match sizes of b/wy/r/n.");
  }
  for (int t = 0; t < nt; ++t) {
    if ((int)b[t].size() != m || (int)wy[t].size() != m || (int)r[t].size() != m) {
      Rcpp::stop("computeG: m does not match marker dimension for some trait.");
    }
    if (n[t] <= 0) {
      Rcpp::stop("computeG: n[t] must be positive.");
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

// Requires rinvwishart(df, S, gen)
// B_eff is m x nt matrix of effective marker effects
void sampleE_exact_sameX(
  int nt,
  int m,
  int nue,
  arma::mat& E,                   // residual covariance
  const arma::mat& B_eff,         // m x nt effective effects
  const arma::mat& WY,            // m x nt, X'y by trait
  const arma::mat& XX,            // m x m, X'X
  const arma::mat& YY,            // nt x nt, y'y cross-trait
  const arma::mat& S0,            // prior scale
  int N,                          // sample size
  std::mt19937& gen)
{
 arma::mat Se = YY - WY.t() * B_eff - B_eff.t() * WY + B_eff.t() * XX * B_eff;
 Se = 0.5 * (Se + Se.t());

 arma::mat S_post = S0 + Se;
 int df_post = nue + N;

 E = rinvwishart(df_post, S_post, gen);
 E = 0.5 * (E + E.t());
}

// Compute genetic covariance matrix G
// void computeG(
//   int nt,
//   int m,
//   const std::vector<std::vector<double>>& b,
//   const std::vector<std::vector<double>>& wy,
//   const std::vector<std::vector<double>>& r,
//   const std::vector<int>& n,
//   arma::mat& G
// ) {
//  for (int t1 = 0; t1 < nt; t1++) {
//   for (int t2 = t1; t2 < nt; t2++) {
//    double ssg = 0.0; // Initialize the sum of squared G values
//
//    if (t1 == t2) {
//     // Calculate the sum of squared G values for the diagonal elements
//     for (int i = 0; i < m; i++) {
//      ssg += b[t1][i] * (wy[t1][i] - r[t1][i]);
//     }
//
//     // Store the G value in the matrix, normalized by the square root of n[t1] and n[t2]
//     G(t1, t2) = ssg / (std::sqrt(static_cast<double>(n[t1])) * std::sqrt(static_cast<double>(n[t2])));
//    }
//
//    if (t1 != t2) {
//     ssg = 0.0; // Reset the sum of squared G values
//
//     // Calculate the sum of squared G values for off-diagonal elements
//     for (int i = 0; i < m; i++) {
//      ssg += b[t1][i] * (wy[t1][i] - r[t1][i]);
//      ssg += b[t2][i] * (wy[t2][i] - r[t2][i]);
//     }
//
//     // Adjust ssg and store the G value in the matrix, normalized by the square root of n[t1] and n[t2]
//     ssg = ssg / 2.0;
//     G(t1, t2) = ssg / (std::sqrt(static_cast<double>(n[t1])) * std::sqrt(static_cast<double>(n[t2])));
//
//     // Since G is symmetric, set the corresponding off-diagonal element
//     G(t2, t1) = G(t1, t2);
//    }
//   }
//  }
// }

// Sample marker effects using multiple trait BayesC
// Corrected: evaluates each model in the active-trait subspace only.
void sampleBetaCMt(int i,
                   int nt,
                   int nmodels,
                   const std::vector<std::vector<int>>& models,
                   std::vector<double>& cmodel,
                   const std::vector<double>& pi,
                   const arma::mat& Ei,
                   const arma::mat& Bi,
                   const std::vector<std::vector<double>>& ww,
                   std::vector<std::vector<double>>& r,
                   std::vector<std::vector<double>>& b,
                   std::vector<std::vector<int>>& d,
                   const std::vector<std::vector<int>>& XXindices,
                   const std::vector<std::vector<std::vector<double>>>& XXvalues,
                   std::mt19937& gen) {

 // Full RHS, still using diagonal Ei as in your current code
 arma::vec rhs_full(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  rhs_full(t) = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
 }

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------
 // 1) Compute log posterior weight for each model
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  std::vector<int> active;
  active.reserve(nt);
  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) active.push_back(t);
  }

  const int q = static_cast<int>(active.size());

  // Null model
  if (q == 0) {
   loglik[k] = std::log(pi[k]);
   continue;
  }

  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) idx(a) = static_cast<arma::uword>(active[a]);

  // Subspace objects
  arma::vec rhs = rhs_full.elem(idx);
  arma::mat C   = Bi.submat(idx, idx);

  // Add likelihood precision on included traits only
  for (int a = 0; a < q; ++a) {
   int t = active[a];
   C(a, a) += ww[t][i] * Ei(t, t);
  }

  // Robust Cholesky
  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   // small jitter escalation
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCMt: C is not SPD.");
  }

  // log |C|
  double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  // C^{-1} rhs via solves, no explicit inverse needed
  arma::vec y = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  // rhs' C^{-1} rhs
  double quad = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // ------------------------------------------------------------
 // 2) Stabilized normalization
 // ------------------------------------------------------------
 double max_loglik = *std::max_element(loglik.begin(), loglik.end());
 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  denom += std::exp(loglik[k] - max_loglik);
 }
 for (int k = 0; k < nmodels; ++k) {
  pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
 }

 // ------------------------------------------------------------
 // 3) Sample model indicator
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);
 double cumprobc = 0.0;
 int mselect = 0;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 4) Sample marker effects conditional on selected model
 // ------------------------------------------------------------
 std::vector<int> active;
 active.reserve(nt);
 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) active.push_back(t);
 }

 arma::rowvec mub(nt, arma::fill::zeros);

 const int q = static_cast<int>(active.size());
 if (q > 0) {
  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) idx(a) = static_cast<arma::uword>(active[a]);

  arma::vec rhs = rhs_full.elem(idx);
  arma::mat C   = Bi.submat(idx, idx);

  for (int a = 0; a < q; ++a) {
   int t = active[a];
   C(a, a) += ww[t][i] * Ei(t, t);
  }

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCMt: selected-model C is not SPD.");
  }

  // Posterior covariance = C^{-1}
  arma::mat I = arma::eye(q, q);
  arma::mat Ci = arma::solve(arma::trimatu(L.t()),
                             arma::solve(arma::trimatl(L), I));

  // Posterior mean = rhs' C^{-1}
  arma::rowvec mean = rhs.t() * Ci;

  // Draw from N(mean, Ci)
  arma::mat sample = mvrnormARMA(Ci);  // assumed 1 x q
  arma::rowvec mub_sub = sample + mean;

  for (int a = 0; a < q; ++a) {
   mub(active[a]) = mub_sub(a);
  }
 }

 // ------------------------------------------------------------
 // 5) Residual update using sparse XX
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}

// Sample marker effects using multiple-trait BayesC
// Faster version:
//  - no explicit inverse
//  - reuse Cholesky of selected model
//  - solve-based Gaussian sampling
//  - less repeated allocation/work
void sampleBetaCMt_fast(
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
  std::vector<std::vector<double>>& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 0) Build full RHS
 // NOTE:
 // This is correct if r[t][i] corresponds to the partial residual
 // with marker i removed, i.e. r = y - Xb + x_i b_i
 // ------------------------------------------------------------
 arma::vec rhs_full(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  rhs_full(t) = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
 }

 // ------------------------------------------------------------
 // 1) Precompute active-trait indices for each model
 // ------------------------------------------------------------
 std::vector<std::vector<int>> active_traits(nmodels);
 for (int k = 0; k < nmodels; ++k) {
  active_traits[k].reserve(nt);
  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) active_traits[k].push_back(t);
  }
 }

 // model log posterior weights
 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // store Cholesky factors for reuse after model selection
 std::vector<arma::mat> chol_store(nmodels);
 std::vector<bool> chol_ok_store(nmodels, false);

 // ------------------------------------------------------------
 // 2) Compute log posterior weight for each model
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  const std::vector<int>& active = active_traits[k];
  const int q = static_cast<int>(active.size());

  // prior weight must be positive
  if (pi[k] <= 0.0) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  // Null model
  if (q == 0) {
   loglik[k] = std::log(pi[k]);
   continue;
  }

  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) {
   idx(a) = static_cast<arma::uword>(active[a]);
  }

  arma::vec rhs = rhs_full.elem(idx);
  arma::mat C   = Bi.submat(idx, idx);

  // add likelihood precision contribution on active traits
  for (int a = 0; a < q; ++a) {
   int t = active[a];
   C(a, a) += ww[t][i] * Ei(t, t);
  }

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCMt_fast: C is not SPD.");
  }

  // cache for possible later reuse
  chol_store[k] = L;
  chol_ok_store[k] = true;

  // log |C|
  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  // mean_col = C^{-1} rhs using triangular solves
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  // quad = rhs' C^{-1} rhs
  const double quad = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }


 // ------------------------------------------------------------
 // 3) Stabilized normalization
 // ------------------------------------------------------------
 const double max_loglik = *std::max_element(loglik.begin(), loglik.end());

 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (denom <= 0.0 || !std::isfinite(denom)) {
  throw std::runtime_error("sampleBetaCMt_fast: invalid model probability normalization.");
 }

 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  } else {
   pmodel[k] = 0.0;
  }
 }

 // ------------------------------------------------------------
 // 4) Sample model indicator
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);

 double cumprobc = 0.0;
 int mselect = nmodels - 1;  // fallback to last model
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 5) Sample marker effects conditional on selected model
 // ------------------------------------------------------------
 arma::vec mub(nt, arma::fill::zeros);

 const std::vector<int>& active = active_traits[mselect];
 const int q = static_cast<int>(active.size());

 if (q > 0) {
  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) {
   idx(a) = static_cast<arma::uword>(active[a]);
  }

  arma::vec rhs = rhs_full.elem(idx);

  // reuse stored Cholesky
  if (!chol_ok_store[mselect]) {
   throw std::runtime_error("sampleBetaCMt_fast: missing Cholesky for selected model.");
  }
  const arma::mat& L = chol_store[mselect];

  // posterior mean: C^{-1} rhs
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  // sample from N(mean, C^{-1})
  // If z ~ N(0, I), then beta = mean + L^{-T} z
  arma::vec z(q);
  for (int a = 0; a < q; ++a) {
   z(a) = std::normal_distribution<double>(0.0, 1.0)(gen);
  }
  arma::vec draw = mean_col + arma::solve(arma::trimatu(L.t()), z);

  for (int a = 0; a < q; ++a) {
   mub(active[a]) = draw(a);
  }
 }

 // ------------------------------------------------------------
 // 6) Residual update using sparse XX
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  const double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   const size_t nnz = XXindices[i].size();
   for (size_t j = 0; j < nnz; ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}

// void sampleBetaCMt(int i,
//                    int nt,
//                    int nmodels,
//                    const std::vector<std::vector<int>>& models,
//                    std::vector<double>& cmodel,
//                    const std::vector<double>& pi,
//                    const arma::mat& Ei,
//                    const arma::mat& Bi,
//                    const std::vector<std::vector<double>>& ww,
//                    std::vector<std::vector<double>>& r,
//                    std::vector<std::vector<double>>& b,
//                    std::vector<std::vector<int>>& d,
//                    const std::vector<std::vector<int>>& XXindices,
//                    const std::vector<std::vector<std::vector<double>>>& XXvalues,
//                    std::mt19937& gen) {
//
//  std::vector<double> rhs(nt);
//  std::vector<double> loglik(nmodels);
//  std::vector<double> pmodel(nmodels);
//
//  // Compute RHS
//  for (int t = 0; t < nt; t++) {
//   rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
//  }
//
//  // Loop over models to compute log-likelihoods
//  for (int k = 0; k < nmodels; k++) {
//   arma::mat C = Bi;
//   for (int t1 = 0; t1 < nt; t1++) {
//    if (models[k][t1] == 1) {
//     C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//    }
//   }
//
//   // Use log_det for numerical stability
//   double logdet, sign;
//   arma::log_det(logdet, sign, C);
//   arma::mat Ci = arma::inv(C);
//
//   // allow more sparsity
//   // int size = 0;
//   // for(int t = 0; t < nt; t++) size += models[k][t];
//   // double lambda = 0.5;   // sharing penalty
//   // double prior = pi[k] * std::pow(lambda, size);
//   // loglik[k] = -0.5 * logdet + std::log(prior);
//   //
//
//   loglik[k] = -0.5 * logdet + std::log(pi[k]);
//
//   for (int t1 = 0; t1 < nt; t1++) {
//    for (int t2 = t1; t2 < nt; t2++) {
//     if (models[k][t1] == 1 && models[k][t2] == 1) {
//      loglik[k] += 0.5 * rhs[t1] * rhs[t2] * (Ci(t1, t2) + (t1 != t2 ? Ci(t2, t1) : 0.0));
//     }
//    }
//   }
//  }
//
//  // Compute numerically stable posterior model probabilities
//  double max_loglik = *std::max_element(loglik.begin(), loglik.end());
//  double denom = 0.0;
//  for (int k = 0; k < nmodels; k++) {
//   denom += std::exp(loglik[k] - max_loglik);
//  }
//  for (int k = 0; k < nmodels; k++) {
//   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
//  }
//
//  // Sample model indicator
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  double u = runif(gen);
//  double cumprobc = 0.0;
//  int mselect = 0;
//  for (int k = 0; k < nmodels; k++) {
//   cumprobc += pmodel[k];
//   if (u < cumprobc) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel[mselect] += 1.0;
//  for (int t = 0; t < nt; t++) {
//   d[t][i] = models[mselect][t];
//  }
//
//  // Sample marker effect from posterior
//  arma::mat C = Bi;
//  for (int t1 = 0; t1 < nt; t1++) {
//   if (models[mselect][t1] == 1) {
//    C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//   }
//  }
//
//  arma::mat Ci = arma::inv(C);
//  arma::mat sample = mvrnormARMA(Ci);  // 1 × nt
//
//  // Compute posterior mean
//  arma::vec rhs_vec(rhs);
//  arma::rowvec mean = (rhs_vec.t() * Ci);  // 1 × nt
//  arma::rowvec mub = sample + mean;
//
//  for (int t = 0; t < nt; t++) {
//   if (models[mselect][t] == 0) {
//    mub(t) = 0.0;
//   }
//  }
//
//  // Update residuals and marker effect
//  for (int t = 0; t < nt; t++) {
//   double diff = mub(t) - b[t][i];
//   if (diff != 0.0) {
//    for (size_t j = 0; j < XXindices[i].size(); j++) {
//     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//    }
//   }
//   b[t][i] = mub(t);
//  }
// }

// Sample marker effects using multi-trait BayesC
// FULL multivariate residual precision version.
// This is consistent when Ei is full and models allow partial trait subsets.
void sampleBetaCMt_fullE(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,   // full residual precision
  const arma::mat& Bi,   // full prior precision for marker effects
  const std::vector<std::vector<double>>& ww,
  std::vector<std::vector<double>>& r,
  std::vector<std::vector<double>>& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------------
 // Build full working response:
 //   z_t = r[t][i] / ww[t][i] + b[t][i]
 //
 // because r[t][i] is the current score-like residual quantity and
 // ww[t][i] is the diagonal marker precision contribution for trait t.
 // ------------------------------------------------------------------
 arma::vec z_full(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  if (ww[t][i] > 0.0) {
   z_full(t) = r[t][i] / ww[t][i] + b[t][i];
  } else {
   z_full(t) = 0.0;
  }
 }

 // ------------------------------------------------------------------
 // 1) Compute posterior model weights
 // ------------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  std::vector<int> active;
  active.reserve(nt);
  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) active.push_back(t);
  }

  const int q = static_cast<int>(active.size());

  // Null model
  if (q == 0) {
   loglik[k] = std::log(pi[k]);
   continue;
  }

  arma::uvec idx(q);
  for (int a = 0; a < q; ++a)
   idx(a) = static_cast<arma::uword>(active[a]);

  // Subspace quantities
  arma::mat Ei_sub = Ei.submat(idx, idx);
  arma::mat Bi_sub = Bi.submat(idx, idx);
  arma::vec z_sub  = z_full.elem(idx);

  // D = diag(sqrt(ww_i))
  arma::vec sqrtw(q, arma::fill::zeros);
  for (int a = 0; a < q; ++a) {
   sqrtw(a) = std::sqrt(std::max(ww[active[a]][i], 0.0));
  }
  arma::mat D = arma::diagmat(sqrtw);

  // Posterior precision:
  //   C = Bi_sub + D * Ei_sub * D
  arma::mat C = Bi_sub + D * Ei_sub * D;

  // rhs:
  //   rhs = D * Ei_sub * D * z_sub
  arma::vec rhs = D * Ei_sub * D * z_sub;

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCMt_fullE: C is not SPD.");
  }

  double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  arma::vec y = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  double quad = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // Stabilized normalization
 double max_loglik = *std::max_element(loglik.begin(), loglik.end());
 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  denom += std::exp(loglik[k] - max_loglik);
 }
 for (int k = 0; k < nmodels; ++k) {
  pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
 }

 // ------------------------------------------------------------------
 // 2) Sample model
 // ------------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);
 double cumprobc = 0.0;
 int mselect = 0;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------------
 // 3) Sample effects under selected model
 // ------------------------------------------------------------------
 std::vector<int> active;
 active.reserve(nt);
 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) active.push_back(t);
 }

 arma::rowvec mub(nt, arma::fill::zeros);

 const int q = static_cast<int>(active.size());
 if (q > 0) {
  arma::uvec idx(q);
  for (int a = 0; a < q; ++a)
   idx(a) = static_cast<arma::uword>(active[a]);

  arma::mat Ei_sub = Ei.submat(idx, idx);
  arma::mat Bi_sub = Bi.submat(idx, idx);
  arma::vec z_sub  = z_full.elem(idx);

  arma::vec sqrtw(q, arma::fill::zeros);
  for (int a = 0; a < q; ++a) {
   sqrtw(a) = std::sqrt(std::max(ww[active[a]][i], 0.0));
  }
  arma::mat D = arma::diagmat(sqrtw);

  arma::mat C   = Bi_sub + D * Ei_sub * D;
  arma::vec rhs = D * Ei_sub * D * z_sub;

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCMt_fullE: selected-model C is not SPD.");
  }

  arma::mat I = arma::eye(q, q);
  arma::mat Ci = arma::solve(arma::trimatu(L.t()),
                             arma::solve(arma::trimatl(L), I));

  arma::rowvec mean = rhs.t() * Ci;
  arma::mat sample  = mvrnormARMA(Ci);   // assumed 1 x q
  arma::rowvec mub_sub = sample + mean;

  for (int a = 0; a < q; ++a) {
   mub(active[a]) = mub_sub(a);
  }
 }

 // ------------------------------------------------------------------
 // 4) Residual update using sparse XX
 // ------------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}

// // Sample marker effects using multiple trait BayesC
// void sampleBetaCStMt(int i,
//                    int nt,
//                    int nmodels,
//                    const std::vector<std::vector<int>>& models,
//                    std::vector<double>& cmodel,
//                    const std::vector<double>& pi,
//                    const arma::mat& Ei,
//                    const arma::mat& Bi,
//                    const std::vector<std::vector<double>>& ww,
//                    std::vector<std::vector<double>>& r,
//                    std::vector<std::vector<double>>& b,
//                    std::vector<std::vector<int>>& d,
//                    const std::vector<std::vector<int>>& XXindices,
//                    const std::vector<std::vector<std::vector<double>>>& XXvalues,
//                    std::mt19937& gen) {
//
//  std::vector<double> rhs(nt);
//  std::vector<double> loglik(nmodels);
//  std::vector<double> pmodel(nmodels);
//
//  // Compute RHS
//  for (int t = 0; t < nt; t++) {
//   rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
//  }
//
//  // Loop over models to compute log-likelihoods
//  for (int k = 0; k < nmodels; k++) {
//   arma::mat C = Bi;
//   for (int t1 = 0; t1 < nt; t1++) {
//    if (models[k][t1] == 1) {
//     C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//    }
//   }
//
//   // Use log_det for numerical stability
//   double logdet, sign;
//   arma::log_det(logdet, sign, C);
//   arma::mat Ci = arma::inv(C);
//
//   // allow more sparsity
//   // int size = 0;
//   // for(int t = 0; t < nt; t++) size += models[k][t];
//   // double lambda = 0.5;   // sharing penalty
//   // double prior = pi[k] * std::pow(lambda, size);
//   // loglik[k] = -0.5 * logdet + std::log(prior);
//   //
//
//   loglik[k] = -0.5 * logdet + std::log(pi[k]);
//
//   for (int t1 = 0; t1 < nt; t1++) {
//    for (int t2 = t1; t2 < nt; t2++) {
//     if (models[k][t1] == 1 && models[k][t2] == 1) {
//      loglik[k] += 0.5 * rhs[t1] * rhs[t2] * (Ci(t1, t2) + (t1 != t2 ? Ci(t2, t1) : 0.0));
//     }
//    }
//   }
//  }
//
//  // Compute numerically stable posterior model probabilities
//  double max_loglik = *std::max_element(loglik.begin(), loglik.end());
//  double denom = 0.0;
//  for (int k = 0; k < nmodels; k++) {
//   denom += std::exp(loglik[k] - max_loglik);
//  }
//  for (int k = 0; k < nmodels; k++) {
//   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
//  }
//
//  // Sample model indicator
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  double u = runif(gen);
//  double cumprobc = 0.0;
//  int mselect = 0;
//  for (int k = 0; k < nmodels; k++) {
//   cumprobc += pmodel[k];
//   if (u < cumprobc) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel[mselect] += 1.0;
//  for (int t = 0; t < nt; t++) {
//   d[t][i] = models[mselect][t];
//  }
//
//  // Sample marker effect from posterior
//  arma::mat C = Bi;
//  for (int t1 = 0; t1 < nt; t1++) {
//   if (models[mselect][t1] == 1) {
//    C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//   }
//  }
//
//  arma::mat Ci = arma::inv(C);
//  arma::mat sample = mvrnormARMA(Ci);  // 1 × nt
//
//  // Compute posterior mean
//  arma::vec rhs_vec(rhs);
//  arma::rowvec mean = (rhs_vec.t() * Ci);  // 1 × nt
//  arma::rowvec mub = sample + mean;
//
//  for (int t = 0; t < nt; t++) {
//   if (models[mselect][t] == 0) {
//    mub(t) = 0.0;
//   }
//  }
//
//  // Update residuals and marker effect
//  for (int t = 0; t < nt; t++) {
//   double diff = mub(t) - b[t][i];
//   if (diff != 0.0) {
//    for (size_t j = 0; j < XXindices[i].size(); j++) {
//     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//    }
//   }
//   b[t][i] = mub(t);
//  }
// }

// Sample marker effects using multiple-trait BayesC
// Correct version: evaluates each model in the active-trait subspace
void sampleBetaCStMt(
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
  std::vector<std::vector<double>>& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 1) Full RHS
 // Assumes r[t][i] is partial residual with marker i removed:
 // r = y - Xb + x_i b_i
 // ------------------------------------------------------------
 arma::vec rhs_full(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  rhs_full(t) = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
 }

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------
 // 2) Compute log posterior weight for each model
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  std::vector<int> active;
  active.reserve(nt);
  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) active.push_back(t);
  }

  const int q = static_cast<int>(active.size());

  // Null model
  if (q == 0) {
   loglik[k] = std::log(pi[k]);
   continue;
  }

  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) {
   idx(a) = static_cast<arma::uword>(active[a]);
  }

  // Active-subspace RHS and precision matrix
  arma::vec rhs = rhs_full.elem(idx);
  arma::mat C   = Bi.submat(idx, idx);

  // Add likelihood precision on active traits only
  for (int a = 0; a < q; ++a) {
   int t = active[a];
   C(a, a) += ww[t][i] * Ei(t, t);
  }

  // Cholesky factorization for SPD check and logdet/solve
  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  // log |C|
  double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  // mean_col = C^{-1} rhs
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  // rhs' C^{-1} rhs
  double quad = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // ------------------------------------------------------------
 // 3) Stabilized normalization
 // ------------------------------------------------------------
 double max_loglik = *std::max_element(loglik.begin(), loglik.end());

 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (denom <= 0.0 || !std::isfinite(denom)) {
  throw std::runtime_error("sampleBetaCStMt: invalid model probability normalization.");
 }

 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  } else {
   pmodel[k] = 0.0;
  }
 }

 // ------------------------------------------------------------
 // 4) Sample model indicator
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);

 double cumprobc = 0.0;
 int mselect = nmodels - 1;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 5) Sample marker effects conditional on selected model
 // ------------------------------------------------------------
 arma::vec mub(nt, arma::fill::zeros);

 std::vector<int> active;
 active.reserve(nt);
 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) active.push_back(t);
 }

 const int q = static_cast<int>(active.size());

 if (q > 0) {
  arma::uvec idx(q);
  for (int a = 0; a < q; ++a) {
   idx(a) = static_cast<arma::uword>(active[a]);
  }

  arma::vec rhs = rhs_full.elem(idx);
  arma::mat C   = Bi.submat(idx, idx);

  for (int a = 0; a < q; ++a) {
   int t = active[a];
   C(a, a) += ww[t][i] * Ei(t, t);
  }

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   throw std::runtime_error("sampleBetaCStMt: selected-model C is not SPD.");
  }

  // posterior mean = C^{-1} rhs
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  // sample from N(mean, C^{-1})
  // if z ~ N(0, I), then draw = mean + L^{-T} z
  arma::vec z(q);
  for (int a = 0; a < q; ++a) {
   z(a) = std::normal_distribution<double>(0.0, 1.0)(gen);
  }
  arma::vec draw = mean_col + arma::solve(arma::trimatu(L.t()), z);

  for (int a = 0; a < q; ++a) {
   mub(active[a]) = draw(a);
  }
 }

 // ------------------------------------------------------------
 // 6) Residual update using sparse XX
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}

// Sample marker effects using multiple-trait BayesCP-G style logic
// Sufficient-statistics version under diagonal Ei = R^{-1}
// Here b[t][i] stores the EFFECTIVE marker effect D_j * beta_j
void sampleBetaCPG_Mt(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,   // assumed diagonal residual precision
  const arma::mat& Bi,   // prior precision = G^{-1}
  const std::vector<std::vector<double>>& ww,   // sum x_ij^2 by trait
  std::vector<std::vector<double>>& r,          // partial residuals
  std::vector<std::vector<double>>& b,          // stores effective effects D_j * beta_j
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 1) Build base RHS contribution
 //
 // Assumes r[t][i] is the partial residual with the current
 // EFFECTIVE marker effect already added back:
 //
 //   r = y - Xb + x_i * b_i_effective
 //
 // Then for a given model, the paper's r_j vector is obtained by
 // zeroing out inactive traits through D_j.
 // ------------------------------------------------------------
 arma::vec rhs_base(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  rhs_base(t) = Ei(t, t) * (r[t][i] + ww[t][i] * b[t][i]);
 }

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------
 // 2) Compute log posterior weight for each model
 //
 // For model k:
 //   C = G^{-1} + D' R^{-1} D * sum(x_ij^2)
 //   rhs = D * base_rhs
 //
 // marginal weight proportional to
 //   |C^{-1}|^{1/2} exp( 1/2 rhs' C^{-1} rhs ) p(model)
 //
 // equivalently:
 //   -1/2 log|C| + 1/2 rhs' C^{-1} rhs + log pi[k]
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  arma::vec rhs(nt, arma::fill::zeros);
  arma::mat C = Bi;

  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) {
    rhs(t) = rhs_base(t);
    C(t, t) += ww[t][i] * Ei(t, t);
   }
  }

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }
  if (!ok) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  // solve C^{-1} rhs
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

  const double quad = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // ------------------------------------------------------------
 // 3) Stabilized normalization
 // ------------------------------------------------------------
 const double max_loglik = *std::max_element(loglik.begin(), loglik.end());

 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (denom <= 0.0 || !std::isfinite(denom)) {
  throw std::runtime_error("sampleBetaCPG_Mt: invalid model probability normalization.");
 }

 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  } else {
   pmodel[k] = 0.0;
  }
 }

 // ------------------------------------------------------------
 // 4) Sample model indicator
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);

 double cumprobc = 0.0;
 int mselect = nmodels - 1;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 5) Sample latent beta_j | d_j, ...
 //
 // Full latent draw:
 //   beta_latent ~ N(C^{-1} rhs, C^{-1})
 //
 // Effective marker effect is:
 //   beta_eff = D * beta_latent
 // ------------------------------------------------------------
 arma::vec rhs(nt, arma::fill::zeros);
 arma::mat C = Bi;

 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) {
   rhs(t) = rhs_base(t);
   C(t, t) += ww[t][i] * Ei(t, t);
  }
 }

 arma::mat L;
 bool ok = arma::chol(L, C, "lower");
 if (!ok) {
  double jitter = 1e-8;
  for (int tries = 0; tries < 7 && !ok; ++tries) {
   C.diag() += jitter;
   ok = arma::chol(L, C, "lower");
   jitter *= 10.0;
  }
 }
 if (!ok) {
  throw std::runtime_error("sampleBetaCPG_Mt: selected-model C is not SPD.");
 }

 // posterior mean = C^{-1} rhs
 arma::vec y        = arma::solve(arma::trimatl(L), rhs);
 arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

 // sample latent beta: mean + L^{-T} z
 arma::vec z(nt);
 for (int t = 0; t < nt; ++t) {
  z(t) = std::normal_distribution<double>(0.0, 1.0)(gen);
 }
 arma::vec beta_latent = mean_col + arma::solve(arma::trimatu(L.t()), z);

 // effective marker effect = D * beta_latent
 arma::vec beta_eff(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) {
   beta_eff(t) = beta_latent(t);
  }
 }

 // ------------------------------------------------------------
 // 6) Residual update using EFFECTIVE effect
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  const double diff = beta_eff(t) - b[t][i];
  if (diff != 0.0) {
   const size_t nnz = XXindices[i].size();
   for (size_t j = 0; j < nnz; ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = beta_eff(t);
 }
}
void sampleBetaCMtMaskedFast(
  int i,
  int region,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  const std::vector<std::vector<int>>& trait_allowed,
  bool warmup_st_only,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const std::vector<ModelCache>& cache,
  double ww_const,
  std::vector<std::vector<double>>& r,
  std::vector<std::vector<double>>& b,
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 const double mt_improve_threshold = 2.0;

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);
 std::vector<int> model_size(nmodels, 0);
 std::vector<double> rhs(nt);

 // RHS using ww_const to match cached matrices
 for (int t = 0; t < nt; t++) {
  rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww_const * b[t][i];
 }

 // Evaluate candidate models
 for (int k = 0; k < nmodels; k++) {
  bool valid = true;
  int size = 0;

  for (int t = 0; t < nt; t++) {
   if (models[k][t] == 1) {
    size++;
    if (trait_allowed[region][t] == 0) {
     valid = false;
     break;
    }
   }
  }

  model_size[k] = size;

  if (warmup_st_only && size > 1) valid = false;
  if (!valid) continue;

  const ModelCache& mc = cache[k];
  if (mc.U.n_elem == 0) continue;

  double quad = quadratic_form(mc.U, rhs);

  loglik[k] =
   -0.5 * mc.logdet +
   std::log(std::max(pi[k], 1e-300)) +
   0.5 * quad;
 }

 // ST baseline: MT must improve over best ST
 double best_st_loglik = -std::numeric_limits<double>::infinity();
 for (int k = 0; k < nmodels; k++) {
  if (model_size[k] <= 1 && std::isfinite(loglik[k])) {
   best_st_loglik = std::max(best_st_loglik, loglik[k]);
  }
 }

 for (int k = 0; k < nmodels; k++) {
  if (model_size[k] > 1 && std::isfinite(loglik[k])) {
   if (loglik[k] < best_st_loglik + mt_improve_threshold) {
    loglik[k] = -std::numeric_limits<double>::infinity();
   }
  }
 }

 // Normalize probabilities
 double max_loglik = -std::numeric_limits<double>::infinity();
 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k])) {
   max_loglik = std::max(max_loglik, loglik[k]);
  }
 }

 // Fallback to null model
 if (!std::isfinite(max_loglik)) {
  int null_model = 0;
  for (int k = 0; k < nmodels; k++) {
   if (model_size[k] == 0) {
    null_model = k;
    break;
   }
  }

  cmodel[null_model] += 1.0;
  for (int t = 0; t < nt; t++) d[t][i] = 0;

  for (int t = 0; t < nt; t++) {
   double diff = -b[t][i];
   if (diff != 0.0) {
    for (size_t j = 0; j < XXindices[i].size(); j++) {
     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
    }
   }
   b[t][i] = 0.0;
  }
  return;
 }

 double denom = 0.0;
 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (!(denom > 0.0) || !std::isfinite(denom)) return;

 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  }
 }

 // Sample model
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);
 double cum = 0.0;
 int mselect = 0;

 for (int k = 0; k < nmodels; k++) {
  cum += pmodel[k];
  if (u < cum) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; t++) d[t][i] = models[mselect][t];

 // Posterior mean and draw
 const ModelCache& mc = cache[mselect];
 if (mc.U.n_elem == 0) return;

 arma::vec mean = posterior_mean(mc.U, rhs);
 arma::vec draw = sample_posterior(mc.U, gen);
 arma::vec mub  = mean + draw;

 for (int t = 0; t < nt; t++) {
  if (models[mselect][t] == 0) mub(t) = 0.0;
 }

 // Update residuals
 for (int t = 0; t < nt; t++) {
  double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); j++) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}


// Sample marker effects using multiple trait BayesC
// with region-specific trait masks, optional ST warmup,
// and ST baseline (MT only if improvement)
// Optimized version: no explicit matrix inversion
void sampleBetaCMtMasked(int i,
                         int region,
                         int nt,
                         int nmodels,
                         const std::vector<std::vector<int>>& models,
                         const std::vector<std::vector<int>>& trait_allowed,
                         bool warmup_st_only,
                         std::vector<double>& cmodel,
                         const std::vector<double>& pi,
                         const arma::mat& Ei,
                         const arma::mat& Bi,
                         const std::vector<std::vector<double>>& ww,
                         std::vector<std::vector<double>>& r,
                         std::vector<std::vector<double>>& b,
                         std::vector<std::vector<int>>& d,
                         const std::vector<std::vector<int>>& XXindices,
                         const std::vector<std::vector<std::vector<double>>>& XXvalues,
                         std::mt19937& gen) {

 const double mt_improve_threshold = 2.0;

 std::vector<double> rhs_std(nt);
 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);
 std::vector<int> model_size(nmodels, 0);

 arma::vec rhs(nt);

 // Compute RHS
 for (int t = 0; t < nt; t++) {
  rhs_std[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
  rhs(t) = rhs_std[t];
 }

 // ============================
 // Evaluate all allowed models
 // ============================
 for (int k = 0; k < nmodels; k++) {

  bool valid = true;
  int size = 0;

  for (int t = 0; t < nt; t++) {
   if (models[k][t] == 1) {
    size++;
    if (trait_allowed[region][t] == 0) {
     valid = false;
     break;
    }
   }
  }

  model_size[k] = size;

  if (warmup_st_only && size > 1) {
   valid = false;
  }

  if (!valid) {
   continue;
  }

  arma::mat C = Bi;

  for (int t = 0; t < nt; t++) {
   if (models[k][t] == 1) {
    C(t, t) += ww[t][i] * Ei(t, t);
   }
  }

  arma::mat U;
  bool chol_ok = arma::chol(U, C);   // U is upper triangular: C = U.t() * U

  if (!chol_ok) {
   continue;
  }

  // log |C| = 2 * sum(log(diag(U)))
  double logdet = 2.0 * arma::sum(arma::log(U.diag()));

  // Solve C^{-1} rhs using Cholesky
  // First solve U.t() y = rhs, then U x = y
  arma::vec y = arma::solve(arma::trimatl(U.t()), rhs);
  arma::vec x = arma::solve(arma::trimatu(U), y);

  double quad = arma::dot(rhs, x);

  loglik[k] = -0.5 * logdet + std::log(std::max(pi[k], 1e-300)) + 0.5 * quad;
 }

 // ==========================================
 // Best ST model (null + single-trait models)
 // ==========================================
 double best_st_loglik = -std::numeric_limits<double>::infinity();

 for (int k = 0; k < nmodels; k++) {
  if (model_size[k] <= 1 && std::isfinite(loglik[k])) {
   best_st_loglik = std::max(best_st_loglik, loglik[k]);
  }
 }

 // ===========================
 // MT must improve over ST
 // ===========================
 for (int k = 0; k < nmodels; k++) {
  if (model_size[k] > 1 && std::isfinite(loglik[k])) {
   if (loglik[k] < best_st_loglik + mt_improve_threshold) {
    loglik[k] = -std::numeric_limits<double>::infinity();
   }
  }
 }

 // ===========================
 // Normalize probabilities
 // ===========================
 double max_loglik = -std::numeric_limits<double>::infinity();

 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k]) && loglik[k] > max_loglik) {
   max_loglik = loglik[k];
  }
 }

 // Fallback: all invalid -> null model
 if (!std::isfinite(max_loglik)) {

  int null_model = 0;
  for (int k = 0; k < nmodels; k++) {
   if (model_size[k] == 0) {
    null_model = k;
    break;
   }
  }

  cmodel[null_model] += 1.0;

  for (int t = 0; t < nt; t++) {
   d[t][i] = 0;
  }

  for (int t = 0; t < nt; t++) {
   double diff = -b[t][i];
   if (diff != 0.0) {
    for (size_t j = 0; j < XXindices[i].size(); j++) {
     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
    }
   }
   b[t][i] = 0.0;
  }

  return;
 }

 double denom = 0.0;
 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 for (int k = 0; k < nmodels; k++) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  }
 }

 // ===========================
 // Sample model indicator
 // ===========================
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);
 double cumprobc = 0.0;
 int mselect = 0;

 for (int k = 0; k < nmodels; k++) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;

 for (int t = 0; t < nt; t++) {
  d[t][i] = models[mselect][t];
 }

 // ===========================
 // Sample marker effects
 // ===========================
 arma::mat C = Bi;

 for (int t = 0; t < nt; t++) {
  if (models[mselect][t] == 1) {
   C(t, t) += ww[t][i] * Ei(t, t);
  }
 }

 arma::mat U;
 bool chol_ok = arma::chol(U, C);

 if (!chol_ok) {
  // conservative fallback: keep old effects
  return;
 }

 // Posterior mean: C^{-1} rhs
 arma::vec y = arma::solve(arma::trimatl(U.t()), rhs);
 arma::vec mean = arma::solve(arma::trimatu(U), y);

 // Sample from N(0, C^{-1})
 // If z ~ N(0,I), then solve(U, z) ~ N(0, (U'U)^-1 ) = N(0, C^{-1})
 arma::vec z = arma::randn(nt);
 arma::vec sample = arma::solve(arma::trimatu(U), z);

 arma::vec mub = mean + sample;

 for (int t = 0; t < nt; t++) {
  if (models[mselect][t] == 0) {
   mub(t) = 0.0;
  }
 }

 // ===========================
 // Update residuals and effects
 // ===========================
 for (int t = 0; t < nt; t++) {

  double diff = mub(t) - b[t][i];

  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); j++) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }

  b[t][i] = mub(t);
 }
}

// Sample marker effects using multiple-trait BayesCP-G
// Fully faithful latent-state version:
//
//   beta[t][i] = latent marker effect b_j
//   b[t][i]    = effective marker effect D_j * b_j
//
// Assumes:
//   - Ei is diagonal residual precision R^{-1}
//   - Bi is prior precision G^{-1}
//   - r contains residuals with current EFFECTIVE effect added back:
//         r = y - X b_eff + x_i * b_eff_i
//
template <class DataView>
void sampleBetaCPG_Mt_latent(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,
  const arma::mat& Bi,
  const DataView& data,
  std::vector<std::vector<double>>& r,
  std::vector<std::vector<double>>& beta,   // latent effects
  std::vector<std::vector<double>>& b,      // effective effects = D * beta
  std::vector<std::vector<int>>& d,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 1) Base RHS using current EFFECTIVE effect
 // ------------------------------------------------------------
 arma::vec rhs_base(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  rhs_base(t) = Ei(t, t) *
   (r[t][i] + sblr::mt::mt_diagonal(data, t, i) * b[t][i]);
 }

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------
 // 2) Compute log posterior weight for each model
 //
 // For model k:
 //   rhs_k = D_k * rhs_base
 //   C_k   = G^{-1} + D_k' R^{-1} D_k * sum(x_ij^2)
 //
 // log weight:
 //   log pi_k - 1/2 log|C_k| + 1/2 rhs_k' C_k^{-1} rhs_k
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  arma::vec rhs(nt, arma::fill::zeros);
  arma::mat C = Bi;

  for (int t = 0; t < nt; ++t) {
   if (models[k][t] == 1) {
    rhs(t) = rhs_base(t);
    C(t, t) += sblr::mt::mt_diagonal(data, t, i) * Ei(t, t);
   }
  }

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }

  if (!ok) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));

  // C^{-1} rhs via triangular solves
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);
  const double quad  = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // ------------------------------------------------------------
 // 3) Stabilized normalization
 // ------------------------------------------------------------
 const double max_loglik = *std::max_element(loglik.begin(), loglik.end());

 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (denom <= 0.0 || !std::isfinite(denom)) {
  throw std::runtime_error("sampleBetaCPG_Mt_latent: invalid model probability normalization.");
 }

 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
  } else {
   pmodel[k] = 0.0;
  }
 }

 // ------------------------------------------------------------
 // 4) Sample model indicator
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);

 double cumprobc = 0.0;
 int mselect = nmodels - 1;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 5) Sample latent beta_j | d_j, ...
 //
 // beta_latent ~ N(C^{-1} rhs, C^{-1})
 // ------------------------------------------------------------
 arma::vec rhs(nt, arma::fill::zeros);
 arma::mat C = Bi;

 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) {
   rhs(t) = rhs_base(t);
   C(t, t) += sblr::mt::mt_diagonal(data, t, i) * Ei(t, t);
  }
 }

 arma::mat L;
 bool ok = arma::chol(L, C, "lower");
 if (!ok) {
  double jitter = 1e-8;
  for (int tries = 0; tries < 7 && !ok; ++tries) {
   C.diag() += jitter;
   ok = arma::chol(L, C, "lower");
   jitter *= 10.0;
  }
 }

 if (!ok) {
  throw std::runtime_error("sampleBetaCPG_Mt_latent: selected-model C is not SPD.");
 }

 // posterior mean = C^{-1} rhs
 arma::vec y        = arma::solve(arma::trimatl(L), rhs);
 arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

 // draw latent beta = mean + L^{-T} z
 arma::vec z(nt);
 std::normal_distribution<double> rnorm(0.0, 1.0);
 for (int t = 0; t < nt; ++t) {
  z(t) = rnorm(gen);
 }

 arma::vec beta_new = mean_col + arma::solve(arma::trimatu(L.t()), z);

 // effective effect = D * beta
 arma::vec b_new(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  if (models[mselect][t] == 1) {
   b_new(t) = beta_new(t);
  }
 }

 // ------------------------------------------------------------
 // 6) Residual update using EFFECTIVE effect only
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  const double diff = b_new(t) - b[t][i];
  if (diff != 0.0) {
   sblr::mt::mt_apply_marker_difference(data, t, i, diff, r[t]);
  }
 }

 // ------------------------------------------------------------
 // 7) Store latent and effective effects
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  beta[t][i] = beta_new(t);  // latent
  b[t][i]    = b_new(t);     // effective
 }
}

void sampleBetaCPG_Mt_latent_fullR(
  int i,
  int nt,
  int nmodels,
  const std::vector<std::vector<int>>& models,
  std::vector<double>& cmodel,
  const std::vector<double>& pi,
  const arma::mat& Ei,   // full residual precision R^{-1}
  const arma::mat& Bi,   // prior precision G^{-1}
  const std::vector<std::vector<double>>& ww,  // trait-specific marker information
  std::vector<std::vector<double>>& r,         // score residuals
  std::vector<std::vector<double>>& beta,      // latent effects
  std::vector<std::vector<double>>& b,         // effective effects
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 1) Build score vector s and sqrt-information matrix
 // ------------------------------------------------------------
 arma::vec s(nt, arma::fill::zeros);
 arma::vec wdiag(nt, arma::fill::zeros);

 for (int t = 0; t < nt; ++t) {
  wdiag(t) = ww[t][i];
  s(t) = r[t][i] + ww[t][i] * b[t][i];
 }

 arma::mat Wsqrt = arma::diagmat(arma::sqrt(wdiag));

 std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
 std::vector<double> pmodel(nmodels, 0.0);

 // ------------------------------------------------------------
 // 2) Model probabilities
 // ------------------------------------------------------------
 for (int k = 0; k < nmodels; ++k) {

  if (pi[k] <= 0.0) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  arma::vec dvec(nt, arma::fill::zeros);
  for (int t = 0; t < nt; ++t) {
   dvec(t) = static_cast<double>(models[k][t]);
  }
  arma::mat D = arma::diagmat(dvec);

  arma::vec rhs = D * (Ei * s);
  arma::mat C   = Bi + D * Wsqrt * Ei * Wsqrt * D;

  arma::mat L;
  bool ok = arma::chol(L, C, "lower");
  if (!ok) {
   double jitter = 1e-8;
   for (int tries = 0; tries < 7 && !ok; ++tries) {
    C.diag() += jitter;
    ok = arma::chol(L, C, "lower");
    jitter *= 10.0;
   }
  }

  if (!ok) {
   loglik[k] = -std::numeric_limits<double>::infinity();
   continue;
  }

  const double logdet = 2.0 * arma::sum(arma::log(L.diag()));
  arma::vec y        = arma::solve(arma::trimatl(L), rhs);
  arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);
  const double quad  = arma::dot(rhs, mean_col);

  loglik[k] = std::log(pi[k]) - 0.5 * logdet + 0.5 * quad;
 }

 // ------------------------------------------------------------
 // 3) Normalize
 // ------------------------------------------------------------
 const double max_loglik = *std::max_element(loglik.begin(), loglik.end());

 double denom = 0.0;
 for (int k = 0; k < nmodels; ++k) {
  if (std::isfinite(loglik[k])) {
   denom += std::exp(loglik[k] - max_loglik);
  }
 }

 if (denom <= 0.0 || !std::isfinite(denom)) {
  throw std::runtime_error("sampleBetaCPG_Mt_latent_fullR: invalid model probability normalization.");
 }

 for (int k = 0; k < nmodels; ++k) {
  pmodel[k] = std::isfinite(loglik[k]) ? std::exp(loglik[k] - max_loglik) / denom : 0.0;
 }

 // ------------------------------------------------------------
 // 4) Sample model
 // ------------------------------------------------------------
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 const double u = runif(gen);

 double cumprobc = 0.0;
 int mselect = nmodels - 1;
 for (int k = 0; k < nmodels; ++k) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; ++t) {
  d[t][i] = models[mselect][t];
 }

 // ------------------------------------------------------------
 // 5) Sample latent beta for selected model
 // ------------------------------------------------------------
 arma::vec dvec(nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  dvec(t) = static_cast<double>(models[mselect][t]);
 }
 arma::mat D = arma::diagmat(dvec);

 arma::vec rhs = D * (Ei * s);
 arma::mat C   = Bi + D * Wsqrt * Ei * Wsqrt * D;

 arma::mat L;
 bool ok = arma::chol(L, C, "lower");
 if (!ok) {
  double jitter = 1e-8;
  for (int tries = 0; tries < 7 && !ok; ++tries) {
   C.diag() += jitter;
   ok = arma::chol(L, C, "lower");
   jitter *= 10.0;
  }
 }

 if (!ok) {
  throw std::runtime_error("sampleBetaCPG_Mt_latent_fullR: selected-model C is not SPD.");
 }

 arma::vec y        = arma::solve(arma::trimatl(L), rhs);
 arma::vec mean_col = arma::solve(arma::trimatu(L.t()), y);

 arma::vec z(nt);
 std::normal_distribution<double> rnorm(0.0, 1.0);
 for (int t = 0; t < nt; ++t) {
  z(t) = rnorm(gen);
 }

 arma::vec beta_new = mean_col + arma::solve(arma::trimatu(L.t()), z);
 arma::vec b_new    = D * beta_new;

 // ------------------------------------------------------------
 // 6) Residual update using effective effects
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  const double diff = b_new(t) - b[t][i];
  if (diff != 0.0) {
   const size_t nnz = XXindices[i].size();
   for (size_t j = 0; j < nnz; ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
 }

 // ------------------------------------------------------------
 // 7) Store
 // ------------------------------------------------------------
 for (int t = 0; t < nt; ++t) {
  beta[t][i] = beta_new(t);
  b[t][i]    = b_new(t);
 }
}

// Sample marker effects using multiple trait BayesC
// with region-specific trait masks, optional ST warmup,
// and ST baseline (MT only if improvement)
// void sampleBetaCMtMasked(int i,
//                          int region,
//                          int nt,
//                          int nmodels,
//                          const std::vector<std::vector<int>>& models,
//                          const std::vector<std::vector<int>>& trait_allowed,
//                          bool warmup_st_only,
//                          std::vector<double>& cmodel,
//                          const std::vector<double>& pi,
//                          const arma::mat& Ei,
//                          const arma::mat& Bi,
//                          const std::vector<std::vector<double>>& ww,
//                          std::vector<std::vector<double>>& r,
//                          std::vector<std::vector<double>>& b,
//                          std::vector<std::vector<int>>& d,
//                          const std::vector<std::vector<int>>& XXindices,
//                          const std::vector<std::vector<std::vector<double>>>& XXvalues,
//                          std::mt19937& gen) {
//
//  const double mt_improve_threshold = 2.0;
//
//  std::vector<double> rhs(nt);
//  std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
//  std::vector<double> pmodel(nmodels, 0.0);
//  std::vector<int> model_size(nmodels, 0);
//
//  // Compute RHS
//  for (int t = 0; t < nt; t++) {
//   rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
//  }
//
//  // Loop over models
//  for (int k = 0; k < nmodels; k++) {
//
//   bool valid = true;
//   int size = 0;
//
//   for (int t = 0; t < nt; t++) {
//
//    if (models[k][t] == 1) {
//
//     size++;
//
//     if (trait_allowed[region][t] == 0) {
//      valid = false;
//      break;
//     }
//    }
//   }
//
//   model_size[k] = size;
//
//   if (warmup_st_only && size > 1) {
//    valid = false;
//   }
//
//   if (!valid) {
//    continue;
//   }
//
//   arma::mat C = Bi;
//
//   for (int t = 0; t < nt; t++) {
//    if (models[k][t] == 1) {
//     C(t, t) += ww[t][i] * Ei(t, t);
//    }
//   }
//
//   double logdet, sign;
//   arma::log_det(logdet, sign, C);
//
//   if (!(sign > 0) || !std::isfinite(logdet)) {
//    continue;
//   }
//
//   arma::mat Ci = arma::inv(C);
//
//   double lk = -0.5 * logdet + std::log(std::max(pi[k], 1e-300));
//
//   for (int t1 = 0; t1 < nt; t1++) {
//    for (int t2 = t1; t2 < nt; t2++) {
//
//     if (models[k][t1] == 1 && models[k][t2] == 1) {
//
//      lk += 0.5 * rhs[t1] * rhs[t2] *
//       (Ci(t1, t2) + (t1 != t2 ? Ci(t2, t1) : 0.0));
//     }
//    }
//   }
//
//   loglik[k] = lk;
//  }
//
//  // ------------------------------------------------
//  // Find best ST model (null + single-trait models)
//  // ------------------------------------------------
//
//  double best_st_loglik = -std::numeric_limits<double>::infinity();
//
//  for (int k = 0; k < nmodels; k++) {
//
//   if (model_size[k] <= 1 && std::isfinite(loglik[k])) {
//
//    best_st_loglik = std::max(best_st_loglik, loglik[k]);
//   }
//  }
//
//  // ------------------------------------------------
//  // Apply MT improvement rule
//  // ------------------------------------------------
//
//  for (int k = 0; k < nmodels; k++) {
//
//   if (model_size[k] > 1 && std::isfinite(loglik[k])) {
//
//    if (loglik[k] < best_st_loglik + mt_improve_threshold) {
//
//     loglik[k] = -std::numeric_limits<double>::infinity();
//    }
//   }
//  }
//
//  // ------------------------------------------------
//  // Compute posterior probabilities
//  // ------------------------------------------------
//
//  double max_loglik = -std::numeric_limits<double>::infinity();
//
//  for (int k = 0; k < nmodels; k++) {
//   if (std::isfinite(loglik[k]) && loglik[k] > max_loglik) {
//    max_loglik = loglik[k];
//   }
//  }
//
//  // Fallback if everything invalid
//  if (!std::isfinite(max_loglik)) {
//
//   int null_model = 0;
//
//   for (int k = 0; k < nmodels; k++) {
//    if (model_size[k] == 0) {
//     null_model = k;
//     break;
//    }
//   }
//
//   cmodel[null_model] += 1.0;
//
//   for (int t = 0; t < nt; t++) d[t][i] = 0;
//
//   for (int t = 0; t < nt; t++) {
//
//    double diff = -b[t][i];
//
//    if (diff != 0.0) {
//     for (size_t j = 0; j < XXindices[i].size(); j++) {
//      r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//     }
//    }
//
//    b[t][i] = 0.0;
//   }
//
//   return;
//  }
//
//  double denom = 0.0;
//
//  for (int k = 0; k < nmodels; k++) {
//   if (std::isfinite(loglik[k])) {
//    denom += std::exp(loglik[k] - max_loglik);
//   }
//  }
//
//  for (int k = 0; k < nmodels; k++) {
//
//   if (std::isfinite(loglik[k])) {
//    pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
//   }
//  }
//
//  // ------------------------------------------------
//  // Sample model
//  // ------------------------------------------------
//
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  double u = runif(gen);
//  double cumprobc = 0.0;
//  int mselect = 0;
//
//  for (int k = 0; k < nmodels; k++) {
//
//   cumprobc += pmodel[k];
//
//   if (u < cumprobc) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel[mselect] += 1.0;
//
//  for (int t = 0; t < nt; t++) {
//   d[t][i] = models[mselect][t];
//  }
//
//  // ------------------------------------------------
//  // Sample marker effects
//  // ------------------------------------------------
//
//  arma::mat C = Bi;
//
//  for (int t = 0; t < nt; t++) {
//   if (models[mselect][t] == 1) {
//    C(t, t) += ww[t][i] * Ei(t, t);
//   }
//  }
//
//  arma::mat Ci = arma::inv(C);
//  arma::mat sample = mvrnormARMA(Ci);
//
//  arma::vec rhs_vec(rhs);
//  arma::rowvec mean = rhs_vec.t() * Ci;
//  arma::rowvec mub = sample + mean;
//
//  for (int t = 0; t < nt; t++) {
//   if (models[mselect][t] == 0) {
//    mub(t) = 0.0;
//   }
//  }
//
//  // Update residuals
//
//  for (int t = 0; t < nt; t++) {
//
//   double diff = mub(t) - b[t][i];
//
//   if (diff != 0.0) {
//
//    for (size_t j = 0; j < XXindices[i].size(); j++) {
//
//     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//    }
//   }
//
//   b[t][i] = mub(t);
//  }
// }
//
// Sample marker effects using multiple trait BayesC
// with region-specific trait masks and optional ST-only warmup
// void sampleBetaCMtMasked(int i,
//                          int region,
//                          int nt,
//                          int nmodels,
//                          const std::vector<std::vector<int>>& models,
//                          const std::vector<std::vector<int>>& trait_allowed,
//                          bool warmup_st_only,
//                          std::vector<double>& cmodel,
//                          const std::vector<double>& pi,
//                          const arma::mat& Ei,
//                          const arma::mat& Bi,
//                          const std::vector<std::vector<double>>& ww,
//                          std::vector<std::vector<double>>& r,
//                          std::vector<std::vector<double>>& b,
//                          std::vector<std::vector<int>>& d,
//                          const std::vector<std::vector<int>>& XXindices,
//                          const std::vector<std::vector<std::vector<double>>>& XXvalues,
//                          std::mt19937& gen) {
//
//  std::vector<double> rhs(nt);
//  std::vector<double> loglik(nmodels, -std::numeric_limits<double>::infinity());
//  std::vector<double> pmodel(nmodels, 0.0);
//
//  // Compute RHS
//  for (int t = 0; t < nt; t++) {
//   rhs[t] = Ei(t, t) * r[t][i] + Ei(t, t) * ww[t][i] * b[t][i];
//  }
//
//  // Loop over models to compute log-likelihoods
//  for (int k = 0; k < nmodels; k++) {
//
//   // Check whether model k is allowed for this region
//   bool valid = true;
//   int model_size = 0;
//
//   for (int t = 0; t < nt; t++) {
//    if (models[k][t] == 1) {
//     model_size++;
//
//     // Trait not allowed in this region
//     if (trait_allowed[region][t] == 0) {
//      valid = false;
//      break;
//     }
//    }
//   }
//
//   // Warmup phase: allow only ST models (null + single-trait models)
//   if (warmup_st_only && model_size > 1) {
//    valid = false;
//   }
//
//   if (!valid) {
//    continue;
//   }
//
//   arma::mat C = Bi;
//   for (int t1 = 0; t1 < nt; t1++) {
//    if (models[k][t1] == 1) {
//     C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//    }
//   }
//
//   // Use log_det for numerical stability
//   double logdet, sign;
//   arma::log_det(logdet, sign, C);
//
//   // Guard against numerical issues
//   if (!(sign > 0) || !std::isfinite(logdet)) {
//    continue;
//   }
//
//   arma::mat Ci = arma::inv(C);
//
//   loglik[k] = -0.5 * logdet + std::log(std::max(pi[k], 1e-300));
//
//   for (int t1 = 0; t1 < nt; t1++) {
//    for (int t2 = t1; t2 < nt; t2++) {
//     if (models[k][t1] == 1 && models[k][t2] == 1) {
//      loglik[k] += 0.5 * rhs[t1] * rhs[t2] *
//       (Ci(t1, t2) + (t1 != t2 ? Ci(t2, t1) : 0.0));
//     }
//    }
//   }
//  }
//
//  // Compute numerically stable posterior model probabilities
//  double max_loglik = -std::numeric_limits<double>::infinity();
//  for (int k = 0; k < nmodels; k++) {
//   if (std::isfinite(loglik[k]) && loglik[k] > max_loglik) {
//    max_loglik = loglik[k];
//   }
//  }
//
//  // Fallback: if all models were excluded, use null model if available
//  if (!std::isfinite(max_loglik)) {
//   int null_model = 0;
//   for (int k = 0; k < nmodels; k++) {
//    int size = 0;
//    for (int t = 0; t < nt; t++) size += models[k][t];
//    if (size == 0) {
//     null_model = k;
//     break;
//    }
//   }
//
//   cmodel[null_model] += 1.0;
//   for (int t = 0; t < nt; t++) d[t][i] = 0;
//
//   for (int t = 0; t < nt; t++) {
//    double diff = -b[t][i];
//    if (diff != 0.0) {
//     for (size_t j = 0; j < XXindices[i].size(); j++) {
//      r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//     }
//    }
//    b[t][i] = 0.0;
//   }
//   return;
//  }
//
//  double denom = 0.0;
//  for (int k = 0; k < nmodels; k++) {
//   if (std::isfinite(loglik[k])) {
//    denom += std::exp(loglik[k] - max_loglik);
//   }
//  }
//
//  for (int k = 0; k < nmodels; k++) {
//   if (std::isfinite(loglik[k])) {
//    pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
//   } else {
//    pmodel[k] = 0.0;
//   }
//  }
//
//  // Sample model indicator
//  std::uniform_real_distribution<double> runif(0.0, 1.0);
//  double u = runif(gen);
//  double cumprobc = 0.0;
//  int mselect = 0;
//
//  for (int k = 0; k < nmodels; k++) {
//   cumprobc += pmodel[k];
//   if (u < cumprobc) {
//    mselect = k;
//    break;
//   }
//  }
//
//  cmodel[mselect] += 1.0;
//  for (int t = 0; t < nt; t++) {
//   d[t][i] = models[mselect][t];
//  }
//
//  // Sample marker effect from posterior
//  arma::mat C = Bi;
//  for (int t1 = 0; t1 < nt; t1++) {
//   if (models[mselect][t1] == 1) {
//    C(t1, t1) += ww[t1][i] * Ei(t1, t1);
//   }
//  }
//
//  arma::mat Ci = arma::inv(C);
//  arma::mat sample = mvrnormARMA(Ci);  // 1 × nt
//
//  // Compute posterior mean
//  arma::vec rhs_vec(rhs);
//  arma::rowvec mean = (rhs_vec.t() * Ci);  // 1 × nt
//  arma::rowvec mub = sample + mean;
//
//  for (int t = 0; t < nt; t++) {
//   if (models[mselect][t] == 0) {
//    mub(t) = 0.0;
//   }
//  }
//
//  // Update residuals and marker effect
//  for (int t = 0; t < nt; t++) {
//   double diff = mub(t) - b[t][i];
//   if (diff != 0.0) {
//    for (size_t j = 0; j < XXindices[i].size(); j++) {
//     r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
//    }
//   }
//   b[t][i] = mub(t);
//  }
// }

void sampleBetaCMt_eigen(int i,
                   int nt,
                   int nmodels,
                   const std::vector<std::vector<int>>& models,
                   std::vector<double>& cmodel,
                   const std::vector<double>& pi,
                   const arma::mat& Ei,
                   const arma::mat& Bi,
                   const std::vector<std::vector<double>>& ww,
                   std::vector<std::vector<double>>& r,
                   std::vector<std::vector<double>>& b,
                   std::vector<std::vector<int>>& d,
                   const std::vector<std::vector<int>>& XXindices,
                   const std::vector<std::vector<std::vector<double>>>& XXvalues,
                   std::mt19937& gen) {

 std::vector<double> rhs(nt);
 std::vector<double> loglik(nmodels);
 std::vector<double> pmodel(nmodels);
 double rhsi;

 // Compute RHS
 for (int t = 0; t < nt; t++) {
  rhsi =0.0;
  for (size_t j = 0; j < XXvalues[t][i].size(); j++) {
   rhsi +=r[t][j]*XXvalues[t][i][j];
  }
  //rhs[t] = rhsi*Ei(t, t) + b[t][i]*Ei(t, t);
  rhs[t] = rhsi*Ei(t, t) + Ei(t, t) * ww[t][i] * b[t][i];
 }


 // Loop over models to compute log-likelihoods
 for (int k = 0; k < nmodels; k++) {
  arma::mat C = Bi;
  for (int t1 = 0; t1 < nt; t1++) {
   if (models[k][t1] == 1) {
    //C(t1, t1) += Ei(t1, t1);
    C(t1, t1) += ww[t1][i] * Ei(t1, t1);
   }
  }


  // Use log_det for numerical stability
  double logdet, sign;
  arma::log_det(logdet, sign, C);
  arma::mat Ci = arma::inv(C);

  loglik[k] = -0.5 * logdet + std::log(pi[k]);

  for (int t1 = 0; t1 < nt; t1++) {
   for (int t2 = t1; t2 < nt; t2++) {
    if (models[k][t1] == 1 && models[k][t2] == 1) {
     loglik[k] += 0.5 * rhs[t1] * rhs[t2] * (Ci(t1, t2) + (t1 != t2 ? Ci(t2, t1) : 0.0));
    }
   }
  }
 }

 // Compute numerically stable posterior model probabilities
 double max_loglik = *std::max_element(loglik.begin(), loglik.end());
 double denom = 0.0;
 for (int k = 0; k < nmodels; k++) {
  denom += std::exp(loglik[k] - max_loglik);
 }
 for (int k = 0; k < nmodels; k++) {
  pmodel[k] = std::exp(loglik[k] - max_loglik) / denom;
 }

 // Sample model indicator
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 double u = runif(gen);
 double cumprobc = 0.0;
 int mselect = 0;
 for (int k = 0; k < nmodels; k++) {
  cumprobc += pmodel[k];
  if (u < cumprobc) {
   mselect = k;
   break;
  }
 }

 cmodel[mselect] += 1.0;
 for (int t = 0; t < nt; t++) {
  d[t][i] = models[mselect][t];
 }


 // Sample marker effect from posterior
 arma::mat C = Bi;
 for (int t1 = 0; t1 < nt; t1++) {
  if (models[mselect][t1] == 1) {
   //C(t1, t1) += Ei(t1, t1);
   C(t1, t1) += ww[t1][i] * Ei(t1, t1);
  }
 }
 arma::mat Ci = arma::inv(C);
 arma::mat sample = mvrnormARMA(Ci);  // 1 × nt

 // Compute posterior mean
 arma::vec rhs_vec(rhs);
 arma::rowvec mean = (rhs_vec.t() * Ci);  // 1 × nt
 arma::rowvec mub = sample + mean;

 for (int t = 0; t < nt; t++) {
  if (models[mselect][t] == 0) {
   mub(t) = 0.0;
  }
 }

 // Update residuals and marker effect
 for (int t = 0; t < nt; t++) {
  double diff = mub(t) - b[t][i];
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); j++) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = mub(t);
 }
}

// Sample pi across trait
void samplePiMt(int nt,
                std::vector<double>& pimarker,
                const std::vector<int>& dmarker,
                std::mt19937& gen) {

 std::vector<double> mc(2);

 // Sample marker specific pi0
 std::fill(mc.begin(), mc.end(), 0.0);
 for (size_t i = 0; i<dmarker.size() ; i++) {
  mc[dmarker[i]] = mc[dmarker[i]] + 1.0;
 }
 double pisum=0.0;
 for (int j = 0; j<2 ; j++) {
  std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
  double rg = rgamma(gen);
  pimarker[j] = rg/dmarker.size();
  pisum = pisum + pimarker[j];
 }
 for (int j = 0; j<2 ; j++) {
  pimarker[j] = pimarker[j]/pisum;
 }

}

// Sample pi within trait BayesC
void samplePiC(int nt,
               std::vector<std::vector<double>>& pitrait,
               const std::vector<std::vector<int>>& d,
               std::mt19937& gen) {

 std::vector<double> mc(2);

 // Sample trait specific pi
 for (int t = 0; t < nt; t++) {
  std::fill(mc.begin(), mc.end(), 0.0);
  for (size_t i = 0; i<d[t].size() ; i++) {
   mc[d[t][i]] = mc[d[t][i]] + 1.0;
  }
  double pisum=0.0;
  for (int j = 0; j<2 ; j++) {
   std::gamma_distribution<double> rgamma(mc[j]+1.0,1.0);
   double rg = rgamma(gen);
   pitrait[t][j] = rg/d[t].size();
   pisum = pisum + pitrait[t][j];
  }
  for (int j = 0; j<2 ; j++) {
   pitrait[t][j] = pitrait[t][j]/pisum;
  }
 }
}


// Sample marker effects based on within and across pi BayesC
void sampleBetaCSt(int i,
                   int nt,
                   std::vector<int>& dmarker,
                   std::vector<double>& pimarker,
                   std::vector<std::vector<double>>& pitrait,
                   const arma::mat& E,
                   const arma::mat& B,
                   const std::vector<std::vector<double>>& ww,
                   std::vector<std::vector<double>>& r,
                   std::vector<std::vector<double>>& b,
                   std::vector<std::vector<int>>& d,
                   const std::vector<std::vector<int>>& XXindices,
                   const std::vector<std::vector<std::vector<double>>>& XXvalues,
                   std::mt19937& gen) {

 std::vector<double> rhs(nt), loglik0(nt), loglik1(nt), bn(nt);
 double lik0t, lik1t, v0, v1, p0, rhs1, lhs1, diff;
 double u;

 // Compute rhs
 lik0t = 0.0;
 lik1t = 0.0;
 for (int t = 0; t < nt; t++) {
  rhs[t] = r[t][i] + ww[t][i] * b[t][i];
  v0 = ww[t][i]*E(t,t);
  v1 = ww[t][i]*E(t,t) + ww[t][i]*ww[t][i]*B(t,t);
  loglik0[t] = -0.5*std::log(v0) -0.5*((rhs[t]*rhs[t])/v0) + std::log(pitrait[t][0]);
  loglik1[t] = -0.5*std::log(v1) -0.5*((rhs[t]*rhs[t])/v1) + std::log(pitrait[t][1]);
  //lik0t = lik0t + loglik0[t];
  lik0t = lik0t - 0.5*std::log(v0) - 0.5*((rhs[t]*rhs[t])/v0);
  lik1t = lik1t + std::log(std::exp(loglik0[t])+std::exp(loglik1[t]));
 }
 lik0t = lik0t + std::log(pimarker[0]);
 lik1t = lik1t + std::log(pimarker[1]);
 p0 = 1.0/(std::exp(lik1t - lik0t) + 1.0);
 dmarker[i]=0;
 std::uniform_real_distribution<double> runif(0.0, 1.0);
 u = runif(gen);
 if(u>p0) dmarker[i]=1;

 for (int t = 0; t<nt ; t++) {
  d[t][i]=0;
  bn[t]=0.0;
 }
 if(dmarker[i]==1) {
  for (int t = 0; t<nt ; t++) {
   p0 = 1.0/(std::exp(loglik1[t] - loglik0[t]) + 1.0);
   d[t][i]=0;
   std::uniform_real_distribution<double> runif(0.0, 1.0);
   u = runif(gen);
   if(u>p0) d[t][i]=1;
   if(d[t][i]==1) {
    rhs1 = r[t][i] + ww[t][i]*b[t][i];
    lhs1 = ww[t][i] + E(t,t)/B(t,t);
    std::normal_distribution<double> rnorm(rhs1/lhs1, sqrt(E(t,t)/lhs1));
    bn[t] = rnorm(gen);
   }
  }
 }
 // Adjust residuals based on sample marker effects
 for (int t = 0; t < nt; t++) {
  diff = (bn[t] - b[t][i]);
  if (diff != 0.0) {
   for (size_t j = 0; j < XXindices[i].size(); j++) {
    r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j] * diff;
   }
  }
  b[t][i] = bn[t];
 }
}

#include "blr_mt_default_core_impl.h"
#include "blr_mt_bed_core_impl.h"
#include "blr_mt_bed_chains_execution_impl.h"
#include "blr_mt_bed_chains_aggregate_impl.h"
#include "blr_mt_bed_convergence_trace_impl.h"
#include "blr_mt_default_finalize_impl.h"
#include "blr_mt_default_legacy_adapter.h"

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtblr(   std::vector<std::vector<double>> wy,
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

 sblr::mt::MtDefaultDataView data{
  wy, ww, yy, XXvalues, XXindices, n
 };
 sblr::mt::MtDefaultModelSpec model{models, sets, method};
 sblr::mt::MtDefaultCovariancePriorView prior{
  ssb_prior, sse_prior, nub, nue
 };
 sblr::mt::MtDefaultExecutionSpec execution{
  updateB, updateE, updatePi, nit, nburn, nthin, seed
 };
 sblr::mt::MtDefaultInitialState initial_state{
  std::move(b), std::move(B), std::move(E), std::move(pi)
 };
 sblr::mt::MtDefaultCoreResult core_result=sblr::mt::run_mt_default_core(
  data, model, prior, execution, std::move(initial_state)
 );
 sblr::mt::MtDefaultFinalResult final_result=
  sblr::mt::finalize_mt_default_result(std::move(core_result));

 return sblr::mt::make_mt_default_legacy_result(
  final_result, wy, nit, nburn);
}

// Internal canonical trait-specific CSR execution. This native maintenance
// route is deliberately not selected by sblr() and is not namespace-exported.
// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtblr_csr_internal(
 std::vector<std::vector<double>> wy,
 std::vector<std::vector<double>> ww,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 std::vector<std::string> ld_prefixes,
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
 int method
) {
 const std::size_t nt=wy.size();
 if (nt==0 || ww.size()!=nt || yy.size()!=nt || n.size()!=nt) {
  throw std::invalid_argument("mtblr_csr_internal: inconsistent trait dimensions");
 }
 const std::size_t m=wy[0].size();
 if (m==0 || (ld_prefixes.size()!=1 && ld_prefixes.size()!=nt)) {
  throw std::invalid_argument(
   "mtblr_csr_internal: ld_prefixes must have length one or trait count");
 }
 for (std::size_t trait=0; trait<nt; ++trait) {
  if (wy[trait].size()!=m || ww[trait].size()!=m) {
   throw std::invalid_argument("mtblr_csr_internal: inconsistent marker dimensions");
  }
 }
 const bool shared=ld_prefixes.size()==1;
 if (shared) {
  for (std::size_t trait=1; trait<nt; ++trait) {
   if (ww[trait]!=ww[0]) {
    throw std::invalid_argument(
     "mtblr_csr_internal: shared LD requires identical trait diagonals");
   }
  }
 }

 std::vector<sblr::core::SparseLdCsrStorage> storage_owners;
 std::vector<arma::rowvec> diagonal_owners;
 const std::size_t owner_count=shared ? 1 : nt;
 storage_owners.reserve(owner_count);
 diagonal_owners.reserve(owner_count);
 for (std::size_t owner=0; owner<owner_count; ++owner) {
  const std::size_t trait=shared ? 0 : owner;
  storage_owners.push_back(read_and_build_st_ld_csr(
   ld_prefixes[shared ? 0 : trait], static_cast<int>(m), ww[trait]));
  arma::rowvec diagonal(m);
  for (std::size_t marker=0; marker<m; ++marker) diagonal(marker)=ww[trait][marker];
  diagonal_owners.push_back(std::move(diagonal));
 }

 std::vector<sblr::core::SparseLdCsrView> trait_views;
 trait_views.reserve(nt);
 for (std::size_t trait=0; trait<nt; ++trait) {
  const std::size_t owner=shared ? 0 : trait;
  const auto& storage=storage_owners[owner];
  sblr::core::SparseLdCsrView view;
  view.marker_count=m;
  view.row_ptr=storage.ptr.data();
  view.row_ptr_size=storage.ptr.size();
  view.column_index=storage.idx.empty() ? nullptr : storage.idx.data();
  view.offdiag_xij=storage.xij.empty() ? nullptr : storage.xij.data();
  view.nonzero_count=storage.idx.size();
  view.diagonal=&diagonal_owners[owner];
  trait_views.push_back(view);
 }

 sblr::mt::MtCsrDataView data{
  wy, yy, n, sblr::mt::MtSparseLdBundleView{m, std::move(trait_views)}
 };
 sblr::mt::MtDefaultModelSpec model{models, sets, method};
 sblr::mt::MtDefaultCovariancePriorView prior{
  ssb_prior, sse_prior, nub, nue
 };
 sblr::mt::MtDefaultExecutionSpec execution{
  updateB, updateE, updatePi, nit, nburn, nthin, seed
 };
 sblr::mt::MtDefaultInitialState initial_state{
  std::move(b), std::move(B), std::move(E), std::move(pi)
 };
 auto core_result=sblr::mt::run_mt_csr_core(
  data, model, prior, execution, std::move(initial_state));
 auto final_result=sblr::mt::finalize_mt_default_result(std::move(core_result));
 return sblr::mt::make_mt_default_legacy_result(
 final_result, wy, nit, nburn);
}

namespace {

struct MtBlockEigenDescriptor {
 std::vector<std::string> bed_files;
 int n_bed=0;
 std::vector<std::vector<int>> cls;
 std::vector<int> rows0;
 std::vector<double> af;
 std::vector<int> block_start;
 EigenFilterMode filter=EigenFilterMode::hard_truncate;
 double tau=0.01;
 double eta=0.0;
};

MtBlockEigenDescriptor parse_mt_block_eigen_descriptor(
 const Rcpp::List& descriptor, std::size_t marker_count) {
 const char* required[]={"bed_files", "n_bed", "cls", "af", "block_start",
                         "eigen_filter", "eigen_tau", "eigen_eta"};
 for (const char* field : required) {
  if (!descriptor.containsElementNamed(field))
   throw std::invalid_argument(std::string("mt block-eigen descriptor lacks '")+field+"'");
 }
 MtBlockEigenDescriptor out;
 out.bed_files=Rcpp::as<std::vector<std::string>>(descriptor["bed_files"]);
 out.n_bed=Rcpp::as<int>(descriptor["n_bed"]);
 out.cls=Rcpp::as<std::vector<std::vector<int>>>(descriptor["cls"]);
 out.af=Rcpp::as<std::vector<double>>(descriptor["af"]);
 out.block_start=Rcpp::as<std::vector<int>>(descriptor["block_start"]);
 out.filter=parse_block_eigen_filter_mode(Rcpp::as<std::string>(descriptor["eigen_filter"]));
 out.tau=Rcpp::as<double>(descriptor["eigen_tau"]);
 out.eta=Rcpp::as<double>(descriptor["eigen_eta"]);
 if (out.bed_files.empty() || out.n_bed<=0 || out.cls.size()!=out.bed_files.size())
  throw std::invalid_argument("mt block-eigen descriptor BED dimensions are invalid");
 std::size_t selected=0;
 for (const auto& file_cls : out.cls) {
  for (int marker : file_cls) {
   if (marker<=0) throw std::invalid_argument("mt block-eigen cls values must be positive and 1-based");
  }
  selected += file_cls.size();
 }
 if (selected!=marker_count || out.af.size()!=marker_count)
  throw std::invalid_argument("mt block-eigen descriptor marker count does not match wy");
 for (double frequency : out.af) {
  if (!std::isfinite(frequency) || frequency<=0.0 || frequency>=1.0)
   throw std::invalid_argument("mt block-eigen allele frequencies must be finite and in (0, 1)");
 }
 if (!std::isfinite(out.tau) || out.tau<0.0 ||
     !std::isfinite(out.eta) || out.eta<0.0)
  throw std::invalid_argument("mt block-eigen tau and eta must be finite and nonnegative");
 if (out.block_start.empty() || out.block_start[0]!=0)
  throw std::invalid_argument("mt block-eigen block_start must begin at zero");
 for (std::size_t i=0; i<out.block_start.size(); ++i) {
  if (out.block_start[i]<0 || static_cast<std::size_t>(out.block_start[i])>=marker_count ||
      (i>0 && out.block_start[i]<=out.block_start[i-1]))
   throw std::invalid_argument("mt block-eigen block_start must be strictly ascending and in range");
 }
 if (descriptor.containsElementNamed("rows") && !Rf_isNull(descriptor["rows"])) {
  out.rows0=Rcpp::as<std::vector<int>>(descriptor["rows"]);
  if (out.rows0.empty()) throw std::invalid_argument("mt block-eigen rows must be nonempty when supplied");
  for (int& row : out.rows0) {
   if (row<=0 || row>out.n_bed)
    throw std::invalid_argument("mt block-eigen rows must be positive, 1-based, and within n_bed");
   --row;
  }
 }
 return out;
}

}  // namespace

namespace {

struct MtBlockEigenAdapterResult {
 sblr::mt::MtDefaultLegacyResult legacy;
 std::vector<std::vector<BlockEigenDiag>> owner_diagnostics;
 std::vector<int> trait_owner;
 std::vector<std::vector<double>> transformed_wy;
 std::vector<MtBlockEigenDescriptor> descriptors;
};

MtBlockEigenAdapterResult run_mt_block_eigen_adapter(
 std::vector<std::vector<double>> wy,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 Rcpp::List operator_descriptors,
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
 int method
) {
 if (method!=4) throw std::invalid_argument("mtblr_block_eigen_internal supports method = 4 only");
 const std::size_t nt=wy.size();
 if (nt==0 || yy.size()!=nt || n.size()!=nt || b.size()!=nt)
  throw std::invalid_argument("mtblr_block_eigen_internal: inconsistent trait dimensions");
 const std::size_t m=wy[0].size();
 if (m==0 || (operator_descriptors.size()!=1 &&
              static_cast<std::size_t>(operator_descriptors.size())!=nt))
  throw std::invalid_argument("mtblr_block_eigen_internal: operator_descriptors must have length one or trait count");
 for (std::size_t trait=0; trait<nt; ++trait) {
  if (wy[trait].size()!=m || b[trait].size()!=m)
   throw std::invalid_argument("mtblr_block_eigen_internal: inconsistent marker dimensions");
 }

 const bool shared=operator_descriptors.size()==1;
 const std::size_t owner_count=shared ? 1 : nt;
 std::vector<MtBlockEigenDescriptor> descriptors;
 descriptors.reserve(owner_count);
 for (std::size_t owner=0; owner<owner_count; ++owner)
  descriptors.push_back(parse_mt_block_eigen_descriptor(
   Rcpp::as<Rcpp::List>(operator_descriptors[static_cast<R_xlen_t>(owner)]), m));

 std::vector<sblr::core::BlockEigenStorage> storage_owners;
 storage_owners.reserve(owner_count);
 std::vector<std::vector<BlockEigenDiag>> owner_diagnostics(owner_count);
 std::vector<std::vector<double>> transformed_wy=wy;
 for (std::size_t owner=0; owner<owner_count; ++owner) {
  const auto& descriptor=descriptors[owner];
  const int* rows=descriptor.rows0.empty() ? nullptr : descriptor.rows0.data();
  PackedBedMatrix packed=read_bedfiles_to_packed_matrix(
   descriptor.bed_files, descriptor.n_bed, rows,
   static_cast<int>(descriptor.rows0.size()), descriptor.cls);
  const std::size_t rows_in_build=shared ? nt : 1;
  arma::mat wy_matrix(rows_in_build, m);
  for (std::size_t row=0; row<rows_in_build; ++row) {
   const std::size_t trait=shared ? row : owner;
   for (std::size_t marker=0; marker<m; ++marker)
    wy_matrix(static_cast<arma::uword>(row), static_cast<arma::uword>(marker))=wy[trait][marker];
  }
  storage_owners.push_back(build_block_eigen(
   packed, descriptor.af, descriptor.block_start, descriptor.filter,
   descriptor.tau, descriptor.eta, wy_matrix, 1,
   &owner_diagnostics[owner]));
  for (std::size_t row=0; row<rows_in_build; ++row) {
   const std::size_t trait=shared ? row : owner;
   for (std::size_t marker=0; marker<m; ++marker)
    transformed_wy[trait][marker]=wy_matrix(
     static_cast<arma::uword>(row), static_cast<arma::uword>(marker));
  }
 }

 std::vector<sblr::core::BlockEigenView> trait_views;
 trait_views.reserve(nt);
 for (std::size_t trait=0; trait<nt; ++trait)
  trait_views.push_back(storage_owners[shared ? 0 : trait].view());

 sblr::mt::MtBlockEigenDataView data{
  transformed_wy, yy, n,
  sblr::mt::MtBlockEigenBundleView{m, std::move(trait_views)}
 };
 sblr::mt::MtDefaultModelSpec model{models, sets, method};
 sblr::mt::MtDefaultCovariancePriorView prior{ssb_prior, sse_prior, nub, nue};
 sblr::mt::MtDefaultExecutionSpec execution{
  updateB, updateE, updatePi, nit, nburn, nthin, seed};
 sblr::mt::MtDefaultInitialState initial_state{
  std::move(b), std::move(B), std::move(E), std::move(pi)};
 auto core_result=sblr::mt::run_mt_block_eigen_core(
  data, model, prior, execution, std::move(initial_state));
 auto final_result=sblr::mt::finalize_mt_default_result(std::move(core_result));
 MtBlockEigenAdapterResult result;
 result.legacy=sblr::mt::make_mt_default_legacy_result(
  final_result, transformed_wy, nit, nburn);
 result.owner_diagnostics=std::move(owner_diagnostics);
 result.trait_owner.resize(nt);
 for (std::size_t trait=0; trait<nt; ++trait)
  result.trait_owner[trait]=static_cast<int>(shared ? 0 : trait);
 result.transformed_wy=std::move(transformed_wy);
 result.descriptors=std::move(descriptors);
 return result;
}

}  // namespace

// Internal canonical trait-specific block-eigen execution. Operator building
// and any matching wy projection complete before the shared MT core creates
// its fit-local RNG. This maintenance route is not namespace-exported.
// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>> mtblr_block_eigen_internal(
 std::vector<std::vector<double>> wy,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 Rcpp::List operator_descriptors,
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
 int method=4
) {
 return run_mt_block_eigen_adapter(
  std::move(wy), std::move(yy), std::move(b), operator_descriptors, sets,
  std::move(B), std::move(E), std::move(ssb_prior), std::move(sse_prior),
  std::move(models), std::move(pi), nub, nue, updateB, updateE, updatePi,
  std::move(n), nit, nburn, nthin, seed, method).legacy;
}

namespace {

Rcpp::NumericMatrix mt_legacy_trait_matrix(
 const std::vector<std::vector<std::vector<double>>>& value,
 int field
) {
 const int nt=static_cast<int>(value[field].size());
 const int nr=nt==0 ? 0 : static_cast<int>(value[field][0].size());
 Rcpp::NumericMatrix out(nr, nt);
 for (int trait=0; trait<nt; ++trait)
  for (int row=0; row<nr; ++row) out(row, trait)=value[field][trait][row];
 return out;
}

Rcpp::IntegerMatrix mt_legacy_state_matrix(
 const std::vector<std::vector<std::vector<double>>>& value,
 int field
) {
 const int nt=static_cast<int>(value[field].size());
 const int nr=nt==0 ? 0 : static_cast<int>(value[field][0].size());
 Rcpp::IntegerMatrix out(nr, nt);
 for (int trait=0; trait<nt; ++trait)
  for (int row=0; row<nr; ++row)
   out(row, trait)=static_cast<int>(value[field][trait][row]);
 return out;
}

Rcpp::NumericVector mt_legacy_vector(
 const std::vector<std::vector<std::vector<double>>>& value,
 int field
) {
 if (value[field].empty()) return Rcpp::NumericVector();
 return Rcpp::NumericVector(Rcpp::wrap(value[field][0]));
}

Rcpp::List mtblr_legacy_to_raw(
 const sblr::mt::MtDefaultLegacyResult& legacy,
 const std::vector<std::vector<int>>& models,
 const std::string& backend,
 const std::string& data_level,
 int nit,
 int nburn,
 int nthin,
 bool updateB,
 bool updateE,
 bool updatePi,
 Rcpp::Nullable<Rcpp::List> extra_diagnostics=R_NilValue
) {
 const int nt=legacy.empty() ? 0 : static_cast<int>(legacy[0].size());
 const int m=(nt==0 || legacy[0][0].empty()) ? 0 :
  static_cast<int>(legacy[0][0].size());
 const int nmodels=static_cast<int>(models.size());
 Rcpp::IntegerMatrix patterns(nmodels, nt);
 for (int model=0; model<nmodels; ++model)
  for (int trait=0; trait<nt; ++trait) patterns(model, trait)=models[model][trait];
 Rcpp::IntegerVector order(m);
 for (int marker=0; marker<m; ++marker)
  order[marker]=static_cast<int>(legacy[6][0][marker]);
 const int marker_count=(nit+nthin-1)/nthin;
 Rcpp::List diagnostics=Rcpp::List::create(
  Rcpp::_["marker"]=marker_count,
  Rcpp::_["covb"]=updateB ? nit : 0,
  Rcpp::_["covg"]=nit,
  Rcpp::_["cove"]=updateE ? nit : 0,
  Rcpp::_["pi"]=updatePi ? nit : 0);
 if (extra_diagnostics.isNotNull()) {
  Rcpp::List extra(extra_diagnostics);
  Rcpp::CharacterVector extra_names=extra.names();
  for (R_xlen_t i=0; i<extra.size(); ++i)
   diagnostics[static_cast<std::string>(extra_names[i])]=extra[i];
 }
 return Rcpp::List::create(
  Rcpp::_["schema"]=Rcpp::List::create(
   Rcpp::_["class"]="mtblr_raw", Rcpp::_["version"]=1),
  Rcpp::_["meta"]=Rcpp::List::create(
   Rcpp::_["model"]="bayesc", Rcpp::_["backend"]=backend,
   Rcpp::_["data_level"]=data_level, Rcpp::_["m"]=m,
   Rcpp::_["nt"]=nt, Rcpp::_["n_trace"]=nit+nburn,
   Rcpp::_["nit"]=nit, Rcpp::_["nburn"]=nburn,
   Rcpp::_["nthin"]=nthin, Rcpp::_["nmodels"]=nmodels),
  Rcpp::_["marker"]=Rcpp::List::create(
   Rcpp::_["bm"]=mt_legacy_trait_matrix(legacy,0),
   Rcpp::_["dm"]=mt_legacy_trait_matrix(legacy,1),
   Rcpp::_["wy"]=mt_legacy_trait_matrix(legacy,2),
   Rcpp::_["r"]=mt_legacy_trait_matrix(legacy,3),
   Rcpp::_["b"]=mt_legacy_trait_matrix(legacy,4),
   Rcpp::_["state"]=mt_legacy_state_matrix(legacy,5),
   Rcpp::_["order"]=order),
  Rcpp::_["trace"]=Rcpp::List::create(
   Rcpp::_["vbs"]=mt_legacy_trait_matrix(legacy,7),
   Rcpp::_["vgs"]=mt_legacy_trait_matrix(legacy,8),
   Rcpp::_["ves"]=mt_legacy_trait_matrix(legacy,9),
   Rcpp::_["vle"]=mt_legacy_trait_matrix(legacy,20),
   Rcpp::_["vld"]=mt_legacy_trait_matrix(legacy,21)),
  Rcpp::_["variance"]=Rcpp::List::create(
   Rcpp::_["covb"]=mt_legacy_trait_matrix(legacy,10),
   Rcpp::_["covg"]=mt_legacy_trait_matrix(legacy,11),
   Rcpp::_["cove"]=mt_legacy_trait_matrix(legacy,12),
   Rcpp::_["vb"]=mt_legacy_trait_matrix(legacy,13),
   Rcpp::_["vg"]=mt_legacy_trait_matrix(legacy,14),
   Rcpp::_["ve"]=mt_legacy_trait_matrix(legacy,15)),
  Rcpp::_["pi"]=Rcpp::List::create(
   Rcpp::_["final"]=mt_legacy_vector(legacy,16),
   Rcpp::_["mean"]=mt_legacy_vector(legacy,17)),
  Rcpp::_["model"]=Rcpp::List::create(Rcpp::_["patterns"]=patterns),
  Rcpp::_["diagnostics"]=diagnostics,
  Rcpp::_["data"]=Rcpp::List::create(),
  Rcpp::_["alignment"]=Rcpp::List::create());
}

std::string mt_block_eigen_filter_name(EigenFilterMode mode) {
 if (mode==EigenFilterMode::hard_truncate) return "hard_truncate";
 if (mode==EigenFilterMode::ridge_fixed) return "ridge_fixed";
 return "ridge_lw";
}

}  // namespace

// Named schema adapter for the public R mtblr_csr() boundary. Numerical
// execution remains exclusively owned by mtblr_csr_internal() and its shared
// Phase 17I core; this function only names and shapes finalized values.
// [[Rcpp::export]]
Rcpp::List mtblr_csr_raw_internal(
 std::vector<std::vector<double>> wy,
 std::vector<std::vector<double>> ww,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 std::vector<std::string> ld_prefixes,
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
 int method
) {
 auto legacy=mtblr_csr_internal(
  std::move(wy), std::move(ww), std::move(yy), std::move(b),
  std::move(ld_prefixes), sets, std::move(B), std::move(E),
  std::move(ssb_prior), std::move(sse_prior), models, std::move(pi),
  nub, nue, updateB, updateE, updatePi, std::move(n), nit, nburn,
  nthin, seed, method);
 return mtblr_legacy_to_raw(
  legacy, models, "mt_csr_bayesc", "summary", nit, nburn, nthin,
 updateB, updateE, updatePi);
}

namespace {

std::vector<int> mt_summary_resolve_chain_seeds(
 int seed, int nchains, const std::vector<int>& requested) {
 if (nchains<1) throw std::invalid_argument("nchains must be positive");
 if (!requested.empty() && static_cast<int>(requested.size())!=nchains)
  throw std::invalid_argument("chain_seeds must have length nchains");
 if (!requested.empty()) return requested;
 std::vector<int> resolved(static_cast<std::size_t>(nchains));
 const std::uint32_t base=static_cast<std::uint32_t>(seed);
 for (int chain=0;chain<nchains;++chain) {
  const std::uint32_t value=base+static_cast<std::uint32_t>(9176u)*
   static_cast<std::uint32_t>(chain);
  resolved[static_cast<std::size_t>(chain)]=static_cast<std::int32_t>(value);
 }
 return resolved;
}

struct MtSummaryChainResult {
 sblr::mt::MtDefaultLegacyResult legacy;
 std::vector<int> component;
 std::vector<std::vector<double>> component_probabilities;
 std::vector<std::vector<double>> pi_trace;
 arma::mat annotation_alpha_final;
 arma::mat annotation_alpha_mean;
 arma::vec annotation_sigma_final;
 arma::vec annotation_sigma_mean;
 std::vector<double> pattern_pi_final;
 std::vector<double> pattern_pi_mean;
 std::vector<std::vector<double>> pattern_pi_trace;
 arma::mat prior_component_probabilities;
 int annotation_updates_attempted=0;
 int annotation_updates_completed=0;
 sblr::mt::MtExtendedTraceResult convergence;
};

Rcpp::List mt_extended_trace_to_list(
 const sblr::mt::MtExtendedTraceResult& trace) {
 auto numeric=[](const std::vector<std::vector<double>>& x) {
  const std::size_t nq=x.size(), n=nq==0 ? 0 : x[0].size();
  Rcpp::NumericMatrix out(static_cast<int>(n),static_cast<int>(nq));
  for (std::size_t q=0;q<nq;++q) for (std::size_t i=0;i<n;++i)
   out(static_cast<int>(i),static_cast<int>(q))=x[q][i];
  return out;
 };
 auto integer=[](const std::vector<std::vector<int>>& x) {
  const std::size_t nq=x.size(), n=nq==0 ? 0 : x[0].size();
  Rcpp::IntegerMatrix out(static_cast<int>(n),static_cast<int>(nq));
  for (std::size_t q=0;q<nq;++q) for (std::size_t i=0;i<n;++i)
   out(static_cast<int>(i),static_cast<int>(q))=x[q][i];
  return out;
 };
 return Rcpp::List::create(
  Rcpp::_["cov_b"]=numeric(trace.cov_b),
  Rcpp::_["cov_g"]=numeric(trace.cov_g),
  Rcpp::_["cov_e"]=numeric(trace.cov_e),
  Rcpp::_["component_pi"]=numeric(trace.component_pi),
  Rcpp::_["pattern_pi"]=numeric(trace.pattern_pi),
  Rcpp::_["joint_pi"]=numeric(trace.joint_pi),
  Rcpp::_["alpha"]=numeric(trace.annotation_alpha),
  Rcpp::_["sigmaSqAlpha"]=numeric(trace.annotation_sigma),
  Rcpp::_["b"]=numeric(trace.selected_b),
  Rcpp::_["d"]=integer(trace.selected_d),
  Rcpp::_["component"]=integer(trace.selected_component));
}

template <class Runner>
std::vector<MtSummaryChainResult> mt_summary_dispatch_chains(
 int nchains, int ncores, const Runner& runner,
 std::vector<std::string>& errors, int& used_workers) {
 if (ncores<1) throw std::invalid_argument("ncores must be positive");
 std::vector<MtSummaryChainResult> results(
  static_cast<std::size_t>(nchains));
 std::vector<int> failed(static_cast<std::size_t>(nchains),0);
 errors.assign(static_cast<std::size_t>(nchains),std::string());
 used_workers=1;
#ifdef _OPENMP
 used_workers=std::min(ncores,nchains);
#pragma omp parallel for num_threads(used_workers) schedule(static)
#endif
 for (int chain=0;chain<nchains;++chain) {
  try {
   results[static_cast<std::size_t>(chain)]=runner(chain);
  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(chain)]=1;
   errors[static_cast<std::size_t>(chain)]=e.what();
  } catch (...) {
   failed[static_cast<std::size_t>(chain)]=1;
   errors[static_cast<std::size_t>(chain)]="unknown error";
  }
 }
 for (int chain=0;chain<nchains;++chain) if (failed[static_cast<std::size_t>(chain)])
  throw std::runtime_error("MT summary chain "+std::to_string(chain+1)+
   " failed: "+errors[static_cast<std::size_t>(chain)]);
 return results;
}

Rcpp::List mt_summary_raw_list(
 const std::vector<MtSummaryChainResult>& chains,
 const std::vector<std::vector<int>>& models, const std::string& backend,
 int nit, int nburn, int nthin, bool updateB, bool updateE, bool updatePi,
 const std::vector<int>& seeds, int requested_workers, int used_workers,
 Rcpp::Nullable<Rcpp::List> extra=R_NilValue) {
 Rcpp::List raws(static_cast<R_xlen_t>(chains.size()));
 for (std::size_t chain=0;chain<chains.size();++chain)
  {
   Rcpp::List raw=mtblr_legacy_to_raw(
    chains[chain].legacy,models,backend,"summary",nit,nburn,nthin,
    updateB,updateE,updatePi,extra);
   if (!chains[chain].component.empty()) {
    Rcpp::List marker=raw["marker"];
    marker["component_final"]=Rcpp::wrap(chains[chain].component);
    const std::size_t m=chains[chain].component_probabilities.size();
    const std::size_t k=m==0 ? 0 : chains[chain].component_probabilities[0].size();
    Rcpp::NumericMatrix probabilities(static_cast<int>(m),static_cast<int>(k));
    for (std::size_t row=0;row<m;++row)
     for (std::size_t col=0;col<k;++col)
      probabilities(static_cast<int>(row),static_cast<int>(col))=
       chains[chain].component_probabilities[row][col];
    marker["component_probabilities"]=probabilities;
    raw["marker"]=marker;
    Rcpp::List pi_namespace=raw["pi"];
    const std::size_t states=chains[chain].pi_trace.size();
    const std::size_t iterations=states==0 ? 0 : chains[chain].pi_trace[0].size();
    Rcpp::NumericMatrix trace(static_cast<int>(iterations),static_cast<int>(states));
    for (std::size_t state=0;state<states;++state)
     for (std::size_t iteration=0;iteration<iterations;++iteration)
      trace(static_cast<int>(iteration),static_cast<int>(state))=
       chains[chain].pi_trace[state][iteration];
    pi_namespace["trace"]=trace; raw["pi"]=pi_namespace;
   }
   if (chains[chain].annotation_alpha_final.n_elem>0) {
    const std::size_t patterns=chains[chain].pattern_pi_trace.size();
    const std::size_t iterations=patterns==0 ? 0 :
     chains[chain].pattern_pi_trace[0].size();
    Rcpp::NumericMatrix pattern_trace(
     static_cast<int>(iterations),static_cast<int>(patterns));
    for (std::size_t p=0;p<patterns;++p)
     for (std::size_t iteration=0;iteration<iterations;++iteration)
      pattern_trace(static_cast<int>(iteration),static_cast<int>(p))=
       chains[chain].pattern_pi_trace[p][iteration];
    raw["annotations"]=Rcpp::List::create(
     Rcpp::_ ["annotation_coefficients_final"]=chains[chain].annotation_alpha_final,
     Rcpp::_ ["annotation_coefficients_mean"]=chains[chain].annotation_alpha_mean,
     Rcpp::_ ["annotation_variances_final"]=chains[chain].annotation_sigma_final,
     Rcpp::_ ["annotation_variances_mean"]=chains[chain].annotation_sigma_mean,
     Rcpp::_ ["pattern_pi_final"]=chains[chain].pattern_pi_final,
     Rcpp::_ ["pattern_pi_mean"]=chains[chain].pattern_pi_mean,
     Rcpp::_ ["pattern_pi_trace"]=pattern_trace,
     Rcpp::_ ["prior_component_probabilities"]=chains[chain].prior_component_probabilities,
     Rcpp::_ ["annotation_updates_attempted"]=chains[chain].annotation_updates_attempted,
     Rcpp::_ ["annotation_updates_completed"]=chains[chain].annotation_updates_completed);
   }
   raw["convergence_capture"]=mt_extended_trace_to_list(
    chains[chain].convergence);
   raws[static_cast<R_xlen_t>(chain)]=raw;
  }
 return Rcpp::List::create(
  Rcpp::_ ["raws"]=raws,
  Rcpp::_ ["chain_seeds"]=Rcpp::wrap(seeds),
  Rcpp::_ ["requested_workers"]=requested_workers,
  Rcpp::_ ["used_workers"]=used_workers,
  Rcpp::_ ["operator_preparations"]=1);
}

}  // namespace

// Shared-preparation deterministic logical-chain execution for MT CSR.
// [[Rcpp::export]]
Rcpp::List mtblr_csr_chains_raw_internal(
 std::vector<std::vector<double>> wy,
 std::vector<std::vector<double>> ww,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 std::vector<std::string> ld_prefixes,
 const std::vector<std::vector<int>>& sets,
 arma::mat B, arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<std::vector<int>> models,
 std::vector<double> pi, double nub, double nue,
 bool updateB, bool updateE, bool updatePi,
 std::vector<int> n, int nit, int nburn, int nthin, int seed, int method,
  int nchains, int ncores, std::vector<int> chain_seeds,
  std::vector<int> joint_component,
  std::vector<double> joint_multiplier,
  std::vector<std::string> joint_names,
  int component_count,
  std::vector<double> marker_scale,
  std::vector<int> component_init,
  std::vector<double> pi_prior,
  std::vector<std::vector<double>> beta_init,
  std::vector<std::vector<int>> state_init,
  arma::mat annotations, arma::mat alpha_init,
  std::vector<double> sigma_alpha_init,
  std::vector<double> pattern_pi_init,
  std::vector<double> pattern_pi_prior,
  bool updateAlpha, bool intercept_flat,
  double sigma_alpha_a, double sigma_alpha_b,
  double pi_floor, int alpha_update_every,
  bool convergence_covariance=false, bool convergence_probability=false,
  bool convergence_annotations=false,
  bool convergence_full_probability=false,
  SEXP convergence_markers=R_NilValue,
  bool convergence_b=false, bool convergence_d=false,
  bool convergence_component=false) {
 const std::size_t nt=wy.size();
 if (nt==0 || ww.size()!=nt || yy.size()!=nt || n.size()!=nt)
  throw std::invalid_argument("mtblr_csr_chains_raw_internal: inconsistent trait dimensions");
 const std::size_t m=wy[0].size();
 if (m==0 || (ld_prefixes.size()!=1 && ld_prefixes.size()!=nt))
  throw std::invalid_argument("mtblr_csr_chains_raw_internal: invalid LD ownership");
 for (std::size_t trait=0;trait<nt;++trait)
  if (wy[trait].size()!=m || ww[trait].size()!=m)
   throw std::invalid_argument("mtblr_csr_chains_raw_internal: inconsistent marker dimensions");
 const bool shared=ld_prefixes.size()==1;
 if (shared) for (std::size_t trait=1;trait<nt;++trait)
  if (ww[trait]!=ww[0]) throw std::invalid_argument(
   "mtblr_csr_chains_raw_internal: shared LD requires identical diagonals");

 std::vector<sblr::core::SparseLdCsrStorage> storage_owners;
 std::vector<arma::rowvec> diagonal_owners;
 const std::size_t owner_count=shared ? 1 : nt;
 storage_owners.reserve(owner_count); diagonal_owners.reserve(owner_count);
 for (std::size_t owner=0;owner<owner_count;++owner) {
  const std::size_t trait=shared ? 0 : owner;
  storage_owners.push_back(read_and_build_st_ld_csr(
   ld_prefixes[shared ? 0 : trait],static_cast<int>(m),ww[trait]));
  arma::rowvec diagonal(m);
  for (std::size_t marker=0;marker<m;++marker) diagonal(marker)=ww[trait][marker];
  diagonal_owners.push_back(std::move(diagonal));
 }
 std::vector<sblr::core::SparseLdCsrView> views;
 views.reserve(nt);
 for (std::size_t trait=0;trait<nt;++trait) {
  const std::size_t owner=shared ? 0 : trait;
  const auto& storage=storage_owners[owner];
  sblr::core::SparseLdCsrView view;
  view.marker_count=m; view.row_ptr=storage.ptr.data();
  view.row_ptr_size=storage.ptr.size();
  view.column_index=storage.idx.empty()?nullptr:storage.idx.data();
  view.offdiag_xij=storage.xij.empty()?nullptr:storage.xij.data();
  view.nonzero_count=storage.idx.size(); view.diagonal=&diagonal_owners[owner];
  views.push_back(view);
 }
 const sblr::mt::MtCsrDataView data{
  wy,yy,n,sblr::mt::MtSparseLdBundleView{m,std::move(views)}};
 sblr::mt::MtJointStateSpec joint;
 if (method==5 || method==6) {
  joint.patterns=models; joint.component=std::move(joint_component);
  joint.multiplier=std::move(joint_multiplier); joint.names=std::move(joint_names);
  joint.component_count=component_count; joint.scaled=true;
 }
 sblr::mt::MtBayesRCSpec bayesrc;
 if (method==6) {
  bayesrc.annotations=&annotations;
  bayesrc.pattern_prior=pattern_pi_prior;
  bayesrc.update_alpha=updateAlpha;
  bayesrc.intercept_flat=intercept_flat;
  bayesrc.sigma_alpha_a=sigma_alpha_a;
  bayesrc.sigma_alpha_b=sigma_alpha_b;
  bayesrc.pi_floor=pi_floor;
  bayesrc.alpha_update_every=alpha_update_every;
 }
 sblr::mt::MtExtendedTraceSpec convergence_spec;
 convergence_spec.covariance=convergence_covariance;
 convergence_spec.probability=convergence_probability;
 convergence_spec.annotations=convergence_annotations;
 convergence_spec.full_probability_states=convergence_full_probability;
 convergence_spec.selected_markers=Rf_isNull(convergence_markers) ?
  std::vector<int>() : Rcpp::as<std::vector<int>>(convergence_markers);
 convergence_spec.selected_b=convergence_b;
 convergence_spec.selected_d=convergence_d;
 convergence_spec.selected_component=convergence_component;
 const bool any_extended=convergence_spec.covariance ||
  convergence_spec.probability || convergence_spec.annotations ||
  !convergence_spec.selected_markers.empty();
 const sblr::mt::MtDefaultModelSpec model{
  models,sets,method,(method==5||method==6)?&joint:nullptr,
  (method==5||method==6)?&marker_scale:nullptr,
  method==5?&pi_prior:nullptr,method==6?&bayesrc:nullptr,
  any_extended?&convergence_spec:nullptr};
 const sblr::mt::MtDefaultCovariancePriorView prior{
  ssb_prior,sse_prior,nub,nue};
 const std::vector<int> seeds=mt_summary_resolve_chain_seeds(
  seed,nchains,chain_seeds);
 std::vector<std::string> errors; int used_workers=1;
 auto runner=[&](int chain) {
  const sblr::mt::MtDefaultExecutionSpec execution{
   updateB,updateE,updatePi,nit,nburn,nthin,
   seeds[static_cast<std::size_t>(chain)]};
  sblr::mt::MtDefaultInitialState initial{
   b,B,E,pi,component_init,beta_init,state_init,alpha_init,
   arma::vec(sigma_alpha_init),pattern_pi_init};
  auto core=sblr::mt::run_mt_csr_core(data,model,prior,execution,std::move(initial));
  auto final=sblr::mt::finalize_mt_default_result(std::move(core));
  MtSummaryChainResult out;
  out.component=final.component;
  out.component_probabilities=final.component_probabilities;
  out.pi_trace=final.pi_trace;
  out.annotation_alpha_final=final.annotation_alpha_final;
  out.annotation_alpha_mean=final.annotation_alpha_mean;
  out.annotation_sigma_final=final.annotation_sigma_final;
  out.annotation_sigma_mean=final.annotation_sigma_mean;
  out.pattern_pi_final=final.pattern_pi_final;
  out.pattern_pi_mean=final.pattern_pi_mean;
  out.pattern_pi_trace=final.pattern_pi_trace;
  out.prior_component_probabilities=final.prior_component_probabilities;
  out.annotation_updates_attempted=final.annotation_updates_attempted;
  out.annotation_updates_completed=final.annotation_updates_completed;
  out.convergence=std::move(final.convergence);
  out.legacy=sblr::mt::make_mt_default_legacy_result(final,wy,nit,nburn);
  return out;
 };
 const auto results=mt_summary_dispatch_chains(
  nchains,ncores,runner,errors,used_workers);
  return mt_summary_raw_list(results,models,
   method==6?"mt_csr_bayesrc":method==5?"mt_csr_bayesr":"mt_csr_bayesc",
   nit,nburn,nthin,
  updateB,updateE,updatePi,seeds,ncores,used_workers);
}

// Named schema adapter for the public R mtblr_block_eigen() boundary. The
// adapter builds operators and executes the Phase 17L core exactly once.
// [[Rcpp::export]]
Rcpp::List mtblr_block_eigen_raw_internal(
 std::vector<std::vector<double>> wy,
 std::vector<double> yy,
 std::vector<std::vector<double>> b,
 Rcpp::List operator_descriptors,
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
 int method=4
) {
 MtBlockEigenAdapterResult adapter=run_mt_block_eigen_adapter(
  std::move(wy), std::move(yy), std::move(b), operator_descriptors, sets,
  std::move(B), std::move(E), std::move(ssb_prior), std::move(sse_prior),
  models, std::move(pi), nub, nue, updateB, updateE, updatePi,
  std::move(n), nit, nburn, nthin, seed, method);
 const int owner_count=static_cast<int>(adapter.descriptors.size());
 Rcpp::List owners(owner_count);
 for (int owner=0; owner<owner_count; ++owner) {
  const auto& descriptor=adapter.descriptors[static_cast<std::size_t>(owner)];
  owners[owner]=Rcpp::List::create(
   Rcpp::_["eigen_filter"]=mt_block_eigen_filter_name(descriptor.filter),
   Rcpp::_["eigen_tau"]=descriptor.tau,
   Rcpp::_["eigen_eta"]=descriptor.eta,
   Rcpp::_["mu_floor"]=0.01,
   Rcpp::_["reference_sample_size"]=descriptor.n_bed,
   Rcpp::_["selected_row_count"]=descriptor.rows0.empty() ?
    descriptor.n_bed : static_cast<int>(descriptor.rows0.size()),
   Rcpp::_["marker_count"]=static_cast<int>(descriptor.af.size()),
   Rcpp::_["blocks"]=block_eigen_diagnostics_to_data_frame(
    adapter.owner_diagnostics[static_cast<std::size_t>(owner)]));
 }
 Rcpp::IntegerVector trait_owner(adapter.trait_owner.size());
 for (std::size_t trait=0; trait<adapter.trait_owner.size(); ++trait)
  trait_owner[static_cast<R_xlen_t>(trait)]=adapter.trait_owner[trait]+1;
 Rcpp::List block_diagnostics=Rcpp::List::create(
  Rcpp::_["owner_count"]=owner_count,
  Rcpp::_["trait_owner"]=trait_owner,
  Rcpp::_["sharing_mode"]=owner_count==1 ?
   "fully_shared_operator" : "trait_specific_operator",
  Rcpp::_["owners"]=owners);
 return mtblr_legacy_to_raw(
  adapter.legacy, models, "mt_block_eigen_bayesc", "summary",
  nit, nburn, nthin,
  updateB, updateE, updatePi,
  Rcpp::List::create(Rcpp::_["block_eigen"]=block_diagnostics));
}

// Shared-preparation deterministic logical-chain execution for MT block eigen.
// [[Rcpp::export]]
Rcpp::List mtblr_block_eigen_chains_raw_internal(
 std::vector<std::vector<double>> wy, std::vector<double> yy,
 std::vector<std::vector<double>> b, Rcpp::List operator_descriptors,
 const std::vector<std::vector<int>>& sets, arma::mat B, arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<std::vector<int>> models, std::vector<double> pi,
 double nub, double nue, bool updateB, bool updateE, bool updatePi,
 std::vector<int> n, int nit, int nburn, int nthin, int seed, int method,
  int nchains, int ncores, std::vector<int> chain_seeds,
  std::vector<int> joint_component,
  std::vector<double> joint_multiplier,
  std::vector<std::string> joint_names,
  int component_count,
  std::vector<double> marker_scale,
  std::vector<int> component_init,
  std::vector<double> pi_prior,
  std::vector<std::vector<double>> beta_init,
  std::vector<std::vector<int>> state_init,
  arma::mat annotations, arma::mat alpha_init,
  std::vector<double> sigma_alpha_init,
  std::vector<double> pattern_pi_init,
  std::vector<double> pattern_pi_prior,
  bool updateAlpha, bool intercept_flat,
  double sigma_alpha_a, double sigma_alpha_b,
  double pi_floor, int alpha_update_every,
  bool convergence_covariance=false, bool convergence_probability=false,
  bool convergence_annotations=false,
  bool convergence_full_probability=false,
  SEXP convergence_markers=R_NilValue,
  bool convergence_b=false, bool convergence_d=false,
  bool convergence_component=false) {
 if (method!=4 && method!=5 && method!=6) throw std::invalid_argument(
  "mtblr_block_eigen_chains_raw_internal supports methods 4, 5, and 6 only");
 const std::size_t nt=wy.size();
 if (nt==0 || yy.size()!=nt || n.size()!=nt || b.size()!=nt)
  throw std::invalid_argument("MT block-eigen chain trait dimensions are inconsistent");
 const std::size_t m=wy[0].size();
 if (m==0 || (operator_descriptors.size()!=1 &&
     static_cast<std::size_t>(operator_descriptors.size())!=nt))
  throw std::invalid_argument("MT block-eigen operator ownership is invalid");
 for (std::size_t trait=0;trait<nt;++trait)
  if (wy[trait].size()!=m || b[trait].size()!=m)
   throw std::invalid_argument("MT block-eigen marker dimensions are inconsistent");

 const bool shared=operator_descriptors.size()==1;
 const std::size_t owner_count=shared ? 1 : nt;
 std::vector<MtBlockEigenDescriptor> descriptors;
 descriptors.reserve(owner_count);
 for (std::size_t owner=0;owner<owner_count;++owner)
  descriptors.push_back(parse_mt_block_eigen_descriptor(
   Rcpp::as<Rcpp::List>(operator_descriptors[static_cast<R_xlen_t>(owner)]),m));
 std::vector<sblr::core::BlockEigenStorage> storage_owners;
 storage_owners.reserve(owner_count);
 std::vector<std::vector<BlockEigenDiag>> owner_diagnostics(owner_count);
 std::vector<std::vector<double>> transformed_wy=wy;
 for (std::size_t owner=0;owner<owner_count;++owner) {
  const auto& descriptor=descriptors[owner];
  const int* rows=descriptor.rows0.empty()?nullptr:descriptor.rows0.data();
  PackedBedMatrix packed=read_bedfiles_to_packed_matrix(
   descriptor.bed_files,descriptor.n_bed,rows,
   static_cast<int>(descriptor.rows0.size()),descriptor.cls);
  const std::size_t rows_in_build=shared ? nt : 1;
  arma::mat wy_matrix(rows_in_build,m);
  for (std::size_t row=0;row<rows_in_build;++row) {
   const std::size_t trait=shared ? row : owner;
   for (std::size_t marker=0;marker<m;++marker)
    wy_matrix(static_cast<arma::uword>(row),static_cast<arma::uword>(marker))=
     wy[trait][marker];
  }
  storage_owners.push_back(build_block_eigen(
   packed,descriptor.af,descriptor.block_start,descriptor.filter,
   descriptor.tau,descriptor.eta,wy_matrix,1,&owner_diagnostics[owner]));
  for (std::size_t row=0;row<rows_in_build;++row) {
   const std::size_t trait=shared ? row : owner;
   for (std::size_t marker=0;marker<m;++marker)
    transformed_wy[trait][marker]=wy_matrix(
     static_cast<arma::uword>(row),static_cast<arma::uword>(marker));
  }
 }
 std::vector<sblr::core::BlockEigenView> views;
 views.reserve(nt);
 for (std::size_t trait=0;trait<nt;++trait)
  views.push_back(storage_owners[shared?0:trait].view());
 const sblr::mt::MtBlockEigenDataView data{
  transformed_wy,yy,n,sblr::mt::MtBlockEigenBundleView{m,std::move(views)}};
 sblr::mt::MtJointStateSpec joint;
 if (method==5 || method==6) {
  joint.patterns=models; joint.component=std::move(joint_component);
  joint.multiplier=std::move(joint_multiplier); joint.names=std::move(joint_names);
  joint.component_count=component_count; joint.scaled=true;
 }
 sblr::mt::MtBayesRCSpec bayesrc;
 if (method==6) {
  bayesrc.annotations=&annotations;
  bayesrc.pattern_prior=pattern_pi_prior;
  bayesrc.update_alpha=updateAlpha;
  bayesrc.intercept_flat=intercept_flat;
  bayesrc.sigma_alpha_a=sigma_alpha_a;
  bayesrc.sigma_alpha_b=sigma_alpha_b;
  bayesrc.pi_floor=pi_floor;
  bayesrc.alpha_update_every=alpha_update_every;
 }
 sblr::mt::MtExtendedTraceSpec convergence_spec;
 convergence_spec.covariance=convergence_covariance;
 convergence_spec.probability=convergence_probability;
 convergence_spec.annotations=convergence_annotations;
 convergence_spec.full_probability_states=convergence_full_probability;
 convergence_spec.selected_markers=Rf_isNull(convergence_markers) ?
  std::vector<int>() : Rcpp::as<std::vector<int>>(convergence_markers);
 convergence_spec.selected_b=convergence_b;
 convergence_spec.selected_d=convergence_d;
 convergence_spec.selected_component=convergence_component;
 const bool any_extended=convergence_spec.covariance ||
  convergence_spec.probability || convergence_spec.annotations ||
  !convergence_spec.selected_markers.empty();
 const sblr::mt::MtDefaultModelSpec model{
  models,sets,method,(method==5||method==6)?&joint:nullptr,
  (method==5||method==6)?&marker_scale:nullptr,
  method==5?&pi_prior:nullptr,method==6?&bayesrc:nullptr,
  any_extended?&convergence_spec:nullptr};
 const sblr::mt::MtDefaultCovariancePriorView prior{ssb_prior,sse_prior,nub,nue};
 const std::vector<int> seeds=mt_summary_resolve_chain_seeds(seed,nchains,chain_seeds);
 std::vector<std::string> errors; int used_workers=1;
 auto runner=[&](int chain) {
  const sblr::mt::MtDefaultExecutionSpec execution{
   updateB,updateE,updatePi,nit,nburn,nthin,seeds[static_cast<std::size_t>(chain)]};
  sblr::mt::MtDefaultInitialState initial{
   b,B,E,pi,component_init,beta_init,state_init,alpha_init,
   arma::vec(sigma_alpha_init),pattern_pi_init};
  auto core=sblr::mt::run_mt_block_eigen_core(
   data,model,prior,execution,std::move(initial));
  auto final=sblr::mt::finalize_mt_default_result(std::move(core));
  MtSummaryChainResult out;
  out.component=final.component;
  out.component_probabilities=final.component_probabilities;
  out.pi_trace=final.pi_trace;
  out.annotation_alpha_final=final.annotation_alpha_final;
  out.annotation_alpha_mean=final.annotation_alpha_mean;
  out.annotation_sigma_final=final.annotation_sigma_final;
  out.annotation_sigma_mean=final.annotation_sigma_mean;
  out.pattern_pi_final=final.pattern_pi_final;
  out.pattern_pi_mean=final.pattern_pi_mean;
  out.pattern_pi_trace=final.pattern_pi_trace;
  out.prior_component_probabilities=final.prior_component_probabilities;
  out.annotation_updates_attempted=final.annotation_updates_attempted;
  out.annotation_updates_completed=final.annotation_updates_completed;
  out.convergence=std::move(final.convergence);
  out.legacy=sblr::mt::make_mt_default_legacy_result(
   final,transformed_wy,nit,nburn);
  return out;
 };
 const auto results=mt_summary_dispatch_chains(
  nchains,ncores,runner,errors,used_workers);
 Rcpp::List owners(static_cast<R_xlen_t>(owner_count));
 for (std::size_t owner=0;owner<owner_count;++owner) {
  const auto& descriptor=descriptors[owner];
  owners[static_cast<R_xlen_t>(owner)]=Rcpp::List::create(
   Rcpp::_["eigen_filter"]=mt_block_eigen_filter_name(descriptor.filter),
   Rcpp::_["eigen_tau"]=descriptor.tau,Rcpp::_["eigen_eta"]=descriptor.eta,
   Rcpp::_["mu_floor"]=0.01,Rcpp::_["reference_sample_size"]=descriptor.n_bed,
   Rcpp::_["selected_row_count"]=descriptor.rows0.empty()?descriptor.n_bed:
    static_cast<int>(descriptor.rows0.size()),
   Rcpp::_["marker_count"]=static_cast<int>(descriptor.af.size()),
   Rcpp::_["blocks"]=block_eigen_diagnostics_to_data_frame(owner_diagnostics[owner]));
 }
 Rcpp::IntegerVector trait_owner(static_cast<R_xlen_t>(nt));
 for (std::size_t trait=0;trait<nt;++trait)
  trait_owner[static_cast<R_xlen_t>(trait)]=static_cast<int>(shared?1:trait+1);
 Rcpp::List block_diagnostics=Rcpp::List::create(
  Rcpp::_["owner_count"]=static_cast<int>(owner_count),
  Rcpp::_["trait_owner"]=trait_owner,
  Rcpp::_["sharing_mode"]=shared?"fully_shared_operator":"trait_specific_operator",
  Rcpp::_["owners"]=owners);
  return mt_summary_raw_list(results,models,
   method==6?"mt_block_eigen_bayesrc":method==5?"mt_block_eigen_bayesr":"mt_block_eigen_bayesc",
  nit,nburn,nthin,updateB,updateE,updatePi,seeds,ncores,used_workers,
  Rcpp::List::create(Rcpp::_["block_eigen"]=block_diagnostics));
}

namespace {

std::vector<std::string> mt_bed_character_vector(
 const Rcpp::CharacterVector& values,
 const char* label
) {
 if (values.size()==0) {
  throw std::invalid_argument(std::string(label)+" must be nonempty");
 }
 std::vector<std::string> out;
 out.reserve(values.size());
 for (R_xlen_t i=0; i<values.size(); ++i) {
  if (Rcpp::CharacterVector::is_na(values[i])) {
   throw std::invalid_argument(std::string(label)+" must not contain NA");
  }
  const std::string value=Rcpp::as<std::string>(values[i]);
  if (value.empty()) {
   throw std::invalid_argument(std::string(label)+" must not contain empty values");
  }
  out.push_back(value);
 }
 return out;
}

void validate_mt_bed_trait_names(const Rcpp::NumericMatrix& phenotype) {
 Rcpp::RObject dimnames=phenotype.attr("dimnames");
 if (dimnames.isNULL()) {
  throw std::invalid_argument(
   "Y must have nonempty unique trait column names");
 }
 Rcpp::List names(dimnames);
 if (names.size()!=2 || Rf_isNull(names[1])) {
  throw std::invalid_argument(
   "Y must have nonempty unique trait column names");
 }
 Rcpp::CharacterVector traits(names[1]);
 if (traits.size()!=phenotype.ncol()) {
  throw std::invalid_argument(
   "Y trait-name count must equal its column count");
 }
 std::set<std::string> seen;
 for (R_xlen_t trait=0; trait<traits.size(); ++trait) {
  if (Rcpp::CharacterVector::is_na(traits[trait])) {
   throw std::invalid_argument(
    "Y trait names must be nonmissing, nonempty, and unique");
  }
  const std::string name=Rcpp::as<std::string>(traits[trait]);
  if (name.empty() || !seen.insert(name).second) {
   throw std::invalid_argument(
    "Y trait names must be nonmissing, nonempty, and unique");
  }
 }
}

struct MtBedPreparedAdapter {
 std::unique_ptr<PackedBedMatrix> owner;
 std::vector<sblr::mt::MtBedMarkerMap> marker_maps;
 arma::mat phenotype;
 arma::mat marker_wy;
 std::vector<int> marker_order;
 std::vector<std::vector<double>> wy_trait_major;

 sblr::core::BedPackedGenotypeView<PackedBedMatrix> genotype_view() const {
  return sblr::core::BedPackedGenotypeView<PackedBedMatrix>{
   *owner,
   owner->data,
   static_cast<std::size_t>(owner->m)*owner->stride,
   static_cast<std::size_t>(owner->m),
   static_cast<std::size_t>(owner->n),
   owner->nbytes,
   owner->stride
  };
 }
};

MtBedPreparedAdapter prepare_mt_bed_adapter(
 Rcpp::CharacterVector bed_files,
 int n_bed,
 Rcpp::List cls,
 Rcpp::Nullable<Rcpp::IntegerVector> rows,
 Rcpp::NumericVector af,
 Rcpp::NumericMatrix Y,
 std::vector<double>& pi,
 const std::string& residual_covariance,
 int nit,
 int nburn,
 int nthin,
 int method
) {
 const std::vector<std::string> bed_paths=
  mt_bed_character_vector(bed_files, "bed_files");
 if (n_bed<=1) throw std::invalid_argument("n_bed must be greater than one");
 if (cls.size()!=bed_files.size()) {
  throw std::invalid_argument("bed_files length must equal cls length");
 }
 std::vector<std::vector<int>> selected_columns(cls.size());
 std::size_t marker_count=0;
 for (R_xlen_t file=0; file<cls.size(); ++file) {
  Rcpp::IntegerVector file_columns=cls[file];
  if (file_columns.size()==0) {
   throw std::invalid_argument("every cls entry must be nonempty");
  }
  selected_columns[file].reserve(file_columns.size());
  for (int value : file_columns) {
   if (value==NA_INTEGER || value<=0) {
    throw std::invalid_argument(
     "cls values must be positive one-based integers");
   }
   selected_columns[file].push_back(value);
   ++marker_count;
  }
 }
 if (af.size()!=static_cast<R_xlen_t>(marker_count)) {
  throw std::invalid_argument("af length must equal the selected marker count");
 }
 std::vector<double> frequencies(af.size());
 for (R_xlen_t marker=0; marker<af.size(); ++marker) {
  const double value=af[marker];
  if (!std::isfinite(value) || value<=0.0 || value>=1.0) {
   throw std::invalid_argument("af values must be finite and in (0, 1)");
  }
  frequencies[marker]=value;
 }

 std::vector<int> rows0;
 const int* row_pointer=nullptr;
 int selected_sample_count=n_bed;
 if (rows.isNotNull()) {
  Rcpp::IntegerVector selected_rows(rows);
  if (selected_rows.size()==0) {
   throw std::invalid_argument("rows must be nonempty when supplied");
  }
  std::set<int> seen;
  rows0.reserve(selected_rows.size());
  for (int value : selected_rows) {
   if (value==NA_INTEGER || value<=0 || value>n_bed) {
    throw std::invalid_argument("rows must be one-based integers in [1, n_bed]");
   }
   if (!seen.insert(value).second) {
    throw std::invalid_argument("rows must not contain duplicates");
   }
   rows0.push_back(value-1);
  }
  row_pointer=rows0.data();
  selected_sample_count=static_cast<int>(rows0.size());
 }
 if (Y.nrow()!=selected_sample_count || Y.ncol()==0) {
  throw std::invalid_argument(
   "Y must have selected sample count rows and at least one trait");
 }
 validate_mt_bed_trait_names(Y);
 for (double value : Y) {
  if (!std::isfinite(value)) {
   throw std::invalid_argument("Y must contain only finite values");
  }
 }
 if ((method!=4 && method!=5 && method!=6) || nit<=0 || nburn<0 || nthin<=0 ||
     (residual_covariance!="full" && residual_covariance!="diagonal")) {
  throw std::invalid_argument(
   "invalid method, residual covariance, or MCMC controls");
 }

 MtBedPreparedAdapter prepared;
 // Phase 17O architecture guard retained: PackedBedMatrix owner=
 // The unique_ptr below constructs that sole owner directly at its stationary
 // address so a prepared-adapter move cannot invalidate borrowed views.
 prepared.owner.reset(new PackedBedMatrix(read_bedfiles_to_packed_matrix(
  bed_paths, n_bed, row_pointer, static_cast<int>(rows0.size()),
  selected_columns)));
 const auto genotype=prepared.genotype_view();
 arma::vec workspace(prepared.owner->n, arma::fill::zeros);
 prepared.marker_maps=sblr::mt::build_mt_bed_marker_maps(
  genotype, frequencies);
 prepared.phenotype=Rcpp::as<arma::mat>(Y);
 prepared.marker_wy=sblr::mt::compute_mt_bed_marker_wy(
  genotype, prepared.marker_maps, prepared.phenotype, workspace);
 prepared.marker_order=sblr::mt::compute_mt_bed_marker_order(
  prepared.marker_wy, prepared.marker_maps);
 prepared.wy_trait_major=sblr::mt::mt_bed_trait_major(prepared.marker_wy);

 double pi_total=0.0;
 for (double value : pi) {
  if (!std::isfinite(value) || value<0.0) {
   throw std::invalid_argument("pi must be finite and nonnegative");
  }
  pi_total+=value;
 }
 if (!std::isfinite(pi_total) || pi_total<=0.0) {
  throw std::invalid_argument("pi must have positive total mass");
 }
 for (double& value : pi) value/=pi_total;
 return prepared;
}

Rcpp::List mt_bed_marker_kernel_to_list(
 const sblr::mt::MtBedMarkerKernelResult& kernel
) {
 const int nmodels=static_cast<int>(kernel.models.size());
 Rcpp::List precision(nmodels);
 Rcpp::List rhs(nmodels);
 Rcpp::List mean(nmodels);
 Rcpp::List covariance(nmodels);
 for (int model=0; model<nmodels; ++model) {
  precision[model]=kernel.models[model].precision;
  rhs[model]=Rcpp::NumericVector(
   kernel.models[model].rhs.begin(), kernel.models[model].rhs.end());
  mean[model]=Rcpp::NumericVector(
   kernel.models[model].mean.begin(), kernel.models[model].mean.end());
  covariance[model]=kernel.models[model].covariance;
 }
 return Rcpp::List::create(
  Rcpp::_["C"]=precision,
  Rcpp::_["rhs"]=rhs,
  Rcpp::_["mean"]=mean,
  Rcpp::_["covariance"]=covariance,
  Rcpp::_["log_weight"]=kernel.log_weight,
  Rcpp::_["probability"]=kernel.probability);
}

}  // namespace

// Deterministic maintenance seam for the production MT BED marker kernel.
// It performs no random draw and does not execute an MCMC iteration.
// [[Rcpp::export]]
Rcpp::List mtblr_bed_marker_contract_internal(
 Rcpp::NumericVector score,
 double xx,
 arma::mat B,
 arma::mat E,
 std::vector<std::vector<int>> models,
 std::vector<double> pi
) {
 if (score.size()==0 || score.size()>64 || models.size()>4096) {
  throw std::invalid_argument(
   "marker inspection requires 1-64 traits and at most 4096 models");
 }
 const std::size_t nt=static_cast<std::size_t>(score.size());
 sblr::mt::validate_mt_bed_spd(B, nt, "B");
 sblr::mt::validate_mt_bed_spd(E, nt, "E");
 if (models.empty() || models.size()!=pi.size()) {
  throw std::invalid_argument(
   "models and pi must be nonempty and have matching lengths");
 }
 double total=0.0;
 for (double value : pi) {
  if (!std::isfinite(value) || value<0.0) {
   throw std::invalid_argument(
    "pi must be finite and nonnegative");
  }
  total+=value;
 }
 if (!std::isfinite(total) || total<=0.0) {
  throw std::invalid_argument("pi must have positive total mass");
 }
 for (double& value : pi) value/=total;
 arma::vec score_vector=Rcpp::as<arma::vec>(score);
 return mt_bed_marker_kernel_to_list(
  sblr::mt::mt_bed_marker_kernel(
   score_vector, xx, arma::inv(B), arma::inv(E),
   models, pi));
}

// Deterministic development oracle for the shared MT BayesR joint-state
// conditional. No sampling or operator access occurs here.
// [[Rcpp::export]]
Rcpp::List mtblr_bayesr_marker_contract_internal(
 Rcpp::NumericVector score, Rcpp::NumericVector diagonal,
 arma::mat B, arma::mat E, std::vector<std::vector<int>> patterns,
 std::vector<int> component, std::vector<double> multiplier,
 std::vector<std::string> state_names, int component_count,
 std::vector<double> pi, double marker_scale
) {
 sblr::mt::MtJointStateSpec spec{
  std::move(patterns),std::move(component),std::move(multiplier),
  std::move(state_names),component_count,true};
 const auto kernel=sblr::mt::mt_joint_marker_kernel(
  Rcpp::as<arma::vec>(score),Rcpp::as<arma::vec>(diagonal),
  arma::inv(B),arma::inv(E),spec,pi,marker_scale);
 Rcpp::List means(kernel.states.size()),covariances(kernel.states.size()),
  precisions(kernel.states.size());
 for (std::size_t state=1;state<kernel.states.size();++state) {
  means[state]=Rcpp::wrap(kernel.states[state].mean);
  precisions[state]=Rcpp::wrap(kernel.states[state].precision);
  covariances[state]=Rcpp::wrap(arma::inv(kernel.states[state].precision));
 }
 return Rcpp::List::create(
  Rcpp::_["probability"]=kernel.probability,
  Rcpp::_["mean"]=means,Rcpp::_["covariance"]=covariances,
  Rcpp::_["precision"]=precisions);
}

// Internal serial one-chain individual-level multivariate BayesC route.
// This maintenance interface is deliberately not namespace-exported.
// [[Rcpp::export]]
Rcpp::List mtblr_bed_internal(
 Rcpp::CharacterVector bed_files,
 int n_bed,
 Rcpp::List cls,
 Rcpp::Nullable<Rcpp::IntegerVector> rows,
 Rcpp::NumericVector af,
 Rcpp::NumericMatrix Y,
 std::vector<std::vector<double>> beta_init,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<int>> state_init,
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
 std::string residual_covariance,
 int nit,
 int nburn,
 int nthin,
 int seed,
 int method=4
) {
 MtBedPreparedAdapter prepared=prepare_mt_bed_adapter(
  bed_files, n_bed, cls, rows, af, Y, pi, residual_covariance,
  nit, nburn, nthin, method);
 const auto genotype=prepared.genotype_view();

 sblr::mt::MtBedDataView<PackedBedMatrix> data{
  genotype, prepared.marker_maps, prepared.phenotype,
  prepared.marker_wy, prepared.marker_order
 };
 sblr::mt::MtBedInitialState initial{
  std::move(beta_init), std::move(b_init), std::move(state_init),
  {}, std::move(B), std::move(E), std::move(pi)
 };
 sblr::mt::MtBedExecutionSpec execution{
  updateB, updateE, updatePi, residual_covariance,
  nit, nburn, nthin, static_cast<std::uint32_t>(seed), method
 };
 sblr::mt::MtBedCoreResult core=sblr::mt::run_mt_bed_bayesc_core(
  data, initial, sets, ssb_prior, sse_prior, models, nub, nue,
  execution);
 const sblr::mt::MtBedCoreDiagnostics diagnostics=core.diagnostics;
 sblr::mt::MtDefaultFinalResult final_result=
  sblr::mt::finalize_mt_default_result(std::move(core.posterior));
 sblr::mt::MtDefaultLegacyResult legacy=
  sblr::mt::make_mt_default_legacy_result(
   final_result, prepared.wy_trait_major, nit, nburn);

 Rcpp::List mt_bed_diagnostics=Rcpp::List::create(
  Rcpp::_["residual_covariance"]=residual_covariance,
  Rcpp::_["sample_count"]=prepared.owner->n,
  Rcpp::_["marker_count"]=prepared.owner->m,
  Rcpp::_["trait_count"]=Y.ncol(),
  Rcpp::_["owner"]="PackedBedMatrix",
  Rcpp::_["view"]="BedPackedGenotypeView",
  Rcpp::_["genotype_scale"]="standardized_genotype",
  Rcpp::_["marker_workspace"]="double",
  Rcpp::_["marker_cholesky_jitter_attempts"]=
   static_cast<double>(diagnostics.marker_cholesky_jitter_attempts),
  Rcpp::_["marker_cholesky_max_increment"]=
   diagnostics.marker_cholesky_max_increment,
  Rcpp::_["full_e_updates"]=
   static_cast<double>(diagnostics.full_e_updates),
  Rcpp::_["diagonal_e_updates"]=
   static_cast<double>(diagnostics.diagonal_e_updates),
  Rcpp::_["sample_residual_returned"]=false,
  Rcpp::_["genetic_values_returned"]=false,
  Rcpp::_["cpo"]="unsupported",
  Rcpp::_["le_ld"]="trait_diagonal_decomposition");
 return mtblr_legacy_to_raw(
  legacy, models, "mt_bed_bayesc", "individual",
  nit, nburn, nthin, updateB, updateE, updatePi,
 Rcpp::List::create(Rcpp::_["mt_bed"]=mt_bed_diagnostics));
}

namespace {

Rcpp::NumericMatrix mt_bed_trait_matrix(
 const std::vector<std::vector<double>>& value
) {
 const int nt=static_cast<int>(value.size());
 const int nr=nt==0 ? 0 : static_cast<int>(value[0].size());
 Rcpp::NumericMatrix out(nr, nt);
 for (int trait=0; trait<nt; ++trait)
  for (int row=0; row<nr; ++row) out(row, trait)=value[trait][row];
 return out;
}

Rcpp::IntegerMatrix mt_bed_state_matrix(
 const std::vector<std::vector<int>>& value
) {
 const int nt=static_cast<int>(value.size());
 const int nr=nt==0 ? 0 : static_cast<int>(value[0].size());
 Rcpp::IntegerMatrix out(nr, nt);
 for (int trait=0; trait<nt; ++trait)
  for (int row=0; row<nr; ++row) out(row, trait)=value[trait][row];
 return out;
}

Rcpp::NumericMatrix mt_bed_row_matrix(
 const std::vector<std::vector<double>>& value
) {
 const int nr=static_cast<int>(value.size());
 const int nc=nr==0 ? 0 : static_cast<int>(value[0].size());
 Rcpp::NumericMatrix out(nr,nc);
 for (int row=0;row<nr;++row)
  for (int col=0;col<nc;++col) out(row,col)=value[row][col];
 return out;
}

Rcpp::NumericVector mt_bed_uint32_vector(
 const std::vector<sblr::mt::MtBedChainTask>& tasks
) {
 Rcpp::NumericVector out(tasks.size());
 for (std::size_t chain=0; chain<tasks.size(); ++chain)
  out[static_cast<R_xlen_t>(chain)]=static_cast<double>(tasks[chain].seed);
 return out;
}

Rcpp::List mt_bed_compact_chain_record(
 const sblr::mt::MtBedChainSummary& chain
) {
 Rcpp::List marker=Rcpp::List::create(
  Rcpp::_["bm"]=mt_bed_trait_matrix(chain.bm),
  Rcpp::_["dm"]=mt_bed_trait_matrix(chain.dm),
  Rcpp::_["b"]=mt_bed_trait_matrix(chain.b),
  Rcpp::_["state"]=mt_bed_state_matrix(chain.state));
 Rcpp::List pi=Rcpp::List::create(
  Rcpp::_["final"]=chain.pi_final,
  Rcpp::_["mean"]=chain.pi_mean);
 if (!chain.component_probabilities.empty() &&
     !chain.component_probabilities.front().empty()) {
  marker["component_final"]=Rcpp::wrap(chain.component_final);
  marker["component_probabilities"]=mt_bed_row_matrix(
   chain.component_probabilities);
  pi["trace"]=mt_bed_trait_matrix(chain.pi_trace);
 }
 const bool bayesrc=chain.annotation_alpha_final.n_elem>0;
 if (bayesrc) {
  pi["final"]=R_NilValue;
  pi["mean"]=R_NilValue;
  pi["trace"]=R_NilValue;
 }
 Rcpp::List out=Rcpp::List::create(
  Rcpp::_["chain"]=chain.chain+1,
  Rcpp::_["seed"]=static_cast<double>(chain.seed),
  Rcpp::_["marker"]=marker,
  Rcpp::_["trace"]=Rcpp::List::create(
   Rcpp::_["vbs"]=mt_bed_trait_matrix(chain.vbs),
   Rcpp::_["vgs"]=mt_bed_trait_matrix(chain.vgs),
   Rcpp::_["ves"]=mt_bed_trait_matrix(chain.ves),
   Rcpp::_["vle"]=mt_bed_trait_matrix(chain.vle),
   Rcpp::_["vld"]=mt_bed_trait_matrix(chain.vld)),
  Rcpp::_["variance"]=Rcpp::List::create(
   Rcpp::_["covb"]=mt_bed_trait_matrix(chain.covb),
   Rcpp::_["covg"]=mt_bed_trait_matrix(chain.covg),
   Rcpp::_["cove"]=mt_bed_trait_matrix(chain.cove),
   Rcpp::_["vb"]=chain.B,
   Rcpp::_["vg"]=chain.G,
   Rcpp::_["ve"]=chain.E),
  Rcpp::_["pi"]=pi,
  Rcpp::_["diagnostics"]=Rcpp::List::create(
   Rcpp::_["marker_cholesky_jitter_attempts"]=static_cast<double>(
    chain.diagnostics.marker_cholesky_jitter_attempts),
   Rcpp::_["marker_cholesky_max_increment"]=
    chain.diagnostics.marker_cholesky_max_increment,
   Rcpp::_["full_e_updates"]=static_cast<double>(
    chain.diagnostics.full_e_updates),
   Rcpp::_["diagonal_e_updates"]=static_cast<double>(
    chain.diagnostics.diagonal_e_updates),
   Rcpp::_["seconds"]=chain.seconds));
 if (bayesrc) out["model_parameters"]=Rcpp::List::create(
  Rcpp::_["annotation_coefficients_final"]=chain.annotation_alpha_final,
  Rcpp::_["annotation_coefficients_mean"]=chain.annotation_alpha_mean,
  Rcpp::_["annotation_variances_final"]=chain.annotation_sigma_final,
  Rcpp::_["annotation_variances_mean"]=chain.annotation_sigma_mean,
  Rcpp::_["pattern_pi_final"]=chain.pattern_pi_final,
  Rcpp::_["pattern_pi_mean"]=chain.pattern_pi_mean,
  Rcpp::_["prior_component_probabilities"]=chain.prior_component_probabilities);
 return out;
}

Rcpp::List mt_bed_convergence_trace_bundle_to_list(
 const sblr::mt::MtBedConvergenceTraceBundle& bundle
) {
 const R_xlen_t quantity_count=
  static_cast<R_xlen_t>(bundle.quantities.size());
 Rcpp::IntegerVector quantity_index(quantity_count);
 Rcpp::CharacterVector group(quantity_count);
 Rcpp::IntegerVector trait_index(quantity_count),trait2_index(quantity_count),
  marker_index(quantity_count),component_index(quantity_count),
  pattern_index(quantity_count),annotation_index(quantity_count),
  stick_index(quantity_count),tier(quantity_count);
 Rcpp::LogicalVector updated(quantity_count),derived(quantity_count),
  structural(quantity_count);
 bool extended=false;
 for (R_xlen_t quantity=0; quantity<quantity_count; ++quantity) {
  const sblr::mt::MtBedConvergenceQuantity& descriptor=
   bundle.quantities[static_cast<std::size_t>(quantity)];
  quantity_index[quantity]=static_cast<int>(quantity)+1;
  group[quantity]=descriptor.group;
  trait_index[quantity]=descriptor.trait+1;
  trait2_index[quantity]=descriptor.trait2+1;
  marker_index[quantity]=descriptor.marker+1;
  component_index[quantity]=descriptor.component+1;
  pattern_index[quantity]=descriptor.pattern+1;
  annotation_index[quantity]=descriptor.annotation+1;
  stick_index[quantity]=descriptor.stick+1;
  tier[quantity]=descriptor.tier;
  extended=extended || descriptor.tier>1;
  updated[quantity]=descriptor.updated;
  derived[quantity]=descriptor.derived;
  structural[quantity]=descriptor.structural;
 }
 Rcpp::NumericVector values(
  bundle.values.begin(), bundle.values.end());
 values.attr("dim")=Rcpp::IntegerVector::create(
  bundle.postburn_draws, bundle.nchains,
  static_cast<int>(quantity_count));
 return Rcpp::List::create(
  Rcpp::_["schema"]=Rcpp::List::create(
   Rcpp::_["class"]="mtblr_convergence_trace_bundle",
   Rcpp::_["version"]=1),
  Rcpp::_["scope"]=extended ? "extended" : "core",
  Rcpp::_["nchains"]=bundle.nchains,
  Rcpp::_["postburn_draws_per_chain"]=bundle.postburn_draws,
  Rcpp::_["quantities"]=Rcpp::DataFrame::create(
   Rcpp::_["quantity_index"]=quantity_index,
   Rcpp::_["group"]=group,
   Rcpp::_["trait_index"]=trait_index,
   Rcpp::_["trait2_index"]=trait2_index,
   Rcpp::_["marker_index"]=marker_index,
   Rcpp::_["component_index"]=component_index,
   Rcpp::_["pattern_index"]=pattern_index,
   Rcpp::_["annotation_index"]=annotation_index,
   Rcpp::_["stick_index"]=stick_index,
   Rcpp::_["tier"]=tier,
   Rcpp::_["updated"]=updated,
   Rcpp::_["derived"]=derived,
   Rcpp::_["structural"]=structural),
  Rcpp::_["values"]=values);
}

}  // namespace

Rcpp::List mtblr_bed_chains_binding_impl(
 Rcpp::CharacterVector bed_files,
 int n_bed,
 Rcpp::List cls,
 Rcpp::Nullable<Rcpp::IntegerVector> rows,
 Rcpp::NumericVector af,
 Rcpp::NumericMatrix Y,
 std::vector<std::vector<double>> beta_init,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<int>> state_init,
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
 std::string residual_covariance,
 int nit,
 int nburn,
 int nthin,
 int seed,
 int method,
 int nchains,
 int ncores,
 std::vector<int> chain_seeds,
 bool keep_chains,
 bool capture_convergence_traces,
 std::vector<int> joint_component,
 std::vector<double> joint_multiplier,
 std::vector<std::string> joint_names,
 int component_count,
 std::vector<double> marker_scale,
 std::vector<double> pi_prior,
 std::vector<int> component_init,
 arma::mat annotations, arma::mat alpha_init,
 std::vector<double> sigma_alpha_init,
 std::vector<double> pattern_pi_init,
 std::vector<double> pattern_pi_prior,
 bool updateAlpha, bool intercept_flat,
 double sigma_alpha_a, double sigma_alpha_b,
 double pi_floor, int alpha_update_every,
 bool convergence_covariance, bool convergence_probability,
 bool convergence_annotations, bool convergence_full_probability,
 std::vector<int> convergence_markers,
 bool convergence_b, bool convergence_d, bool convergence_component
) {
 if (nchains<=0 || ncores<=0) {
  throw std::invalid_argument("nchains and ncores must be positive integers");
 }
 if (!chain_seeds.empty() &&
     chain_seeds.size()!=static_cast<std::size_t>(nchains)) {
  throw std::invalid_argument("chain_seeds must be empty or have length nchains");
 }
 MtBedPreparedAdapter prepared=prepare_mt_bed_adapter(
  bed_files, n_bed, cls, rows, af, Y, pi, residual_covariance,
  nit, nburn, nthin, method);
 const auto genotype=prepared.genotype_view();
 sblr::mt::MtBedDataView<PackedBedMatrix> data{
  genotype, prepared.marker_maps, prepared.phenotype,
  prepared.marker_wy, prepared.marker_order
 };
 const sblr::mt::MtBedInitialState initial{
  std::move(beta_init), std::move(b_init), std::move(state_init),
  std::move(component_init), std::move(B), std::move(E), std::move(pi),
  alpha_init,arma::vec(sigma_alpha_init),pattern_pi_init
 };
 sblr::mt::MtJointStateSpec joint;
 if (method==5 || method==6) {
  joint.patterns=models;
  joint.component=std::move(joint_component);
  joint.multiplier=std::move(joint_multiplier);
  joint.names=std::move(joint_names);
  joint.component_count=component_count;
  joint.scaled=true;
 }
 sblr::mt::MtBayesRCSpec bayesrc;
 if (method==6) {
  bayesrc.annotations=&annotations;
  bayesrc.pattern_prior=pattern_pi_prior;
  bayesrc.update_alpha=updateAlpha;
  bayesrc.intercept_flat=intercept_flat;
  bayesrc.sigma_alpha_a=sigma_alpha_a;
  bayesrc.sigma_alpha_b=sigma_alpha_b;
  bayesrc.pi_floor=pi_floor;
  bayesrc.alpha_update_every=alpha_update_every;
 }
 sblr::mt::MtExtendedTraceSpec convergence_spec;
 convergence_spec.covariance=convergence_covariance;
 convergence_spec.probability=convergence_probability;
 convergence_spec.annotations=convergence_annotations;
 convergence_spec.full_probability_states=convergence_full_probability;
 convergence_spec.selected_markers=std::move(convergence_markers);
 convergence_spec.selected_b=convergence_b;
 convergence_spec.selected_d=convergence_d;
 convergence_spec.selected_component=convergence_component;
 const bool any_extended=convergence_spec.covariance ||
  convergence_spec.probability || convergence_spec.annotations ||
  !convergence_spec.selected_markers.empty();
 sblr::mt::MtBedExecutionSpec execution{
  updateB, updateE, updatePi, residual_covariance,
  nit, nburn, nthin, static_cast<std::uint32_t>(seed), method
 };
 std::vector<sblr::mt::MtBedChainTask> tasks(
  static_cast<std::size_t>(nchains));
 for (int chain=0; chain<nchains; ++chain) {
  tasks[static_cast<std::size_t>(chain)].chain=chain;
  tasks[static_cast<std::size_t>(chain)].seed=chain_seeds.empty() ?
   sblr::mt::resolve_mt_bed_chain_seed(seed, chain) :
   sblr::mt::resolve_mt_bed_explicit_chain_seed(
    chain_seeds[static_cast<std::size_t>(chain)]);
 }
#ifdef _OPENMP
 const bool openmp_available=true;
 const int used_workers=std::min(ncores, nchains);
#else
 const bool openmp_available=false;
 const int used_workers=1;
 if (ncores>1) {
  Rcpp::warning(
   "OpenMP chain parallelism is unavailable; using serial execution");
 }
#endif
 const auto dispatch_start=std::chrono::steady_clock::now();
 std::vector<sblr::mt::MtBedChainExecutionResult> results=
  sblr::mt::dispatch_mt_bed_chain_tasks(
   tasks, used_workers, data, initial, sets, ssb_prior, sse_prior,
   models, nub, nue, execution, (method==5||method==6)?&joint:nullptr,
   (method==5||method==6)?&marker_scale:nullptr,
   method==5?&pi_prior:nullptr,method==6?&bayesrc:nullptr,
   any_extended?&convergence_spec:nullptr);
 const double dispatch_seconds=std::chrono::duration<double>(
  std::chrono::steady_clock::now()-dispatch_start).count();

 std::ostringstream failures;
 bool any_failure=false;
 for (std::size_t chain=0; chain<results.size(); ++chain) {
  if (results[chain].failed) {
   if (!any_failure) failures << (capture_convergence_traces ?
    "mtblr_bed_convergence_trace_internal failed:" :
    "mtblr_bed_chains_internal failed:");
   failures << "\nchain " << chain+1 << ": " << results[chain].error;
   any_failure=true;
  }
 }
 if (any_failure) throw std::runtime_error(failures.str());

 sblr::mt::MtBedConvergenceTraceBundle convergence_bundle;
 if (capture_convergence_traces) {
  convergence_bundle=sblr::mt::build_mt_bed_convergence_trace_bundle(
   results, nburn, nit, updateB, updateE, updatePi, updateAlpha, method,
   residual_covariance, any_extended?&convergence_spec:nullptr);
 }

 sblr::mt::MtBedChainsAggregateResult aggregate=
  sblr::mt::aggregate_mt_bed_chains(results, keep_chains);
 const double marker_count=aggregate.pooled.marker_retained_count;
 const double covb_count=aggregate.pooled.covb_retained_count;
 const double covg_count=aggregate.pooled.covg_retained_count;
 const double cove_count=aggregate.pooled.cove_retained_count;
 const double pi_count=aggregate.pooled.pi_retained_count;
 sblr::mt::MtDefaultFinalResult final_result=
  sblr::mt::finalize_mt_bed_chains_result(std::move(aggregate.pooled));
 const std::vector<int> component_final=final_result.component;
 const std::vector<std::vector<double>> component_probabilities=
  final_result.component_probabilities;
 const std::vector<std::vector<double>> pi_trace=final_result.pi_trace;
 sblr::mt::MtDefaultLegacyResult legacy=
  sblr::mt::make_mt_default_legacy_result(
   final_result, prepared.wy_trait_major, nit, nburn);

 Rcpp::NumericVector chain_seconds(results.size());
 double seconds_sum=0.0;
 double seconds_max=0.0;
 for (std::size_t chain=0; chain<results.size(); ++chain) {
  chain_seconds[static_cast<R_xlen_t>(chain)]=results[chain].seconds;
  seconds_sum+=results[chain].seconds;
  seconds_max=std::max(seconds_max, results[chain].seconds);
 }
 const sblr::mt::MtBedChainsDiagnostics& diagnostics=aggregate.diagnostics;
 Rcpp::List mt_bed_diagnostics=Rcpp::List::create(
  Rcpp::_["residual_covariance"]=residual_covariance,
  Rcpp::_["sample_count"]=prepared.owner->n,
  Rcpp::_["marker_count"]=prepared.owner->m,
  Rcpp::_["trait_count"]=Y.ncol(),
  Rcpp::_["owner"]="PackedBedMatrix",
  Rcpp::_["view"]="BedPackedGenotypeView",
  Rcpp::_["genotype_scale"]="standardized_genotype",
  Rcpp::_["marker_workspace"]="double",
  Rcpp::_["marker_cholesky_jitter_attempts"]=static_cast<double>(
   diagnostics.marker_cholesky_jitter_attempts),
  Rcpp::_["marker_cholesky_max_increment"]=
   diagnostics.marker_cholesky_max_increment,
  Rcpp::_["full_e_updates"]=static_cast<double>(diagnostics.full_e_updates),
  Rcpp::_["diagonal_e_updates"]=static_cast<double>(
   diagnostics.diagonal_e_updates),
  Rcpp::_["sample_residual_returned"]=false,
  Rcpp::_["genetic_values_returned"]=false,
  Rcpp::_["cpo"]="unsupported",
  Rcpp::_["le_ld"]="trait_diagonal_decomposition",
  Rcpp::_["requested_cores"]=ncores,
  Rcpp::_["used_workers"]=used_workers,
  Rcpp::_["openmp_available"]=openmp_available,
  Rcpp::_["chain_seeds"]=mt_bed_uint32_vector(tasks),
  Rcpp::_["chain_seconds"]=chain_seconds,
  Rcpp::_["seconds_mean"]=seconds_sum/static_cast<double>(nchains),
  Rcpp::_["seconds_max"]=seconds_max,
  Rcpp::_["dispatch_seconds"]=dispatch_seconds,
  Rcpp::_["primary_chain"]=1,
  Rcpp::_["final_state_policy"]="primary_chain",
  Rcpp::_["posterior_summary_policy"]="pooled_retained_samples",
  Rcpp::_["trace_policy"]="iterationwise_chain_mean",
  Rcpp::_["chain_marker_cholesky_jitter_attempts"]=
   diagnostics.chain_marker_cholesky_jitter_attempts,
  Rcpp::_["chain_marker_cholesky_max_increment"]=
   diagnostics.chain_marker_cholesky_max_increment,
  Rcpp::_["chain_full_e_updates"]=diagnostics.chain_full_e_updates,
  Rcpp::_["chain_diagonal_e_updates"]=
   diagnostics.chain_diagonal_e_updates);
 Rcpp::List raw=mtblr_legacy_to_raw(
  legacy, models, method==6?"mt_bed_bayesrc":method==5?"mt_bed_bayesr":"mt_bed_bayesc", "individual",
  nit, nburn, nthin, updateB, updateE, updatePi,
  Rcpp::List::create(Rcpp::_["mt_bed"]=mt_bed_diagnostics));
 Rcpp::List meta=raw["meta"];
 meta["nchains"]=nchains;
 meta["keep_chains"]=keep_chains;
 raw["meta"]=meta;
 Rcpp::List marker=raw["marker"];
 marker["bm_sd"]=mt_bed_trait_matrix(aggregate.bm_sd);
 marker["bm_min"]=mt_bed_trait_matrix(aggregate.bm_min);
 marker["bm_max"]=mt_bed_trait_matrix(aggregate.bm_max);
 marker["dm_sd"]=mt_bed_trait_matrix(aggregate.dm_sd);
 marker["dm_min"]=mt_bed_trait_matrix(aggregate.dm_min);
 marker["dm_max"]=mt_bed_trait_matrix(aggregate.dm_max);
 if (method==5 || method==6) {
  marker["component_final"]=Rcpp::wrap(component_final);
  marker["component_probabilities"]=mt_bed_row_matrix(
   component_probabilities);
  Rcpp::List pi_namespace=raw["pi"];
  pi_namespace["trace"]=mt_bed_trait_matrix(pi_trace);
  raw["pi"]=pi_namespace;
 }
 if (method==6) {
  raw["annotations"]=Rcpp::List::create(
   Rcpp::_ ["annotation_coefficients_final"]=final_result.annotation_alpha_final,
   Rcpp::_ ["annotation_coefficients_mean"]=final_result.annotation_alpha_mean,
   Rcpp::_ ["annotation_variances_final"]=final_result.annotation_sigma_final,
   Rcpp::_ ["annotation_variances_mean"]=final_result.annotation_sigma_mean,
   Rcpp::_ ["pattern_pi_final"]=final_result.pattern_pi_final,
   Rcpp::_ ["pattern_pi_mean"]=final_result.pattern_pi_mean,
   Rcpp::_ ["pattern_pi_trace"]=mt_bed_trait_matrix(final_result.pattern_pi_trace),
   Rcpp::_ ["prior_component_probabilities"]=final_result.prior_component_probabilities,
   Rcpp::_ ["annotation_updates_attempted"]=final_result.annotation_updates_attempted,
   Rcpp::_ ["annotation_updates_completed"]=final_result.annotation_updates_completed);
 }
 raw["marker"]=marker;
 Rcpp::List raw_diagnostics=raw["diagnostics"];
 raw_diagnostics["marker"]=static_cast<int>(marker_count);
 raw_diagnostics["covb"]=static_cast<int>(covb_count);
 raw_diagnostics["covg"]=static_cast<int>(covg_count);
 raw_diagnostics["cove"]=static_cast<int>(cove_count);
 raw_diagnostics["pi"]=static_cast<int>(pi_count);
 raw["diagnostics"]=raw_diagnostics;
 if (keep_chains) {
  Rcpp::List chains(aggregate.chains.size());
  Rcpp::CharacterVector names(aggregate.chains.size());
  for (std::size_t chain=0; chain<aggregate.chains.size(); ++chain) {
   chains[static_cast<R_xlen_t>(chain)]=
    mt_bed_compact_chain_record(aggregate.chains[chain]);
   names[static_cast<R_xlen_t>(chain)]="chain"+std::to_string(chain+1);
  }
  chains.attr("names")=names;
  raw["chains"]=chains;
 } else {
  raw["chains"]=R_NilValue;
 }
 if (capture_convergence_traces) {
  return Rcpp::List::create(
   Rcpp::_["raw"]=raw,
   Rcpp::_["trace_bundle"]=
    mt_bed_convergence_trace_bundle_to_list(convergence_bundle));
 }
 return raw;
}

// Internal multichain individual-level multivariate BayesC route. It prepares
// one immutable packed-BED dataset and dispatches one unchanged core per chain.
// [[Rcpp::export]]
Rcpp::List mtblr_bed_chains_internal(
 Rcpp::CharacterVector bed_files, int n_bed, Rcpp::List cls,
 Rcpp::Nullable<Rcpp::IntegerVector> rows, Rcpp::NumericVector af,
 Rcpp::NumericMatrix Y, std::vector<std::vector<double>> beta_init,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<int>> state_init,
 const std::vector<std::vector<int>>& sets, arma::mat B, arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<std::vector<int>> models, std::vector<double> pi,
 double nub, double nue, bool updateB, bool updateE, bool updatePi,
 std::string residual_covariance, int nit, int nburn, int nthin,
 int seed, int method, int nchains, int ncores,
 std::vector<int> chain_seeds, bool keep_chains,
 std::vector<int> joint_component,
 std::vector<double> joint_multiplier,
 std::vector<std::string> joint_names,
 int component_count,
 std::vector<double> marker_scale,
 std::vector<double> pi_prior,
 std::vector<int> component_init,
 arma::mat annotations, arma::mat alpha_init,
 std::vector<double> sigma_alpha_init,
 std::vector<double> pattern_pi_init,
 std::vector<double> pattern_pi_prior,
 bool updateAlpha, bool intercept_flat,
 double sigma_alpha_a, double sigma_alpha_b,
 double pi_floor, int alpha_update_every,
 bool convergence_covariance=false, bool convergence_probability=false,
 bool convergence_annotations=false,
 bool convergence_full_probability=false,
 SEXP convergence_markers=R_NilValue,
 bool convergence_b=false, bool convergence_d=false,
 bool convergence_component=false
) {
 return mtblr_bed_chains_binding_impl(
  bed_files, n_bed, cls, rows, af, Y, std::move(beta_init),
  std::move(b_init), std::move(state_init), sets, std::move(B),
  std::move(E), std::move(ssb_prior), std::move(sse_prior),
  std::move(models), std::move(pi), nub, nue, updateB, updateE,
  updatePi, residual_covariance, nit, nburn, nthin, seed, method,
  nchains, ncores, std::move(chain_seeds), keep_chains, false,
  std::move(joint_component),std::move(joint_multiplier),
  std::move(joint_names),component_count,std::move(marker_scale),
  std::move(pi_prior),std::move(component_init),std::move(annotations),
  std::move(alpha_init),std::move(sigma_alpha_init),
  std::move(pattern_pi_init),std::move(pattern_pi_prior),updateAlpha,
  intercept_flat,sigma_alpha_a,sigma_alpha_b,pi_floor,alpha_update_every,
  convergence_covariance,convergence_probability,convergence_annotations,
  convergence_full_probability,Rf_isNull(convergence_markers) ?
   std::vector<int>() : Rcpp::as<std::vector<int>>(convergence_markers),
  convergence_b,convergence_d,convergence_component);
}

// Internal MT BED multichain route with a post-burn Tier 1 trace bundle.
// It is deliberately not namespace-exported and computes no diagnostics.
// [[Rcpp::export]]
Rcpp::List mtblr_bed_convergence_trace_internal(
 Rcpp::CharacterVector bed_files, int n_bed, Rcpp::List cls,
 Rcpp::Nullable<Rcpp::IntegerVector> rows, Rcpp::NumericVector af,
 Rcpp::NumericMatrix Y, std::vector<std::vector<double>> beta_init,
 std::vector<std::vector<double>> b_init,
 std::vector<std::vector<int>> state_init,
 const std::vector<std::vector<int>>& sets, arma::mat B, arma::mat E,
 std::vector<std::vector<double>> ssb_prior,
 std::vector<std::vector<double>> sse_prior,
 std::vector<std::vector<int>> models, std::vector<double> pi,
 double nub, double nue, bool updateB, bool updateE, bool updatePi,
 std::string residual_covariance, int nit, int nburn, int nthin,
 int seed, int method, int nchains, int ncores,
 std::vector<int> chain_seeds, bool keep_chains,
 std::vector<int> joint_component,
 std::vector<double> joint_multiplier,
 std::vector<std::string> joint_names,
 int component_count,
 std::vector<double> marker_scale,
 std::vector<double> pi_prior,
 std::vector<int> component_init,
 arma::mat annotations, arma::mat alpha_init,
 std::vector<double> sigma_alpha_init,
 std::vector<double> pattern_pi_init,
 std::vector<double> pattern_pi_prior,
 bool updateAlpha, bool intercept_flat,
 double sigma_alpha_a, double sigma_alpha_b,
 double pi_floor, int alpha_update_every,
 bool convergence_covariance=false, bool convergence_probability=false,
 bool convergence_annotations=false,
 bool convergence_full_probability=false,
 SEXP convergence_markers=R_NilValue,
 bool convergence_b=false, bool convergence_d=false,
 bool convergence_component=false
) {
 return mtblr_bed_chains_binding_impl(
  bed_files, n_bed, cls, rows, af, Y, std::move(beta_init),
  std::move(b_init), std::move(state_init), sets, std::move(B),
  std::move(E), std::move(ssb_prior), std::move(sse_prior),
  std::move(models), std::move(pi), nub, nue, updateB, updateE,
  updatePi, residual_covariance, nit, nburn, nthin, seed, method,
  nchains, ncores, std::move(chain_seeds), keep_chains, true,
  std::move(joint_component),std::move(joint_multiplier),
  std::move(joint_names),component_count,std::move(marker_scale),
  std::move(pi_prior),std::move(component_init),std::move(annotations),
  std::move(alpha_init),std::move(sigma_alpha_init),
  std::move(pattern_pi_init),std::move(pattern_pi_prior),updateAlpha,
  intercept_flat,sigma_alpha_a,sigma_alpha_b,pi_floor,alpha_update_every,
  convergence_covariance,convergence_probability,convergence_annotations,
  convergence_full_probability,Rf_isNull(convergence_markers) ?
   std::vector<int>() : Rcpp::as<std::vector<int>>(convergence_markers),
  convergence_b,convergence_d,convergence_component);
}

// INTERNAL RESEARCH ONLY: not publicly routed or supported. Retained until the
// shared scalar/MT block-eigen representation and per-trait operators exist.
// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtblr_eigen(   std::vector<std::vector<double>> wy,
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
 int q = wy[0].size();
 int m = b[0].size();
 int nmodels = models.size();
 double nsamples=0.0;

 //double logliksum, detC, diff, cumprobc;
 //int mselect;

 std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));

 std::vector<double> mu(nt), rhs(nt), conv(nt);
 std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> ves(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit+nburn, 0.0));
 std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> mus(nt, std::vector<double>(nit+nburn, 0.0));

 std::vector<double> x2t(m);
 std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
 std::vector<int> order(m);
 std::vector<double> gamma(4, 0.0);

 std::vector<std::vector<double>> pitrait(nt, std::vector<double>(4, 0.0));
 std::vector<std::vector<double>> pistrait(nt, std::vector<double>(4, 0.0));
 std::vector<double> pimarker(2, 0.0);
 std::vector<double> pismarker(2, 0.0);
 std::vector<int> dmarker(m, 0);

 gamma[0] = 0.0;
 gamma[1] = 0.01;
 gamma[2] = 0.1;
 gamma[3] = 1.0;

 for (int t = 0; t < nt; t++) {
  pitrait[t][0] = 1.0-pi[0];
  pitrait[t][1] = pi[0];
 }
 if(method==5){
  for (int t = 0; t < nt; t++) {
   pitrait[t][0] = 0.95;
   pitrait[t][1] = 0.02;
   pitrait[t][2] = 0.02;
   pitrait[t][3] = 0.01;
  }
 }

 pimarker[0] = 1.0-pi[0];
 pimarker[1] = pi[0];

 // Initialize variables
 //for (int i = 0; i < m; i++) {
 for (int i = 0; i < q; i++) {
  for (int t = 0; t < nt; t++) {
   r[t][i] = wy[t][i];
  }
 }

 for (int i = 0; i < q; i++) {
  for (int t = 0; t < nt; t++) {
   double rhsi =0.0;
   for (size_t j = 0; j < XXvalues[t][i].size(); j++) {
    rhsi +=r[t][j]*XXvalues[t][i][j];
   }
   x2[t][i] = rhsi*rhsi;
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


 // Initialize (co)variance matrices
 arma::mat C(nt,nt, fill::zeros);
 arma::mat Bi = arma::inv(B);
 arma::mat Ei = arma::inv(E);
 //arma::mat G(nt,nt, fill::zeros);
 arma::mat G = arma::eye<arma::mat>(nt, nt);
 arma::rowvec probs(nmodels, fill::zeros);

 // Start Gibbs sampler
 std::random_device rd;
 std::mt19937 gen(seed);

 for ( int it = 0; it < nit+nburn; it++) {

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

  //Sample marker effects (BayesC)
  if (method==4) {
   for (size_t s = 0; s < sets.size(); ++s) {
    const std::vector<int>& set = sets[s];

    // Step 1: Sample B using only markers in the current set
    //sampleBset(nt, m, nub, B, d, b, ssb_prior, set, gen);

    // Step 2: Invert B to get Bi
    //arma::mat Bi = arma::inv(B);

    // Step 3: Sample marker effects for markers in this set
    // Convert current set to a std::unordered_set for fast lookup
    std::unordered_set<int> current_set(set.begin(), set.end());

    for (int isort = 0; isort < m; isort++) {
     int i = order[isort];

     // Only update if marker i is in the current set
     if (current_set.find(i) != current_set.end()) {
      sampleBetaCMt_eigen(i, nt,
                          nmodels, models, cmodel, pi,
                          Ei, Bi,
                          ww, r, b, d,
                          XXindices, XXvalues,
                          gen);
     }
    }
   }
  }

  arma::mat G(nt, nt, arma::fill::zeros);
  for (int t1 = 0; t1 < nt; t1++) {
   for (int t2 = t1; t2 < nt; t2++) {
    double ssg = 0.0;
    //for (int i = 0; i < m; i++) {
    for (int i = 0; i < q; i++) {
     double e1 = wy[t1][i] - r[t1][i];
     double e2 = wy[t2][i] - r[t2][i];
     ssg += e1 * e2;
    }
    G(t1, t2) = ssg;
    if (t1 != t2) G(t2, t1) = ssg; // Symmetric assignment
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
   samplePi(cmodel, pi, gen);
   for (int k = 0; k<nmodels ; k++) {
    if(it>nburn) pis[k] = pis[k] + pi[k];
   }
  }

  // Sample marker variance
  if(updateB && method==4) {
   sampleB(nt, m, nub, B, d, b, ssb_prior, gen);
   Bi = arma::inv(B);
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
  //computeG(nt, m, b, wy, r, n, G);
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
   sampleE_eigen(nt, m, nue, E, b, wy, r, sse_prior, yy, n, gen);
   for (int t = 0; t < nt; t++) {
    ves[t][it] = E(t,t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if(it>nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
    }
   }
   // Compute inverse (optional)
   arma::mat Ei = arma::inv(E);  // use or return this if needed

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
   result[18][t][i] = pistrait[t][i]/nsamples;
  }
  for (int i=0; i < 2; i++) {
   result[19][t][i] = pismarker[i]/nsamples;
  }
 }
 return result;
}

