// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

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
 const arma::uword p = V.n_rows;

 if (V.n_rows != V.n_cols) {
  throw std::runtime_error("rwishart: V must be square.");
 }
 if (df <= p - 1) {
  throw std::runtime_error("rwishart: df must be > p - 1.");
 }

 arma::mat C;
 bool ok = arma::chol(C, V, "lower");
 if (!ok) {
  throw std::runtime_error("rwishart: V must be SPD.");
 }

 arma::mat A(p, p, arma::fill::zeros);

 for (arma::uword i = 0; i < p; ++i) {
  std::chi_squared_distribution<double> rchisq(df - i);
  A(i, i) = std::sqrt(rchisq(gen));

  std::normal_distribution<double> rnorm(0.0, 1.0);
  for (arma::uword j = 0; j < i; ++j) {
   A(i, j) = rnorm(gen);
  }
 }

 arma::mat W = C * A * A.t() * C.t();
 W = 0.5 * (W + W.t());
 return W;
}

// Draw Sigma ~ InvWishart(df, S)
// using: if W ~ Wishart(df, S^{-1}), then Sigma = W^{-1}
arma::mat rinvwishart(unsigned int df,
                      const arma::mat& S,
                      std::mt19937& gen) {
 if (S.n_rows != S.n_cols) {
  throw std::runtime_error("rinvwishart: S must be square.");
 }

 arma::mat S_inv = arma::inv_sympd(S);
 arma::mat W = rwishart(df, S_inv, gen);
 arma::mat Sigma = arma::inv_sympd(W);

 Sigma = 0.5 * (Sigma + Sigma.t());
 return Sigma;
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
void sampleBetaCPG_Mt_latent(
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
  std::vector<std::vector<double>>& beta,   // latent effects
  std::vector<std::vector<double>>& b,      // effective effects = D * beta
  std::vector<std::vector<int>>& d,
  const std::vector<std::vector<int>>& XXindices,
  const std::vector<std::vector<std::vector<double>>>& XXvalues,
  std::mt19937& gen)
{
 // ------------------------------------------------------------
 // 1) Base RHS using current EFFECTIVE effect
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
   const size_t nnz = XXindices[i].size();
   for (size_t j = 0; j < nnz; ++j) {
    r[t][XXindices[i][j]] -= XXvalues[t][i][j] * diff;
   }
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

 const int nt=core_result.nt;
 const int m=core_result.m;
 const int nmodels=core_result.nmodels;
 const double marker_retained_count=core_result.marker_retained_count;
 const double covb_retained_count=core_result.covb_retained_count;
 const double covg_retained_count=core_result.covg_retained_count;
 const double cove_retained_count=core_result.cove_retained_count;
 const double pi_retained_count=core_result.pi_retained_count;
 const auto& bm=core_result.bm;
 const auto& dm=core_result.dm;
 const auto& r=core_result.r;
 const auto& final_b=core_result.b;
 const auto& d=core_result.d;
 const auto& order=core_result.order;
 const auto& vbs=core_result.vbs;
 const auto& vgs=core_result.vgs;
 const auto& ves=core_result.ves;
 const auto& cvbm=core_result.cvbm;
 const auto& cvgm=core_result.cvgm;
 const auto& cvem=core_result.cvem;
 const auto& final_B=core_result.B;
 const auto& G=core_result.G;
 const auto& final_E=core_result.E;
 const auto& final_pi=core_result.pi;
 const auto& pis=core_result.pis;
 const auto& pistrait=core_result.pistrait;
 const auto& pismarker=core_result.pismarker;

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
   result[0][t][i] = bm[t][i]/marker_retained_count;
   result[1][t][i] = dm[t][i]/marker_retained_count;
   result[2][t][i] = wy[t][i];
   result[3][t][i] = r[t][i];
   result[4][t][i] = final_b[t][i];
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
   result[10][t1][t2] = covb_retained_count > 0.0 ?
    cvbm[t1][t2] / covb_retained_count : 0.0;
   result[11][t1][t2] = covg_retained_count > 0.0 ?
    cvgm[t1][t2] / covg_retained_count : 0.0;
   result[12][t1][t2] = cove_retained_count > 0.0 ?
    cvem[t1][t2] / cove_retained_count : 0.0;
   result[13][t1][t2] = final_B(t1,t2);
   result[14][t1][t2] = G(t1,t2);
   result[15][t1][t2] = final_E(t1,t2);
  }
 }
 for (int t=0; t < nt; t++) {
  for (int i=0; i < nmodels; i++) {
   result[16][t][i] = final_pi[i];
   result[17][t][i] = pi_retained_count > 0.0 ?
    pis[i] / pi_retained_count : 0.0;
  }
 }
 for (int t=0; t < nt; t++) {
  for (int i=0; i < 4; i++) {
   result[18][t][i] = pistrait[t][i]/nit;
  }
  for (int i=0; i < 2; i++) {
   result[19][t][i] = pismarker[i]/nit;
  }
 }
 return result;
}

// [[Rcpp::export]]
std::vector<std::vector<std::vector<double>>>  mtblr_hybrid(
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
  int method) {

 // Define local variables
 int nt = wy.size();
 int m = wy[0].size();
 int nmodels = models.size();
 int nsets = sets.size();
 double nsamples = 0.0;

 // Hybrid settings
 int warmup_iter = 100;

 std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));

 std::vector<double> mu(nt), rhs(nt), conv(nt);
 std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> ves(nt, std::vector<double>(nit + nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(nit + nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(nit + nburn, 0.0));
 std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> mus(nt, std::vector<double>(nit + nburn, 0.0));

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

 // region x trait mask
 std::vector<std::vector<int>> trait_allowed(nsets, std::vector<int>(nt, 1));

 gamma[0] = 0.0;
 gamma[1] = 0.01;
 gamma[2] = 0.1;
 gamma[3] = 1.0;

 for (int t = 0; t < nt; t++) {
  pitrait[t][0] = 1.0 - pi[0];
  pitrait[t][1] = pi[0];
 }
 if (method == 5) {
  for (int t = 0; t < nt; t++) {
   pitrait[t][0] = 0.95;
   pitrait[t][1] = 0.02;
   pitrait[t][2] = 0.02;
   pitrait[t][3] = 0.01;
  }
 }

 pimarker[0] = 1.0 - pi[0];
 pimarker[1] = pi[0];

 // ------------------------------------------
 // Check whether ww is approximately constant
 // ------------------------------------------
 double ww_sum = 0.0;
 long long ww_count = 0;

 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   ww_sum += ww[t][i];
   ww_count++;
  }
 }

 double ww_const = ww_sum / std::max(1LL, ww_count);

 double max_rel_dev = 0.0;
 for (int t = 0; t < nt; ++t) {
  for (int i = 0; i < m; ++i) {
   double rel = std::abs(ww[t][i] - ww_const) / std::max(1e-12, std::abs(ww_const));
   if (rel > max_rel_dev) max_rel_dev = rel;
  }
 }

 // Use cache only if ww is effectively constant
 //bool use_cached_models = (max_rel_dev < 1e-6);
 bool use_cached_models = true;

 // Initialize variables
 for (int i = 0; i < m; i++) {
  for (int t = 0; t < nt; t++) {
   r[t][i] = wy[t][i];
   double wwi = std::max(std::abs(ww[t][i]), 1e-12);
   x2[t][i] = (wy[t][i] / wwi) * (wy[t][i] / wwi);
   //x2[t][i] = (wy[t][i] / ww[t][i]) * (wy[t][i] / ww[t][i]);
  }
 }

 // Wy - W'Wb
 for (int i = 0; i < m; i++) {
  for (int t = 0; t < nt; t++) {
   if (b[t][i] != 0.0) {
    for (size_t j = 0; j < XXindices[i].size(); j++) {
     r[t][XXindices[i][j]] = r[t][XXindices[i][j]] - XXvalues[t][i][j] * b[t][i];
    }
   }
  }
 }

 for (int i = 0; i < m; i++) {
  x2t[i] = 0.0;
  for (int t = 0; t < nt; t++) {
   if (x2[t][i] > x2t[i]) {
    x2t[i] = x2[t][i];
   }
  }
 }

 // Establish order of markers as they are entered into the model
 std::iota(order.begin(), order.end(), 0);
 std::sort(std::begin(order),
           std::end(order),
           [&](int i1, int i2) { return x2t[i1] > x2t[i2]; });

 // Initialize (co)variance matrices
 arma::mat C(nt, nt, arma::fill::zeros);
 arma::mat Bi = arma::inv(B);
 arma::mat Ei = arma::inv(E);
 arma::mat G(nt, nt, arma::fill::zeros);
 arma::rowvec probs(nmodels, arma::fill::zeros);

 // Start Gibbs sampler
 std::mt19937 gen(seed);

 for (int it = 0; it < nit + nburn; it++) {

  bool warmup_st_only = (it < warmup_iter);

  std::fill(cmodel.begin(), cmodel.end(), 0.0);

  // Sample marker effects (BayesC / hybrid)
  if (method == 4) {
   arma::mat Bi_current = arma::inv(B);

   for (size_t s = 0; s < sets.size(); ++s) {
    const std::vector<int>& set = sets[s];

    std::vector<ModelCache> cache;
    if (use_cached_models) {
     precomputeModelCache(
      nt,
      nmodels,
      models,
      Ei,
      Bi_current,
      ww_const,
      trait_allowed[s],
                   cache
     );
    }

    std::unordered_set<int> current_set(set.begin(), set.end());

    for (int isort = 0; isort < m; isort++) {
     int i = order[isort];

     if (current_set.find(i) != current_set.end()) {
      if (use_cached_models) {
       sampleBetaCMtMaskedFast(
        i,
        static_cast<int>(s),
        nt,
        nmodels,
        models,
        trait_allowed,
        warmup_st_only,
        cmodel,
        pi,
        Ei,
        cache,
        ww_const,
        r,
        b,
        d,
        XXindices,
        XXvalues,
        gen);
      }
     }
    }
   }
  }

  // After warmup: detect relevant traits in each region
  if (it == warmup_iter) {

   int k_keep = std::min(3, nt);

   for (size_t s = 0; s < sets.size(); ++s) {

    std::vector<double> score(nt, 0.0);

    for (int idx : sets[s]) {
     for (int t = 0; t < nt; ++t) {
      score[t] += std::abs(b[t][idx]);
     }
    }

    // rank traits
    std::vector<int> order_trait(nt);
    std::iota(order_trait.begin(), order_trait.end(), 0);

    std::sort(order_trait.begin(), order_trait.end(),
              [&](int a, int b){ return score[a] > score[b]; });

    // reset mask
    for (int t = 0; t < nt; ++t)
     trait_allowed[s][t] = 0;

    // keep top traits
    for (int k = 0; k < k_keep; ++k)
     trait_allowed[s][order_trait[k]] = 1;

    // random exploration
    std::uniform_real_distribution<double> runif(0.0, 1.0);
    for (int t = 0; t < nt; ++t) {
     if (trait_allowed[s][t] == 0) {
      if (runif(gen) < 0.05)
       trait_allowed[s][t] = 1;
     }
    }
   }
  }

  // Store values
  for (int t = 0; t < nt; t++) {
   for (int i = 0; i < m; i++) {
    if (d[t][i] > 0) {
     if ((it > nburn) && (it % nthin == 0)) {
      dm[t][i] = dm[t][i] + 1.0;
      bm[t][i] = bm[t][i] + b[t][i];
     }
    }
   }
  }

  // Sample pi for Bayes C
  if (updatePi && method == 4) {
   samplePi(cmodel, pi, gen);
   for (int k = 0; k < nmodels; k++) {
    if (it > nburn) pis[k] = pis[k] + pi[k];
   }
  }

  // Sample marker variance
  if (updateB && method == 4) {
   sampleB(nt, m, nub, B, d, b, ssb_prior, gen);
   for (int t = 0; t < nt; t++) {
    vbs[t][it] = B(t, t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if (it > nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1, t2);
    }
   }
  }

  // Update genetic variance
  computeG(nt, m, b, wy, r, n, G);
  for (int t = 0; t < nt; t++) {
   vgs[t][it] = G(t, t);
  }
  for (int t1 = 0; t1 < nt; t1++) {
   for (int t2 = 0; t2 < nt; t2++) {
    if (it > nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1, t2);
   }
  }

  // Sample residual variance
  if (updateE) {
   sampleE(nt, m, nue, E, b, wy, r, sse_prior, yy, n, gen);
   Ei = arma::inv(E);
   for (int t = 0; t < nt; t++) {
    ves[t][it] = E(t, t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if (it > nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1, t2);
    }
   }
  }

  if ((it > nburn) && (it % nthin == 0)) {
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

 for (int t = 0; t < nt; t++) {
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

 for (int t = 0; t < nt; t++) {
  for (int i = 0; i < m; i++) {
   result[0][t][i] = bm[t][i] / nsamples;
   result[1][t][i] = dm[t][i] / nsamples;
   result[2][t][i] = wy[t][i];
   result[3][t][i] = r[t][i];
   result[4][t][i] = b[t][i];
   result[5][t][i] = d[t][i];
   result[6][t][i] = order[i];
  }
 }

 for (int t = 0; t < nt; t++) {
  for (int i = 0; i < nit + nburn; i++) {
   result[7][t][i] = vbs[t][i];
   result[8][t][i] = vgs[t][i];
   result[9][t][i] = ves[t][i];
  }
 }

 for (int t1 = 0; t1 < nt; t1++) {
  for (int t2 = 0; t2 < nt; t2++) {
   result[10][t1][t2] = cvbm[t1][t2] / nsamples;
   result[11][t1][t2] = cvgm[t1][t2] / nsamples;
   result[12][t1][t2] = cvem[t1][t2] / nsamples;
   result[13][t1][t2] = B(t1, t2);
   result[14][t1][t2] = G(t1, t2);
   result[15][t1][t2] = E(t1, t2);
  }
 }

 for (int t = 0; t < nt; t++) {
  for (int i = 0; i < nmodels; i++) {
   result[16][t][i] = pi[i];
   result[17][t][i] = pis[i] / nsamples;
  }
 }

 for (int t = 0; t < nt; t++) {
  for (int i = 0; i < 4; i++) {
   result[18][t][i] = pistrait[t][i] / nit;
  }
  for (int i = 0; i < 2; i++) {
   result[19][t][i] = pismarker[i] / nit;
  }
 }

 return result;
}

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

