// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "distributions.h"

#include <random>
#include <stdexcept>
#include <cmath>
#include <numeric>

// g++ -O3 -std=c++17 -I/path/to/armadillo/include -c distributions.cpp

// ======================================================
// Scalar RNG helpers
// ======================================================

double rnorm_cpp(double mean, double sd, std::mt19937& gen) {
 if (sd < 0.0) {
  throw std::runtime_error("rnorm_cpp: sd must be non-negative.");
 }
 std::normal_distribution<double> dist(mean, sd);
 return dist(gen);
}

double rgamma_cpp(double shape, double rate, std::mt19937& gen) {
 if (shape <= 0.0 || rate <= 0.0) {
  throw std::runtime_error("rgamma_cpp: shape and rate must be positive.");
 }
 std::gamma_distribution<double> dist(shape, 1.0 / rate); // scale = 1/rate
 return dist(gen);
}

double rinvgamma_cpp(double shape, double rate, std::mt19937& gen) {
 return 1.0 / rgamma_cpp(shape, rate, gen);
}

double runif_cpp(std::mt19937& gen) {
 std::uniform_real_distribution<double> dist(0.0, 1.0);
 return dist(gen);
}

int rbernoulli_cpp(double p, std::mt19937& gen) {
 if (p < 0.0 || p > 1.0) {
  throw std::runtime_error("rbernoulli_cpp: p must be in [0, 1].");
 }
 std::bernoulli_distribution dist(p);
 return static_cast<int>(dist(gen));
}

// ======================================================
// Gaussian
// ======================================================

arma::vec rnorm_vec_cpp(int n, double mean, double sd, std::mt19937& gen) {
 if (n < 0) {
  throw std::runtime_error("rnorm_vec_cpp: n must be non-negative.");
 }

 arma::vec x(n);
 for (int i = 0; i < n; ++i) {
  x(i) = rnorm_cpp(mean, sd, gen);
 }
 return x;
}

arma::vec rmvnorm0_cpp(const arma::mat& Sigma, std::mt19937& gen) {
 if (Sigma.n_rows != Sigma.n_cols) {
  throw std::runtime_error("rmvnorm0_cpp: Sigma must be square.");
 }

 arma::mat C;
 if (!arma::chol(C, Sigma)) {
  throw std::runtime_error("rmvnorm0_cpp: Sigma must be SPD.");
 }

 const arma::uword p = Sigma.n_cols;
 arma::vec z(p);

 for (arma::uword i = 0; i < p; ++i) {
  z(i) = rnorm_cpp(0.0, 1.0, gen);
 }

 return C.t() * z;
}

arma::vec rmvnorm_cpp(const arma::vec& mu,
                      const arma::mat& Sigma,
                      std::mt19937& gen) {
 if (Sigma.n_rows != Sigma.n_cols) {
  throw std::runtime_error("rmvnorm_cpp: Sigma must be square.");
 }
 if (mu.n_elem != Sigma.n_rows) {
  throw std::runtime_error("rmvnorm_cpp: dimension mismatch.");
 }

 return mu + rmvnorm0_cpp(Sigma, gen);
}

// ======================================================
// Covariance
// ======================================================

arma::mat rwishart_cpp(unsigned int df,
                       const arma::mat& V,
                       std::mt19937& gen) {
 const arma::uword p = V.n_rows;

 if (V.n_rows != V.n_cols) {
  throw std::runtime_error("rwishart_cpp: V must be square.");
 }
 if (df <= p - 1) {
  throw std::runtime_error("rwishart_cpp: df must be > p - 1.");
 }

 arma::mat C;
 if (!arma::chol(C, V, "lower")) {
  throw std::runtime_error("rwishart_cpp: V must be SPD.");
 }

 arma::mat A(p, p, arma::fill::zeros);

 for (arma::uword i = 0; i < p; ++i) {
  std::chi_squared_distribution<double> chisq(df - i);
  A(i, i) = std::sqrt(chisq(gen));

  std::normal_distribution<double> norm(0.0, 1.0);
  for (arma::uword j = 0; j < i; ++j) {
   A(i, j) = norm(gen);
  }
 }

 arma::mat W = C * A * A.t() * C.t();
 return 0.5 * (W + W.t());
}

arma::mat rinvwishart_cpp(unsigned int df,
                          const arma::mat& S,
                          std::mt19937& gen) {
 if (S.n_rows != S.n_cols) {
  throw std::runtime_error("rinvwishart_cpp: S must be square.");
 }

 arma::mat S_inv;
 if (!arma::inv_sympd(S_inv, S)) {
  throw std::runtime_error("rinvwishart_cpp: S must be SPD.");
 }

 arma::mat W = rwishart_cpp(df, S_inv, gen);

 arma::mat Sigma;
 if (!arma::inv_sympd(Sigma, W)) {
  throw std::runtime_error("rinvwishart_cpp: Wishart draw inversion failed.");
 }

 return 0.5 * (Sigma + Sigma.t());
}

// ======================================================
// Discrete
// ======================================================

int rcategorical_cpp(const arma::vec& probs, std::mt19937& gen) {
 if (probs.n_elem == 0) {
  throw std::runtime_error("rcategorical_cpp: probs is empty.");
 }

 double sum_probs = arma::accu(probs);
 if (!std::isfinite(sum_probs) || sum_probs <= 0.0) {
  throw std::runtime_error("rcategorical_cpp: probabilities must sum to a positive finite value.");
 }

 double u = runif_cpp(gen);
 double cum = 0.0;

 for (arma::uword i = 0; i < probs.n_elem; ++i) {
  double p = probs(i) / sum_probs;
  if (p < 0.0 || !std::isfinite(p)) {
   throw std::runtime_error("rcategorical_cpp: probabilities must be non-negative and finite.");
  }
  cum += p;
  if (u <= cum) {
   return static_cast<int>(i);
  }
 }

 return static_cast<int>(probs.n_elem - 1);
}

// ======================================================
// Dirichlet
// ======================================================

arma::vec rdirichlet_cpp(const arma::vec& alpha, std::mt19937& gen) {
 if (alpha.n_elem == 0) {
  throw std::runtime_error("rdirichlet_cpp: alpha is empty.");
 }

 arma::vec x(alpha.n_elem);
 double sum_x = 0.0;

 for (arma::uword i = 0; i < alpha.n_elem; ++i) {
  if (alpha(i) <= 0.0 || !std::isfinite(alpha(i))) {
   throw std::runtime_error("rdirichlet_cpp: alpha entries must be positive and finite.");
  }
  x(i) = rgamma_cpp(alpha(i), 1.0, gen);
  sum_x += x(i);
 }

 if (!std::isfinite(sum_x) || sum_x <= 0.0) {
  throw std::runtime_error("rdirichlet_cpp: invalid gamma sum.");
 }

 return x / sum_x;
}
