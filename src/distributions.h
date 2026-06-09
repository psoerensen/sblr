#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <random>

// ======================================================
// Distribution sampling utilities (pure C++)
//
// Uses:
//  - std::mt19937 RNG
//  - Armadillo for linear algebra
//
// No R / Rcpp dependencies
// ======================================================


// ======================================================
// Scalar RNG
// ======================================================

double rnorm_cpp(double mean, double sd, std::mt19937& gen);
double rgamma_cpp(double shape, double rate, std::mt19937& gen);
double rinvgamma_cpp(double shape, double rate, std::mt19937& gen);
double runif_cpp(std::mt19937& gen);


// ======================================================
// Gaussian
// ======================================================

arma::vec rnorm_vec_cpp(int n, double mean, double sd, std::mt19937& gen);

arma::vec rmvnorm0_cpp(const arma::mat& Sigma,
                       std::mt19937& gen);

arma::vec rmvnorm_cpp(const arma::vec& mu,
                      const arma::mat& Sigma,
                      std::mt19937& gen);


// ======================================================
// Covariance (Wishart)
// ======================================================

arma::mat rwishart_cpp(unsigned int df,
                       const arma::mat& V,
                       std::mt19937& gen);

arma::mat rinvwishart_cpp(unsigned int df,
                          const arma::mat& S,
                          std::mt19937& gen);


// ======================================================
// Discrete
// ======================================================

int rbernoulli_cpp(double p, std::mt19937& gen);

int rcategorical_cpp(const arma::vec& probs,
                     std::mt19937& gen);


// ======================================================
// Dirichlet
// ======================================================

arma::vec rdirichlet_cpp(const arma::vec& alpha,
                         std::mt19937& gen);
