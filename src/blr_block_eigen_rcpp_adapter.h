#ifndef SBLR_BLR_BLOCK_EIGEN_RCPP_ADAPTER_H
#define SBLR_BLR_BLOCK_EIGEN_RCPP_ADAPTER_H

#include <RcppArmadillo.h>

#include "blr_csr_bayesc_policy.h"
#include "blr_csr_bayesr_policy.h"
#include "blr_phase3_execution.h"

Rcpp::List stblr_cpg_omp_csr_block_eigen_with_policy(
  std::vector<std::vector<double>>, std::vector<std::vector<double>>,
  std::vector<double>, std::vector<std::vector<double>>,
  std::vector<std::vector<double>>, bool,
  std::vector<std::vector<double>>, bool, bool, std::string,
  arma::mat, arma::mat, std::vector<std::vector<double>>,
  std::vector<std::vector<double>>, std::vector<double>, double, double,
  bool, bool, bool, double, std::vector<int>, int, int, int, double, double,
  int, int, int, bool, std::vector<int>, bool, double, double, int, int,
  Rcpp::Nullable<Rcpp::NumericVector>, bool, double, Rcpp::NumericVector,
  double, Rcpp::Nullable<Rcpp::NumericVector>, Rcpp::IntegerVector,
  bool, bool, Rcpp::CharacterVector, int, Rcpp::List,
  Rcpp::Nullable<Rcpp::IntegerVector>, Rcpp::NumericVector,
  Rcpp::IntegerVector, std::string, double, double, std::string, double, int,
  CsrBayesCPolicyFactory*, const BlrPhase3ExecutionContract&
);

Rcpp::List stblr_cpg_omp_csr_bayesr_block_eigen_with_policy(
  std::vector<std::vector<double>>, std::vector<std::vector<double>>,
  std::vector<double>, std::vector<std::vector<double>>,
  std::vector<std::vector<double>>, bool,
  std::vector<std::vector<double>>, bool, bool, std::string,
  arma::mat, arma::mat, std::vector<std::vector<double>>,
  std::vector<std::vector<double>>, std::vector<double>,
  std::vector<double>, std::vector<double>, double, double, bool, bool, bool,
  double, std::vector<int>, int, int, int, int, int, int, bool,
  std::vector<int>, int, int, bool, double, double, int, int,
  Rcpp::Nullable<Rcpp::NumericVector>, bool, double, Rcpp::NumericVector,
  double, Rcpp::Nullable<Rcpp::NumericVector>, Rcpp::IntegerVector,
  bool, bool, bool, bool, Rcpp::CharacterVector, int, Rcpp::List,
  Rcpp::Nullable<Rcpp::IntegerVector>, Rcpp::NumericVector,
  Rcpp::IntegerVector, std::string, double, double, std::string, double, int,
  Rcpp::List, CsrBayesRPolicyFactory*, const BlrPhase3ExecutionContract&
);

#endif
