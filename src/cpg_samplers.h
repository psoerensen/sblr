#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
#include <random>

void sampleB_cpg_arma(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const arma::mat& beta,
  const arma::mat& ssb_prior,
  std::mt19937& gen);

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
  std::mt19937& gen);

void computeG_cpg_arma(
  int nt,
  int m,
  const arma::mat& b,
  const arma::mat& wy,
  const arma::mat& r,
  const std::vector<int>& n,
  arma::mat& G);

void samplePi_cpg(
  std::vector<double>& cmodel,
  std::vector<double>& pi,
  std::mt19937& gen);

void sampleB_cpg(
  int nt,
  int m,
  int nub,
  arma::mat& B,
  const std::vector<std::vector<double>>& beta,
  const std::vector<std::vector<double>>& ssb_prior,
  std::mt19937& gen);

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
  std::mt19937& gen);

void computeG_cpg(
  int nt,
  int m,
  const std::vector<std::vector<double>>& b,
  const std::vector<std::vector<double>>& wy,
  const std::vector<std::vector<double>>& r,
  const std::vector<int>& n,
  arma::mat& G);
