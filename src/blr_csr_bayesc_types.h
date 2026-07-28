#ifndef SBLR_CORE_BLR_CSR_BAYESC_TYPES_H
#define SBLR_CORE_BLR_CSR_BAYESC_TYPES_H

// Match the package's established Armadillo index ABI while keeping this
// header independent of any binding library.
#if !defined(ARMA_USE_LAPACK)
#define ARMA_USE_LAPACK
#endif
#if !defined(ARMA_USE_BLAS)
#define ARMA_USE_BLAS
#endif
#define ARMA_HAVE_STD_ISFINITE
#define ARMA_HAVE_STD_ISINF
#define ARMA_HAVE_STD_ISNAN
#if defined(_WIN32) && !defined(ARMA_USE_OPENMP)
#define ARMA_USE_OPENMP
#define ARMA_DONT_PRINT_OPENMP_WARNING 1
#endif
#if !defined(ARMA_64BIT_WORD) && !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD 1
#endif
#include <armadillo>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "blr_result.h"
#include "blr_sparse_ld_csr.h"
#include "blr_spec.h"

namespace sblr {
namespace core {

// Borrowed immutable view of the one shared, pre-scaled symmetric CSR object.
// The binding owns every referenced buffer and must keep it alive until
// run_csr_bayesc() and all trait-chain tasks have returned. No member exposes
// mutable CSR access, and no chain result or chain state owns a CSR payload.
struct CsrBayesCDataView {
  std::size_t marker_count = 0;
  std::size_t trait_count = 0;
  SparseLdCsrView ld;
  const arma::mat* wy = nullptr;  // traits x markers
  const arma::vec* yy = nullptr;  // traits
  const std::vector<int>* sample_size = nullptr;
};

struct CsrBayesCPriors {
  const arma::mat* marker_variance = nullptr;
  const arma::mat* residual_variance = nullptr;
  const arma::mat* marker_scale_prior = nullptr;
  const arma::mat* residual_scale_prior = nullptr;
  std::vector<double> inclusion_probability;
  double marker_degrees_freedom = 0.0;
  double residual_degrees_freedom = 0.0;
  double inclusion_prior_active = 0.0;
  double inclusion_prior_null = 0.0;
};

struct CsrBayesCInitialState {
  const arma::mat* effects = nullptr;  // traits x markers
  const std::vector<std::vector<double>>* inclusion = nullptr;
  const std::vector<std::vector<double>>* residual = nullptr;
  bool use_inclusion = false;
  bool use_residual = false;
};

struct CsrBayesCControls {
  int nit = 0;
  int nburn = 0;
  int nthin = 1;
  int nchains = 1;
  int ncores = 1;
  int seed = 1;
  std::vector<int> chain_seeds;
  bool keep_chains = false;
  bool update_marker_variance = true;
  bool update_residual_variance = true;
  bool update_inclusion_probability = true;
  bool rebuild_residual_before_update = false;
  double residual_adjustment = 0.0;
  bool update_ld_swap = false;
  double ld_swap_probability = 0.0;
  int ld_swap_moves = 0;
  bool use_fixed_maf_effect_scale = false;
  const arma::rowvec* fixed_maf_effect_scale = nullptr;
  bool estimate_maf_effect_s = false;
  double maf_effect_s_initial = 0.0;
  double maf_effect_s_prior_lower = -3.0;
  double maf_effect_s_prior_upper = 2.0;
  double maf_effect_s_proposal_sd = 0.35;
  const arma::rowvec* maf_effect_s_log_h = nullptr;
  std::vector<int> convergence_markers;
  bool convergence_b = false;
  bool convergence_d = false;
};

struct CsrBayesCOutputSpec {
  bool keep_chains = false;
};

struct CsrBayesCLdFriendsView {
  const std::uint64_t* row_ptr = nullptr;
  std::size_t row_ptr_size = 0;
  const int* index = nullptr;
  std::size_t friend_count = 0;
};

struct CsrBayesCExecutionInput {
  ResolvedSpec specification;
  CsrBayesCDataView data;
  CsrBayesCPriors priors;
  CsrBayesCInitialState initial;
  CsrBayesCControls controls;
  CsrBayesCOutputSpec output;
  CsrBayesCLdFriendsView ld_friends;
  const std::vector<int>* marker_order = nullptr;
};

// Chain-owned mutable output/state vocabulary. Effects, residuals, inclusion
// state, parameter traces, accumulators, RNG state, and workspace are local to
// one task in the core. Deliberately contains no CSR values or indices.
struct CsrBayesCChainResult {
  arma::rowvec marker_mean;
  arma::rowvec marker_pip;
  arma::rowvec final_effect;
  arma::rowvec final_residual;
  arma::Row<int> final_state;
  arma::rowvec marker_variance_trace;
  arma::rowvec genetic_variance_trace;
  arma::rowvec residual_variance_trace;
  arma::rowvec inclusion_trace;
  arma::rowvec le_variance_trace;
  arma::rowvec ld_variance_trace;
  arma::rowvec maf_effect_s_trace;
  arma::mat convergence_b;
  arma::imat convergence_d;
  double final_marker_variance = 0.0;
  double final_genetic_variance = 0.0;
  double final_residual_variance = 0.0;
  double final_le_variance = 0.0;
  double final_ld_variance = 0.0;
  double final_inclusion_probability = 0.0;
  double retained_samples = 0.0;
  double ld_swap_attempted = 0.0;
  double ld_swap_accepted = 0.0;
  double maf_effect_s_attempted = 0.0;
  double maf_effect_s_accepted = 0.0;
  int thread_used = 0;
  double seconds = 0.0;
};

// BayesC extension of the canonical typed result vocabulary. Matrices retain
// the existing trait-major native orientation; the binding alone transposes
// them into the canonical markers x traits and samples x traits R shapes.
struct CsrBayesCResult : BlrResult {
  std::size_t marker_count = 0;
  std::size_t trait_count = 0;
  std::size_t trace_count = 0;
  int chain_count = 0;
  arma::mat marker_mean;
  arma::mat marker_pip;
  arma::mat marker_score;
  arma::mat final_residual;
  arma::mat final_effect;
  arma::mat final_state;
  arma::mat marker_mean_sd;
  arma::mat marker_mean_min;
  arma::mat marker_mean_max;
  arma::mat marker_pip_sd;
  arma::mat marker_pip_min;
  arma::mat marker_pip_max;
  arma::mat marker_variance_trace;
  arma::mat genetic_variance_trace;
  arma::mat residual_variance_trace;
  arma::mat inclusion_trace;
  arma::mat le_variance_trace;
  arma::mat ld_variance_trace;
  arma::mat maf_effect_s_trace;
  arma::vec final_marker_variance;
  arma::vec final_genetic_variance;
  arma::vec final_residual_variance;
  arma::vec final_le_variance;
  arma::vec final_ld_variance;
  arma::vec final_inclusion_probability;
  arma::vec retained_samples;
  arma::vec ld_swap_attempted;
  arma::vec ld_swap_accepted;
  arma::vec maf_effect_s_attempted;
  arma::vec maf_effect_s_accepted;
  std::vector<CsrBayesCChainResult> chains;
};

void validate_csr_bayesc_execution_input(
  const CsrBayesCExecutionInput& input
);

CsrBayesCResult run_csr_bayesc(const CsrBayesCExecutionInput& input);

}  // namespace core
}  // namespace sblr

#endif
