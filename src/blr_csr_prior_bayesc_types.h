#ifndef SBLR_BLR_CSR_PRIOR_BAYESC_TYPES_H
#define SBLR_BLR_CSR_PRIOR_BAYESC_TYPES_H

#include <armadillo>

#include "blr_csr_annotation_bayesc_types.h"
#include "blr_phase3_execution.h"
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace sblr { namespace core {

struct CsrPriorBayesCLdFriendsView {
 std::vector<std::uint64_t> ptr;
 std::vector<int> idx;
 std::vector<double> r2;
};

// All storage is prepared and owned by the binding function. The callable core
// borrows it for the duration of execution; mutable sampler state is identified
// explicitly and no binding-backed object is constructed in worker regions.
struct CsrPriorBayesCExecutionContext {
 int marker_count=0, trait_count=0;
 int iterations=0, burnin=0, thinning=1, cores=1, seed=0;
 int chain_index=0, chain_count=1;
 bool use_initial_inclusion=false, use_initial_residual=false;
 bool rebuild_residual_before_update=false;
 bool use_marker_probability=false, use_marker_multiplier=false;
 bool update_marker_variance=true, update_residual_variance=true;
 bool update_global_probability=true, update_ld_swap=false;
 double marker_degrees_freedom=0.0, residual_degrees_freedom=0.0;
 double residual_adjustment=0.0;
 double inclusion_prior_active=0.0, inclusion_prior_null=0.0;
 double ld_swap_probability=0.0;
 int ld_swap_moves=0;

 arma::mat* marker_score=nullptr;
 const arma::mat* marker_diagonal=nullptr;
 arma::mat* effect=nullptr;
 arma::mat* residual=nullptr;
 arma::Mat<int>* inclusion=nullptr;
 const arma::mat* marker_probability=nullptr;
 const arma::mat* marker_multiplier=nullptr;
 const arma::vec* phenotype_sum_squares=nullptr;
 const arma::mat* marker_scale_prior=nullptr;
 const arma::mat* residual_scale_prior=nullptr;
 const arma::mat* marker_variance_initial=nullptr;
 const arma::mat* residual_variance_initial=nullptr;
 const std::vector<double>* global_probability=nullptr;
 const std::vector<int>* sample_size=nullptr;
 const std::vector<std::vector<double>>* initial_inclusion=nullptr;
 const std::vector<std::vector<double>>* initial_residual=nullptr;
 const void* ld_storage=nullptr;
 std::size_t ld_row_ptr_count=0;
 const CsrPriorBayesCLdFriendsView* ld_friends=nullptr;
 const std::vector<int>* marker_order=nullptr;
 const std::vector<int>* convergence_markers=nullptr;
 bool convergence_b=false, convergence_d=false;
 const BlrPhase3ExecutionContract* execution_contract=nullptr;
};

struct CsrPriorBayesCExecutionResult {
 std::vector<std::vector<std::vector<double>>> raw;
 std::vector<arma::mat> convergence_b;
 std::vector<arma::imat> convergence_d;
 int requested_thread_count=1, configured_thread_count=1;
 int actual_team_size=1;
 std::vector<int> trait_worker_id;
};

inline void validate_csr_prior_bayesc_execution_context(
 const CsrPriorBayesCExecutionContext& x
) {
 AnnotationBayesCDataView data;
 data.marker_count=static_cast<std::size_t>(x.marker_count);
 data.trait_count=static_cast<std::size_t>(x.trait_count);
 data.row_ptr_count=x.ld_row_ptr_count;
 data.diagonal_count=x.marker_diagonal ? x.marker_diagonal->n_cols : 0;
 data.sample_size_count=x.sample_size ? x.sample_size->size() : 0;
 AnnotationBayesCControls controls;
 controls.iterations=x.iterations; controls.burnin=x.burnin;
 controls.thinning=x.thinning; controls.cores=x.cores;
 controls.update_ld_swap=x.update_ld_swap;
 controls.ld_swap_probability=x.ld_swap_probability;
 controls.ld_swap_moves=x.ld_swap_moves;
 validate_annotation_bayesc_common(data, controls);

 FixedPriorBayesCPolicyView policy;
 policy.marker_probability=x.marker_probability ? x.marker_probability->memptr() : nullptr;
 policy.probability_count=x.use_marker_probability && x.marker_probability ? x.marker_probability->n_elem : 0;
 policy.marker_multiplier=x.marker_multiplier ? x.marker_multiplier->memptr() : nullptr;
 policy.multiplier_count=x.use_marker_multiplier && x.marker_multiplier ? x.marker_multiplier->n_elem : 0;
 policy.use_marker_probability=x.use_marker_probability;
 policy.use_marker_multiplier=x.use_marker_multiplier;
 validate_fixed_prior_policy(policy, data.marker_count, data.trait_count);

 if (!x.marker_score || !x.effect || !x.residual || !x.inclusion ||
     !x.phenotype_sum_squares || !x.marker_scale_prior ||
     !x.residual_scale_prior || !x.marker_variance_initial ||
     !x.residual_variance_initial || !x.global_probability || !x.ld_storage ||
     !x.initial_inclusion || !x.initial_residual || !x.ld_friends ||
     !x.marker_order || !x.convergence_markers)
  throw std::invalid_argument("fixed-prior execution context has a null dependency");
 if (x.marker_order->size()!=data.marker_count)
  throw std::invalid_argument("fixed-prior marker order length mismatch");
 for (const int marker : *x.convergence_markers)
  if (marker < 0 || marker >= x.marker_count)
   throw std::invalid_argument("fixed-prior convergence marker index is out of range");
}

CsrPriorBayesCExecutionResult run_csr_prior_bayesc(
 const CsrPriorBayesCExecutionContext& context
);

} }

#endif
