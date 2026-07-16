#ifndef SBLR_BLR_CSR_GROUP_BAYESC_TYPES_H
#define SBLR_BLR_CSR_GROUP_BAYESC_TYPES_H

#include <armadillo>

#include "blr_csr_annotation_bayesc_types.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sblr { namespace core {

// Prepared storage is owned by the native adapter and must outlive the call.
// The core borrows all shared inputs; trait-local sampler, group, RNG, and
// accumulator state is created and owned inside run_csr_group_bayesc().
struct CsrGroupBayesCExecutionContext {
 int marker_count=0, trait_count=0, group_count=0;
 int iterations=0, burnin=0, thinning=1, cores=1, seed=0;
 bool use_initial_inclusion=false, use_initial_residual=false;
 bool rebuild_residual_before_update=false;
 bool update_group_multiplier=true, normalize_group_multiplier=true;
 bool update_marker_variance=true, update_residual_variance=true;
 bool update_group_probability=true, update_ld_swap=false;
 double group_multiplier_prior_df=0.0, group_multiplier_prior_scale=0.0;
 double marker_degrees_freedom=0.0, residual_degrees_freedom=0.0;
 double residual_adjustment=0.0, ld_swap_probability=0.0;
 int ld_swap_moves=0;

 arma::mat* marker_score=nullptr;
 const arma::mat* marker_diagonal=nullptr;
 arma::mat* effect=nullptr;
 arma::mat* residual=nullptr;
 arma::Mat<int>* inclusion=nullptr;
 const arma::vec* phenotype_sum_squares=nullptr;
 const arma::mat* marker_scale_prior=nullptr;
 const arma::mat* residual_scale_prior=nullptr;
 const arma::mat* marker_variance_initial=nullptr;
 const arma::mat* residual_variance_initial=nullptr;
 const std::vector<double>* global_probability=nullptr;
 const std::vector<int>* sample_size=nullptr;
 const std::vector<std::vector<double>>* initial_inclusion=nullptr;
 const std::vector<std::vector<double>>* initial_residual=nullptr;
 const std::vector<std::vector<double>>* initial_group_probability=nullptr;
 const std::vector<std::vector<double>>* initial_group_multiplier=nullptr;
 const arma::rowvec* probability_prior_a=nullptr;
 const arma::rowvec* probability_prior_b=nullptr;
 const arma::Row<int>* marker_group=nullptr;
 const arma::rowvec* group_size=nullptr;
 const std::vector<int>* marker_order=nullptr;
 const void* ld_storage=nullptr;
 std::size_t ld_row_ptr_count=0;
 const void* ld_friends_storage=nullptr;
 GroupBayesCPolicyView group_policy;
};

struct CsrGroupBayesCExecutionResult {
 std::vector<std::vector<std::vector<double>>> raw;
};

inline void validate_csr_group_bayesc_execution_context(
 const CsrGroupBayesCExecutionContext& x
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

 if (!x.marker_score || !x.effect || !x.residual || !x.inclusion ||
     !x.phenotype_sum_squares || !x.marker_scale_prior ||
     !x.residual_scale_prior || !x.marker_variance_initial ||
     !x.residual_variance_initial || !x.global_probability ||
     !x.initial_inclusion || !x.initial_residual ||
     !x.initial_group_probability || !x.initial_group_multiplier ||
     !x.probability_prior_a || !x.probability_prior_b ||
     !x.marker_group || !x.group_size || !x.marker_order ||
     !x.ld_storage || !x.ld_friends_storage)
  throw std::invalid_argument("group execution context has a null dependency");
 if (x.marker_order->size()!=data.marker_count)
  throw std::invalid_argument("group marker order length mismatch");
 if (x.initial_group_probability->size()!=data.trait_count ||
     x.initial_group_multiplier->size()!=data.trait_count)
  throw std::invalid_argument("group initial state trait count mismatch");

 for (std::size_t t=0; t<data.trait_count; ++t) {
  GroupBayesCPolicyView policy=x.group_policy;
  policy.initial_probability=(*x.initial_group_probability)[t].data();
  policy.probability_count=policy.group_count;
  policy.initial_multiplier=(*x.initial_group_multiplier)[t].data();
  policy.multiplier_count=policy.group_count;
  validate_group_policy(policy, data.marker_count, 1);
 }
}

CsrGroupBayesCExecutionResult run_csr_group_bayesc(
 const CsrGroupBayesCExecutionContext& context
);

} }

#endif
