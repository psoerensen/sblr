#ifndef SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_TYPES_H
#define SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_TYPES_H

#include <armadillo>

#include "blr_csr_annotation_bayesc_types.h"

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace sblr { namespace core {

// Prepared storage is owned by the native adapter and must outlive the call.
// Shared CSR, annotation, ordering, priors, and initial coefficients are
// borrowed. Trait-local sampler, coefficient, RNG, and accumulator state is
// created and owned inside run_csr_learned_annotation_bayesc().
struct CsrLearnedAnnotationBayesCExecutionContext {
 int marker_count=0, trait_count=0, annotation_count=0;
 int iterations=0, burnin=0, thinning=1, cores=1, seed=0;
 bool use_initial_inclusion=false, use_initial_residual=false;
 bool rebuild_residual_before_update=false;
 bool learn_probability=true, learn_multiplier=false;
 bool update_marker_variance=true, update_residual_variance=true;
 bool update_global_probability=true, update_ld_swap=false;
 int annotation_update_every=1, ld_swap_moves=0;
 double probability_prior_sd=1.0, multiplier_prior_sd=1.0;
 double probability_proposal_sd=0.0, multiplier_proposal_sd=0.0;
 double probability_min=0.0, probability_max=1.0;
 double multiplier_min=0.0, multiplier_max=1.0;
 double marker_degrees_freedom=0.0, residual_degrees_freedom=0.0;
 double residual_adjustment=0.0, global_probability_prior_a=0.0;
 double global_probability_prior_b=0.0, ld_swap_probability=0.0;

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
 const arma::mat* annotation=nullptr;
 const arma::mat* initial_probability_coefficient=nullptr;
 const arma::mat* initial_multiplier_coefficient=nullptr;
 const std::vector<int>* marker_order=nullptr;
 const void* ld_storage=nullptr;
 std::size_t ld_row_ptr_count=0;
 const void* ld_friends_storage=nullptr;
 LearnedAnnotationBayesCPolicyView annotation_policy;
};

struct CsrLearnedAnnotationBayesCExecutionResult {
 std::vector<std::vector<std::vector<double>>> raw;
};

inline void validate_csr_learned_annotation_bayesc_execution_context(
 const CsrLearnedAnnotationBayesCExecutionContext& x
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
     !x.sample_size || !x.initial_inclusion || !x.initial_residual ||
     !x.annotation || !x.initial_probability_coefficient ||
     !x.initial_multiplier_coefficient || !x.marker_order ||
     !x.ld_storage || !x.ld_friends_storage)
  throw std::invalid_argument("learned-annotation execution context has a null dependency");
 if (x.marker_order->size()!=data.marker_count)
  throw std::invalid_argument("learned-annotation marker order length mismatch");
 if (x.annotation->n_rows!=data.marker_count ||
     x.annotation->n_cols!=static_cast<arma::uword>(x.annotation_count))
  throw std::invalid_argument("learned-annotation matrix dimensions mismatch");
 validate_learned_annotation_policy(x.annotation_policy, data.marker_count, data.trait_count);
}

CsrLearnedAnnotationBayesCExecutionResult
run_csr_learned_annotation_bayesc(
 const CsrLearnedAnnotationBayesCExecutionContext& context
);

} }

#endif
