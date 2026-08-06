#ifndef SBLR_BLR_BED_BAYESRC_TYPES_H
#define SBLR_BLR_BED_BAYESRC_TYPES_H

#include "blr_bed_family_types.h"

#include <armadillo>
#include "blr_aggregate_component_trace.h"
#include "st_bayesrc_annotation_prior.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

struct BedBayesRCComponentSpec {
 const std::vector<double>& scales;
 std::size_t null_component;
 std::size_t residual_component;
};

template <class AnnotationMatrix>
struct BedBayesRCAnnotationSpec {
 const AnnotationMatrix& matrix;
 std::size_t marker_count;
 std::size_t annotation_count;
 std::size_t intercept_column;
};

struct BedBayesRCCoefficientPriorSpec {
 const arma::mat& initial_alpha;
 const arma::vec& initial_step_variances;
 StBayesRCInterceptPrior intercept_prior;
 double inverse_chisq_df;
 double inverse_chisq_scale;
 bool update_coefficients;
 int update_every;
};

template <class PackedGenotype>
struct BedBayesRCPackedGenotypeView {
 const PackedGenotype& storage;
 std::size_t marker_count;
 std::size_t sample_count;
 std::size_t bytes_per_marker;
};

struct BedBayesRCChainExecutionResult {
 arma::rowvec bm, dm, component_mean, b, state;
 arma::rowvec vbs, vgs, ves, vles, vlds, pis;
 arma::mat comp_prob;
 arma::mat annot_alpha_mean, annot_alpha_final;
 arma::vec annot_sigma_mean, annot_sigma_final;
 arma::mat convergence_alpha, convergence_sigma, convergence_b;
 arma::imat convergence_d, convergence_component;
 AggregateComponentTrace convergence_aggregate;
 arma::imat coupling_replica_identity;
 arma::mat coupling_active_count, coupling_expected_active, coupling_swap;
 double coupling_transition_seconds=0.0, coupling_swap_seconds=0.0;
 bool coupling_tempering=false;
 arma::rowvec mean_prior;
 arma::vec residual;
 double final_vb=0.0, final_vg=0.0, final_ve=0.0;
 double log_cpo=std::numeric_limits<double>::quiet_NaN();
 double mean_log_cpo=std::numeric_limits<double>::quiet_NaN();
 double nsamples=0.0;
 int failed=0;
 std::string error;
};

struct BedBayesRCAggregationContext {
 const arma::mat& annotation;
 const std::vector<double>& component_scales;
 std::size_t marker_count;
 std::size_t annotation_count;
 std::size_t component_count;
 std::size_t trait_count;
 std::size_t chain_count;
 std::size_t trace_length;
 double pi_floor;
 bool keep_chains;
};

struct BedBayesRCExecutionResult {
 arma::mat bm, dm, b, state, component_mean;
 arma::mat vbs, vgs, ves, vle, vld, pis;
 arma::mat final_prior, mean_prior;
 std::vector<arma::mat> comp_prob, marker_prior_final;
 std::vector<arma::mat> alpha_mean, alpha_final;
 arma::mat sigma_mean, sigma_final;
 arma::vec final_vb, final_vg, final_ve;
 arma::vec log_cpo, mean_log_cpo, nsamples;
 arma::mat component_counts;
 arma::mat wy, residual_marker_score;
 std::vector<BedBayesRCChainExecutionResult> retained_chains;
 std::vector<std::string> failures;
 std::size_t marker_count=0, annotation_count=0, component_count=0;
 std::size_t trait_count=0, chain_count=0, trace_length=0;
};

template <class PackedGenotype, class AnnotationMatrix, class MarkerMap>
struct BedBayesRCChainExecutionContext {
 BedBayesRCPackedGenotypeView<PackedGenotype> genotype;
 BedBayesRCAnnotationSpec<AnnotationMatrix> annotation;
 BedBayesRCComponentSpec components;
 BedBayesRCCoefficientPriorSpec coefficient_prior;
 const std::vector<MarkerMap>& marker_maps;
 const std::vector<int>& marker_order;
 const arma::mat& phenotype;
 const std::vector<std::vector<double>>& initial_effects;
 const arma::mat& initial_B;
 const arma::mat& initial_E;
 const arma::mat& ssb_prior;
 const arma::mat& sse_prior;
 double pi_floor, nub, nue, adjE;
 bool update_marker_variance, update_residual_variance;
 int iterations, burnin, thinning, rebuild_every;
 std::uint64_t chain_seed;
 int trait_index, chain_index;
 const std::vector<int>& convergence_markers;
 bool convergence_annotations, convergence_b, convergence_d,
  convergence_component;
};

template <class PackedGenotype, class AnnotationMatrix, class MarkerMap>
inline void validate_bed_bayesrc_chain_context(
 const BedBayesRCChainExecutionContext<PackedGenotype,AnnotationMatrix,MarkerMap>& x
) {
 const auto& g=x.genotype;
 if (g.marker_count==0 || g.sample_count==0 ||
     g.marker_count!=static_cast<std::size_t>(g.storage.m) ||
     g.sample_count!=static_cast<std::size_t>(g.storage.n) || g.bytes_per_marker==0)
  throw std::invalid_argument("BayesRC packed genotype view is invalid");
 if (x.marker_maps.size()!=g.marker_count || x.marker_order.size()!=g.marker_count)
  throw std::invalid_argument("BayesRC marker metadata dimensions are inconsistent");
 if (x.phenotype.n_rows!=g.sample_count || x.trait_index<0 ||
     static_cast<arma::uword>(x.trait_index)>=x.phenotype.n_cols)
  throw std::invalid_argument("BayesRC phenotype dimensions are inconsistent");
 const auto& a=x.annotation;
 if (a.marker_count!=g.marker_count || a.annotation_count==0 ||
     a.matrix.n_rows!=a.marker_count || a.matrix.n_cols!=a.annotation_count ||
     a.intercept_column!=0)
  throw std::invalid_argument("BayesRC annotation dimensions or intercept are invalid");
 const std::size_t k=x.components.scales.size();
 if (k<2 || x.components.null_component!=0 ||
     x.components.residual_component+1!=k || x.components.scales[0]!=0.0)
  throw std::invalid_argument("BayesRC component dimensions or ordering are invalid");
 for (std::size_t i=1;i<k;++i)
  if (!std::isfinite(x.components.scales[i]) || x.components.scales[i]<=0.0)
   throw std::invalid_argument("BayesRC active component scales must be positive finite");
 const auto& p=x.coefficient_prior;
 if (p.initial_alpha.n_rows!=a.annotation_count || p.initial_alpha.n_cols!=k-1 ||
     p.initial_step_variances.n_elem!=k-1 || !p.initial_alpha.is_finite() ||
     !p.initial_step_variances.is_finite() || arma::any(p.initial_step_variances<=0.0))
  throw std::invalid_argument("BayesRC coefficient-prior dimensions or values are invalid");
 if (!std::isfinite(p.inverse_chisq_df) || p.inverse_chisq_df<=0.0 ||
     !std::isfinite(p.inverse_chisq_scale) || p.inverse_chisq_scale<=0.0 ||
     p.update_every<=0)
  throw std::invalid_argument("BayesRC coefficient-prior controls are invalid");
 if (x.trait_index>=static_cast<int>(x.initial_effects.size()) ||
     x.initial_effects[static_cast<std::size_t>(x.trait_index)].size()!=g.marker_count)
  throw std::invalid_argument("BayesRC initial-effect dimensions are inconsistent");
 const arma::uword nt=x.phenotype.n_cols;
 if (x.initial_B.n_rows!=nt || x.initial_B.n_cols!=nt ||
     x.initial_E.n_rows!=nt || x.initial_E.n_cols!=nt ||
     x.ssb_prior.n_rows!=nt || x.ssb_prior.n_cols!=nt ||
     x.sse_prior.n_rows!=nt || x.sse_prior.n_cols!=nt)
  throw std::invalid_argument("BayesRC variance/prior dimensions are inconsistent");
 if (x.iterations<=0 || x.burnin<0 || x.thinning<=0 || x.rebuild_every<0 ||
     x.chain_index<0 || !std::isfinite(x.pi_floor) || x.pi_floor<=0.0 || x.pi_floor>=1.0 ||
     !std::isfinite(x.nub) || !std::isfinite(x.nue) || !std::isfinite(x.adjE))
  throw std::invalid_argument("BayesRC execution controls are invalid");
 for (int marker: x.convergence_markers)
  if (marker<0 || static_cast<std::size_t>(marker)>=g.marker_count)
   throw std::invalid_argument("BayesRC convergence marker index is out of range");
}

} }

#endif
