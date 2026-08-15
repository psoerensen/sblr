#ifndef SBLR_BLR_BED_BAYESR_TYPES_H
#define SBLR_BLR_BED_BAYESR_TYPES_H

#include "blr_scheduled_execution_types.h"
#include "blr_bed_family_types.h"

#include <armadillo>
#include "blr_aggregate_component_trace.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include "blr_phase3_execution.h"
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

template <class PackedGenotype>
using BedBayesRPackedGenotypeView=BedPackedGenotypeView<PackedGenotype>;

struct BedBayesRComponentSpec {
 std::size_t null_component;
 const std::vector<double>& scales;
 const std::vector<double>& initial_probabilities;
 const std::vector<double>& dirichlet_prior;
};

struct BedBayesRSchedulerControl {
 ScheduledSweepControl sweep;
 NullSkipControl skip;
 CandidateControl candidate;
};

struct BedBayesRProgressEvent {
 int iteration=0;
 int total_iterations=0;
 int trait_index=0;
 int chain_index=0;
 double vb=0.0, ve=0.0, vg=0.0, vle=0.0, vld=0.0, vei=0.0;
 double pi_nonnull=0.0, included=0.0;
 std::size_t component_count=0, active_count=0, candidate_count=0;
};

struct BedBayesRChainExecutionResult {
 arma::rowvec bm, dm, component_mean, b, d_as_double;
 arma::rowvec vbs, vgs, ves, vles, vlds;
 arma::mat pip_k;
 double final_vb=0.0, final_vg=0.0, final_ve=0.0;
 std::vector<double> final_pi;
 double final_vle=0.0, final_vld=0.0;
 double log_cpo=std::numeric_limits<double>::quiet_NaN();
 double mean_log_cpo=std::numeric_limits<double>::quiet_NaN();
 std::vector<double> mean_pi;
 arma::mat convergence_pi, convergence_b;
 arma::imat convergence_d, convergence_component;
 AggregateComponentTrace convergence_aggregate;
 double nsamples=0.0, seconds=0.0;
 int failed=0;
 std::string error;
 std::vector<BedBayesRProgressEvent> progress_events;
};

struct BedBayesRAggregationContext {
 std::size_t marker_count=0, trait_count=0, chain_count=0, trace_length=0;
 std::size_t component_count=0, null_component=0;
};

struct BedBayesRExecutionResult {
 std::size_t marker_count=0, trait_count=0, chain_count=0, trace_length=0;
 std::size_t component_count=0;
 arma::mat bm, dm, component_mean, final_effects, final_states;
 arma::mat bm_sd, bm_min, bm_max, dm_sd, dm_min, dm_max;
 arma::mat wy, residual_score;
 arma::mat vbs, vgs, ves, vles, vlds;
 std::vector<arma::mat> component_probability;
 arma::vec final_vb, final_vg, final_ve, final_vle, final_vld;
 arma::mat final_pi, mean_pi;
 arma::vec log_cpo, mean_log_cpo, retained_samples;
 arma::vec seconds_mean, seconds_max;
 int failed_tasks=0;
};

template <class PackedGenotype, class MarkerMap>
struct BedBayesRChainExecutionContext {
 BedBayesRPackedGenotypeView<PackedGenotype> genotype;
 const std::vector<MarkerMap>& marker_maps;
 const std::vector<int>& marker_order;
 const arma::mat& phenotype;
 const std::vector<std::vector<double>>& initial_effects;
 const arma::mat& initial_B;
 const arma::mat& initial_E;
 const arma::mat& ssb_prior;
 const arma::mat& sse_prior;
 BedBayesRComponentSpec components;
 BedBayesRSchedulerControl scheduler;
 double nub, nue, adjE;
 int iterations, burnin, thinning, rebuild_every, progress_every;
 std::uint64_t chain_seed;
 int trait_index, chain_index;
 bool updateB, updateE, updatePi;
 const std::vector<int>& convergence_markers;
 bool convergence_probability, convergence_b, convergence_d;
 bool convergence_component;
 const BlrPhase3ExecutionContract* execution_contract = nullptr;
};

template <class PackedGenotype, class MarkerMap>
inline void validate_bed_bayesr_chain_context(
 const BedBayesRChainExecutionContext<PackedGenotype,MarkerMap>& x
) {
 const auto& g=x.genotype;
 if (g.packed_markers==nullptr || g.marker_count==0 || g.sample_count==0)
  throw std::invalid_argument("BayesR packed genotype view is invalid");
 if (g.bytes_per_marker!=(g.sample_count+3u)/4u || g.stride<g.bytes_per_marker ||
     g.packed_size<g.marker_count*g.stride)
  throw std::invalid_argument("BayesR packed genotype dimensions are inconsistent");
 if (x.marker_maps.size()!=g.marker_count || x.marker_order.size()!=g.marker_count)
  throw std::invalid_argument("BayesR marker metadata dimensions are inconsistent");
 if (x.trait_index<0 || static_cast<arma::uword>(x.trait_index)>=x.phenotype.n_cols ||
     x.phenotype.n_rows!=g.sample_count)
  throw std::invalid_argument("BayesR phenotype dimensions are inconsistent");
 if (x.trait_index>=static_cast<int>(x.initial_effects.size()) ||
     x.initial_effects[static_cast<std::size_t>(x.trait_index)].size()!=g.marker_count)
  throw std::invalid_argument("BayesR initial-effect dimensions are inconsistent");
 const arma::uword nt=x.phenotype.n_cols;
 if (x.initial_B.n_rows!=nt || x.initial_B.n_cols!=nt ||
     x.initial_E.n_rows!=nt || x.initial_E.n_cols!=nt ||
     x.ssb_prior.n_rows!=nt || x.ssb_prior.n_cols!=nt ||
     x.sse_prior.n_rows!=nt || x.sse_prior.n_cols!=nt)
  throw std::invalid_argument("BayesR variance/prior dimensions are inconsistent");
 const std::size_t k=x.components.scales.size();
 if (k<2 || x.components.null_component!=0 ||
     x.components.initial_probabilities.size()!=k || x.components.dirichlet_prior.size()!=k)
  throw std::invalid_argument("BayesR component dimensions or null index are invalid");
 if (x.components.scales[0]!=0.0)
  throw std::invalid_argument("BayesR null component scale must be zero");
 double psum=0.0;
 for (std::size_t i=0;i<k;++i) {
  if (!std::isfinite(x.components.scales[i]) || (i>0 && x.components.scales[i]<=0.0) ||
      !std::isfinite(x.components.initial_probabilities[i]) || x.components.initial_probabilities[i]<0.0 ||
      !std::isfinite(x.components.dirichlet_prior[i]) || x.components.dirichlet_prior[i]<=0.0)
   throw std::invalid_argument("BayesR component values are invalid");
  psum+=x.components.initial_probabilities[i];
 }
 if (!std::isfinite(psum) || psum<=0.0)
  throw std::invalid_argument("BayesR initial probabilities have invalid sum");
 if (x.iterations<=0 || x.burnin<0 || x.thinning<=0 || x.rebuild_every<0 || x.progress_every<0)
  throw std::invalid_argument("BayesR iteration controls are invalid");
 if (x.scheduler.sweep.full_sweep_every<0 || x.scheduler.skip.base_interval<=0 ||
     x.scheduler.skip.maximum_interval<0 ||
     !std::isfinite(x.scheduler.candidate.probability_threshold) ||
     x.scheduler.candidate.probability_threshold<0.0 ||
     x.scheduler.candidate.probability_threshold>1.0 ||
     x.scheduler.candidate.lifetime<0)
  throw std::invalid_argument("BayesR scheduler controls are invalid");
 if (x.scheduler.skip.growth_rule!="probability_adaptive")
  throw std::invalid_argument("BayesR null-skip growth rule is unsupported");
 if (x.chain_index<0 || !std::isfinite(x.nub) || !std::isfinite(x.nue))
  throw std::invalid_argument("BayesR chain or variance controls are invalid");
 for (int marker : x.convergence_markers) if (marker<0 ||
     static_cast<std::size_t>(marker)>=g.marker_count)
  throw std::invalid_argument("BayesR convergence marker index is out of range");
}

} }

#endif
