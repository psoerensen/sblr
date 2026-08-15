#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_TYPES_H
#define SBLR_BLR_BED_SCHEDULED_BAYESC_TYPES_H

#include "blr_scheduled_execution_types.h"
#include "blr_bed_family_types.h"

#include <armadillo>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include "blr_phase3_execution.h"
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

struct BedScheduledBayesCMarkerMap {
 double val[4];
 double xx;
};

struct BedScheduledBayesCChainExecutionResult {
 arma::rowvec bm;
 arma::rowvec dm;
 arma::rowvec b;
 arma::rowvec d_as_double;
 arma::rowvec vbs;
 arma::rowvec vgs;
 arma::rowvec ves;
 arma::rowvec pis;
 arma::rowvec vles;
 arma::rowvec vlds;
 arma::mat convergence_b;
 arma::imat convergence_d;
 double final_vb=0.0;
 double final_vg=0.0;
 double final_ve=0.0;
 double final_pi=0.0;
 double final_vle=0.0;
 double final_vld=0.0;
 double log_cpo=std::numeric_limits<double>::quiet_NaN();
 double mean_log_cpo=std::numeric_limits<double>::quiet_NaN();
 double mean_pi=0.0;
 double nsamples=0.0;
 double seconds=0.0;
 int failed=0;
 std::string error;
};

struct BedScheduledBayesCExecutionResult {
 arma::mat bm;
 arma::mat dm;
 arma::mat bm_sd;
 arma::mat dm_sd;
 arma::mat bm_min;
 arma::mat dm_min;
 arma::mat bm_max;
 arma::mat dm_max;
 arma::mat final_effects;
 arma::mat final_states;
 arma::mat wy;
 arma::mat residual_scores;
 arma::mat marker_variance_trace;
 arma::mat genetic_variance_trace;
 arma::mat residual_variance_trace;
 arma::mat inclusion_trace;
 arma::mat vle_trace;
 arma::mat vld_trace;
 arma::vec final_marker_variance;
 arma::vec final_genetic_variance;
 arma::vec final_residual_variance;
 arma::vec final_inclusion_probability;
 arma::vec final_vle;
 arma::vec final_vld;
 arma::vec mean_inclusion_probability;
 arma::vec mean_total_log_cpo;
 arma::vec mean_log_cpo;
 arma::vec mean_retained_samples;
 arma::vec mean_seconds;
 arma::vec max_seconds;
 int marker_count=0;
 int sample_count=0;
 int trait_count=0;
 int chain_count=0;
 int task_count=0;
 int trace_length=0;
};

template <class PackedGenotype>
struct BedScheduledBayesCAggregationContext {
 BedPackedGenotypeView<PackedGenotype> genotype;
 const std::vector<BedScheduledBayesCMarkerMap>& marker_maps;
 const std::vector<int>& marker_order;
 const arma::mat& phenotype;
 int trait_count;
 int chain_count;
 int trace_length;
 bool return_wy;
 bool return_residual_scores;
};

template <class PackedGenotype>
struct BedScheduledBayesCChainExecutionContext {
 BedPackedGenotypeView<PackedGenotype> genotype;
 const std::vector<BedScheduledBayesCMarkerMap>& marker_maps;
 const std::vector<int>& marker_order;
 const arma::mat& phenotype;
 const std::vector<std::vector<double>>& initial_effects;
 const arma::mat& initial_B;
 const arma::mat& initial_E;
 const arma::mat& ssb_prior;
 const arma::mat& sse_prior;
 const std::vector<double>& initial_pi;
 ScheduledSweepControl sweep;
 NullSkipControl skip;
 CandidateControl candidate;
 double nub;
 double nue;
 double adjE;
 double pi_prior_a;
 double pi_prior_b;
 int iterations;
 int burnin;
 int thinning;
 int rebuild_every;
 int progress_every;
 std::uint64_t chain_seed;
 int trait_index;
 int chain_index;
 bool updateB;
 bool updateE;
 bool updatePi;
 const std::vector<int>& convergence_markers;
 bool convergence_b;
 bool convergence_d;
 const BlrPhase3ExecutionContract* execution_contract = nullptr;
};

template <class PackedGenotype>
inline void validate_bed_scheduled_bayesc_chain_context(
 const BedScheduledBayesCChainExecutionContext<PackedGenotype>& x
) {
 const auto& g=x.genotype;
 if (g.packed_markers==nullptr) throw std::invalid_argument("BED chain packed genotype storage is null");
 if (g.marker_count==0) throw std::invalid_argument("BED chain marker_count must be positive");
 if (g.sample_count==0) throw std::invalid_argument("BED chain sample_count must be positive");
 if (g.bytes_per_marker!=(g.sample_count+3u)/4u)
  throw std::invalid_argument("BED chain bytes_per_marker is inconsistent with sample_count");
 if (g.stride<g.bytes_per_marker)
  throw std::invalid_argument("BED chain genotype stride is smaller than bytes_per_marker");
 if (g.packed_size<g.marker_count*g.stride)
  throw std::invalid_argument("BED chain packed genotype size is inconsistent with dimensions");
 if (x.marker_maps.size()!=g.marker_count)
  throw std::invalid_argument("BED chain marker-map dimension mismatch");
 if (x.marker_order.size()!=g.marker_count)
  throw std::invalid_argument("BED chain marker-order dimension mismatch");
 if (x.trait_index<0 || static_cast<arma::uword>(x.trait_index)>=x.phenotype.n_cols ||
     x.phenotype.n_rows!=g.sample_count)
  throw std::invalid_argument("BED chain phenotype dimensions are inconsistent");
 if (x.trait_index>=static_cast<int>(x.initial_effects.size()) ||
     x.initial_effects[static_cast<std::size_t>(x.trait_index)].size()!=g.marker_count)
  throw std::invalid_argument("BED chain initial-effect dimensions are inconsistent");
 const arma::uword nt=x.phenotype.n_cols;
 if (x.initial_B.n_rows!=nt || x.initial_B.n_cols!=nt ||
     x.initial_E.n_rows!=nt || x.initial_E.n_cols!=nt ||
     x.ssb_prior.n_rows!=nt || x.ssb_prior.n_cols!=nt ||
     x.sse_prior.n_rows!=nt || x.sse_prior.n_cols!=nt)
  throw std::invalid_argument("BED chain prior dimensions are inconsistent");
 if (x.initial_pi.size()!=2)
  throw std::invalid_argument("BED chain initial_pi must have length two");
 if (x.iterations<=0) throw std::invalid_argument("BED chain iterations must be positive");
 if (x.burnin<0) throw std::invalid_argument("BED chain burnin must be non-negative");
 if (x.thinning<=0) throw std::invalid_argument("BED chain thinning must be positive");
 if (x.rebuild_every<0) throw std::invalid_argument("BED chain rebuild_every must be non-negative");
 if (x.progress_every<0) throw std::invalid_argument("BED chain progress_every must be non-negative");
 if (x.sweep.full_sweep_every<0)
  throw std::invalid_argument("BED chain full_sweep_every must be non-negative");
 if (x.skip.base_interval<=0) throw std::invalid_argument("BED chain null_skip_base must be positive");
 if (x.skip.maximum_interval<0) throw std::invalid_argument("BED chain null_skip_max must be non-negative");
 if (x.skip.growth_rule!="probability_adaptive")
  throw std::invalid_argument("BED chain null-skip growth rule is unsupported");
 if (!std::isfinite(x.candidate.probability_threshold) ||
     x.candidate.probability_threshold<0.0 || x.candidate.probability_threshold>1.0)
  throw std::invalid_argument("BED chain candidate_threshold must be in [0,1]");
 if (x.candidate.lifetime<0)
  throw std::invalid_argument("BED chain candidate_lifetime must be non-negative");
 if (x.chain_index<0) throw std::invalid_argument("BED chain index must be non-negative");
 if (!std::isfinite(x.nub) || !std::isfinite(x.nue))
  throw std::invalid_argument("BED chain variance degrees of freedom must be finite");
 if (!std::isfinite(x.pi_prior_a) || x.pi_prior_a<=0.0 ||
     !std::isfinite(x.pi_prior_b) || x.pi_prior_b<=0.0)
  throw std::invalid_argument("BED chain pi prior parameters must be finite and positive");
 for (int marker : x.convergence_markers) if (marker<0 ||
     static_cast<std::size_t>(marker)>=g.marker_count)
  throw std::invalid_argument("BED chain convergence marker index is out of range");
}

} }

#endif
