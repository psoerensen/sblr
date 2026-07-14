#ifndef SBLR_CORE_BLR_CSR_SBAYESRC_TYPES_H
#define SBLR_CORE_BLR_CSR_SBAYESRC_TYPES_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

struct CsrSBayesRCDataView {
  std::size_t marker_count=0, trait_count=0, nonzero_count=0;
  const std::uint64_t* row_ptr=nullptr; std::size_t row_ptr_count=0;
  const std::uint32_t* column_index=nullptr; const float* value=nullptr;
  const double* diagonal=nullptr; std::size_t diagonal_count=0;
  const int* sample_size=nullptr; std::size_t sample_size_count=0;
  bool shared_read_only=true, per_chain_payload=false, storage_outlives_execution=true;
};

// Borrowed column-major marker-by-annotation design prepared by the R binding.
struct SBayesRCAnnotationDesignView {
  std::size_t marker_count=0, annotation_count=0;
  const double* values=nullptr; std::size_t value_count=0;
  std::vector<std::string> annotation_order;
  std::string layout="column_major";
  bool includes_intercept=true, standardized=true, centered_binary=false;
  bool shared_read_only=true, per_chain_payload=false, storage_outlives_execution=true;
};

struct SBayesRCComponentSpec {
  std::vector<double> scales;
  std::size_t null_component=0;
  std::string scale_interpretation="variance_multiplier";
};

struct SBayesRCAlphaSpec {
  std::size_t annotation_count=0, step_count=0;
  std::vector<double> initial_values;       // column-major annotation x step
  std::vector<double> initial_variance;     // one value per step
  bool intercept_flat=true, update=true;
  double variance_prior_a=2.0, variance_prior_b=2.0;
  int update_every=10;
};

struct SBayesRCProbabilityPolicy {
  std::string transformation="probit_stick_breaking";
  std::size_t reference_component=0;
  std::vector<std::size_t> stick_order;
  double probability_floor=1e-12;
  bool floor_then_normalize=true;
};

struct SBayesRCPriors { double marker_df=4.0, residual_df=4.0; };

struct SBayesRCControls {
  int iterations=0, burnin=0, thinning=1, chains=1, cores=1, seed=1;
  std::vector<int> chain_seeds;
  bool keep_chains=false, update_marker_variance=true, update_residual_variance=true;
  bool rebuild_residual_before_update=false, update_ld_swap=false;
  double ld_swap_probability=.05, ld_swap_r2=.8;
  int ld_swap_max_friends=50, ld_swap_moves=1;
};

struct SBayesRCOutputControl { bool keep_chains=false, diagnostics=true; };

struct CsrSBayesRCExecutionInput {
  CsrSBayesRCDataView data;
  SBayesRCAnnotationDesignView annotation;
  SBayesRCComponentSpec components;
  SBayesRCAlphaSpec alpha;
  SBayesRCProbabilityPolicy probability;
  SBayesRCPriors priors;
  SBayesRCControls controls;
  SBayesRCOutputControl output;
  std::vector<std::string> marker_order, trait_order;
};

struct SBayesRCArray { std::vector<double> values; std::vector<std::size_t> dimensions; };
struct SBayesRCChainResult {
  SBayesRCArray marker_mean, marker_pip, component_probability, component_count;
  SBayesRCArray alpha, alpha_variance, variance_trace, diagnostics;
  bool failed=false; std::string failure_message; double elapsed_seconds=0.0;
};
struct CsrSBayesRCResult {
  std::size_t marker_count=0, trait_count=0, annotation_count=0, component_dimension=0;
  SBayesRCArray marker_mean, marker_pip, component_probability, component_assignment;
  SBayesRCArray component_count, component_probability_trace, final_component_probability;
  SBayesRCArray alpha_trace, alpha_mean, alpha_variance, alpha_diagnostics;
  SBayesRCArray marker_variance_trace, genetic_variance_trace, residual_variance_trace;
  SBayesRCArray genic_variance_trace, ld_variance_trace, chain_timing, diagnostics;
  std::vector<SBayesRCChainResult> chains;
  std::vector<std::string> marker_order, trait_order, annotation_order, component_names;
  bool chains_retained=false;
};

inline void validate_csr_sbayesrc_execution_input(const CsrSBayesRCExecutionInput& x) {
  const std::size_t m=x.data.marker_count, nt=x.data.trait_count;
  if(m==0) throw std::invalid_argument("data.marker_count must be positive.");
  if(nt==0) throw std::invalid_argument("data.trait_count must be positive.");
  if(!x.data.shared_read_only||x.data.per_chain_payload||!x.data.storage_outlives_execution)
    throw std::invalid_argument("CSR data must be borrowed shared read-only storage that outlives execution, without per-chain payload.");
  if(!x.data.row_ptr||x.data.row_ptr_count!=m+1) throw std::invalid_argument("data.row_ptr length mismatch.");
  if(x.data.nonzero_count && (!x.data.column_index||!x.data.value)) throw std::invalid_argument("CSR indices and values are required.");
  if(!x.data.diagonal||x.data.diagonal_count!=m) throw std::invalid_argument("data.diagonal length mismatch.");
  if(!x.data.sample_size||x.data.sample_size_count!=nt) throw std::invalid_argument("data.sample_size length mismatch.");
  if(x.marker_order.size()!=m||x.trait_order.size()!=nt) throw std::invalid_argument("marker or trait order length mismatch.");
  const auto& a=x.annotation;
  if(a.marker_count!=m||a.annotation_count==0||!a.values||a.value_count!=m*a.annotation_count)
    throw std::invalid_argument("annotation dimensions must be marker_count by positive annotation_count.");
  if(a.annotation_order.size()!=a.annotation_count) throw std::invalid_argument("annotation_order length mismatch.");
  if(a.layout!="column_major") throw std::invalid_argument("annotation layout must be column_major.");
  if(!a.shared_read_only||a.per_chain_payload||!a.storage_outlives_execution)
    throw std::invalid_argument("annotation data must be borrowed shared read-only storage that outlives execution, without per-chain payload.");
  for(std::size_t i=0;i<a.value_count;++i) if(!std::isfinite(a.values[i])) throw std::invalid_argument("annotation values must be finite.");
  const std::size_t k=x.components.scales.size();
  if(k<2||x.components.null_component!=0||x.components.scales[0]!=0.0) throw std::invalid_argument("components require null component zero at index 0.");
  for(std::size_t i=0;i<k;++i) if(!std::isfinite(x.components.scales[i])||(i&&x.components.scales[i]<=0.0)) throw std::invalid_argument("component scales are invalid.");
  if(x.alpha.annotation_count!=a.annotation_count||x.alpha.step_count!=k-1||x.alpha.initial_values.size()!=a.annotation_count*(k-1))
    throw std::invalid_argument("alpha dimensions must be annotation_count by component_count - 1.");
  if(x.alpha.initial_variance.size()!=k-1) throw std::invalid_argument("alpha variance length mismatch.");
  for(double v:x.alpha.initial_values) if(!std::isfinite(v)) throw std::invalid_argument("alpha initial values must be finite.");
  for(double v:x.alpha.initial_variance) if(!std::isfinite(v)||v<=0) throw std::invalid_argument("alpha variances must be positive finite.");
  if(!std::isfinite(x.alpha.variance_prior_a)||x.alpha.variance_prior_a<=0||!std::isfinite(x.alpha.variance_prior_b)||x.alpha.variance_prior_b<=0)
    throw std::invalid_argument("alpha variance priors must be positive finite.");
  if(x.alpha.update&&x.alpha.update_every<=0) throw std::invalid_argument("alpha.update_every must be positive when alpha is learned.");
  if(x.probability.transformation!="probit_stick_breaking"||x.probability.reference_component!=0||x.probability.stick_order.size()!=k-1)
    throw std::invalid_argument("probability policy must use the current ordered probit stick contract.");
  for(std::size_t j=0;j<k-1;++j) if(x.probability.stick_order[j]!=j+1) throw std::invalid_argument("probability stick_order is invalid.");
  if(!std::isfinite(x.probability.probability_floor)||x.probability.probability_floor<=0||x.probability.probability_floor>=1)
    throw std::invalid_argument("probability_floor must be in (0, 1).");
  if(x.controls.iterations<=0||x.controls.burnin<0||x.controls.thinning<=0) throw std::invalid_argument("invalid MCMC controls.");
  if(x.controls.chains<=0||x.controls.cores<=0) throw std::invalid_argument("chain and core counts must be positive.");
  if(!x.controls.chain_seeds.empty()&&x.controls.chain_seeds.size()!=static_cast<std::size_t>(x.controls.chains)) throw std::invalid_argument("chain_seeds length mismatch.");
  if(!std::isfinite(x.controls.ld_swap_probability)||x.controls.ld_swap_probability<0||x.controls.ld_swap_probability>1||
     !std::isfinite(x.controls.ld_swap_r2)||x.controls.ld_swap_r2<0||x.controls.ld_swap_r2>1||x.controls.ld_swap_max_friends<=0||x.controls.ld_swap_moves<0)
    throw std::invalid_argument("invalid LD-swap controls.");
  if(x.output.keep_chains!=x.controls.keep_chains) throw std::invalid_argument("output.keep_chains must match controls.keep_chains.");
}

} }
#endif
