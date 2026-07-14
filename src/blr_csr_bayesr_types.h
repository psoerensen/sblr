#ifndef SBLR_CORE_BLR_CSR_BAYESR_TYPES_H
#define SBLR_CORE_BLR_CSR_BAYESR_TYPES_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

// Borrowed immutable ordinary-CSR storage. The binding owns these buffers and
// must keep them alive until all chain tasks return. No chain owns CSR data.
struct CsrBayesRDataView {
  std::size_t marker_count = 0;
  std::size_t trait_count = 0;
  const std::uint64_t* row_ptr = nullptr;
  std::size_t row_ptr_count = 0;
  const std::uint32_t* column_index = nullptr;
  const float* value = nullptr;
  std::size_t nonzero_count = 0;
  const double* diagonal = nullptr;
  std::size_t diagonal_count = 0;
  const int* sample_size = nullptr;
  std::size_t sample_size_count = 0;
  bool shared_read_only = true;
  bool per_chain_payload = false;
  bool storage_outlives_execution = true;
};

struct BayesRComponentSpec {
  std::vector<double> scales;
  std::vector<double> initial_probability;
  std::vector<double> dirichlet_prior;
  std::size_t null_component = 0;
  bool update_probability = true;
  std::string scale_interpretation = "variance_multiplier";
};

struct CsrBayesRPriors {
  double marker_degrees_freedom = 0.0;
  double residual_degrees_freedom = 0.0;
};

struct CsrBayesRControls {
  int iterations = 0;
  int burnin = 0;
  int thinning = 1;
  int chains = 1;
  int cores = 1;
  int seed = 1;
  std::vector<int> chain_seeds;
  bool keep_chains = false;
  bool update_marker_variance = true;
  bool update_residual_variance = true;
  bool rebuild_residual_before_update = false;
  int residual_update_start = 0;
  int residual_update_every = 1;
  bool update_ld_swap = false;
  double ld_swap_probability = 0.05;
  double ld_swap_r2 = 0.8;
  int ld_swap_max_friends = 50;
  int ld_swap_moves = 1;
};

struct CsrBayesROutputControl { bool keep_chains = false; };

struct CsrBayesRExecutionInput {
  CsrBayesRDataView data;
  BayesRComponentSpec components;
  CsrBayesRPriors priors;
  CsrBayesRControls controls;
  CsrBayesROutputControl output;
  std::vector<std::string> marker_order;
  std::vector<std::string> trait_order;
};

struct CsrBayesRArray {
  std::vector<double> values;
  std::vector<std::size_t> dimensions;
};

struct CsrBayesRChainResult {
  CsrBayesRArray marker_mean, marker_pip, component_probability;
  CsrBayesRArray component_mean, final_effect, final_component;
  CsrBayesRArray marker_variance_trace, genetic_variance_trace;
  CsrBayesRArray residual_variance_trace, genic_variance_trace;
  CsrBayesRArray ld_variance_trace, component_probability_trace;
  CsrBayesRArray diagnostics;
  bool failed = false;
  std::string failure_message;
  double elapsed_seconds = 0.0;
  double retained_samples = 0.0;
};

// Stable binding-neutral result vocabulary for CSR BayesR contracts.
struct CsrBayesRResult {
  std::size_t marker_count = 0, trait_count = 0, component_dimension = 0;
  CsrBayesRArray marker_mean, marker_pip, component_probability;
  CsrBayesRArray component_mean, final_effect, final_residual, final_component;
  CsrBayesRArray component_count, final_component_probability;
  CsrBayesRArray mean_component_probability, component_probability_trace;
  CsrBayesRArray marker_variance_trace, genetic_variance_trace;
  CsrBayesRArray residual_variance_trace, genic_variance_trace, ld_variance_trace;
  CsrBayesRArray chain_timing, diagnostics;
  std::vector<CsrBayesRChainResult> chains;
  std::vector<std::string> marker_order, trait_order, component_names;
  bool chains_retained = false;
};

inline void validate_csr_bayesr_execution_input(const CsrBayesRExecutionInput& x) {
  if (x.data.marker_count == 0) throw std::invalid_argument("data.marker_count must be positive.");
  if (x.data.trait_count == 0) throw std::invalid_argument("data.trait_count must be positive.");
  if (!x.data.shared_read_only) throw std::invalid_argument("data must be borrowed shared read-only storage.");
  if (x.data.per_chain_payload) throw std::invalid_argument("data must not contain per-chain CSR payload.");
  if (!x.data.storage_outlives_execution) throw std::invalid_argument("data storage must outlive execution.");
  if (x.data.row_ptr == nullptr || x.data.row_ptr_count != x.data.marker_count + 1)
    throw std::invalid_argument("data.row_ptr must contain marker_count + 1 entries.");
  if (x.data.nonzero_count > 0 && (x.data.column_index == nullptr || x.data.value == nullptr))
    throw std::invalid_argument("data CSR indices and values are required for nonzero storage.");
  if (x.data.diagonal == nullptr || x.data.diagonal_count != x.data.marker_count)
    throw std::invalid_argument("data.diagonal must match marker_count.");
  if (x.data.sample_size == nullptr || x.data.sample_size_count != x.data.trait_count)
    throw std::invalid_argument("data.sample_size must match trait_count.");
  if (x.marker_order.size() != x.data.marker_count) throw std::invalid_argument("marker_order length mismatch.");
  if (x.trait_order.size() != x.data.trait_count) throw std::invalid_argument("trait_order length mismatch.");
  const std::size_t k = x.components.scales.size();
  if (k < 2) throw std::invalid_argument("components must contain at least two scales.");
  if (x.components.null_component >= k) throw std::invalid_argument("components.null_component is out of range.");
  if (x.components.scales[x.components.null_component] != 0.0)
    throw std::invalid_argument("the null component scale must equal zero.");
  for (std::size_t i=0;i<k;++i) {
    const double s=x.components.scales[i];
    if (!std::isfinite(s) || (i != x.components.null_component && s <= 0.0))
      throw std::invalid_argument("component scales must be finite, with positive non-null scales.");
  }
  if (x.components.initial_probability.size()!=k || x.components.dirichlet_prior.size()!=k)
    throw std::invalid_argument("component probability and prior dimensions must match scales.");
  double psum=0.0;
  for(double p:x.components.initial_probability){if(!std::isfinite(p)||p<0.0)throw std::invalid_argument("component probabilities must be finite and non-negative."); psum+=p;}
  if(!std::isfinite(psum)||psum<=0.0)throw std::invalid_argument("component probabilities must have positive sum.");
  for(double a:x.components.dirichlet_prior)if(!std::isfinite(a)||a<=0.0)throw std::invalid_argument("Dirichlet priors must be finite and positive.");
  if (x.controls.iterations<=0 || x.controls.burnin<0 || x.controls.thinning<=0)
    throw std::invalid_argument("invalid iterations, burnin, or thinning.");
  if (x.controls.chains<=0) throw std::invalid_argument("controls.chains must be positive.");
  if (x.controls.cores<=0) throw std::invalid_argument("controls.cores must be positive.");
  if (!x.controls.chain_seeds.empty() && x.controls.chain_seeds.size()!=static_cast<std::size_t>(x.controls.chains))
    throw std::invalid_argument("controls.chain_seeds must have length controls.chains.");
  if (x.controls.residual_update_start<0 || x.controls.residual_update_every<=0)
    throw std::invalid_argument("invalid residual update timing.");
  if (!std::isfinite(x.controls.ld_swap_probability)||x.controls.ld_swap_probability<0.0||x.controls.ld_swap_probability>1.0)
    throw std::invalid_argument("ld_swap_probability must be in [0, 1].");
  if (!std::isfinite(x.controls.ld_swap_r2)||x.controls.ld_swap_r2<0.0||x.controls.ld_swap_r2>1.0)
    throw std::invalid_argument("ld_swap_r2 must be in [0, 1].");
  if (x.controls.ld_swap_max_friends<=0 || x.controls.ld_swap_moves<0)
    throw std::invalid_argument("invalid LD-swap counts.");
  if (x.output.keep_chains != x.controls.keep_chains)
    throw std::invalid_argument("output.keep_chains must match controls.keep_chains.");
}

} }
#endif
