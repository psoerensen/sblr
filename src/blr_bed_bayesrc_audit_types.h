#ifndef SBLR_BLR_BED_BAYESRC_AUDIT_TYPES_H
#define SBLR_BLR_BED_BAYESRC_AUDIT_TYPES_H

// Phase 14A binding-neutral vocabulary only. It does not invoke the sampler.
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sblr { namespace audit {

enum class BedBayesRCAnnotationOrientation { marker_by_annotation };

struct BedBayesRCComponentSpec {
 std::size_t component_count = 0;
 std::size_t null_component = 0;
 std::size_t stick_count = 0;
 std::vector<double> scales;
 void validate() const {
  if (component_count < 2) throw std::invalid_argument("component_count must be at least two");
  if (null_component != 0) throw std::invalid_argument("null_component must be zero");
  if (stick_count + 1 != component_count) throw std::invalid_argument("stick_count must equal component_count minus one");
  if (scales.size() != component_count) throw std::invalid_argument("scale vector length mismatch");
  if (scales[0] != 0.0) throw std::invalid_argument("null component scale must be zero");
  for (std::size_t k=1;k<scales.size();++k)
   if (!std::isfinite(scales[k]) || scales[k] <= 0.0) throw std::invalid_argument("active scales must be positive finite");
 }
};

struct BedBayesRCAnnotationSpec {
 std::size_t marker_count = 0;
 std::size_t annotation_count = 0;
 std::size_t coefficient_rows = 0;
 std::size_t coefficient_columns = 0;
 bool intercept_present = false;
 bool immutable_fit_owned = true;
 BedBayesRCAnnotationOrientation orientation = BedBayesRCAnnotationOrientation::marker_by_annotation;
 void validate(const BedBayesRCComponentSpec& component) const {
  if (!marker_count || !annotation_count) throw std::invalid_argument("annotation dimensions must be positive");
  if (coefficient_rows != annotation_count || coefficient_columns != component.stick_count)
   throw std::invalid_argument("annotation coefficient dimensions mismatch");
  if (!immutable_fit_owned) throw std::invalid_argument("annotation storage must be immutable and fit-owned");
 }
};

struct BedBayesRCCoefficientPriorSpec {
 std::vector<double> initial_variances;
 double inverse_chisq_df = 2.0;
 double inverse_chisq_scale = 2.0;
 bool intercept_flat = true;
 bool update_coefficients = true;
 std::size_t update_every = 1;
 void validate(std::size_t sticks) const {
  if (initial_variances.size()!=sticks) throw std::invalid_argument("coefficient-prior dimension mismatch");
  for (double x:initial_variances) if (!std::isfinite(x)||x<=0.0) throw std::invalid_argument("coefficient variances must be positive finite");
  if (!std::isfinite(inverse_chisq_df)||inverse_chisq_df<=0.0||!std::isfinite(inverse_chisq_scale)||inverse_chisq_scale<=0.0)
   throw std::invalid_argument("coefficient hyperparameters must be positive finite");
  if (!update_every) throw std::invalid_argument("coefficient update interval must be positive");
 }
};

struct BedBayesRCExecutionAuditSpec {
 std::size_t iterations=0,burnin=0,thinning=0,chains=0,cores=0;
 unsigned int seed=0;
 bool full_sweep_only=true;
 void validate() const {
  if (!iterations||!thinning||!chains||!cores) throw std::invalid_argument("invalid execution controls");
  if (!full_sweep_only) throw std::invalid_argument("packed-BED BayesRC is full-sweep only");
 }
};

struct BedBayesRCOwnershipAuditSpec {
 bool genotype_fit_owned=true,annotation_fit_owned=true,chain_state_private=true;
 bool rng_chain_owned=true,no_per_chain_genotype_copy=true,no_mcmc_disk_access=true;
 void validate() const {
  if (!(genotype_fit_owned&&annotation_fit_owned&&chain_state_private&&rng_chain_owned&&no_per_chain_genotype_copy&&no_mcmc_disk_access))
   throw std::invalid_argument("invalid BayesRC ownership contract");
 }
};

struct BedBayesRCChainResultVocabulary {
 bool marker_posterior=true,marker_component_probability=true,annotation_coefficients=true;
 bool variance_traces=true,cpo=true,retained_count=true,failure=true;
};
struct BedBayesRCExecutionResultVocabulary {
 bool cross_chain_marker_summary=true,cross_chain_annotation_summary=true;
 bool component_and_probability_summary=true,diagnostics=true;
};

} }
#endif
