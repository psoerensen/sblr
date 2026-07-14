#ifndef SBLR_BLR_CSR_ANNOTATION_BAYESC_TYPES_H
#define SBLR_BLR_CSR_ANNOTATION_BAYESC_TYPES_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

enum class AnnotationBayesCProbabilityPolicy { fixed_marker, group, learned_logistic };
enum class AnnotationBayesCScalePolicy { fixed_marker, group, learned_exponential };

struct AnnotationBayesCDataView {
 std::size_t marker_count=0, trait_count=0;
 const std::uint64_t* row_ptr=nullptr; std::size_t row_ptr_count=0;
 const double* diagonal=nullptr; std::size_t diagonal_count=0;
 const int* sample_size=nullptr; std::size_t sample_size_count=0;
 bool shared_read_only=true, per_chain_payload=false, storage_outlives_execution=true;
};

struct AnnotationBayesCControls {
 int iterations=0, burnin=0, thinning=1, chains=1, cores=1, seed=0;
 std::vector<int> chain_seeds;
 bool keep_chains=false, update_marker_variance=true, update_residual_variance=true;
 bool update_global_probability=true, update_ld_swap=false;
 double ld_swap_probability=0.05, ld_swap_r2=0.8;
 int ld_swap_max_friends=50, ld_swap_moves=1;
};

struct AnnotationBayesCOutputControl { bool keep_chains=false, diagnostics=true; };

struct FixedPriorBayesCPolicyView {
 const double* marker_probability=nullptr; std::size_t probability_count=0;
 const double* marker_multiplier=nullptr; std::size_t multiplier_count=0;
 bool use_marker_probability=false, use_marker_multiplier=false;
 bool shared_read_only=true, per_chain_payload=false, storage_outlives_execution=true;
 std::string global_probability_interaction="marker_probability_replaces_global_in_marker_draw";
 std::string global_variance_interaction="marker_multiplier_scales_global_marker_variance";
};

struct GroupBayesCPolicyView {
 const int* marker_group=nullptr; std::size_t marker_group_count=0;
 std::size_t group_count=0; std::vector<std::string> group_order;
 const double* initial_probability=nullptr; std::size_t probability_count=0;
 const double* probability_prior_a=nullptr; const double* probability_prior_b=nullptr;
 std::size_t probability_prior_count=0;
 const double* initial_multiplier=nullptr; std::size_t multiplier_count=0;
 bool update_probability=true, update_multiplier=true, normalize_multiplier=true;
 double multiplier_prior_df=4.0, multiplier_prior_scale=1.0;
 bool zero_based_index=true, shared_read_only=true, per_chain_payload=false;
 bool storage_outlives_execution=true;
};

struct LearnedAnnotationBayesCPolicyView {
 std::size_t marker_count=0, annotation_count=0, trait_count=0;
 const double* annotation=nullptr; std::size_t annotation_value_count=0;
 std::vector<std::string> annotation_order;
 std::string layout="column_major_marker_by_annotation";
 bool includes_intercept=false;
 const double* eta_probability_init=nullptr; std::size_t eta_probability_count=0;
 const double* eta_multiplier_init=nullptr; std::size_t eta_multiplier_count=0;
 bool learn_probability=true, learn_multiplier=false;
 double probability_prior_sd=1.0, multiplier_prior_sd=1.0;
 double probability_proposal_sd=0.05, multiplier_proposal_sd=0.05;
 int update_every=10;
 double probability_min=1e-6, probability_max=1.0-1e-6;
 double multiplier_min=1e-3, multiplier_max=1e3;
 std::string probability_link="centered_logit_offset";
 std::string multiplier_link="centered_exponential";
 bool shared_read_only=true, per_chain_payload=false, storage_outlives_execution=true;
};

template<class Policy> struct CsrAnnotationBayesCExecutionInput {
 AnnotationBayesCDataView data;
 std::vector<std::string> marker_order, trait_order;
 Policy policy;
 AnnotationBayesCControls controls;
 AnnotationBayesCOutputControl output;
};

struct AnnotationBayesCSharedResult {
 std::vector<double> marker_mean, marker_pip, variance_marker, variance_genetic,
  variance_residual, vle, vld, task_seconds;
 std::vector<int> failed, thread_used;
 std::vector<std::string> errors;
};
struct FixedPriorBayesCResult { AnnotationBayesCSharedResult shared; std::vector<double> effective_probability, effective_multiplier; };
struct GroupBayesCResult { AnnotationBayesCSharedResult shared; std::vector<double> probability_mean, multiplier_mean, included_mean; std::vector<std::string> group_order; };
struct LearnedAnnotationBayesCResult { AnnotationBayesCSharedResult shared; std::vector<double> eta_probability_mean, eta_multiplier_mean, effective_probability, effective_multiplier, acceptance; };

inline void validate_annotation_bayesc_common(const AnnotationBayesCDataView& d,
 const AnnotationBayesCControls& c) {
 if (!d.marker_count) throw std::invalid_argument("data.marker_count must be positive");
 if (!d.trait_count) throw std::invalid_argument("data.trait_count must be positive");
 if (d.row_ptr_count != d.marker_count+1) throw std::invalid_argument("data.row_ptr_count must equal marker_count + 1");
 if (d.diagonal_count != d.marker_count) throw std::invalid_argument("data.diagonal_count must equal marker_count");
 if (d.sample_size_count != d.trait_count) throw std::invalid_argument("data.sample_size_count must equal trait_count");
 if (!d.shared_read_only || d.per_chain_payload || !d.storage_outlives_execution) throw std::invalid_argument("data must be borrowed, shared read-only, and outlive execution");
 if (c.iterations<=0 || c.burnin<0 || c.thinning<=0) throw std::invalid_argument("controls MCMC values are invalid");
 if (c.chains<=0 || c.cores<=0) throw std::invalid_argument("controls chains and cores must be positive");
 if (!c.chain_seeds.empty() && c.chain_seeds.size()!=static_cast<std::size_t>(c.chains)) throw std::invalid_argument("controls.chain_seeds length must equal chains");
 if (c.ld_swap_probability<0 || c.ld_swap_probability>1 || c.ld_swap_r2<0 || c.ld_swap_r2>1 || c.ld_swap_max_friends<0 || c.ld_swap_moves<0) throw std::invalid_argument("controls LD-swap values are invalid");
}
inline void validate_fixed_prior_policy(const FixedPriorBayesCPolicyView& p, std::size_t m, std::size_t nt) {
 const std::size_t n=m*nt;
 if (p.use_marker_probability && p.probability_count!=n) throw std::invalid_argument("fixed_prior.marker_probability length must equal markers * traits");
 if (p.use_marker_multiplier && p.multiplier_count!=n) throw std::invalid_argument("fixed_prior.marker_multiplier length must equal markers * traits");
 for (std::size_t i=0;i<p.probability_count;++i) if (!std::isfinite(p.marker_probability[i]) || p.marker_probability[i]<=0 || p.marker_probability[i]>=1) throw std::invalid_argument("fixed_prior.marker_probability must be finite and inside (0,1)");
 for (std::size_t i=0;i<p.multiplier_count;++i) if (!std::isfinite(p.marker_multiplier[i]) || p.marker_multiplier[i]<=0) throw std::invalid_argument("fixed_prior.marker_multiplier must be finite and positive");
 if (!p.shared_read_only || p.per_chain_payload || !p.storage_outlives_execution) throw std::invalid_argument("fixed_prior storage must be borrowed immutable and outlive execution");
}
inline void validate_group_policy(const GroupBayesCPolicyView& p, std::size_t m, std::size_t nt) {
 if (!p.group_count || p.marker_group_count!=m || p.group_order.size()!=p.group_count) throw std::invalid_argument("group dimensions are inconsistent");
 if (!p.zero_based_index) throw std::invalid_argument("group.marker_group must be zero-based at the typed boundary");
 std::vector<bool> seen(p.group_count,false); for(std::size_t i=0;i<m;++i){ if(p.marker_group[i]<0 || static_cast<std::size_t>(p.marker_group[i])>=p.group_count) throw std::invalid_argument("group.marker_group index is invalid"); seen[p.marker_group[i]]=true; }
 for(bool x:seen) if(!x) throw std::invalid_argument("group policy does not permit empty groups");
 const std::size_t n=p.group_count*nt;
 if(p.probability_count!=n || p.multiplier_count!=n || p.probability_prior_count!=p.group_count) throw std::invalid_argument("group value dimensions are inconsistent");
 for(std::size_t i=0;i<n;++i){ if(!std::isfinite(p.initial_probability[i])||p.initial_probability[i]<=0||p.initial_probability[i]>=1) throw std::invalid_argument("group.initial_probability must be inside (0,1)"); if(!std::isfinite(p.initial_multiplier[i])||p.initial_multiplier[i]<=0) throw std::invalid_argument("group.initial_multiplier must be positive"); }
 if(!p.shared_read_only||p.per_chain_payload||!p.storage_outlives_execution) throw std::invalid_argument("group storage must be borrowed immutable and outlive execution");
}
inline void validate_learned_annotation_policy(const LearnedAnnotationBayesCPolicyView& p, std::size_t m, std::size_t nt) {
 if(p.marker_count!=m||p.trait_count!=nt||!p.annotation_count||p.annotation_value_count!=m*p.annotation_count) throw std::invalid_argument("annotation dimensions are inconsistent");
 if(p.annotation_order.size()!=p.annotation_count) throw std::invalid_argument("annotation.annotation_order length is invalid");
 if(p.layout!="column_major_marker_by_annotation") throw std::invalid_argument("annotation.layout is unsupported");
 for(std::size_t i=0;i<p.annotation_value_count;++i) if(!std::isfinite(p.annotation[i])) throw std::invalid_argument("annotation values must be finite");
 const std::size_t n=p.annotation_count*nt;
 if(p.eta_probability_count!=n||p.eta_multiplier_count!=n) throw std::invalid_argument("annotation coefficient dimensions are inconsistent");
 if(p.probability_prior_sd<=0||p.multiplier_prior_sd<=0||p.probability_proposal_sd<0||p.multiplier_proposal_sd<0||p.update_every<=0) throw std::invalid_argument("annotation prior/proposal/update controls are invalid");
 if(!(p.probability_min>0&&p.probability_min<p.probability_max&&p.probability_max<1)) throw std::invalid_argument("annotation probability bounds are invalid");
 if(!(p.multiplier_min>0&&p.multiplier_min<p.multiplier_max)) throw std::invalid_argument("annotation multiplier bounds are invalid");
 if(p.probability_link!="centered_logit_offset"||p.multiplier_link!="centered_exponential") throw std::invalid_argument("annotation link tag is incompatible with production semantics");
 if(!p.shared_read_only||p.per_chain_payload||!p.storage_outlives_execution) throw std::invalid_argument("annotation storage must be borrowed immutable and outlive execution");
}

} }
#endif
