#ifndef SBLR_BLR_MT_BAYESR_KERNEL_IMPL_H
#define SBLR_BLR_MT_BAYESR_KERNEL_IMPL_H

#include "blr_mt_bayesr_types.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace sblr {
namespace mt {

inline void validate_mt_joint_state_spec(
 const MtJointStateSpec& spec, std::size_t nt
) {
 const std::size_t states=spec.patterns.size();
 if (states<2 || states>4096 || spec.component.size()!=states ||
     spec.multiplier.size()!=states || spec.names.size()!=states ||
     spec.component_count<2) {
  throw std::invalid_argument("invalid MT BayesR joint-state specification");
 }
 for (std::size_t state=0;state<states;++state) {
  if (spec.patterns[state].size()!=nt)
   throw std::invalid_argument("joint-state pattern width must equal nt");
  bool any=false;
  for (int value:spec.patterns[state]) {
   if (value!=0 && value!=1)
    throw std::invalid_argument("joint-state patterns must be binary");
   any=any || value==1;
  }
  if (state==0) {
   if (any || spec.component[state]!=0 || spec.multiplier[state]!=0.0 ||
       spec.names[state]!="null")
    throw std::invalid_argument("joint state 1 must be the unique null state");
  } else if (!any || spec.component[state]<=0 ||
             spec.component[state]>=spec.component_count ||
             !std::isfinite(spec.multiplier[state]) ||
             spec.multiplier[state]<=0.0) {
   throw std::invalid_argument("non-null joint states require positive components");
  }
 }
 const int positive_count=spec.component_count-1;
 if (states<2 || (states-1)%static_cast<std::size_t>(positive_count)!=0)
  throw std::invalid_argument("joint-state count must equal 1 + patterns * components");
 std::vector<double> component_multiplier(
  static_cast<std::size_t>(spec.component_count),0.0);
 for (std::size_t state=1;state<states;++state) {
  const int expected=1+static_cast<int>((state-1)%positive_count);
  if (spec.component[state]!=expected)
   throw std::invalid_argument("joint states must use ascending components within pattern");
  if (state<=static_cast<std::size_t>(positive_count))
   component_multiplier[expected]=spec.multiplier[state];
  else if (spec.multiplier[state]!=component_multiplier[expected])
   throw std::invalid_argument("component multipliers must be shared across patterns");
  if (expected>1 && component_multiplier[expected]<=component_multiplier[expected-1])
   throw std::invalid_argument("positive component multipliers must be strictly ordered");
  if (expected>1 && spec.patterns[state]!=spec.patterns[state-1])
   throw std::invalid_argument("trait pattern must remain fixed within a component block");
 }
}

inline arma::mat mt_joint_cholesky(arma::mat precision) {
 arma::mat lower;
 if (arma::chol(lower,precision,"lower")) return lower;
 double jitter=1e-10;
 for (int attempt=0;attempt<8;++attempt) {
  precision.diag()+=jitter;
  if (arma::chol(lower,precision,"lower")) return lower;
  jitter*=10.0;
 }
 throw std::runtime_error("MT BayesR marker precision is not positive definite");
}

inline MtJointMarkerKernelResult mt_joint_marker_kernel(
 const arma::vec& score,
 const arma::vec& diagonal,
 const arma::mat& B_inverse,
 const arma::mat& E_inverse,
 const MtJointStateSpec& spec,
 const std::vector<double>& pi,
 double marker_scale
) {
 const std::size_t nt=score.n_elem;
 validate_mt_joint_state_spec(spec,nt);
 if (diagonal.n_elem!=nt || B_inverse.n_rows!=nt ||
     B_inverse.n_cols!=nt || E_inverse.n_rows!=nt ||
     E_inverse.n_cols!=nt || pi.size()!=spec.patterns.size() ||
     !score.is_finite() || !diagonal.is_finite() ||
     arma::any(diagonal<=0.0) || !B_inverse.is_finite() ||
     !E_inverse.is_finite() || !std::isfinite(marker_scale) ||
     marker_scale<=0.0)
  throw std::invalid_argument("invalid MT BayesR marker-kernel input");

 MtJointMarkerKernelResult out;
 const std::size_t states=spec.patterns.size();
 out.states.resize(states); out.probability.assign(states,0.0);
 std::vector<double> log_weight(states,-std::numeric_limits<double>::infinity());
 if (std::isfinite(pi[0]) && pi[0]>0.0) log_weight[0]=std::log(pi[0]);
 for (std::size_t state=1;state<states;++state) {
  if (!std::isfinite(pi[state]) || pi[state]<=0.0) continue;
  const double scale=spec.multiplier[state]*marker_scale;
  arma::mat D(nt,nt,arma::fill::zeros);
  for (std::size_t trait=0;trait<nt;++trait)
   D(trait,trait)=static_cast<double>(spec.patterns[state][trait]);
  arma::mat prior_precision=B_inverse/scale;
  const arma::mat likelihood_precision=D*E_inverse*D;
  auto& candidate=out.states[state];
  candidate.precision=prior_precision;
  candidate.rhs=D*E_inverse*score;
  for (std::size_t trait=0;trait<nt;++trait)
   candidate.precision(trait,trait)+=
    diagonal[trait]*likelihood_precision(trait,trait);
  // Preserve full residual-precision cross terms when trait diagonals agree.
  const double common=diagonal[0];
  bool same=true;
  for (std::size_t trait=1;trait<nt;++trait) same=same && diagonal[trait]==common;
  if (same) candidate.precision=prior_precision+common*likelihood_precision;
  candidate.lower=mt_joint_cholesky(candidate.precision);
  candidate.mean=arma::solve(arma::trimatu(candidate.lower.t()),
    arma::solve(arma::trimatl(candidate.lower),candidate.rhs));
  arma::mat prior_lower=mt_joint_cholesky(prior_precision);
  const double log_prior_det_half=arma::sum(arma::log(prior_lower.diag()));
  const double log_post_det_half=arma::sum(arma::log(candidate.lower.diag()));
  candidate.log_weight=std::log(pi[state])+log_prior_det_half-
   log_post_det_half+0.5*arma::dot(candidate.rhs,candidate.mean);
  log_weight[state]=candidate.log_weight;
 }
 const double maximum=*std::max_element(log_weight.begin(),log_weight.end());
 if (!std::isfinite(maximum))
  throw std::runtime_error("MT BayesR joint probabilities cannot be normalized");
 double total=0.0;
 for (std::size_t state=0;state<states;++state) {
  out.probability[state]=std::isfinite(log_weight[state]) ?
   std::exp(log_weight[state]-maximum) : 0.0;
  total+=out.probability[state];
 }
 if (!std::isfinite(total) || total<=0.0)
  throw std::runtime_error("MT BayesR joint probabilities cannot be normalized");
 for (double& value:out.probability) value/=total;
 return out;
}

inline std::size_t mt_joint_draw_state(
 const std::vector<double>& probability, std::mt19937& rng
) {
 std::uniform_real_distribution<double> uniform(0.0,1.0);
 const double draw=uniform(rng); double cumulative=0.0;
 for (std::size_t state=0;state<probability.size();++state) {
  cumulative+=probability[state];
  if (draw<cumulative || state+1==probability.size()) return state;
 }
 throw std::runtime_error("failed to draw MT BayesR joint state");
}

inline arma::vec mt_joint_draw_beta(
 const MtJointMarkerKernelState& state, std::size_t nt, std::mt19937& rng
) {
 arma::vec z(nt); std::normal_distribution<double> normal(0.0,1.0);
 for (std::size_t trait=0;trait<nt;++trait) z[trait]=normal(rng);
 return state.mean+arma::solve(arma::trimatu(state.lower.t()),z);
}

inline std::vector<std::vector<double>> mt_bayesr_base_effects(
 const std::vector<std::vector<double>>& effects,
 const std::vector<int>& component,
 const MtJointStateSpec& joint,
 const std::vector<double>& marker_scale
) {
 std::vector<std::vector<double>> out=effects;
 const std::size_t nt=effects.size();
 const std::size_t m=nt==0 ? 0 : effects[0].size();
 for (std::size_t marker=0;marker<m;++marker) {
  const int comp=component[marker];
  double multiplier=0.0;
  for (std::size_t state=1;state<joint.component.size();++state)
   if (joint.component[state]==comp) { multiplier=joint.multiplier[state]; break; }
  const double scale=comp>0 ? multiplier*marker_scale[marker] : 1.0;
  const double root=std::sqrt(scale);
  for (std::size_t trait=0;trait<nt;++trait)
   out[trait][marker]=comp>0 ? effects[trait][marker]/root : 0.0;
 }
 return out;
}

}  // namespace mt
}  // namespace sblr

#endif
