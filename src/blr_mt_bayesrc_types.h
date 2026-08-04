#ifndef SBLR_BLR_MT_BAYESRC_TYPES_H
#define SBLR_BLR_MT_BAYESRC_TYPES_H

#include <armadillo>
#include "st_bayesrc_annotation_prior.h"
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "blr_mt_bayesr_types.h"
#include "st_bayesrc_annotation_prior.h"

namespace sblr {
namespace mt {

struct MtBayesRCSpec {
 const arma::mat* annotations=nullptr;
 std::vector<double> pattern_prior;
 bool update_alpha=true;
 StBayesRCInterceptPrior intercept_prior{false, arma::vec(), arma::vec()};
 double sigma_alpha_a=2.0;
 double sigma_alpha_b=2.0;
 double pi_floor=1e-12;
 int alpha_update_every=1;
};

inline std::size_t mt_bayesrc_positive_component_count(
 const MtJointStateSpec& joint
) {
 if (joint.component_count < 2)
  throw std::invalid_argument("MT BayesRC requires a null and a positive component");
 return static_cast<std::size_t>(joint.component_count-1);
}

inline std::size_t mt_bayesrc_pattern_count(const MtJointStateSpec& joint) {
 const std::size_t positive=mt_bayesrc_positive_component_count(joint);
 if (joint.patterns.size()<2 || (joint.patterns.size()-1)%positive!=0)
  throw std::invalid_argument("MT BayesRC joint-state ordering is inconsistent");
 return (joint.patterns.size()-1)/positive;
}

inline std::size_t mt_bayesrc_pattern_index(
 std::size_t state, const MtJointStateSpec& joint
) {
 if (state==0 || state>=joint.patterns.size())
  throw std::invalid_argument("MT BayesRC active state index is out of range");
 return (state-1)/mt_bayesrc_positive_component_count(joint);
}

inline void validate_mt_bayesrc_spec(
 const MtBayesRCSpec& spec, const MtJointStateSpec& joint,
 std::size_t marker_count
) {
 if (spec.annotations==nullptr || spec.annotations->n_rows!=marker_count ||
     spec.annotations->n_cols==0 || !spec.annotations->is_finite())
  throw std::invalid_argument("MT BayesRC annotations must be finite m by q data");
 const std::size_t patterns=mt_bayesrc_pattern_count(joint);
 if (spec.pattern_prior.size()!=patterns)
  throw std::invalid_argument("MT BayesRC pattern prior length differs from active patterns");
 for (double value:spec.pattern_prior)
  if (!std::isfinite(value) || value<=0.0)
   throw std::invalid_argument("MT BayesRC pattern prior must be finite and positive");
 if (!std::isfinite(spec.sigma_alpha_a) || spec.sigma_alpha_a<=0.0 ||
     !std::isfinite(spec.sigma_alpha_b) || spec.sigma_alpha_b<=0.0 ||
     !std::isfinite(spec.pi_floor) || spec.pi_floor<=0.0 ||
     spec.pi_floor>=1.0 || spec.alpha_update_every<=0)
  throw std::invalid_argument("MT BayesRC annotation controls are invalid");
}

inline std::vector<double> mt_bayesrc_marker_prior(
 const arma::rowvec& theta, const std::vector<double>& omega,
 const MtJointStateSpec& joint
) {
 if (theta.n_elem!=static_cast<arma::uword>(joint.component_count) ||
     omega.size()!=mt_bayesrc_pattern_count(joint))
  throw std::invalid_argument("MT BayesRC probability dimensions differ");
 std::vector<double> prior(joint.patterns.size(),0.0);
 prior[0]=theta[0];
 for (std::size_t state=1;state<prior.size();++state) {
  const int component=joint.component[state];
  prior[state]=theta[static_cast<arma::uword>(component)]*
   omega[mt_bayesrc_pattern_index(state,joint)];
 }
 double total=0.0;
 for (double value:prior) {
  if (!std::isfinite(value) || value<0.0)
   throw std::runtime_error("MT BayesRC joint prior is nonfinite");
  total+=value;
 }
 if (!std::isfinite(total) || total<=0.0)
  throw std::runtime_error("MT BayesRC joint prior cannot be normalized");
 for (double& value:prior) value/=total;
 return prior;
}

}  // namespace mt
}  // namespace sblr

#endif
