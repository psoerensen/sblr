#ifndef SBLR_BLR_NORMAL_PROBABILITY_H
#define SBLR_BLR_NORMAL_PROBABILITY_H

// Binding-neutral interface with the established Rmath-backed implementation.
namespace sblr { namespace core {

struct StandardNormalProbability {
 static inline double cdf(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
 }

 static inline double quantile(double p) {
  return R::qnorm(p, 0.0, 1.0, 1, 0);
 }
};

} }

#endif
