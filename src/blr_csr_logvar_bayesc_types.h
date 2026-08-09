#ifndef SBLR_BLR_CSR_LOGVAR_BAYESC_TYPES_H
#define SBLR_BLR_CSR_LOGVAR_BAYESC_TYPES_H

#include <armadillo>

#include <cstddef>
#include <vector>

#include "st_logvar_annotation_prior.h"

namespace sblr {
namespace logvar {

struct CsrLogvarBayesCChainOutput {
 arma::vec theta;
 arma::rowvec prior_scale;
 arma::mat theta_trace;
 arma::vec theta_sum;
 arma::rowvec prior_scale_sum;
 double retained_samples = 0.0;
 EssUpdateDiagnostics diagnostics;
};

struct CsrLogvarBayesCPolicyInput {
 const arma::mat* annotation = nullptr;
 const arma::mat* theta_initial = nullptr;  // annotations x traits
 double theta_prior_sd = kThetaPriorSdV1;
 bool update_theta = true;
 int trace_count = 0;
 int burnin = 0;
 int thinning = 1;
 std::vector<CsrLogvarBayesCChainOutput>* outputs = nullptr;
};

}  // namespace logvar
}  // namespace sblr

#endif
