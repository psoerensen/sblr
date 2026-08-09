#ifndef SBLR_BLR_CSR_LOGVAR_BAYESR_TYPES_H
#define SBLR_BLR_CSR_LOGVAR_BAYESR_TYPES_H

#include <armadillo>

#include <vector>

#include "st_logvar_annotation_prior.h"

namespace sblr {
namespace logvar {

struct CsrLogvarBayesRChainOutput {
 arma::vec theta;
 arma::rowvec prior_scale;
 arma::mat theta_trace;
 arma::vec theta_sum;
 arma::rowvec prior_scale_sum;
 double retained_samples = 0.0;
 EssUpdateDiagnostics diagnostics;
};

struct CsrLogvarBayesRPolicyInput {
 const arma::mat* annotation = nullptr;
 const arma::mat* theta_initial = nullptr;
 double theta_prior_sd = kThetaPriorSdV1;
 bool update_theta = true;
 int trace_count = 0;
 std::vector<CsrLogvarBayesRChainOutput>* outputs = nullptr;
};

}  // namespace logvar
}  // namespace sblr

#endif
