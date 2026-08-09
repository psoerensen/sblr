#ifndef SBLR_BLR_CSR_BAYESC_RCPP_ADAPTER_H
#define SBLR_BLR_CSR_BAYESC_RCPP_ADAPTER_H

#include <RcppArmadillo.h>

#include "blr_csr_bayesc_types.h"

struct CsrBayesCRawConversionContext {
 int marker_count;
 int trait_count;
 int nit;
 int nburn;
 int nthin;
 int ncores;
 int nchains;
 bool keep_chains;
 double pi_prior_a;
 double pi_prior_b;
 bool update_ld_swap;
 bool use_fixed_maf_effect_scale;
 bool estimate_maf_effect_s;
 const std::vector<int>* sample_size;
 const std::vector<int>* convergence_markers;
};

Rcpp::List stblr_csr_bayesc_result_to_raw(
  const sblr::core::CsrBayesCResult& result,
  const CsrBayesCRawConversionContext& context
);

#endif
