#ifndef SBLR_BLR_MT_DEFAULT_FINALIZE_IMPL_H
#define SBLR_BLR_MT_DEFAULT_FINALIZE_IMPL_H

#include "blr_mt_default_types.h"

namespace sblr {
namespace mt {

inline MtDefaultFinalResult finalize_mt_default_result(
 MtDefaultCoreResult core_result
) {
 MtDefaultFinalResult result;
 result.nt=core_result.nt;
 result.m=core_result.m;
 result.nmodels=core_result.nmodels;
 result.marker_retained_count=core_result.marker_retained_count;
 result.covb_retained_count=core_result.covb_retained_count;
 result.covg_retained_count=core_result.covg_retained_count;
 result.cove_retained_count=core_result.cove_retained_count;
 result.pi_retained_count=core_result.pi_retained_count;

 result.bm=std::move(core_result.bm);
 result.dm=std::move(core_result.dm);
 for (int t=0; t < result.nt; t++) {
  for (int i=0; i < result.m; i++) {
   result.bm[t][i] = result.bm[t][i]/result.marker_retained_count;
   result.dm[t][i] = result.dm[t][i]/result.marker_retained_count;
  }
 }

 result.r=std::move(core_result.r);
 result.b=std::move(core_result.b);
 result.d=std::move(core_result.d);
 result.component=std::move(core_result.component);
 result.component_probabilities=std::move(core_result.component_counts);
 for (auto& marker:result.component_probabilities)
  for (double& value:marker) value/=result.marker_retained_count;
 result.marker_order=std::move(core_result.order);
 result.vbs=std::move(core_result.vbs);
 result.vgs=std::move(core_result.vgs);
 result.ves=std::move(core_result.ves);
 result.vle=std::move(core_result.vle);
 result.vld=std::move(core_result.vld);

 result.covb=std::move(core_result.cvbm);
 result.covg=std::move(core_result.cvgm);
 result.cove=std::move(core_result.cvem);
 for (int t1=0; t1 < result.nt; t1++) {
  for (int t2=0; t2 < result.nt; t2++) {
   result.covb[t1][t2] = result.covb_retained_count > 0.0 ?
    result.covb[t1][t2] / result.covb_retained_count : 0.0;
   result.covg[t1][t2] = result.covg_retained_count > 0.0 ?
    result.covg[t1][t2] / result.covg_retained_count : 0.0;
   result.cove[t1][t2] = result.cove_retained_count > 0.0 ?
    result.cove[t1][t2] / result.cove_retained_count : 0.0;
  }
 }

 result.vb=std::move(core_result.B);
 result.vg=std::move(core_result.G);
 result.ve=std::move(core_result.E);
 result.pi_final=std::move(core_result.pi);
 result.pi_mean.assign(static_cast<std::size_t>(result.nmodels), 0.0);
 for (int i=0; i < result.nmodels; i++) {
  result.pi_mean[i] = result.pi_retained_count > 0.0 ?
   core_result.pis[i] / result.pi_retained_count : 0.0;
 }
 result.pi_trace=std::move(core_result.pi_trace);

 // These legacy fields are unsupported by the authoritative method and their
 // core accumulators remain zero. Moving them preserves their zero semantics.
 result.pitrait=std::move(core_result.pistrait);
 result.pimarker=std::move(core_result.pismarker);
 result.annotation_alpha_final=std::move(core_result.annotation_alpha_final);
 result.annotation_alpha_mean=std::move(core_result.annotation_alpha_mean);
 result.annotation_sigma_final=std::move(core_result.annotation_sigma_final);
 result.annotation_sigma_mean=std::move(core_result.annotation_sigma_mean);
 result.pattern_pi_final=std::move(core_result.pattern_pi_final);
 result.pattern_pi_mean=std::move(core_result.pattern_pi_mean);
 result.pattern_pi_trace=std::move(core_result.pattern_pi_trace);
 result.prior_component_probabilities=
  std::move(core_result.prior_component_probabilities);
 result.annotation_updates_attempted=core_result.annotation_updates_attempted;
 result.annotation_updates_completed=core_result.annotation_updates_completed;
 result.convergence=std::move(core_result.convergence);
 return result;
}

}  // namespace mt
}  // namespace sblr

#endif
