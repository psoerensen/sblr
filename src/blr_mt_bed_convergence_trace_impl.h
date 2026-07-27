#ifndef SBLR_BLR_MT_BED_CONVERGENCE_TRACE_IMPL_H
#define SBLR_BLR_MT_BED_CONVERGENCE_TRACE_IMPL_H

#include "blr_mt_bed_chains_types.h"
#include "blr_mt_bed_convergence_types.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace mt {

inline MtBedConvergenceTraceBundle build_mt_bed_convergence_trace_bundle(
    const std::vector<MtBedChainExecutionResult>& results,
    int nburn,
    int nit,
    bool updateB,
    bool updateE) {
 if (results.empty() || nburn<0 || nit<=0) {
  throw std::invalid_argument("invalid MT BED convergence trace request");
 }
 const MtDefaultCoreResult& reference=results[0].core.posterior;
 const int nt=reference.nt;
 if (nt<=0) {
  throw std::invalid_argument("MT BED convergence traces require traits");
 }
 const std::size_t expected=static_cast<std::size_t>(nburn+nit);
 MtBedConvergenceTraceBundle bundle;
 bundle.nchains=static_cast<int>(results.size());
 bundle.postburn_draws=nit;
 bundle.quantities.reserve(static_cast<std::size_t>(5*nt));
 for (const char* group : {"vbs", "vgs", "ves", "vle", "vld"}) {
  for (int trait=0; trait<nt; ++trait) {
   const bool updated=std::string(group)=="vbs" ? updateB :
    (std::string(group)=="ves" ? updateE : true);
   bundle.quantities.push_back(MtBedConvergenceQuantity{
    group, trait, updated});
  }
 }
 bundle.values.reserve(bundle.quantities.size()*results.size()*
                       static_cast<std::size_t>(nit));
 for (const MtBedConvergenceQuantity& quantity : bundle.quantities) {
  for (std::size_t chain=0; chain<results.size(); ++chain) {
   const MtBedChainExecutionResult& result=results[chain];
   if (result.failed || result.chain!=static_cast<int>(chain) ||
       result.core.posterior.nt!=nt) {
    throw std::invalid_argument(
     "MT BED convergence traces require successful ordered chains");
   }
   const std::vector<std::vector<double>>* traces=nullptr;
   if (quantity.group=="vbs") traces=&result.core.posterior.vbs;
   else if (quantity.group=="vgs") traces=&result.core.posterior.vgs;
   else if (quantity.group=="ves") traces=&result.core.posterior.ves;
   else if (quantity.group=="vle") traces=&result.core.posterior.vle;
   else traces=&result.core.posterior.vld;
   if (traces->size()!=static_cast<std::size_t>(nt) ||
       (*traces)[static_cast<std::size_t>(quantity.trait)].size()!=expected) {
    throw std::invalid_argument(
     "inconsistent MT BED convergence trace dimensions");
   }
   const std::vector<double>& trace=
    (*traces)[static_cast<std::size_t>(quantity.trait)];
   for (int iteration=0; iteration<nit; ++iteration) {
    const double value=trace[static_cast<std::size_t>(nburn+iteration)];
    if (quantity.updated && !std::isfinite(value)) {
     throw std::invalid_argument(
      "updated MT BED convergence traces must be finite");
    }
    bundle.values.push_back(value);
   }
  }
 }
 return bundle;
}

}  // namespace mt
}  // namespace sblr

#endif
