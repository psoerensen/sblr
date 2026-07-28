#ifndef SBLR_BLR_MT_BED_CONVERGENCE_TRACE_IMPL_H
#define SBLR_BLR_MT_BED_CONVERGENCE_TRACE_IMPL_H

#include "blr_mt_bed_chains_types.h"
#include "blr_mt_bed_convergence_types.h"

#include <algorithm>
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
    bool updateE,
    bool updatePi,
    bool updateAlpha,
    int method,
    const std::string& residual_covariance,
    const MtExtendedTraceSpec* extended_spec=nullptr) {
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
   MtBedConvergenceQuantity quantity;
   quantity.group=group; quantity.trait=trait; quantity.updated=updated;
   quantity.derived=std::string(group)=="vgs" ||
    std::string(group)=="vle" || std::string(group)=="vld";
   bundle.quantities.push_back(quantity);
  }
 }
 if (extended_spec!=nullptr) {
  if (extended_spec->covariance) {
   for (const char* group:{"cov_b","cov_g","cov_e"})
    for (int col=0;col<nt;++col) for (int row=col+1;row<nt;++row) {
     MtBedConvergenceQuantity q; q.group=group; q.trait=col; q.trait2=row;
     q.tier=2; q.derived=std::string(group)=="cov_g";
     q.structural=std::string(group)=="cov_e" &&
      residual_covariance=="diagonal";
     q.updated=std::string(group)=="cov_b" ? updateB :
      (std::string(group)=="cov_e" ? updateE : true);
     bundle.quantities.push_back(q);
    }
  }
  if (extended_spec->probability) {
   const auto& ext=reference.convergence;
   for (std::size_t k=0;k<ext.component_pi.size();++k) {
    MtBedConvergenceQuantity q; q.group="component_pi"; q.component=k;
    q.tier=2; q.updated=method==6 ? updateAlpha : updatePi;
    q.derived=method==5; bundle.quantities.push_back(q);
   }
   for (std::size_t p=0;p<ext.pattern_pi.size();++p) {
    MtBedConvergenceQuantity q; q.group="pattern_pi"; q.pattern=p;
    q.tier=2; q.updated=updatePi; q.derived=method==5;
    bundle.quantities.push_back(q);
   }
   for (std::size_t state=0;state<ext.joint_pi.size();++state) {
    MtBedConvergenceQuantity q; q.group="joint_pi"; q.pattern=state;
    q.tier=2; q.updated=updatePi; bundle.quantities.push_back(q);
   }
  }
  if (extended_spec->annotations && method==6) {
   const auto& ext=reference.convergence;
   for (std::size_t qindex=0;qindex<ext.annotation_alpha.size();++qindex) {
    MtBedConvergenceQuantity q; q.group="alpha"; q.annotation=qindex;
    q.tier=2; q.updated=updateAlpha; bundle.quantities.push_back(q);
   }
   for (std::size_t qindex=0;qindex<ext.annotation_sigma.size();++qindex) {
    MtBedConvergenceQuantity q; q.group="sigmaSqAlpha"; q.stick=qindex;
    q.tier=2; q.updated=updateAlpha; bundle.quantities.push_back(q);
   }
  }
  if (extended_spec->selected_b)
   for (int marker:extended_spec->selected_markers)
    for (int trait=0;trait<nt;++trait) {
    MtBedConvergenceQuantity q; q.group="b"; q.marker=marker; q.trait=trait;
    q.tier=3; q.updated=true; bundle.quantities.push_back(q);
   }
  if (extended_spec->selected_d)
   for (int marker:extended_spec->selected_markers)
    for (int trait=0;trait<nt;++trait) {
    MtBedConvergenceQuantity q; q.group="d"; q.marker=marker; q.trait=trait;
    q.tier=3; q.updated=true; q.derived=true; bundle.quantities.push_back(q);
   }
  if (extended_spec->selected_component && (method==5 || method==6))
   for (int marker:extended_spec->selected_markers) {
    MtBedConvergenceQuantity q; q.group="component"; q.marker=marker;
    q.tier=3; q.updated=true; bundle.quantities.push_back(q);
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
   const std::vector<std::vector<int>>* integer_traces=nullptr;
   if (quantity.group=="vbs") traces=&result.core.posterior.vbs;
   else if (quantity.group=="vgs") traces=&result.core.posterior.vgs;
   else if (quantity.group=="ves") traces=&result.core.posterior.ves;
   else if (quantity.group=="vle") traces=&result.core.posterior.vle;
   else if (quantity.group=="vld") traces=&result.core.posterior.vld;
   else {
    const auto& ext=result.core.posterior.convergence;
    int pair=-1;
    if (quantity.trait2>=0) {
     pair=0;
     for (int col=0;col<nt;++col) for (int row=col+1;row<nt;++row,++pair)
      if (col==quantity.trait && row==quantity.trait2) goto pair_found;
pair_found: ;
    }
    if (quantity.group=="cov_b") traces=&ext.cov_b;
    else if (quantity.group=="cov_g") traces=&ext.cov_g;
    else if (quantity.group=="cov_e") traces=&ext.cov_e;
    else if (quantity.group=="component_pi") traces=&ext.component_pi;
    else if (quantity.group=="pattern_pi") traces=&ext.pattern_pi;
    else if (quantity.group=="joint_pi") traces=&ext.joint_pi;
    else if (quantity.group=="alpha") traces=&ext.annotation_alpha;
    else if (quantity.group=="sigmaSqAlpha") traces=&ext.annotation_sigma;
    else if (quantity.group=="b") traces=&ext.selected_b;
    else if (quantity.group=="d") integer_traces=&ext.selected_d;
    else if (quantity.group=="component") integer_traces=&ext.selected_component;
    int index=0;
    if (pair>=0) index=pair;
    else if (quantity.group=="component_pi") index=quantity.component;
    else if (quantity.group=="pattern_pi" || quantity.group=="joint_pi") index=quantity.pattern;
    else if (quantity.group=="alpha") index=quantity.annotation;
    else if (quantity.group=="sigmaSqAlpha") index=quantity.stick;
    else {
     auto found=std::find(extended_spec->selected_markers.begin(),
      extended_spec->selected_markers.end(),quantity.marker);
     const int selected=static_cast<int>(found-extended_spec->selected_markers.begin());
     index=(quantity.group=="component") ? selected : selected*nt+quantity.trait;
    }
    if (quantity.structural) {
     for (int iteration=0;iteration<nit;++iteration) bundle.values.push_back(0.0);
     continue;
    }
    if (traces!=nullptr) {
     if (index<0 || static_cast<std::size_t>(index)>=traces->size())
      throw std::invalid_argument("extended convergence trace index differs");
     const auto& trace=(*traces)[static_cast<std::size_t>(index)];
     if (trace.size()!=expected) throw std::invalid_argument("extended convergence trace length differs");
     for (int iteration=0;iteration<nit;++iteration)
      bundle.values.push_back(trace[static_cast<std::size_t>(nburn+iteration)]);
     continue;
    }
    if (integer_traces!=nullptr) {
     if (index<0 || static_cast<std::size_t>(index)>=integer_traces->size())
      throw std::invalid_argument("extended state trace index differs");
     const auto& trace=(*integer_traces)[static_cast<std::size_t>(index)];
     if (trace.size()!=expected) throw std::invalid_argument("extended state trace length differs");
     for (int iteration=0;iteration<nit;++iteration)
      bundle.values.push_back(static_cast<double>(trace[static_cast<std::size_t>(nburn+iteration)]));
     continue;
    }
   }
   if (traces==nullptr || traces->size()!=static_cast<std::size_t>(nt) ||
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
