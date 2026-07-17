#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_AGGREGATE_IMPL_H
#define SBLR_BLR_BED_SCHEDULED_BAYESC_AGGREGATE_IMPL_H

// Implementation detail: include only from the packed-BED BayesC binding TU,
// after its packed-genotype numerical helpers have been declared.

#include "blr_bed_scheduled_bayesc_types.h"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sblr { namespace core {

template <class PackedGenotype>
BedScheduledBayesCExecutionResult aggregate_bed_scheduled_bayesc_results(
 const std::vector<BedScheduledBayesCChainExecutionResult>& chain_results,
 const BedScheduledBayesCAggregationContext<PackedGenotype>& context
) {
 const int nt=context.trait_count;
 const int nchains=context.chain_count;
 const int m=static_cast<int>(context.genotype.marker_count);
 const int n=static_cast<int>(context.genotype.sample_count);
 const int n_trace=context.trace_length;
 const int njobs=nt*nchains;
 if (nt<=0 || nchains<=0 || m<=0 || n<=0 || n_trace<=0)
  throw std::invalid_argument("BED BayesC aggregation dimensions must be positive");
 if (chain_results.size()!=static_cast<std::size_t>(njobs))
  throw std::invalid_argument("BED BayesC chain-result count is inconsistent");
 if (context.marker_maps.size()!=static_cast<std::size_t>(m) ||
     context.marker_order.size()!=static_cast<std::size_t>(m) ||
     context.phenotype.n_rows!=static_cast<arma::uword>(n) ||
     context.phenotype.n_cols!=static_cast<arma::uword>(nt))
  throw std::invalid_argument("BED BayesC aggregation input dimensions are inconsistent");

 BedScheduledBayesCExecutionResult out;
 out.bm.zeros(nt,m); out.dm.zeros(nt,m);
 out.bm_sd.zeros(nt,m); out.dm_sd.zeros(nt,m);
 out.bm_min.zeros(nt,m); out.dm_min.zeros(nt,m);
 out.bm_max.zeros(nt,m); out.dm_max.zeros(nt,m);
 out.final_effects.zeros(nt,m); out.final_states.zeros(nt,m);
 out.wy.zeros(nt,m); out.residual_scores.zeros(nt,m);
 out.marker_variance_trace.zeros(nt,n_trace);
 out.genetic_variance_trace.zeros(nt,n_trace);
 out.residual_variance_trace.zeros(nt,n_trace);
 out.inclusion_trace.zeros(nt,n_trace);
 out.vle_trace.zeros(nt,n_trace); out.vld_trace.zeros(nt,n_trace);
 out.final_marker_variance.zeros(nt); out.final_genetic_variance.zeros(nt);
 out.final_residual_variance.zeros(nt); out.final_inclusion_probability.zeros(nt);
 out.final_vle.zeros(nt); out.final_vld.zeros(nt);
 out.mean_inclusion_probability.zeros(nt); out.mean_total_log_cpo.zeros(nt);
 out.mean_log_cpo.zeros(nt); out.mean_retained_samples.zeros(nt);
 out.mean_seconds.zeros(nt); out.max_seconds.zeros(nt);
 out.marker_count=m; out.sample_count=n; out.trait_count=nt;
 out.chain_count=nchains; out.task_count=njobs; out.trace_length=n_trace;

 for (int chain=0; chain<nchains; ++chain) {
  for (int t=0; t<nt; ++t) {
   const auto& r=chain_results[static_cast<std::size_t>(chain*nt+t)];
   const arma::uword tu=static_cast<arma::uword>(t);
   if (r.bm.n_elem!=static_cast<arma::uword>(m) ||
       r.dm.n_elem!=static_cast<arma::uword>(m) ||
       r.b.n_elem!=static_cast<arma::uword>(m) ||
       r.d_as_double.n_elem!=static_cast<arma::uword>(m) ||
       r.vbs.n_elem!=static_cast<arma::uword>(n_trace) ||
       r.vgs.n_elem!=static_cast<arma::uword>(n_trace) ||
       r.ves.n_elem!=static_cast<arma::uword>(n_trace) ||
       r.pis.n_elem!=static_cast<arma::uword>(n_trace) ||
       r.vles.n_elem!=static_cast<arma::uword>(n_trace) ||
       r.vlds.n_elem!=static_cast<arma::uword>(n_trace))
    throw std::invalid_argument("BED BayesC chain-result dimensions are inconsistent");
   out.bm.row(tu)+=r.bm; out.dm.row(tu)+=r.dm;
   out.final_effects.row(tu)+=r.b; out.final_states.row(tu)+=r.d_as_double;
   out.marker_variance_trace.row(tu)+=r.vbs;
   out.genetic_variance_trace.row(tu)+=r.vgs;
   out.residual_variance_trace.row(tu)+=r.ves;
   out.inclusion_trace.row(tu)+=r.pis;
   out.vle_trace.row(tu)+=r.vles; out.vld_trace.row(tu)+=r.vlds;
   out.final_marker_variance(tu)+=r.final_vb;
   out.final_genetic_variance(tu)+=r.final_vg;
   out.final_residual_variance(tu)+=r.final_ve;
   out.final_vle(tu)+=r.final_vle; out.final_vld(tu)+=r.final_vld;
   out.final_inclusion_probability(tu)+=r.final_pi;
   out.mean_inclusion_probability(tu)+=r.mean_pi;
   out.mean_total_log_cpo(tu)+=r.log_cpo;
   out.mean_log_cpo(tu)+=r.mean_log_cpo;
   out.mean_retained_samples(tu)+=r.nsamples;
   out.mean_seconds(tu)+=r.seconds;
   out.max_seconds(tu)=std::max(out.max_seconds(tu),r.seconds);
  }
 }

 const double inv_chains=1.0/static_cast<double>(nchains);
 out.bm*=inv_chains; out.dm*=inv_chains;
 out.final_effects*=inv_chains; out.final_states*=inv_chains;
 out.marker_variance_trace*=inv_chains; out.genetic_variance_trace*=inv_chains;
 out.residual_variance_trace*=inv_chains; out.inclusion_trace*=inv_chains;
 out.vle_trace*=inv_chains; out.vld_trace*=inv_chains;
 out.final_marker_variance*=inv_chains; out.final_genetic_variance*=inv_chains;
 out.final_residual_variance*=inv_chains; out.final_vle*=inv_chains;
 out.final_vld*=inv_chains; out.final_inclusion_probability*=inv_chains;
 out.mean_inclusion_probability*=inv_chains; out.mean_total_log_cpo*=inv_chains;
 out.mean_log_cpo*=inv_chains; out.mean_retained_samples*=inv_chains;
 out.mean_seconds*=inv_chains;

 for (int t=0; t<nt; ++t) {
  const arma::uword tu=static_cast<arma::uword>(t);
  out.bm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  out.dm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  out.bm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  out.dm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  for (int chain=0; chain<nchains; ++chain) {
   const auto& r=chain_results[static_cast<std::size_t>(chain*nt+t)];
   for (int j=0; j<m; ++j) {
    const arma::uword ju=static_cast<arma::uword>(j);
    out.bm_min(tu,ju)=std::min(out.bm_min(tu,ju),r.bm(ju));
    out.dm_min(tu,ju)=std::min(out.dm_min(tu,ju),r.dm(ju));
    out.bm_max(tu,ju)=std::max(out.bm_max(tu,ju),r.bm(ju));
    out.dm_max(tu,ju)=std::max(out.dm_max(tu,ju),r.dm(ju));
   }
   if (nchains>1) {
    const arma::rowvec bm_diff=r.bm-out.bm.row(tu);
    const arma::rowvec dm_diff=r.dm-out.dm.row(tu);
    out.bm_sd.row(tu)+=bm_diff%bm_diff;
    out.dm_sd.row(tu)+=dm_diff%dm_diff;
   }
  }
  if (nchains>1) {
   out.bm_sd.row(tu)=arma::sqrt(out.bm_sd.row(tu)/static_cast<double>(nchains-1));
   out.dm_sd.row(tu)=arma::sqrt(out.dm_sd.row(tu)/static_cast<double>(nchains-1));
  }
 }

 if (context.return_wy || context.return_residual_scores) {
  const auto& G=context.genotype.storage;
  for (int t=0; t<nt; ++t) {
   const arma::vec y_t=context.phenotype.col(static_cast<arma::uword>(t));
   const arma::rowvec b_t=out.final_effects.row(static_cast<arma::uword>(t));
   const arma::vec xb_t=bed_xb_from_b_scheduled_chains(
    G,context.marker_maps,context.marker_order,b_t);
   const arma::vec e_t=y_t-xb_t;
   for (int j=0; j<m; ++j) {
    if (context.return_wy)
     out.wy(static_cast<arma::uword>(t),static_cast<arma::uword>(j))=
      bed_marker_dot_residual_scheduled_chains(
       G,j,context.marker_maps[static_cast<std::size_t>(j)],y_t.memptr());
    if (context.return_residual_scores)
     out.residual_scores(static_cast<arma::uword>(t),static_cast<arma::uword>(j))=
      bed_marker_dot_residual_scheduled_chains(
       G,j,context.marker_maps[static_cast<std::size_t>(j)],e_t.memptr());
   }
  }
 }
 return out;
}

} }

#endif
