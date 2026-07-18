#ifndef SBLR_BLR_BED_BAYESRC_AGGREGATE_IMPL_H
#define SBLR_BLR_BED_BAYESRC_AGGREGATE_IMPL_H

#include "blr_bed_bayesrc_types.h"
#include "st_bayesrc_annotation_prior.h"

// Implementation detail: included only by the packed-BED BayesRC binding unit.
namespace sblr { namespace core {

static BedBayesRCExecutionResult aggregate_bed_bayesrc_results(
 const std::vector<BedBayesRCChainExecutionResult>& chain_results,
 const BedBayesRCAggregationContext& context
) {
 const std::size_t m=context.marker_count, na=context.annotation_count;
 const std::size_t k=context.component_count, nt=context.trait_count;
 const std::size_t nchains=context.chain_count, ntrace=context.trace_length;
 if (m==0 || na==0 || k<2 || nt==0 || nchains==0 ||
     chain_results.size()!=nt*nchains || context.annotation.n_rows!=m ||
     context.annotation.n_cols!=na || context.component_scales.size()!=k)
  throw std::invalid_argument("BayesRC aggregation dimensions are inconsistent");

 BedBayesRCExecutionResult out;
 out.marker_count=m; out.annotation_count=na; out.component_count=k;
 out.trait_count=nt; out.chain_count=nchains; out.trace_length=ntrace;
 out.bm.zeros(m,nt); out.dm.zeros(m,nt); out.b.zeros(m,nt);
 out.state.zeros(m,nt); out.component_mean.zeros(m,nt);
 out.vbs.zeros(ntrace,nt); out.vgs.zeros(ntrace,nt); out.ves.zeros(ntrace,nt);
 out.vle.zeros(ntrace,nt); out.vld.zeros(ntrace,nt); out.pis.zeros(ntrace,nt);
 out.final_prior.zeros(nt,k); out.mean_prior.zeros(nt,k);
 out.comp_prob.resize(nt); out.marker_prior_final.resize(nt);
 out.alpha_mean.resize(nt); out.alpha_final.resize(nt);
 out.sigma_mean.zeros(nt,k-1); out.sigma_final.zeros(nt,k-1);
 out.final_vb.zeros(nt); out.final_vg.zeros(nt); out.final_ve.zeros(nt);
 out.log_cpo.zeros(nt); out.mean_log_cpo.zeros(nt); out.nsamples.zeros(nt);
 out.component_counts.zeros(nt,k);
 if (context.keep_chains) out.retained_chains=chain_results;

 for (std::size_t t=0;t<nt;++t) {
  out.comp_prob[t].zeros(m,k); out.marker_prior_final[t].zeros(m,k);
  out.alpha_mean[t].zeros(na,k-1); out.alpha_final[t].zeros(na,k-1);
  for (std::size_t ch=0;ch<nchains;++ch) {
   const BedBayesRCChainExecutionResult& z=chain_results[ch*nt+t];
   if (z.failed) out.failures.push_back(z.error);
   out.bm.col(t)+=z.bm.t(); out.dm.col(t)+=z.dm.t(); out.b.col(t)+=z.b.t();
   out.state.col(t)+=z.state.t(); out.component_mean.col(t)+=z.component_mean.t();
   out.vbs.col(t)+=z.vbs.t(); out.vgs.col(t)+=z.vgs.t(); out.ves.col(t)+=z.ves.t();
   out.vle.col(t)+=z.vles.t(); out.vld.col(t)+=z.vlds.t(); out.pis.col(t)+=z.pis.t();
   out.comp_prob[t]+=z.comp_prob; out.alpha_mean[t]+=z.annot_alpha_mean;
   out.alpha_final[t]+=z.annot_alpha_final;
   out.sigma_mean.row(t)+=z.annot_sigma_mean.t(); out.mean_prior.row(t)+=z.mean_prior;
   out.sigma_final.row(t)+=z.annot_sigma_final.t();
   const arma::mat marker_prior=st_bayesrc_compute_snp_pi(
    context.annotation,z.annot_alpha_final,context.pi_floor
   );
   out.marker_prior_final[t]+=marker_prior;
   out.final_prior.row(t)+=arma::mean(marker_prior,0);
   out.final_vb(t)+=z.final_vb; out.final_vg(t)+=z.final_vg; out.final_ve(t)+=z.final_ve;
   out.log_cpo(t)+=z.log_cpo; out.mean_log_cpo(t)+=z.mean_log_cpo;
   out.nsamples(t)+=z.nsamples;
  }
 }
 const double inv=1.0/static_cast<double>(nchains);
 out.bm*=inv; out.dm*=inv; out.b*=inv; out.state*=inv; out.component_mean*=inv;
 out.vbs*=inv; out.vgs*=inv; out.ves*=inv; out.vle*=inv; out.vld*=inv; out.pis*=inv;
 out.final_prior*=inv; out.mean_prior*=inv; out.sigma_mean*=inv; out.sigma_final*=inv;
 out.final_vb*=inv; out.final_vg*=inv; out.final_ve*=inv;
 out.log_cpo*=inv; out.mean_log_cpo*=inv; out.nsamples*=inv;
 for (std::size_t t=0;t<nt;++t) {
  out.comp_prob[t]*=inv; out.marker_prior_final[t]*=inv;
  out.alpha_mean[t]*=inv; out.alpha_final[t]*=inv;
  out.component_counts.row(t)=arma::sum(out.comp_prob[t],0);
 }
 return out;
}

} }

#endif
