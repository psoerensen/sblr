#ifndef SBLR_BLR_BED_BAYESR_AGGREGATE_IMPL_H
#define SBLR_BLR_BED_BAYESR_AGGREGATE_IMPL_H

// Binding-neutral implementation detail for packed-BED BayesR aggregation.
namespace sblr { namespace core {

inline BedBayesRExecutionResult aggregate_bed_bayesr_results(
 const std::vector<BedBayesRChainExecutionResult>& chain_results,
 const BedBayesRAggregationContext& context
) {
 const int m=static_cast<int>(context.marker_count);
 const int nt=static_cast<int>(context.trait_count);
 const int nchains=static_cast<int>(context.chain_count);
 const int ntrace=static_cast<int>(context.trace_length);
 const int K=static_cast<int>(context.component_count);
 if (m<=0 || nt<=0 || nchains<=0 || ntrace<=0 || K<2 || context.null_component!=0 ||
     chain_results.size()!=context.trait_count*context.chain_count)
  throw std::invalid_argument("BayesR aggregation dimensions are inconsistent");

 BedBayesRExecutionResult out;
 out.marker_count=context.marker_count; out.trait_count=context.trait_count;
 out.chain_count=context.chain_count; out.trace_length=context.trace_length;
 out.component_count=context.component_count;
 out.bm=arma::mat(nt,m,arma::fill::zeros); out.dm=arma::mat(nt,m,arma::fill::zeros);
 out.component_mean=arma::mat(nt,m,arma::fill::zeros);
 out.final_effects=arma::mat(nt,m,arma::fill::zeros);
 out.final_states=arma::mat(nt,m,arma::fill::zeros);
 out.bm_sd=arma::mat(nt,m,arma::fill::zeros); out.dm_sd=arma::mat(nt,m,arma::fill::zeros);
 out.bm_min=arma::mat(nt,m,arma::fill::zeros); out.dm_min=arma::mat(nt,m,arma::fill::zeros);
 out.bm_max=arma::mat(nt,m,arma::fill::zeros); out.dm_max=arma::mat(nt,m,arma::fill::zeros);
 out.wy=arma::mat(nt,m,arma::fill::zeros); out.residual_score=arma::mat(nt,m,arma::fill::zeros);
 out.vbs=arma::mat(nt,ntrace,arma::fill::zeros); out.vgs=arma::mat(nt,ntrace,arma::fill::zeros);
 out.ves=arma::mat(nt,ntrace,arma::fill::zeros); out.vles=arma::mat(nt,ntrace,arma::fill::zeros);
 out.vlds=arma::mat(nt,ntrace,arma::fill::zeros);
 out.final_vb=arma::vec(nt,arma::fill::zeros); out.final_vg=arma::vec(nt,arma::fill::zeros);
 out.final_ve=arma::vec(nt,arma::fill::zeros); out.final_vle=arma::vec(nt,arma::fill::zeros);
 out.final_vld=arma::vec(nt,arma::fill::zeros);
 out.final_pi=arma::mat(K,nt,arma::fill::zeros); out.mean_pi=arma::mat(K,nt,arma::fill::zeros);
 out.log_cpo=arma::vec(nt,arma::fill::zeros); out.mean_log_cpo=arma::vec(nt,arma::fill::zeros);
 out.retained_samples=arma::vec(nt,arma::fill::zeros);
 out.seconds_mean=arma::vec(nt,arma::fill::zeros); out.seconds_max=arma::vec(nt,arma::fill::zeros);
 out.component_probability.assign(static_cast<std::size_t>(nt),
  arma::mat(K,m,arma::fill::zeros));

 for (int ch=0;ch<nchains;++ch) for (int t=0;t<nt;++t) {
  const auto& r=chain_results[static_cast<std::size_t>(ch*nt+t)];
  const arma::uword tu=static_cast<arma::uword>(t);
  out.bm.row(tu)+=r.bm; out.dm.row(tu)+=r.dm;
  out.component_mean.row(tu)+=r.component_mean;
  out.final_effects.row(tu)+=r.b; out.final_states.row(tu)+=r.d_as_double;
  out.vbs.row(tu)+=r.vbs; out.vgs.row(tu)+=r.vgs; out.ves.row(tu)+=r.ves;
  out.vles.row(tu)+=r.vles; out.vlds.row(tu)+=r.vlds;
  out.component_probability[static_cast<std::size_t>(t)]+=r.pip_k;
  out.final_vb(tu)+=r.final_vb; out.final_vg(tu)+=r.final_vg; out.final_ve(tu)+=r.final_ve;
  out.final_vle(tu)+=r.final_vle; out.final_vld(tu)+=r.final_vld;
  out.log_cpo(tu)+=r.log_cpo; out.mean_log_cpo(tu)+=r.mean_log_cpo;
  out.retained_samples(tu)+=r.nsamples; out.seconds_mean(tu)+=r.seconds;
  out.seconds_max(tu)=std::max(out.seconds_max(tu),r.seconds);
  out.failed_tasks+=r.failed;
  for (int k=0;k<K;++k) {
   out.final_pi(k,tu)+=r.final_pi[static_cast<std::size_t>(k)];
   out.mean_pi(k,tu)+=r.mean_pi[static_cast<std::size_t>(k)];
  }
 }
 const double inv=1.0/static_cast<double>(nchains);
 out.bm*=inv; out.dm*=inv; out.component_mean*=inv; out.final_effects*=inv; out.final_states*=inv;
 out.vbs*=inv; out.vgs*=inv; out.ves*=inv; out.vles*=inv; out.vlds*=inv;
 for (auto& x:out.component_probability) x*=inv;
 out.final_vb*=inv; out.final_vg*=inv; out.final_ve*=inv; out.final_vle*=inv; out.final_vld*=inv;
 out.final_pi*=inv; out.mean_pi*=inv; out.log_cpo*=inv; out.mean_log_cpo*=inv;
 out.retained_samples*=inv; out.seconds_mean*=inv;

 for (int t=0;t<nt;++t) {
  const arma::uword tu=static_cast<arma::uword>(t);
  out.bm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  out.dm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  out.bm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  out.dm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  for (int ch=0;ch<nchains;++ch) {
   const auto& r=chain_results[static_cast<std::size_t>(ch*nt+t)];
   for (int j=0;j<m;++j) {
    const arma::uword ju=static_cast<arma::uword>(j);
    out.bm_min(tu,ju)=std::min(out.bm_min(tu,ju),r.bm(ju));
    out.dm_min(tu,ju)=std::min(out.dm_min(tu,ju),r.dm(ju));
    out.bm_max(tu,ju)=std::max(out.bm_max(tu,ju),r.bm(ju));
    out.dm_max(tu,ju)=std::max(out.dm_max(tu,ju),r.dm(ju));
   }
   if (nchains>1) {
    const arma::rowvec bd=r.bm-out.bm.row(tu), dd=r.dm-out.dm.row(tu);
    out.bm_sd.row(tu)+=bd%bd; out.dm_sd.row(tu)+=dd%dd;
   }
  }
  if (nchains>1) {
   out.bm_sd.row(tu)=arma::sqrt(out.bm_sd.row(tu)/static_cast<double>(nchains-1));
   out.dm_sd.row(tu)=arma::sqrt(out.dm_sd.row(tu)/static_cast<double>(nchains-1));
  }
 }
 return out;
}

} }
#endif
