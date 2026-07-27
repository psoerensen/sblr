#ifndef SBLR_BLR_MT_DEFAULT_CORE_IMPL_H
#define SBLR_BLR_MT_DEFAULT_CORE_IMPL_H

#include "blr_mt_default_types.h"
#include "blr_mt_ld_access.h"
#include "blr_mt_bayesr_kernel_impl.h"

namespace sblr {
namespace mt {

template <class DataView>
inline void sample_mt_bayesr_marker(
 int marker, int nt, const DataView& data, const MtJointStateSpec& joint,
 const std::vector<double>& marker_scale, std::vector<double>& counts,
 const std::vector<double>& pi, const arma::mat& Ei, const arma::mat& Bi,
 std::vector<std::vector<double>>& residual,
 std::vector<std::vector<double>>& beta,
 std::vector<std::vector<double>>& b,
 std::vector<std::vector<int>>& d, std::vector<int>& component,
 std::mt19937& rng
) {
 arma::vec score(nt), diagonal(nt);
 for (int trait=0;trait<nt;++trait) {
  diagonal[trait]=mt_diagonal(data,trait,marker);
  score[trait]=residual[trait][marker]+diagonal[trait]*b[trait][marker];
 }
 MtJointMarkerKernelResult kernel=mt_joint_marker_kernel(
  score,diagonal,Bi,Ei,joint,pi,marker_scale[marker]);
 const std::size_t selected=mt_joint_draw_state(kernel.probability,rng);
 counts[selected]+=1.0;
 arma::vec beta_new(nt,arma::fill::zeros);
 if (selected>0) beta_new=mt_joint_draw_beta(kernel.states[selected],nt,rng);
 for (int trait=0;trait<nt;++trait) {
  const double effective=selected>0 ?
   static_cast<double>(joint.patterns[selected][trait])*beta_new[trait] : 0.0;
  const double difference=effective-b[trait][marker];
  if (difference!=0.0)
   mt_apply_marker_difference(data,trait,marker,difference,residual[trait]);
  beta[trait][marker]=selected>0 ? beta_new[trait] : 0.0;
  b[trait][marker]=effective;
  d[trait][marker]=selected>0 ? joint.patterns[selected][trait] : 0;
 }
 component[marker]=selected>0 ? joint.component[selected] : 0;
}

template <class DataView>
inline MtDefaultCoreResult run_mt_bayesc_core_impl(
 const DataView& data,
 const MtDefaultModelSpec& model,
 const MtDefaultCovariancePriorView& prior,
 const MtDefaultExecutionSpec& execution,
 MtDefaultInitialState initial_state
) {
 const auto& wy=data.wy;
 const auto& yy=data.yy;
 const auto& n=data.n;
 const auto& models=model.models;
 const auto& sets=model.sets;
 const int method=model.method;
 const auto& ssb_prior=prior.ssb_prior;
 const auto& sse_prior=prior.sse_prior;
 const double nub=prior.nub;
 const double nue=prior.nue;
 auto b=std::move(initial_state.b);
 auto B=std::move(initial_state.B);
 auto E=std::move(initial_state.E);
 auto pi=std::move(initial_state.pi);
 auto component=std::move(initial_state.component);

 // Define local variables
 int nt = static_cast<int>(mt_trait_count(data));
 int m = static_cast<int>(mt_marker_count(data));
 int nmodels = models.size();
 double marker_retained_count = 0.0;
 double covb_retained_count = 0.0;
 double covg_retained_count = 0.0;
 double cove_retained_count = 0.0;
 double pi_retained_count = 0.0;

 //double logliksum, detC, diff, cumprobc;
 //int mselect;

 std::vector<std::vector<int>> d=initial_state.state.empty() ?
  std::vector<std::vector<int>>(nt,std::vector<int>(m,0)) :
  std::move(initial_state.state);
 std::vector<std::vector<double>> beta=initial_state.beta.empty() ?
  std::vector<std::vector<double>>(nt,std::vector<double>(m,0.0)) :
  std::move(initial_state.beta);
 if (method==5 && component.empty()) component.assign(m,0);
 if (method==5 && component.size()!=static_cast<std::size_t>(m))
  throw std::invalid_argument("MT component initialization must have length m");

 std::vector<double> mu(nt), rhs(nt), conv(nt);
 std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);
 std::vector<std::vector<double>> pi_trace(
  nmodels,std::vector<double>(execution.nit+execution.nburn,0.0));
 int component_count=0;
 if (method==5) {
  if (model.joint==nullptr || model.marker_scale==nullptr)
   throw std::invalid_argument("MT BayesR requires joint states and marker scales");
  validate_mt_joint_state_spec(*model.joint,static_cast<std::size_t>(nt));
  if (model.joint->patterns.size()!=static_cast<std::size_t>(nmodels) ||
      model.marker_scale->size()!=static_cast<std::size_t>(m))
   throw std::invalid_argument("MT BayesR state or marker-scale dimensions differ");
  if (d.size()!=static_cast<std::size_t>(nt) ||
      beta.size()!=static_cast<std::size_t>(nt))
   throw std::invalid_argument("MT BayesR initialization trait dimensions differ");
  for (int trait=0;trait<nt;++trait)
   if (d[trait].size()!=static_cast<std::size_t>(m) ||
       beta[trait].size()!=static_cast<std::size_t>(m))
    throw std::invalid_argument("MT BayesR initialization marker dimensions differ");
  component_count=model.joint->component_count;
  if (model.pi_prior!=nullptr) {
   if (model.pi_prior->size()!=static_cast<std::size_t>(nmodels))
    throw std::invalid_argument("MT BayesR pi prior length differs from joint states");
   for (double value:*model.pi_prior)
    if (!std::isfinite(value) || value<=0.0)
     throw std::invalid_argument("MT BayesR pi prior must be finite and positive");
  }
 }
 std::vector<std::vector<double>> component_counts(
  static_cast<std::size_t>(m),std::vector<double>(component_count,0.0));

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> ves(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vle(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vld(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> cvbm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvem(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> cvgm(nt, std::vector<double>(nt, 0.0));
 std::vector<std::vector<double>> mus(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));

 std::vector<double> x2t(m);
 std::vector<std::vector<double>> x2(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> r(nt, std::vector<double>(m, 0.0));
 std::vector<int> order(m);
 std::vector<double> gamma(4, 0.0);

 std::vector<std::vector<double>> pitrait(nt, std::vector<double>(4, 0.0));
 std::vector<std::vector<double>> pistrait(nt, std::vector<double>(4, 0.0));
 std::vector<double> pimarker(2, 0.0);
 std::vector<double> pismarker(2, 0.0);
 std::vector<int> dmarker(m, 0);

 gamma[0] = 0.0;
 gamma[1] = 0.01;
 gamma[2] = 0.1;
 gamma[3] = 1.0;

 for (int t = 0; t < nt; t++) {
  pitrait[t][0] = 1.0-pi[0];
  pitrait[t][1] = pi[0];
 }
 if(method==5){
  for (int t = 0; t < nt; t++) {
   pitrait[t][0] = 0.95;
   pitrait[t][1] = 0.02;
   pitrait[t][2] = 0.02;
   pitrait[t][3] = 0.01;
  }
 }

 pimarker[0] = 1.0-pi[0];
 pimarker[1] = pi[0];

 // Initialize representation-neutral marker scores.
 for (int i = 0; i < m; i++) {
  for (int t = 0; t < nt; t++) {
   const double diagonal=mt_diagonal(data, t, i);
   x2[t][i] = (wy[t][i]/diagonal)*(wy[t][i]/diagonal);
  }
 }

 // Wy - W'Wb, preserving each representation's stored traversal order.
 mt_rebuild_residuals(data, b, r);

 for ( int i = 0; i < m; i++) {
  x2t[i] = 0.0;
  for ( int t = 0; t < nt; t++) {
   if(x2[t][i]>x2t[i]) {
    x2t[i] = x2[t][i];
   }
  }
 }


 // Establish order of markers as they are entered into the model
 std::iota(order.begin(), order.end(), 0);
 std::sort(  std::begin(order),
             std::end(order),
             [&](int i1, int i2) { return x2t[i1] > x2t[i2]; } );


 // Initialize (co)variance matrices
 arma::mat C(nt,nt, fill::zeros);
 arma::mat Bi = arma::inv(B);
 arma::mat Ei = arma::inv(E);
 arma::mat G(nt,nt, fill::zeros);
 arma::rowvec probs(nmodels, fill::zeros);

 // Start Gibbs sampler
 std::random_device rd;
 std::mt19937 gen(execution.seed);

 for ( int it = 0; it < execution.nit+execution.nburn; it++) {

  if (model.pi_prior==nullptr) std::fill(cmodel.begin(), cmodel.end(), 1.0);
  else cmodel=*model.pi_prior;

  //Sample marker effects (BayesC)
  if (method==4) {
   for (size_t s = 0; s < sets.size(); ++s) {
    const std::vector<int>& set = sets[s];

    // Step 1: Sample B using only markers in the current set
    if (execution.updateB) {
     sampleBset(nt, m, nub, B, d, b, ssb_prior, set, gen);
     //sampleB(nt, m, nub, B, d, b, ssb_prior, gen);
     sampleB_latent(nt, m, nub, B, beta, ssb_prior, gen);
    }

    // Step 2: Invert B to get Bi
    arma::mat Bi = arma::inv(B);

    // Step 3: Sample marker effects for markers in this set
    // Convert current set to a std::unordered_set for fast lookup
    std::unordered_set<int> current_set(set.begin(), set.end());

    for (int isort = 0; isort < m; isort++) {
     int i = order[isort];

     // Only update if marker i is in the current set
     if (current_set.find(i) != current_set.end()) {
      // sampleBetaCMt_fast(i, nt,
      //                 nmodels, models, cmodel, pi,
      //                 Ei, Bi,
      //                 ww, r, b, d,
      //                 XXindices, XXvalues,
      //                 gen);
     //  sampleBetaCStMt(i, nt,
     //                nmodels, models, cmodel, pi,
     //                Ei, Bi,
     //                ww, r, b, d,
     //                XXindices, XXvalues,
     //                gen);

      // sampleBetaCPG_Mt(i, nt,
      //               nmodels, models, cmodel, pi,
      //               Ei, Bi,
      //               ww, r, b, d,
      //               XXindices, XXvalues,
      //               gen);

      sampleBetaCPG_Mt_latent(i, nt,
                       nmodels, models, cmodel, pi,
                       Ei, Bi,
                       data, r, beta, b, d,
                       gen);

     }
    }
   }
  }
  if (method==5) {
   const auto& joint=*model.joint;
   const auto& scales=*model.marker_scale;
   for (size_t s=0;s<sets.size();++s) {
    const auto& set=sets[s];
    if (execution.updateB) {
     const auto base=mt_bayesr_base_effects(b,component,joint,scales);
     const auto base_beta=mt_bayesr_base_effects(beta,component,joint,scales);
     sampleBset(nt,m,nub,B,d,base,ssb_prior,set,gen);
     sampleB_latent(nt,m,nub,B,base_beta,ssb_prior,gen);
    }
    Bi=arma::inv(B);
    std::unordered_set<int> current_set(set.begin(),set.end());
    for (int isort=0;isort<m;++isort) {
     const int marker=order[isort];
     if (current_set.find(marker)!=current_set.end())
      sample_mt_bayesr_marker(marker,nt,data,joint,scales,cmodel,pi,Ei,Bi,
       r,beta,b,d,component,gen);
    }
   }
  }
  // Store values
  for (int t = 0; t < nt; t++) {
   for ( int i = 0; i < m; i++) {
    //if (d[t][i] == 1) {
    if (d[t][i] > 0) {
     if ((it >= execution.nburn) && ((it - execution.nburn) % execution.nthin == 0)) {
      dm[t][i] = dm[t][i] + 1.0;
       bm[t][i] = bm[t][i] + b[t][i];
     }
    }
   }
  }
  if (method==5 && it>=execution.nburn &&
      ((it-execution.nburn)%execution.nthin==0))
   for (int marker=0;marker<m;++marker)
    component_counts[marker][component[marker]]+=1.0;

  // Sample pi for Bayes C
  if(execution.updatePi && (method==4 || method==5)) {
   samplePi(cmodel, pi, gen);
   for (int k = 0; k<nmodels ; k++) {
    if(it >= execution.nburn) pis[k] = pis[k] + pi[k];
   }
   if (it >= execution.nburn) pi_retained_count = pi_retained_count + 1.0;
  }
  if (method==5 && !execution.updatePi && it>=execution.nburn) {
   for (int k=0;k<nmodels;++k) pis[k]+=pi[k];
   pi_retained_count+=1.0;
  }
  for (int state=0;state<nmodels;++state) pi_trace[state][it]=pi[state];

  // Sample marker variance
  if(execution.updateB && (method==4 || method==5)) {
   if (method==4) sampleB(nt,m,nub,B,d,b,ssb_prior,gen);
   else {
    const auto base=mt_bayesr_base_effects(
     b,component,*model.joint,*model.marker_scale);
    sampleB(nt,m,nub,B,d,base,ssb_prior,gen);
   }
   for (int t = 0; t < nt; t++) {
    vbs[t][it] = B(t,t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if(it >= execution.nburn) cvbm[t1][t2] = cvbm[t1][t2] + B(t1,t2);
    }
   }
   if (it >= execution.nburn) covb_retained_count = covb_retained_count + 1.0;
  }

  //Update genetic variance
  computeG(nt, m, b, wy, r, n, G);
  for (int t = 0; t < nt; t++) {
   vgs[t][it] = G(t,t);
   long double diagonal_contribution=0.0L;
   for (int i=0; i<m; ++i) {
    const long double effect=static_cast<long double>(b[t][i]);
    diagonal_contribution+=static_cast<long double>(mt_diagonal(data,t,i))*
     effect*effect;
   }
   vle[t][it]=static_cast<double>(
    diagonal_contribution/static_cast<long double>(n[t]));
   vld[t][it]=vgs[t][it]-vle[t][it];
  }
  for (int t1 = 0; t1 < nt; t1++  ) {
   for (int t2 = 0; t2 < nt; t2++) {
    if(it >= execution.nburn) cvgm[t1][t2] = cvgm[t1][t2] + G(t1,t2);
   }
  }
  if (it >= execution.nburn) covg_retained_count = covg_retained_count + 1.0;

  // Sample residual variance
  if(execution.updateE) {
   sampleE(nt, m, nue, E, b, wy, r, sse_prior, yy, n, gen);
   Ei = arma::inv(E);
   for (int t = 0; t < nt; t++) {
    ves[t][it] = E(t,t);
   }
   for (int t1 = 0; t1 < nt; t1++) {
    for (int t2 = 0; t2 < nt; t2++) {
     if(it >= execution.nburn) cvem[t1][t2] = cvem[t1][t2] + E(t1,t2);
    }
   }
   if (it >= execution.nburn) cove_retained_count = cove_retained_count + 1.0;
  }

  if ((it >= execution.nburn) && ((it - execution.nburn) % execution.nthin == 0)) {
   marker_retained_count = marker_retained_count + 1.0;
  }

 }
 if (marker_retained_count <= 0.0) {
  throw std::runtime_error("mtblr: no retained marker-summary samples.");
 }

 MtDefaultCoreResult result;
 result.nt=nt;
 result.m=m;
 result.nmodels=nmodels;
 result.marker_retained_count=marker_retained_count;
 result.covb_retained_count=covb_retained_count;
 result.covg_retained_count=covg_retained_count;
 result.cove_retained_count=cove_retained_count;
 result.pi_retained_count=pi_retained_count;
 result.bm=std::move(bm);
 result.dm=std::move(dm);
 result.r=std::move(r);
 result.b=std::move(b);
 result.d=std::move(d);
 result.component=std::move(component);
 result.component_counts=std::move(component_counts);
 result.order=std::move(order);
 result.vbs=std::move(vbs);
 result.vgs=std::move(vgs);
 result.ves=std::move(ves);
 result.vle=std::move(vle);
 result.vld=std::move(vld);
 result.cvbm=std::move(cvbm);
 result.cvgm=std::move(cvgm);
 result.cvem=std::move(cvem);
 result.B=std::move(B);
 result.G=std::move(G);
 result.E=std::move(E);
 result.pi=std::move(pi);
 result.pis=std::move(pis);
 result.pi_trace=std::move(pi_trace);
 result.pistrait=std::move(pistrait);
 result.pismarker=std::move(pismarker);
 return result;
}

inline MtDefaultCoreResult run_mt_default_core(
 const MtDefaultDataView& data,
 const MtDefaultModelSpec& model,
 const MtDefaultCovariancePriorView& prior,
 const MtDefaultExecutionSpec& execution,
 MtDefaultInitialState initial_state
) {
 return run_mt_bayesc_core_impl(
  data, model, prior, execution, std::move(initial_state));
}

inline MtDefaultCoreResult run_mt_csr_core(
 const MtCsrDataView& data,
 const MtDefaultModelSpec& model,
 const MtDefaultCovariancePriorView& prior,
 const MtDefaultExecutionSpec& execution,
 MtDefaultInitialState initial_state
) {
 validate_mt_csr_data(data, model);
 return run_mt_bayesc_core_impl(
  data, model, prior, execution, std::move(initial_state));
}

inline MtDefaultCoreResult run_mt_block_eigen_core(
 const MtBlockEigenDataView& data,
 const MtDefaultModelSpec& model,
 const MtDefaultCovariancePriorView& prior,
 const MtDefaultExecutionSpec& execution,
 MtDefaultInitialState initial_state
) {
 validate_mt_block_eigen_data(data, model);
 return run_mt_bayesc_core_impl(
  data, model, prior, execution, std::move(initial_state));
}

}  // namespace mt
}  // namespace sblr

#endif
