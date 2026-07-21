#ifndef SBLR_BLR_MT_DEFAULT_CORE_IMPL_H
#define SBLR_BLR_MT_DEFAULT_CORE_IMPL_H

#include "blr_mt_default_types.h"
#include "blr_mt_ld_access.h"

namespace sblr {
namespace mt {

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

 std::vector<std::vector<int>> d(nt, std::vector<int>(m, 0));
 std::vector<std::vector<double>> beta(nt, std::vector<double>(m, 0));

 std::vector<double> mu(nt), rhs(nt), conv(nt);
 std::vector<double> pmodel(nmodels), pcum(nmodels), loglik(nmodels), cmodel(nmodels);
 std::vector<double> pis(nmodels, 0.0);

 std::vector<std::vector<double>> bm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> dm(nt, std::vector<double>(m, 0.0));
 std::vector<std::vector<double>> ves(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vbs(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
 std::vector<std::vector<double>> vgs(nt, std::vector<double>(execution.nit+execution.nburn, 0.0));
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

  std::fill(cmodel.begin(), cmodel.end(), 1.0);

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

  // Sample pi for Bayes C
  if(execution.updatePi && method==4) {
   samplePi(cmodel, pi, gen);
   for (int k = 0; k<nmodels ; k++) {
    if(it >= execution.nburn) pis[k] = pis[k] + pi[k];
   }
   if (it >= execution.nburn) pi_retained_count = pi_retained_count + 1.0;
  }

  // Sample marker variance
  if(execution.updateB && method==4) {
   sampleB(nt, m, nub, B, d, b, ssb_prior, gen);
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
 result.order=std::move(order);
 result.vbs=std::move(vbs);
 result.vgs=std::move(vgs);
 result.ves=std::move(ves);
 result.cvbm=std::move(cvbm);
 result.cvgm=std::move(cvgm);
 result.cvem=std::move(cvem);
 result.B=std::move(B);
 result.G=std::move(G);
 result.E=std::move(E);
 result.pi=std::move(pi);
 result.pis=std::move(pis);
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
