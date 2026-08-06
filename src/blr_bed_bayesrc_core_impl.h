#ifndef SBLR_BLR_BED_BAYESRC_CORE_IMPL_H
#define SBLR_BLR_BED_BAYESRC_CORE_IMPL_H

#include "blr_bed_bayesrc_types.h"
#include "st_bayesrc_annotation_prior.h"

#include <array>
#include <chrono>

// Implementation detail: included only by the packed-BED BayesRC binding unit.
namespace sblr { namespace core {

static inline void sample_marker_bayesrc(
  const FastPackedBedMatrixBR& G, int marker, const MarkerMapBayesR& map,
  const arma::rowvec& marker_prior, const std::vector<double>& gamma,
  double vb, double vei, arma::vec& residual, double& marker_effect,
  int& component, std::mt19937& gen,
  std::uniform_real_distribution<double>& runif,
  std::normal_distribution<double>& norm01
) {
 const int ncomponent = static_cast<int>(gamma.size());
 const double vei_safe = std::max(vei, 1e-300);
 const double score = br_dot_residual(G, marker, map, residual.memptr()) +
  map.xx * marker_effect;
 std::vector<double> logp(static_cast<std::size_t>(ncomponent));
 logp[0] = std::log(std::max(marker_prior(0), 1e-300));
 for (int k = 1; k < ncomponent; ++k) {
  const double vbk = vb * gamma[static_cast<std::size_t>(k)];
  const double denom = std::max(vei_safe + map.xx * vbk, 1e-300);
  logp[static_cast<std::size_t>(k)] =
   std::log(std::max(marker_prior(static_cast<arma::uword>(k)), 1e-300)) +
   0.5 * std::log(vei_safe / denom) +
   0.5 * score * score * vbk / (vei_safe * denom);
 }
 const double mx = *std::max_element(logp.begin(), logp.end());
 std::vector<double> prob(static_cast<std::size_t>(ncomponent));
 double psum = 0.0;
 for (int k = 0; k < ncomponent; ++k) {
  prob[static_cast<std::size_t>(k)] = std::exp(logp[static_cast<std::size_t>(k)] - mx);
  psum += prob[static_cast<std::size_t>(k)];
 }
 int component_new = ncomponent - 1;
 const double u = runif(gen);
 double cumulative = 0.0;
 for (int k = 0; k < ncomponent - 1; ++k) {
  cumulative += prob[static_cast<std::size_t>(k)] / psum;
  if (u < cumulative) { component_new = k; break; }
 }
 double effect_new = 0.0;
 if (component_new > 0) {
  const double vbk = vb * gamma[static_cast<std::size_t>(component_new)];
  const double lhs = map.xx + vei_safe / vbk;
  effect_new = score / lhs + std::sqrt(vei_safe / lhs) * norm01(gen);
 }
 const double diff = effect_new - marker_effect;
 if (diff != 0.0)
  br_update_residual(G, marker, map, residual.memptr(), diff);
 marker_effect = effect_new;
 component = component_new;
}

struct BedBayesRCTemperedReplicaState {
 arma::rowvec effect;
 arma::Row<int> component;
 arma::vec residual;
 double marker_variance=0.0, genetic_variance=0.0;
 double residual_variance=0.0, adjusted_residual_variance=0.0;
 arma::mat alpha, marker_probability;
 arma::vec sigma_sq_alpha;
 int identity=0;
 std::mt19937 generator;
 std::uniform_real_distribution<double> uniform{0.0,1.0};
 std::normal_distribution<double> normal{0.0,1.0};
};

static inline std::uint32_t bed_bayesrc_tempered_seed(
 std::uint64_t ensemble_seed, std::uint64_t stream
) {
 std::uint64_t value = ensemble_seed + 0x9e3779b97f4a7c15ULL * (stream + 1ULL);
 value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
 value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
 value ^= value >> 31U;
 return static_cast<std::uint32_t>(value);
}

template <class PackedGenotype, class AnnotationMatrix, class MarkerMap>
static BedBayesRCChainExecutionResult run_bed_bayesrc_coupling_tempering_chain(
 const BedBayesRCChainExecutionContext<PackedGenotype,AnnotationMatrix,MarkerMap>& context
) {
 validate_bed_bayesrc_chain_context(context);
 const auto& prior = context.coefficient_prior.intercept_prior;
 if (!prior.coupling_tempering || prior.coupling_swap_every <= 0 ||
     prior.legacy_flat || !context.coefficient_prior.update_coefficients ||
     context.coefficient_prior.update_every != 1 ||
     prior.allocation_updates_per_cycle != 1 ||
     prior.annotation_updates_per_cycle != 1) {
  throw std::runtime_error(
   "BayesRC coupling tempering requires a proper intercept prior, alpha updates every cycle, and a 1/1 kernel schedule.");
 }
 const auto& G=context.genotype.storage;
 const auto& maps=context.marker_maps;
 const auto& marker_order=context.marker_order;
 const auto& annotation=context.annotation.matrix;
 const auto& gamma=context.components.scales;
 const int trait=context.trait_index;
 const int m=G.m, K=static_cast<int>(gamma.size());
 const int total_it=context.iterations+context.burnin;
 const std::array<double,3> coupling{{0.0,0.5,1.0}};
 const arma::vec baseline=prior.mean;
 const arma::vec phenotype=context.phenotype.col(static_cast<arma::uword>(trait));

 BedBayesRCChainExecutionResult out;
 out.coupling_tempering=true;
 out.bm.zeros(m); out.dm.zeros(m); out.component_mean.zeros(m);
 out.b.zeros(m); out.state.zeros(m);
 out.vbs.zeros(total_it); out.vgs.zeros(total_it); out.ves.zeros(total_it);
 out.vles.zeros(total_it); out.vlds.zeros(total_it); out.pis.zeros(total_it);
 out.comp_prob.zeros(m,K); out.mean_prior.zeros(K);
 if (context.convergence_annotations) {
  out.convergence_alpha.zeros(context.iterations,annotation.n_cols*(K-1));
  out.convergence_sigma.zeros(context.iterations,K-1);
 }
 if (context.convergence_b)
  out.convergence_b.zeros(context.iterations,context.convergence_markers.size());
 if (context.convergence_d)
  out.convergence_d.zeros(context.iterations,context.convergence_markers.size());
 if (context.convergence_component)
  out.convergence_component.zeros(context.iterations,context.convergence_markers.size());
 if (context.convergence_component)
  allocate_aggregate_component_trace(
   out.convergence_aggregate,context.iterations,K);
 out.coupling_replica_identity.zeros(context.iterations,3);
 out.coupling_active_count.zeros(context.iterations,3);
 out.coupling_expected_active.zeros(context.iterations,3);
 const int maximum_attempts=total_it/prior.coupling_swap_every;
 out.coupling_swap.zeros(maximum_attempts,9);

 std::vector<BedBayesRCTemperedReplicaState> replica(3);
 for (int slot=0;slot<3;++slot) {
  auto& state=replica[static_cast<std::size_t>(slot)];
  state.effect.zeros(m); state.component.zeros(m);
  for (int marker=0;marker<m;++marker) {
   state.effect(static_cast<arma::uword>(marker))=
    context.initial_effects[static_cast<std::size_t>(trait)][static_cast<std::size_t>(marker)];
   state.component(static_cast<arma::uword>(marker))=
    state.effect(static_cast<arma::uword>(marker))!=0.0 ? 1 : 0;
  }
  state.residual=phenotype-br_xb(G,maps,marker_order,state.effect);
  state.marker_variance=context.initial_B(trait,trait);
  state.residual_variance=context.initial_E(trait,trait);
  state.genetic_variance=br_computeG(phenotype,state.residual);
  state.adjusted_residual_variance=state.residual_variance+
   context.adjE*state.genetic_variance;
  state.alpha=context.coefficient_prior.initial_alpha;
  state.sigma_sq_alpha=context.coefficient_prior.initial_step_variances;
  state.marker_probability=st_bayesrc_compute_tempered_snp_pi(
   annotation,state.alpha,baseline,coupling[static_cast<std::size_t>(slot)],
   context.pi_floor);
  state.identity=slot;
  state.generator.seed(bed_bayesrc_tempered_seed(context.chain_seed,slot));
 }
 std::mt19937 swap_generator(bed_bayesrc_tempered_seed(context.chain_seed,17));
 std::uniform_real_distribution<double> swap_uniform(0.0,1.0);
 arma::mat alpha_acc(annotation.n_cols,K-1,arma::fill::zeros);
 arma::vec sigma_acc(K-1,arma::fill::zeros);
 arma::vec log_inv_cpo(G.n,arma::fill::value(-std::numeric_limits<double>::infinity()));
 double nsamples=0.0;
 int attempt=0;

 for (int it=0;it<total_it;++it) {
  const auto transition_start=std::chrono::steady_clock::now();
  for (int slot=0;slot<3;++slot) {
   auto& state=replica[static_cast<std::size_t>(slot)];
   for (int marker:marker_order) {
    const arma::uword index=static_cast<arma::uword>(marker);
    double effect=state.effect(index);
    int component=state.component(index);
    sample_marker_bayesrc(
     G,marker,maps[static_cast<std::size_t>(marker)],
     state.marker_probability.row(index),gamma,state.marker_variance,
     state.adjusted_residual_variance,state.residual,effect,component,
     state.generator,state.uniform,state.normal);
    state.effect(index)=effect;
    state.component(index)=component;
   }
   if (context.update_residual_variance && context.rebuild_every>0 &&
       ((it+1)%context.rebuild_every==0))
    state.residual=phenotype-br_xb(G,maps,marker_order,state.effect);
   if (context.update_marker_variance)
    sampleB_bayesr(m,gamma,context.nub,state.marker_variance,state.effect,
     state.component,context.ssb_prior(trait,trait),state.generator);
   if (context.update_residual_variance)
    sampleE_bayesr(context.nue,state.residual_variance,state.residual,
     context.sse_prior(trait,trait),state.generator);
   st_bayesrc_update_tempered_annotation_prior(
    annotation,state.component,state.alpha,state.sigma_sq_alpha,prior,
    context.coefficient_prior.inverse_chisq_df,
    context.coefficient_prior.inverse_chisq_scale,baseline,
    coupling[static_cast<std::size_t>(slot)],state.generator);
   state.marker_probability=st_bayesrc_compute_tempered_snp_pi(
    annotation,state.alpha,baseline,coupling[static_cast<std::size_t>(slot)],
    context.pi_floor);
   state.genetic_variance=br_computeG(phenotype,state.residual);
   state.adjusted_residual_variance=state.residual_variance+
    context.adjE*state.genetic_variance;
  }
  out.coupling_transition_seconds+=std::chrono::duration<double>(
   std::chrono::steady_clock::now()-transition_start).count();

  if ((it+1)%prior.coupling_swap_every==0) {
   const auto swap_start=std::chrono::steady_clock::now();
   const int lower=attempt%2;
   const int upper=lower+1;
   const auto& lower_state=replica[static_cast<std::size_t>(lower)];
   const auto& upper_state=replica[static_cast<std::size_t>(upper)];
   const double log_ratio=
    st_bayesrc_tempered_log_allocation_prior(
     annotation,upper_state.alpha,baseline,coupling[static_cast<std::size_t>(lower)],
     context.pi_floor,upper_state.component)+
    st_bayesrc_tempered_log_allocation_prior(
     annotation,lower_state.alpha,baseline,coupling[static_cast<std::size_t>(upper)],
     context.pi_floor,lower_state.component)-
    st_bayesrc_tempered_log_allocation_prior(
     annotation,lower_state.alpha,baseline,coupling[static_cast<std::size_t>(lower)],
     context.pi_floor,lower_state.component)-
    st_bayesrc_tempered_log_allocation_prior(
     annotation,upper_state.alpha,baseline,coupling[static_cast<std::size_t>(upper)],
     context.pi_floor,upper_state.component);
   if (!std::isfinite(log_ratio))
    throw std::runtime_error("BayesRC coupling-tempering swap ratio is non-finite.");
   const double probability=std::exp(std::min(0.0,log_ratio));
   const int lower_identity=lower_state.identity;
   const int upper_identity=upper_state.identity;
   const double target_before=arma::accu(replica[2].component>0);
   const bool accepted=swap_uniform(swap_generator)<probability;
   if (accepted) {
    std::swap(replica[static_cast<std::size_t>(lower)],
              replica[static_cast<std::size_t>(upper)]);
    for (int slot:std::array<int,2>{{lower,upper}}) {
     auto& state=replica[static_cast<std::size_t>(slot)];
     state.marker_probability=st_bayesrc_compute_tempered_snp_pi(
      annotation,state.alpha,baseline,coupling[static_cast<std::size_t>(slot)],
      context.pi_floor);
    }
   }
   out.coupling_swap(attempt,0)=it+1;
   out.coupling_swap(attempt,1)=lower;
   out.coupling_swap(attempt,2)=accepted ? 1.0 : 0.0;
   out.coupling_swap(attempt,3)=probability;
   out.coupling_swap(attempt,4)=log_ratio;
   out.coupling_swap(attempt,5)=lower_identity;
   out.coupling_swap(attempt,6)=upper_identity;
   out.coupling_swap(attempt,7)=target_before;
   out.coupling_swap(attempt,8)=arma::accu(replica[2].component>0);
   ++attempt;
   out.coupling_swap_seconds+=std::chrono::duration<double>(
    std::chrono::steady_clock::now()-swap_start).count();
  }

  auto& target=replica[2];
  target.genetic_variance=br_computeG(phenotype,target.residual);
  target.adjusted_residual_variance=target.residual_variance+
   context.adjE*target.genetic_variance;
  out.vbs(static_cast<arma::uword>(it))=target.marker_variance;
  out.vgs(static_cast<arma::uword>(it))=target.genetic_variance;
  out.ves(static_cast<arma::uword>(it))=target.residual_variance;
  const double vle=br_computeLE(target.effect,maps,G.n);
  out.vles(static_cast<arma::uword>(it))=vle;
  out.vlds(static_cast<arma::uword>(it))=target.genetic_variance-vle;
  out.pis(static_cast<arma::uword>(it))=
   1.0-arma::mean(target.marker_probability.col(0));
  if (it>=context.burnin) {
   const arma::uword draw=static_cast<arma::uword>(it-context.burnin);
   for (int slot=0;slot<3;++slot) {
    const auto& state=replica[static_cast<std::size_t>(slot)];
    out.coupling_replica_identity(draw,slot)=state.identity;
    out.coupling_active_count(draw,slot)=arma::accu(state.component>0);
    out.coupling_expected_active(draw,slot)=arma::accu(
     1.0-state.marker_probability.col(0));
   }
   if (context.convergence_annotations) {
    arma::uword q=0;
    for (int stick=0;stick<K-1;++stick)
     for (arma::uword annotation_index=0;
          annotation_index<annotation.n_cols;++annotation_index)
      out.convergence_alpha(draw,q++)=target.alpha(
       annotation_index,static_cast<arma::uword>(stick));
    out.convergence_sigma.row(draw)=target.sigma_sq_alpha.t();
   }
   for (std::size_t q=0;q<context.convergence_markers.size();++q) {
    const arma::uword marker=static_cast<arma::uword>(context.convergence_markers[q]);
    const int component=target.component(marker);
    if (context.convergence_b) out.convergence_b(draw,q)=target.effect(marker);
    if (context.convergence_d) out.convergence_d(draw,q)=component>0 ? 1 : 0;
    if (context.convergence_component)
     out.convergence_component(draw,q)=component;
   }
   if (context.convergence_component)
    capture_aggregate_component_trace(
     target.component,static_cast<int>(m),K,static_cast<int>(draw),
     out.convergence_aggregate);
  }
  if (it>=context.burnin && ((it-context.burnin)%context.thinning==0)) {
   nsamples+=1.0;
   br_update_log_inv_cpo(target.residual,target.residual_variance,log_inv_cpo);
   for (int marker=0;marker<m;++marker) {
    const arma::uword index=static_cast<arma::uword>(marker);
    out.bm(index)+=target.effect(index);
    out.dm(index)+=target.component(index)>0 ? 1.0 : 0.0;
    out.component_mean(index)+=target.component(index);
    out.comp_prob(index,static_cast<arma::uword>(target.component(index)))+=1.0;
   }
   alpha_acc+=target.alpha;
   sigma_acc+=target.sigma_sq_alpha;
   out.mean_prior+=arma::mean(target.marker_probability,0);
  }
 }
 if (nsamples<=0.0) nsamples=1.0;
 out.bm/=nsamples; out.dm/=nsamples; out.component_mean/=nsamples;
 out.comp_prob/=nsamples; alpha_acc/=nsamples; sigma_acc/=nsamples;
 out.mean_prior/=nsamples;
 const auto& target=replica[2];
 out.b=target.effect;
 for (int marker=0;marker<m;++marker) out.state(marker)=target.component(marker);
 out.annot_alpha_mean=alpha_acc; out.annot_alpha_final=target.alpha;
 out.annot_sigma_mean=sigma_acc; out.annot_sigma_final=target.sigma_sq_alpha;
 out.residual=target.residual; out.final_vb=target.marker_variance;
 out.final_vg=target.genetic_variance; out.final_ve=target.residual_variance;
 out.nsamples=nsamples;
 out.log_cpo=br_compute_total_log_cpo(log_inv_cpo,nsamples);
 out.mean_log_cpo=out.log_cpo/G.n;
 return out;
}

template <class PackedGenotype, class AnnotationMatrix, class MarkerMap>
static BedBayesRCChainExecutionResult run_bed_bayesrc_chain(
 const BedBayesRCChainExecutionContext<PackedGenotype,AnnotationMatrix,MarkerMap>& context
) {
 if (context.coefficient_prior.intercept_prior.coupling_tempering) {
  try {
   return run_bed_bayesrc_coupling_tempering_chain(context);
  } catch (const std::exception& exception) {
   BedBayesRCChainExecutionResult out;
   out.failed=1; out.error=exception.what();
   return out;
  }
 }
 const auto& G=context.genotype.storage;
 const auto& maps=context.marker_maps;
 const auto& marker_order=context.marker_order;
 const auto& y=context.phenotype;
 const auto& b_init=context.initial_effects;
 const auto& B=context.initial_B;
 const auto& E=context.initial_E;
 const auto& ssb_prior=context.ssb_prior;
 const auto& sse_prior=context.sse_prior;
 const auto& annotation=context.annotation.matrix;
 const auto& gamma=context.components.scales;
 const auto& annot_alpha_init=context.coefficient_prior.initial_alpha;
 const auto& annot_sigma_init=context.coefficient_prior.initial_step_variances;
 const auto& intercept_prior=context.coefficient_prior.intercept_prior;
 const double sigma_a=context.coefficient_prior.inverse_chisq_df;
 const double sigma_b=context.coefficient_prior.inverse_chisq_scale;
 const double pi_floor=context.pi_floor;
 const double nub=context.nub, nue=context.nue;
 const bool update_alpha=context.coefficient_prior.update_coefficients;
 const bool update_b=context.update_marker_variance;
 const bool update_e=context.update_residual_variance;
 const int alpha_every=context.coefficient_prior.update_every;
 const double adjE=context.adjE;
 const int nit=context.iterations, nburn=context.burnin, nthin=context.thinning;
 const int rebuild_every=context.rebuild_every;
 const int trait=context.trait_index;
 const int m = G.m, K = static_cast<int>(gamma.size());
 const int total_it = nit + nburn;
 BedBayesRCChainExecutionResult out;
 out.bm.zeros(m); out.dm.zeros(m); out.component_mean.zeros(m);
 out.b.zeros(m); out.state.zeros(m);
 out.vbs.zeros(total_it); out.vgs.zeros(total_it); out.ves.zeros(total_it);
 out.vles.zeros(total_it); out.vlds.zeros(total_it); out.pis.zeros(total_it);
 out.comp_prob.zeros(m, K); out.mean_prior.zeros(K);
 if (context.convergence_annotations) {
  out.convergence_alpha.zeros(nit,annotation.n_cols*(K-1));
  out.convergence_sigma.zeros(nit,K-1);
 }
 if (context.convergence_b)
  out.convergence_b.zeros(nit,context.convergence_markers.size());
 if (context.convergence_d)
  out.convergence_d.zeros(nit,context.convergence_markers.size());
 if (context.convergence_component)
  out.convergence_component.zeros(nit,context.convergence_markers.size());
 if (context.convergence_component)
  allocate_aggregate_component_trace(
   out.convergence_aggregate,nit,K);
 try {
  validate_bed_bayesrc_chain_context(context);
  std::mt19937 gen(static_cast<unsigned int>(context.chain_seed));
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  std::normal_distribution<double> norm01(0.0, 1.0);
  arma::vec y_t = y.col(static_cast<arma::uword>(trait));
  arma::rowvec b_t(m, arma::fill::zeros);
  arma::Row<int> component_t(m, arma::fill::zeros);
  for (int j = 0; j < m; ++j) {
   b_t(static_cast<arma::uword>(j)) =
    b_init[static_cast<std::size_t>(trait)][static_cast<std::size_t>(j)];
   component_t(static_cast<arma::uword>(j)) = b_t(static_cast<arma::uword>(j)) != 0.0 ? 1 : 0;
  }
  arma::vec residual = y_t - br_xb(G, maps, marker_order, b_t);
  double vb = B(trait, trait), ve = E(trait, trait);
  double vg = br_computeG(y_t, residual), vei = ve + adjE * vg;
  arma::mat annot_alpha = annot_alpha_init;
  arma::vec annot_sigma = annot_sigma_init;
  arma::mat snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
  arma::mat alpha_acc(annotation.n_cols, K - 1, arma::fill::zeros);
  arma::vec sigma_acc(K - 1, arma::fill::zeros);
  arma::vec log_inv_cpo(G.n, arma::fill::value(-std::numeric_limits<double>::infinity()));
  double nsamples = 0.0;
  const int allocation_updates = intercept_prior.allocation_updates_per_cycle;
  const int annotation_updates = intercept_prior.annotation_updates_per_cycle;
  const bool legacy_schedule = allocation_updates == 1 && annotation_updates == 1;

  for (int it = 0; it < total_it; ++it) {
   if (legacy_schedule) {
    // Preserve the historical 1/1 transition and RNG order exactly.
    for (int marker : marker_order) {
     const arma::uword ju = static_cast<arma::uword>(marker);
     double effect = b_t(ju);
     int component = component_t(ju);
     sample_marker_bayesrc(
      G, marker, maps[static_cast<std::size_t>(marker)], snp_pi.row(ju), gamma,
      vb, vei, residual, effect, component, gen, runif, norm01
     );
     b_t(ju) = effect;
     component_t(ju) = component;
    }
    if (update_e && rebuild_every > 0 && ((it + 1) % rebuild_every == 0))
     residual = y_t - br_xb(G, maps, marker_order, b_t);
    if (update_b)
     sampleB_bayesr(m, gamma, nub, vb, b_t, component_t, ssb_prior(trait, trait), gen);
    if (update_e)
     sampleE_bayesr(nue, ve, residual, sse_prior(trait, trait), gen);
    if (update_alpha && ((it + 1) % alpha_every == 0)) {
     st_bayesrc_update_annotation_prior(
      annotation, component_t, annot_alpha, annot_sigma, intercept_prior,
      sigma_a, sigma_b, gen
     );
     snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
    }
   } else {
    for (int allocation_rep = 0; allocation_rep < allocation_updates;
         ++allocation_rep) {
     for (int marker : marker_order) {
      const arma::uword ju = static_cast<arma::uword>(marker);
      double effect = b_t(ju);
      int component = component_t(ju);
      sample_marker_bayesrc(
       G, marker, maps[static_cast<std::size_t>(marker)], snp_pi.row(ju), gamma,
       vb, vei, residual, effect, component, gen, runif, norm01
      );
      b_t(ju) = effect;
      component_t(ju) = component;
     }
     const int allocation_index = it * allocation_updates + allocation_rep + 1;
     if (update_e && rebuild_every > 0 && allocation_index % rebuild_every == 0)
      residual = y_t - br_xb(G, maps, marker_order, b_t);
     if (update_b)
      sampleB_bayesr(m, gamma, nub, vb, b_t, component_t, ssb_prior(trait, trait), gen);
     if (update_e)
      sampleE_bayesr(nue, ve, residual, sse_prior(trait, trait), gen);
     vg = br_computeG(y_t, residual);
     vei = ve + adjE * vg;
    }
    if (update_alpha && ((it + 1) % alpha_every == 0)) {
     for (int annotation_rep = 0; annotation_rep < annotation_updates;
          ++annotation_rep) {
      st_bayesrc_update_annotation_prior(
       annotation, component_t, annot_alpha, annot_sigma, intercept_prior,
       sigma_a, sigma_b, gen
      );
      snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
     }
    }
   }
   vg = br_computeG(y_t, residual);
   const double vle = br_computeLE(b_t, maps, G.n);
   vei = ve + adjE * vg;
   out.vbs(static_cast<arma::uword>(it)) = vb;
   out.vgs(static_cast<arma::uword>(it)) = vg;
   out.ves(static_cast<arma::uword>(it)) = ve;
   out.vles(static_cast<arma::uword>(it)) = vle;
   out.vlds(static_cast<arma::uword>(it)) = vg - vle;
   out.pis(static_cast<arma::uword>(it)) = 1.0 - arma::mean(snp_pi.col(0));
   if (it>=nburn) {
    const arma::uword draw=static_cast<arma::uword>(it-nburn);
    if (context.convergence_annotations) {
     arma::uword q=0;
     for (int stick=0; stick<K-1; ++stick)
      for (arma::uword annotation_index=0;
           annotation_index<annotation.n_cols; ++annotation_index)
       out.convergence_alpha(draw,q++)=annot_alpha(
        annotation_index,static_cast<arma::uword>(stick));
     out.convergence_sigma.row(draw)=annot_sigma.t();
    }
    for (std::size_t q=0;q<context.convergence_markers.size();++q) {
     const arma::uword marker=static_cast<arma::uword>(
      context.convergence_markers[q]);
     const int component=component_t(marker);
     if (context.convergence_b) out.convergence_b(draw,q)=b_t(marker);
     if (context.convergence_d) out.convergence_d(draw,q)=component>0 ? 1 : 0;
     if (context.convergence_component)
      out.convergence_component(draw,q)=component;
    }
    if (context.convergence_component)
     capture_aggregate_component_trace(
      component_t,static_cast<int>(m),K,static_cast<int>(draw),
      out.convergence_aggregate);
   }
   if (it >= nburn && ((it - nburn) % nthin == 0)) {
    nsamples += 1.0;
    br_update_log_inv_cpo(residual, ve, log_inv_cpo);
    for (int j = 0; j < m; ++j) {
     const arma::uword ju = static_cast<arma::uword>(j);
     out.bm(ju) += b_t(ju);
     out.dm(ju) += component_t(ju) > 0 ? 1.0 : 0.0;
     out.component_mean(ju) += component_t(ju);
     out.comp_prob(ju, static_cast<arma::uword>(component_t(ju))) += 1.0;
    }
    alpha_acc += annot_alpha;
    sigma_acc += annot_sigma;
    out.mean_prior += arma::mean(snp_pi, 0);
   }
  }
  if (nsamples <= 0.0) nsamples = 1.0;
  out.bm /= nsamples; out.dm /= nsamples; out.component_mean /= nsamples;
  out.comp_prob /= nsamples; alpha_acc /= nsamples; sigma_acc /= nsamples;
  out.mean_prior /= nsamples;
  out.b = b_t;
  for (int j = 0; j < m; ++j) out.state(j) = component_t(j);
  out.annot_alpha_mean = alpha_acc; out.annot_alpha_final = annot_alpha;
  out.annot_sigma_mean = sigma_acc; out.annot_sigma_final = annot_sigma;
  out.residual = residual; out.final_vb = vb; out.final_vg = vg; out.final_ve = ve;
  out.nsamples = nsamples;
  out.log_cpo = br_compute_total_log_cpo(log_inv_cpo, nsamples);
  out.mean_log_cpo = out.log_cpo / G.n;
 } catch (const std::exception& ex) { out.failed = 1; out.error = ex.what(); }
 return out;
}

} }

#endif
