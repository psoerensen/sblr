#ifndef SBLR_BLR_BED_BAYESRC_CORE_IMPL_H
#define SBLR_BLR_BED_BAYESRC_CORE_IMPL_H

#include "blr_bed_bayesrc_types.h"
#include "st_bayesrc_annotation_prior.h"

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

template <class PackedGenotype, class AnnotationMatrix, class MarkerMap>
static BedBayesRCChainExecutionResult run_bed_bayesrc_chain(
 const BedBayesRCChainExecutionContext<PackedGenotype,AnnotationMatrix,MarkerMap>& context
) {
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
 const bool intercept_flat=context.coefficient_prior.intercept_flat;
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

  for (int it = 0; it < total_it; ++it) {
   // Exact full sweep: every marker is visited once in marker_order every iteration.
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
     annotation, component_t, annot_alpha, annot_sigma, intercept_flat,
     sigma_a, sigma_b, gen
    );
    snp_pi = st_bayesrc_compute_snp_pi(annotation, annot_alpha, pi_floor);
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
