#ifndef SBLR_BLR_CSR_SBAYESRC_CORE_IMPL_H
#define SBLR_BLR_CSR_SBAYESRC_CORE_IMPL_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "blr_aggregate_component_trace.h"

// Implementation detail only. This header is included by
// st_sbayesrc_omp_csr.cpp after the package's Armadillo configuration and the
// concrete operator types have been established. All borrowed objects in the
// context must outlive run_csr_sbayesrc(); this is the native realization of
// the canonical storage_outlives_execution ownership contract.

namespace sblr { namespace core {

struct CsrSBayesRCFailureState {
 bool captured = false;
 int trait = -1, chain = -1, iteration = -1, internal_iteration = -1;
 int sample_size = 0, active_marker_count = 0;
 double residual_df = 0.0, yy = 0.0, sse_prior = 0.0;
 double prior_contribution = 0.0, maintained_sse = 0.0;
 double rebuilt_sse = 0.0, quadratic_sse = 0.0;
 double maintained_scale = 0.0, rebuilt_scale = 0.0, quadratic_scale = 0.0;
 double b_wy = 0.0, b_r = 0.0, fitted_quadratic = 0.0;
 double independent_quadratic = 0.0, marker_variance = 0.0;
 double genetic_variance = 0.0, genic_variance = 0.0, ld_variance = 0.0;
 double residual_variance = 0.0, adjusted_residual_variance = 0.0;
 double heritability = 0.0, maximum_absolute_effect = 0.0;
 double effect_norm = 0.0, residual_norm = 0.0, rebuilt_residual_norm = 0.0;
 double maximum_residual_drift = 0.0, relative_residual_drift = 0.0;
 bool effects_finite = false, residual_finite = false, rebuilt_residual_finite = false;
 std::string operator_name;
 arma::rowvec effects, residual, rebuilt_residual, score, diagonal;
 arma::Row<int> component;
 std::vector<double> preceding_marker_variance, preceding_genetic_variance;
 std::vector<double> preceding_residual_variance, preceding_heritability;
 std::vector<double> preceding_maximum_effect;
};

inline void write_csr_sbayesrc_failure_state(
 const CsrSBayesRCFailureState& x,
 const std::string& path
) {
 if (!x.captured || path.empty()) return;
 std::ofstream out(path.c_str(), std::ios::out | std::ios::trunc);
 if (!out) throw std::runtime_error(
  "could not open SBayesRC failure-state diagnostic path: " + path);
 out << std::setprecision(17);
 const auto scalar = [&out](const char* key, const auto& value) {
  out << "meta\t" << key << "\t" << value << '\n';
 };
 scalar("operator", x.operator_name);
 scalar("trait", x.trait); scalar("chain", x.chain);
 scalar("iteration", x.iteration); scalar("internal_iteration", x.internal_iteration);
 scalar("sample_size", x.sample_size); scalar("residual_df", x.residual_df);
 scalar("yy", x.yy); scalar("sse_prior", x.sse_prior);
 scalar("prior_contribution", x.prior_contribution);
 scalar("maintained_sse", x.maintained_sse); scalar("rebuilt_sse", x.rebuilt_sse);
 scalar("quadratic_sse", x.quadratic_sse);
 scalar("maintained_scale", x.maintained_scale); scalar("rebuilt_scale", x.rebuilt_scale);
 scalar("quadratic_scale", x.quadratic_scale);
 scalar("b_wy", x.b_wy); scalar("b_r", x.b_r);
 scalar("fitted_quadratic", x.fitted_quadratic);
 scalar("independent_quadratic", x.independent_quadratic);
 scalar("marker_variance", x.marker_variance);
 scalar("genetic_variance", x.genetic_variance);
 scalar("genic_variance", x.genic_variance); scalar("ld_variance", x.ld_variance);
 scalar("residual_variance", x.residual_variance);
 scalar("adjusted_residual_variance", x.adjusted_residual_variance);
 scalar("heritability", x.heritability);
 scalar("active_marker_count", x.active_marker_count);
 scalar("maximum_absolute_effect", x.maximum_absolute_effect);
 scalar("effect_norm", x.effect_norm); scalar("residual_norm", x.residual_norm);
 scalar("rebuilt_residual_norm", x.rebuilt_residual_norm);
 scalar("maximum_residual_drift", x.maximum_residual_drift);
 scalar("relative_residual_drift", x.relative_residual_drift);
 scalar("effects_finite", static_cast<int>(x.effects_finite));
 scalar("residual_finite", static_cast<int>(x.residual_finite));
 scalar("rebuilt_residual_finite", static_cast<int>(x.rebuilt_residual_finite));
 const int maximum_component = x.component.is_empty() ? -1 : x.component.max();
 for (int k = 0; k <= maximum_component; ++k)
  scalar(("component_count_" + std::to_string(k)).c_str(),
         arma::accu(x.component == k));
 const auto trajectory = [&out](const char* key, const std::vector<double>& value) {
  for (std::size_t i = 0; i < value.size(); ++i)
   out << "trajectory\t" << key << '\t' << i << '\t' << value[i] << '\n';
 };
 trajectory("marker_variance", x.preceding_marker_variance);
 trajectory("genetic_variance", x.preceding_genetic_variance);
 trajectory("residual_variance", x.preceding_residual_variance);
 trajectory("heritability", x.preceding_heritability);
 trajectory("maximum_absolute_effect", x.preceding_maximum_effect);
 out << "marker\tindex\teffect\tresidual\trebuilt_residual\tscore\tdiagonal\tcomponent\n";
 for (arma::uword i = 0; i < x.effects.n_elem; ++i)
  out << "marker\t" << (i + 1) << '\t' << x.effects(i) << '\t' << x.residual(i)
      << '\t' << x.rebuilt_residual(i) << '\t' << x.score(i) << '\t'
      << x.diagonal(i) << '\t' << x.component(i) << '\n';
 if (!out) throw std::runtime_error(
  "failed while writing SBayesRC failure-state diagnostic: " + path);
}

template <class Operator>
struct CsrSBayesRCExecutionContext {
 const Operator& op;
 const SBayesRCLDLDFriends& ld_swap_friends;
 const arma::mat& wy_mat;
 const arma::mat& ww_mat;
 arma::mat& b_mat;
 arma::mat& r_mat;
 arma::mat& comp_mat;
 const arma::vec& yy_vec;
 const arma::mat& ssb_prior_mat;
 const arma::mat& sse_prior_mat;
 const arma::mat& B;
 const arma::mat& E;
 const arma::mat& A;
 const arma::vec& gamma;
 const arma::mat& alpha_init;
 const arma::vec& sigmaSqAlpha_init;
 const std::vector<std::vector<double>>& comp_init;
 const std::vector<std::vector<double>>& r_init;
 const std::vector<int>& sample_size;
 const std::vector<int>& chain_seeds;
 const std::vector<int>& marker_order;
 const arma::rowvec& prior_scale;
 const arma::rowvec& maf_effect_s_log_h;

 // Contract views describe the same borrowed prepared inputs. The
 // numerical body continues to use the native Armadillo/operator references.
 const SBayesRCAnnotationDesignView& annotation_contract;
 const SBayesRCComponentSpec& component_contract;
 const SBayesRCAlphaSpec& alpha_contract;
 const SBayesRCProbabilityPolicy& probability_contract;
 const SBayesRCPriors& prior_contract;
 const SBayesRCControls& control_contract;
 const SBayesRCOutputControl& output_contract;
 const StBayesRCSelectionGenomicConfig& selection_config;

 int marker_count;
 int trait_count;
 int annotation_count;
 int component_count;
 int step_count;
 int iterations;
 int burnin;
 int thinning;
 int cores;
 int seed;
 int chains;
 bool use_component_initialization;
 bool use_residual_initialization;
 bool rebuild_residual_before_update;
 int low_rank_residual_rebuild_every;
 StBayesRCInterceptPrior intercept_prior;
 double alpha_variance_prior_a;
 double alpha_variance_prior_b;
 double probability_floor;
 double marker_variance_df;
 double residual_variance_df;
 bool update_alpha;
 bool update_marker_variance;
 bool update_residual_variance;
 int alpha_update_every;
 double residual_adjustment;
 bool update_ld_swap;
 double ld_swap_probability;
 int ld_swap_moves;
 bool use_selection_prior_scale;
 bool estimate_maf_effect_s;
 double maf_effect_s_initial;
 double maf_effect_s_prior_lower;
 double maf_effect_s_prior_upper;
 double maf_effect_s_proposal_sd;
 const std::vector<int>& convergence_markers;
 const std::vector<double>& phenotype_variance;
 BlockResidualControl block_residual_control;
 bool convergence_annotations;
 bool convergence_b;
 bool convergence_d;
 bool convergence_component;
};

struct CsrSBayesRCExecutionResult {
 arma::mat bm_task, dm_task, b_task, r_task, comp_task_double;
 arma::mat vbs_task, vgs_task, ves_task, pis_task, vles_task, vlds_task;
 arma::mat maf_effect_s_task;
 arma::vec final_vb_task, final_vg_task, final_ve_task;
 arma::vec final_pi_active_task, final_vle_task, final_vld_task;
 arma::mat final_pi_component_task;
 arma::vec nsamples_task, ld_swap_attempted_task, ld_swap_accepted_task;
 arma::vec low_rank_residual_rebuild_count_task;
 arma::vec low_rank_residual_max_abs_drift_task;
 arma::vec maf_effect_s_attempted_task, maf_effect_s_accepted_task;
 std::vector<arma::mat> alpha_mean_task, comp_prob_mean_task;
 std::vector<arma::vec> sigmaSqAlpha_mean_task, ncomp_mean_task;
 std::vector<arma::mat> convergence_alpha_task, convergence_sigma_task;
 std::vector<arma::mat> convergence_b_task;
 std::vector<arma::imat> convergence_d_task, convergence_component_task;
 std::vector<AggregateComponentTrace> convergence_aggregate_task;
 std::vector<arma::vec> selection_pip_task, selection_pi_a_mean_task;
 std::vector<arma::vec> selection_tau2_mean_task, selection_included_mean_task;
 std::vector<arma::vec> selection_switches_task;
 std::vector<arma::mat> selection_alpha_conditional_mean_task;
 std::vector<arma::uvec> selection_delta_final_task;
 std::vector<arma::mat> selection_empty_stick_diagnostics_task;
 std::vector<arma::umat> selection_delta_trace_task;
 std::vector<arma::vec> selection_pi_a_trace_task;
 std::vector<arma::uvec> selection_included_trace_task;

 arma::mat bm_mat, dm_mat, bm_sd_mat, dm_sd_mat;
 arma::mat bm_min_mat, dm_min_mat, bm_max_mat, dm_max_mat;
 arma::mat b_mat, r_mat, comp_mat;
 arma::mat vbs_mat, vgs_mat, ves_mat, pis_mat, vles_mat, vlds_mat;
 arma::mat maf_effect_s_mat;
 arma::vec final_vb, final_vg, final_ve, final_pi_active;
 arma::mat final_pi_component;
 arma::vec final_vle, final_vld, nsamples_vec;
 arma::vec maf_effect_s_attempted_vec, maf_effect_s_accepted_vec;
 std::vector<arma::mat> alpha_mean, comp_prob_mean;
 std::vector<arma::vec> sigmaSqAlpha_mean, ncomp_mean;
 std::vector<int> failed, thread_used;
 std::vector<std::string> errors;
 std::vector<double> task_seconds;
 arma::mat ld_swap_diagnostics, ld_swap_chain_diagnostics;
 arma::mat low_rank_residual_diagnostics, low_rank_residual_chain_diagnostics;
 arma::mat block_ve_posterior_mean_task, block_ve_final_task;
 arma::mat block_ve_resampled_task, block_ve_reset_task;
 arma::mat summary_heritability_task;
 std::vector<arma::mat> block_ve_history_task;
};

template <class Operator>
CsrSBayesRCExecutionResult run_csr_sbayesrc(
 CsrSBayesRCExecutionContext<Operator>& context
) {
 const Operator& op = context.op;
 const SBayesRCLDLDFriends& ld_swap_friends = context.ld_swap_friends;
 const arma::mat& wy_mat = context.wy_mat;
 arma::mat& b_mat = context.b_mat;
 arma::mat& r_mat = context.r_mat;
 arma::mat& comp_mat = context.comp_mat;
 const arma::vec& yy_vec = context.yy_vec;
 const arma::mat& ssb_prior_mat = context.ssb_prior_mat;
 const arma::mat& sse_prior_mat = context.sse_prior_mat;
 const arma::mat& B = context.B;
 const arma::mat& E = context.E;
 const arma::mat& A = context.A;
 const arma::vec& gamma = context.gamma;
 const arma::mat& alpha_init = context.alpha_init;
 const arma::vec& sigmaSqAlpha_init = context.sigmaSqAlpha_init;
 const std::vector<std::vector<double>>& comp_init = context.comp_init;
 const std::vector<std::vector<double>>& r_init = context.r_init;
 const std::vector<int>& n = context.sample_size;
 const std::vector<int>& chain_seeds_vec = context.chain_seeds;
 const std::vector<int>& order = context.marker_order;
 const std::vector<double>& phenotype_variance = context.phenotype_variance;
 const BlockResidualControl& block_residual_control =
  context.block_residual_control;
 const StBayesRCSelectionGenomicConfig& selection_config =
  context.selection_config;
 const arma::rowvec& prior_scale = context.prior_scale;
 const arma::rowvec& maf_effect_s_log_h_row = context.maf_effect_s_log_h;
 const int m = context.marker_count;
 const int nt = context.trait_count;
 const int nAnno = context.annotation_count;
 const int Kgamma = context.component_count;
 const int nstep = context.step_count;
 const int nit = context.iterations;
 const int nburn = context.burnin;
 const int nthin = context.thinning;
 const int ncores = context.cores;
 const int seed = context.seed;
 const int nchains = context.chains;
 const bool use_comp_init = context.use_component_initialization;
 const bool use_r_init = context.use_residual_initialization;
 const bool rebuild_r_before_updateE = context.rebuild_residual_before_update;
 const int low_rank_residual_rebuild_every =
  context.low_rank_residual_rebuild_every;
 const auto& intercept_prior = context.intercept_prior;
 const double sigmaSqAlpha_a = context.alpha_variance_prior_a;
 const double sigmaSqAlpha_b = context.alpha_variance_prior_b;
 const double pi_floor = context.probability_floor;
 const double nub = context.marker_variance_df;
 const double nue = context.residual_variance_df;
 const bool updateAlpha = context.update_alpha;
 const bool updateB = context.update_marker_variance;
 const bool updateE = context.update_residual_variance;
 const int alpha_update_every = context.alpha_update_every;
 const double adjE = context.residual_adjustment;
 const bool updateLDswap = context.update_ld_swap;
 const double ld_swap_prob = context.ld_swap_probability;
 const int ld_swap_moves = context.ld_swap_moves;
 const bool use_maf_effect_s_prior_scale = context.use_selection_prior_scale;
 const bool estimate_maf_effect_s = context.estimate_maf_effect_s;
 const double maf_effect_s_init = context.maf_effect_s_initial;
 const double maf_effect_s_prior_lower = context.maf_effect_s_prior_lower;
 const double maf_effect_s_prior_upper = context.maf_effect_s_prior_upper;
 const double maf_effect_s_proposal_sd = context.maf_effect_s_proposal_sd;
 const std::vector<int>& convergence_markers=context.convergence_markers;
 const int ntasks = stblr_num_chain_tasks(nt, nchains);
 const int block_count = op.block_count();
 if (block_residual_control.uses_block_variance() &&
     (block_count <= 0 || static_cast<int>(phenotype_variance.size()) != nt))
  throw std::runtime_error(
   "block residual policy requires retained blocks and one phenotype variance per trait.");
 const char* failure_path_env = std::getenv("SBLR_SBAYESRC_FAILURE_STATE_PATH");
 const std::string failure_state_path = failure_path_env == nullptr
  ? std::string()
  : std::string(failure_path_env);
 const bool capture_failure_state = !failure_state_path.empty();

 arma::mat bm_task(ntasks, m, arma::fill::zeros);
 arma::mat dm_task(ntasks, m, arma::fill::zeros);
 arma::mat b_task(ntasks, m, arma::fill::zeros);
 arma::mat r_task(ntasks, m, arma::fill::zeros);
 arma::mat comp_task_double(ntasks, m, arma::fill::zeros);

 arma::mat vbs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vgs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat ves_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat pis_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vles_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vlds_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat maf_effect_s_task(ntasks, nit + nburn, arma::fill::zeros);

 arma::vec final_vb_task(ntasks, arma::fill::zeros);
 arma::vec final_vg_task(ntasks, arma::fill::zeros);
 arma::vec final_ve_task(ntasks, arma::fill::zeros);
 arma::vec final_pi_active_task(ntasks, arma::fill::zeros);
 arma::mat final_pi_component_task(ntasks, Kgamma, arma::fill::zeros);
 arma::vec final_vle_task(ntasks, arma::fill::zeros);
 arma::vec final_vld_task(ntasks, arma::fill::zeros);
 arma::vec nsamples_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_attempted_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_accepted_task(ntasks, arma::fill::zeros);
 arma::vec low_rank_residual_rebuild_count_task(ntasks, arma::fill::zeros);
 arma::vec low_rank_residual_max_abs_drift_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_attempted_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_accepted_task(ntasks, arma::fill::zeros);
 arma::mat block_ve_posterior_mean_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_final_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_resampled_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_reset_task(ntasks, block_count, arma::fill::zeros);
 arma::mat summary_heritability_task(ntasks, nit + nburn, arma::fill::zeros);
 std::vector<arma::mat> block_ve_history_task(static_cast<std::size_t>(ntasks));

 std::vector<arma::mat> alpha_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> sigmaSqAlpha_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> comp_prob_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> ncomp_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_alpha_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_sigma_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_b_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_d_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_component_task(static_cast<std::size_t>(ntasks));
 std::vector<AggregateComponentTrace> convergence_aggregate_task(
  static_cast<std::size_t>(ntasks));
 const int selectable_annotation_count = selection_config.enabled ? nAnno - 1 : 0;
 std::vector<arma::vec> selection_pip_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> selection_pi_a_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> selection_tau2_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> selection_included_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> selection_switches_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> selection_alpha_conditional_mean_task(
  static_cast<std::size_t>(ntasks));
 std::vector<arma::uvec> selection_delta_final_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> selection_empty_stick_diagnostics_task(
  static_cast<std::size_t>(ntasks));
 std::vector<arma::umat> selection_delta_trace_task(
  static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> selection_pi_a_trace_task(
  static_cast<std::size_t>(ntasks));
 std::vector<arma::uvec> selection_included_trace_task(
  static_cast<std::size_t>(ntasks));

 if (selection_config.enabled && selectable_annotation_count <= 0)
  throw std::invalid_argument("SBayesRC-S requires an intercept and at least one selectable annotation");

 for (int marker: convergence_markers) if (marker<0 || marker>=m)
  throw std::invalid_argument("SBayesRC convergence marker index is out of range.");

 for (int task = 0; task < ntasks; ++task) {
  alpha_mean_task[static_cast<std::size_t>(task)] = arma::mat(nAnno, nstep, arma::fill::zeros);
  sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)] = arma::vec(nstep, arma::fill::zeros);
  comp_prob_mean_task[static_cast<std::size_t>(task)] = arma::mat(m, Kgamma, arma::fill::zeros);
  ncomp_mean_task[static_cast<std::size_t>(task)] = arma::vec(Kgamma, arma::fill::zeros);
  if (selection_config.enabled) {
   selection_pip_task[static_cast<std::size_t>(task)] =
    arma::vec(selectable_annotation_count, arma::fill::zeros);
   selection_pi_a_mean_task[static_cast<std::size_t>(task)] =
    arma::vec(1u, arma::fill::zeros);
   selection_tau2_mean_task[static_cast<std::size_t>(task)] =
    arma::vec(nstep, arma::fill::zeros);
   selection_included_mean_task[static_cast<std::size_t>(task)] =
    arma::vec(1u, arma::fill::zeros);
   selection_switches_task[static_cast<std::size_t>(task)] =
    arma::vec(selectable_annotation_count, arma::fill::zeros);
   selection_alpha_conditional_mean_task[static_cast<std::size_t>(task)] =
    arma::mat(selectable_annotation_count, nstep, arma::fill::zeros);
   selection_empty_stick_diagnostics_task[static_cast<std::size_t>(task)] =
    arma::mat(nstep, 4u, arma::fill::zeros);
  }
  if (context.convergence_annotations) {
   convergence_alpha_task[static_cast<std::size_t>(task)]=
    arma::mat(nit,nAnno*nstep,arma::fill::zeros);
   convergence_sigma_task[static_cast<std::size_t>(task)]=
    arma::mat(nit,nstep,arma::fill::zeros);
   if (selection_config.enabled) {
    selection_delta_trace_task[static_cast<std::size_t>(task)] =
     arma::umat(nit, selectable_annotation_count, arma::fill::zeros);
    selection_pi_a_trace_task[static_cast<std::size_t>(task)] =
     arma::vec(nit, arma::fill::zeros);
    selection_included_trace_task[static_cast<std::size_t>(task)] =
     arma::uvec(nit, arma::fill::zeros);
   }
  }
  if (context.convergence_b)
   convergence_b_task[static_cast<std::size_t>(task)]=
    arma::mat(nit,convergence_markers.size(),arma::fill::zeros);
  if (context.convergence_d)
   convergence_d_task[static_cast<std::size_t>(task)]=
    arma::imat(nit,convergence_markers.size(),arma::fill::zeros);
  if (context.convergence_component)
   convergence_component_task[static_cast<std::size_t>(task)]=
    arma::imat(nit,convergence_markers.size(),arma::fill::zeros);
  if (context.convergence_component)
   allocate_aggregate_component_trace(
    convergence_aggregate_task[static_cast<std::size_t>(task)],nit,Kgamma);
 }

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat bm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat dm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat bm_min_mat(nt, m, arma::fill::zeros);
 arma::mat dm_min_mat(nt, m, arma::fill::zeros);
 arma::mat bm_max_mat(nt, m, arma::fill::zeros);
 arma::mat dm_max_mat(nt, m, arma::fill::zeros);

 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat maf_effect_s_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi_active(nt, arma::fill::zeros);
 arma::mat final_pi_component(nt, Kgamma, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec maf_effect_s_attempted_vec(nt, arma::fill::zeros);
 arma::vec maf_effect_s_accepted_vec(nt, arma::fill::zeros);

 std::vector<arma::mat> alpha_mean(static_cast<std::size_t>(nt));
 std::vector<arma::vec> sigmaSqAlpha_mean(static_cast<std::size_t>(nt));
 std::vector<arma::mat> comp_prob_mean(static_cast<std::size_t>(nt));
 std::vector<arma::vec> ncomp_mean(static_cast<std::size_t>(nt));

 for (int t = 0; t < nt; ++t) {
  alpha_mean[static_cast<std::size_t>(t)] = arma::mat(nAnno, nstep, arma::fill::zeros);
  sigmaSqAlpha_mean[static_cast<std::size_t>(t)] = arma::vec(nstep, arma::fill::zeros);
  comp_prob_mean[static_cast<std::size_t>(t)] = arma::mat(m, Kgamma, arma::fill::zeros);
  ncomp_mean[static_cast<std::size_t>(t)] = arma::vec(Kgamma, arma::fill::zeros);
 }

 std::vector<int> failed(static_cast<std::size_t>(ntasks), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(ntasks));
 std::vector<int> thread_used(static_cast<std::size_t>(ntasks), 0);
 std::vector<double> task_seconds(static_cast<std::size_t>(ntasks), 0.0);
 std::vector<CsrSBayesRCFailureState> failure_states(
  static_cast<std::size_t>(ntasks));

 int nthreads = 1;

#ifdef _OPENMP
 omp_set_dynamic(0);
 nthreads = stblr_num_threads_for_tasks(ncores, ntasks);
 omp_set_num_threads(nthreads);

#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int task = 0; task < ntasks; ++task) {
  const int t = stblr_task_trait(task, nchains);
  const int chain = stblr_task_chain(task, nchains);
  const arma::uword task_u = static_cast<arma::uword>(task);
#ifdef _OPENMP
  const double wall_start = omp_get_wtime();
  thread_used[static_cast<std::size_t>(task)] = omp_get_thread_num();
#else
  const double wall_start = 0.0;
  thread_used[static_cast<std::size_t>(task)] = 0;
#endif

  try {
   unsigned int task_seed = 0u;
   if (!chain_seeds_vec.empty()) {
    task_seed = stblr_seed_with_chain_base(
     chain_seeds_vec[static_cast<std::size_t>(chain)],
     t
    );
   } else if (nchains == 1) {
    task_seed = stblr_trait_seed(seed, t);
   } else {
    task_seed = stblr_chain_seed(seed, t, chain);
   }
   std::mt19937 gen_t(task_seed);

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   const arma::rowvec& ww_t = op.diag();

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) {
    b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }

   arma::rowvec r_t(op.residual_size(), arma::fill::zeros);
   arma::Row<int> comp_t(m, arma::fill::zeros);

   if (use_comp_init) {
    for (int i = 0; i < m; ++i) {
     const int k = static_cast<int>(std::round(comp_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(i)]));
     if (k < 0 || k >= Kgamma) {
      throw std::runtime_error("comp_init contains component outside 0..Kgamma-1.");
     }
     comp_t(static_cast<arma::uword>(i)) = k;
    }
   } else {
    for (int i = 0; i < m; ++i) {
     comp_t(static_cast<arma::uword>(i)) =
      (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("stblr_cpg_omp_csr_sbayesrc: r_init contains NaN/Inf.");
    }
   } else {
    op.rebuild(t, wy_t, b_t, r_t);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_operator(op, t, b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   arma::mat alpha_t = alpha_init;
   arma::vec sigmaSqAlpha_t = sigmaSqAlpha_init;
   const arma::mat selection_annotation = selection_config.enabled ?
    A.cols(1u, A.n_cols - 1u) : arma::mat();
   StBayesRCSelectionState selection_state;
   arma::uvec selection_previous_empty(nstep, arma::fill::zeros);
   arma::uvec selection_empty_run(nstep, arma::fill::zeros);
   arma::mat selection_empty_diagnostics(nstep, 4u, arma::fill::zeros);
   if (selection_config.enabled) {
    selection_state.delta = selection_config.delta_init;
    selection_state.alpha = alpha_init;
    selection_state.pi_a = selection_config.pi_a_init;
    selection_state.tau2 = selection_config.tau2_init;
    std::vector<arma::uvec> initial_eligible;
    std::vector<arma::ivec> initial_outcome;
    st_bayesrc_selection_build_observed_sticks(
     comp_t, nstep, initial_eligible, initial_outcome);
    for (int stick = 0; stick < nstep; ++stick) {
     const bool empty = initial_eligible[static_cast<std::size_t>(stick)].n_elem == 0u;
     selection_previous_empty(static_cast<arma::uword>(stick)) = empty ? 1u : 0u;
     if (empty) {
      const arma::uword index = static_cast<arma::uword>(stick);
      selection_empty_diagnostics(index, 0u) = 1.0;
      selection_empty_diagnostics(index, 1u) = 1.0;
      selection_empty_diagnostics(index, 3u) = 1.0;
      selection_empty_run(index) = 1u;
     }
    }
    st_bayesrc_selection_validate(
     selection_annotation, initial_eligible, &initial_outcome,
     selection_state, selection_config.hyper,
     selection_config.intercept_mean,
     selection_config.intercept_precision);
    alpha_t = selection_state.alpha;
    sigmaSqAlpha_t = selection_state.tau2;
   }

   for (int j = 0; j < nstep; ++j) {
    if (!std::isfinite(sigmaSqAlpha_t(static_cast<arma::uword>(j))) ||
        sigmaSqAlpha_t(static_cast<arma::uword>(j)) <= 0.0) {
     throw std::runtime_error("sigmaSqAlpha_init contains invalid value.");
    }
   }

   arma::mat snpPi_t = selection_config.enabled ?
    st_bayesrc_selection_compute_pi(selection_annotation, alpha_t, pi_floor) :
    st_bayesrc_compute_snp_pi(A, alpha_t, pi_floor);

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);

   arma::mat alpha_accum(nAnno, nstep, arma::fill::zeros);
   arma::vec sigmaSqAlpha_accum(nstep, arma::fill::zeros);
   arma::mat comp_prob_accum(m, Kgamma, arma::fill::zeros);
   arma::vec ncomp_accum(Kgamma, arma::fill::zeros);
   arma::vec selection_delta_accum(selectable_annotation_count,
                                   arma::fill::zeros);
   arma::mat selection_alpha_included_accum(selectable_annotation_count,
                                             nstep, arma::fill::zeros);
   arma::vec selection_alpha_included_count(selectable_annotation_count,
                                             arma::fill::zeros);
   arma::vec selection_tau2_accum(nstep, arma::fill::zeros);
   arma::vec selection_switches(selectable_annotation_count, arma::fill::zeros);
   double selection_pi_a_accum = 0.0;
   double selection_included_accum = 0.0;

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;
   double low_rank_residual_rebuild_count_t = 0.0;
   double low_rank_residual_max_abs_drift_t = 0.0;
   double maf_effect_s_current = maf_effect_s_init;
   double maf_effect_s_attempted_t = 0.0;
   double maf_effect_s_accepted_t = 0.0;
   arma::rowvec dynamic_prior_scale;
   BlockResidualChainState block_ve_state;
   if (block_residual_control.uses_block_variance())
    block_ve_state = make_block_residual_chain_state(
     block_count, nit + nburn,
     phenotype_variance[static_cast<std::size_t>(t)],
     block_residual_control.keep_history);
   std::vector<double> maximum_effect_history;
   if (capture_failure_state)
    maximum_effect_history.reserve(static_cast<std::size_t>(nit + nburn));
   const int allocation_updates = intercept_prior.allocation_updates_per_cycle;
   const int annotation_updates = intercept_prior.annotation_updates_per_cycle;
   const bool legacy_schedule = allocation_updates == 1 && annotation_updates == 1;
   const auto update_annotation_hierarchy = [&]() {
    if (!selection_config.enabled) {
     st_bayesrc_update_annotation_prior(
      A, comp_t, alpha_t, sigmaSqAlpha_t, intercept_prior,
      sigmaSqAlpha_a, sigmaSqAlpha_b, gen_t);
     snpPi_t = st_bayesrc_compute_snp_pi(A, alpha_t, pi_floor);
     return;
    }
    std::vector<arma::uvec> eligible;
    std::vector<arma::ivec> outcome;
    st_bayesrc_selection_build_observed_sticks(
     comp_t, nstep, eligible, outcome);
    for (int stick = 0; stick < nstep; ++stick) {
     const arma::uword index = static_cast<arma::uword>(stick);
     const bool empty = eligible[static_cast<std::size_t>(stick)].n_elem == 0u;
     const bool was_empty = selection_previous_empty(index) == 1u;
     if (empty) {
      selection_empty_diagnostics(index, 0u) += 1.0;
      if (!was_empty) selection_empty_diagnostics(index, 1u) += 1.0;
      selection_empty_run(index) += 1u;
      selection_empty_diagnostics(index, 3u) = std::max(
       selection_empty_diagnostics(index, 3u),
       static_cast<double>(selection_empty_run(index)));
     } else {
      if (was_empty) selection_empty_diagnostics(index, 2u) += 1.0;
      selection_empty_run(index) = 0u;
     }
     selection_previous_empty(index) = empty ? 1u : 0u;
    }
    const arma::uvec previous_delta = selection_state.delta;
    const std::vector<arma::vec> latent = st_bayesrc_selection_sample_latent(
     selection_annotation, eligible, outcome, selection_state.alpha, gen_t);
    const arma::ivec* fixed = selection_config.fixed_delta ?
     &selection_config.fixed_delta_value : nullptr;
    st_bayesrc_selection_delta_sweep(
     selection_annotation, eligible, latent, selection_state, gen_t, fixed);
    st_bayesrc_selection_blocked_redraw(
     selection_annotation, eligible, latent, selection_state,
     selection_config.intercept_mean,
     selection_config.intercept_precision, gen_t);
    st_bayesrc_selection_update_hyperparameters(
     selection_state, selection_config.hyper, gen_t,
     selection_config.update_pi_a, selection_config.update_tau2);
    selection_switches += arma::conv_to<arma::vec>::from(
     selection_state.delta != previous_delta);
    alpha_t = selection_state.alpha;
    sigmaSqAlpha_t = selection_state.tau2;
    snpPi_t = st_bayesrc_selection_compute_pi(
     selection_annotation, alpha_t, pi_floor);
   };

   for (int it = 0; it < nit + nburn; ++it) {
    const int allocation_repetitions = legacy_schedule ? 1 : allocation_updates;
    for (int allocation_rep = 0; allocation_rep < allocation_repetitions;
         ++allocation_rep) {
    const bool use_scaled_prior = use_maf_effect_s_prior_scale || estimate_maf_effect_s;
    if (estimate_maf_effect_s) {
     fill_maf_effect_s_prior_scale_sbayesrc(
      m,
      maf_effect_s_current,
      maf_effect_s_log_h_row,
      dynamic_prior_scale
     );
    }
    const arma::rowvec& current_prior_scale =
     estimate_maf_effect_s ? dynamic_prior_scale : prior_scale;

    if (use_scaled_prior) {
     for (int isort = 0; isort < m; ++isort) {
      const int i = order[static_cast<std::size_t>(isort)];

      sampleBeta_SBayesRC_ST_csr_scaled(
       i,
       snpPi_t.row(static_cast<arma::uword>(i)),
       gamma,
       vb_t,
       current_prior_scale,
       block_residual_variance_for_marker(
        op, i, vei_t, block_residual_control, block_ve_state),
       ww_t,
       r_t,
       b_t,
       comp_t,
       op,
       gen_t
      );
     }
    } else {
     for (int isort = 0; isort < m; ++isort) {
      const int i = order[static_cast<std::size_t>(isort)];

      sampleBeta_SBayesRC_ST_csr(
       i,
       snpPi_t.row(static_cast<arma::uword>(i)),
       gamma,
       vb_t,
       block_residual_variance_for_marker(
        op, i, vei_t, block_residual_control, block_ve_state),
       ww_t,
       r_t,
       b_t,
       comp_t,
       op,
       gen_t
      );
     }
    }

    if (updateLDswap && ld_swap_moves > 0 && ld_swap_prob > 0.0) {
     std::uniform_real_distribution<double> runif(0.0, 1.0);
     if (runif(gen_t) < ld_swap_prob) {
      for (int move = 0; move < ld_swap_moves; ++move) {
       bool attempted = false;
       const bool accepted = use_scaled_prior ?
        attempt_ld_swap_sbayesrc_ST_csr_scaled(
         m,
         t,
         chain,
         it,
         vei_t,
         vb_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
         gamma,
         snpPi_t,
         current_prior_scale,
         r_t,
         b_t,
         comp_t,
         op,
         ld_swap_friends,
        gen_t,
        attempted
       ) :
        attempt_ld_swap_sbayesrc_ST_csr(
         m,
         t,
         chain,
         it,
         vei_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
         snpPi_t,
         r_t,
         b_t,
         comp_t,
         op,
         ld_swap_friends,
         gen_t,
         attempted
        );
       if (attempted) {
        ld_swap_attempted_t += 1.0;
        if (accepted) ld_swap_accepted_t += 1.0;
       }
      }
     }
    }

    if (legacy_schedule && updateAlpha &&
        ((it + 1) % alpha_update_every == 0)) {
     update_annotation_hierarchy();
    }

    if (updateB) {
     if (use_scaled_prior) {
      sampleB_SBayesRC_ST_csr_scaled(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       gamma,
       current_prior_scale,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     } else {
      sampleB_SBayesRC_ST_csr(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       gamma,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     }

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error(
        "vb became invalid after sampleB. iter=" +
         std::to_string(it) +
         ", vb=" + std::to_string(vb_t)
      );
     }
    }

    if (estimate_maf_effect_s) {
     maf_effect_s_attempted_t += 1.0;
     const bool maf_effect_s_accepted = update_maf_effect_s_sbayesrc(
      maf_effect_s_current,
      b_t,
      comp_t,
      vb_t,
      gamma,
      maf_effect_s_log_h_row,
      maf_effect_s_prior_lower,
      maf_effect_s_prior_upper,
      maf_effect_s_proposal_sd,
      gen_t
     );
     if (maf_effect_s_accepted) maf_effect_s_accepted_t += 1.0;
    }

    if (capture_failure_state && allocation_rep + 1 == allocation_repetitions)
     maximum_effect_history.push_back(arma::abs(b_t).max());

    const int allocation_index = legacy_schedule
     ? it + 1
     : it * allocation_updates + allocation_rep + 1;
    if (op.uses_retained_low_rank() &&
        low_rank_residual_rebuild_every > 0 &&
        (allocation_index % low_rank_residual_rebuild_every == 0)) {
     const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
     low_rank_residual_max_abs_drift_t = std::max(
      low_rank_residual_max_abs_drift_t, drift);
     low_rank_residual_rebuild_count_t += 1.0;
    }

    if (updateE && !block_residual_control.uses_block_variance()) {
     if (rebuild_r_before_updateE) {
      op.rebuild(t, wy_t, b_t, r_t);
     }

     const double yy_t = yy_vec(static_cast<arma::uword>(t));
     const double sse_prior_t =
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
     const double maintained_sse = op.residual_sse(t, yy_t, b_t, wy_t, r_t);
     const double residual_scale = maintained_sse + nue * sse_prior_t;
     if (!std::isfinite(residual_scale) || residual_scale <= 0.0) {
      if (capture_failure_state) {
       CsrSBayesRCFailureState& state =
        failure_states[static_cast<std::size_t>(task)];
       arma::rowvec rebuilt = r_t;
       op.rebuild(t, wy_t, b_t, rebuilt);
       arma::rowvec marker_residual, marker_rebuilt_residual;
       op.materialize_residual(t, r_t, marker_residual);
       op.materialize_residual(t, rebuilt, marker_rebuilt_residual);
       const double b_wy = op.projected_score_dot(t, b_t, wy_t);
       const double independent_quadratic = op.quadratic_form(b_t);
       const double rebuilt_sse = op.residual_sse(
        t, yy_t, b_t, wy_t, rebuilt);
       const double quadratic_sse =
        yy_t - 2.0 * b_wy + independent_quadratic;
       const arma::rowvec residual_difference =
        marker_residual - marker_rebuilt_residual;
       const double maximum_residual_drift = residual_difference.is_empty()
        ? 0.0
        : arma::abs(residual_difference).max();
       const double rebuilt_norm = arma::norm(marker_rebuilt_residual, 2);
       const double independent_vg = independent_quadratic /
        static_cast<double>(n[t]);
       const double le = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
       state.captured = true;
       state.trait = t;
       state.chain = chain;
       state.iteration = it + 1;
       state.internal_iteration = it;
       state.sample_size = n[t];
       state.residual_df = nue;
       state.yy = yy_t;
       state.sse_prior = sse_prior_t;
       state.prior_contribution = nue * sse_prior_t;
       state.maintained_sse = maintained_sse;
       state.rebuilt_sse = rebuilt_sse;
       state.quadratic_sse = quadratic_sse;
       state.maintained_scale = residual_scale;
       state.rebuilt_scale = rebuilt_sse + state.prior_contribution;
       state.quadratic_scale = quadratic_sse + state.prior_contribution;
       state.b_wy = b_wy;
       state.b_r = arma::dot(b_t, marker_residual);
       state.fitted_quadratic = op.fitted_quadratic(t, b_t, wy_t, r_t);
       state.independent_quadratic = independent_quadratic;
       state.marker_variance = vb_t;
       state.genetic_variance = independent_vg;
       state.genic_variance = le;
       state.ld_variance = independent_vg - le;
       state.residual_variance = ve_t;
       state.adjusted_residual_variance = vei_t;
       state.heritability = independent_vg / (independent_vg + ve_t);
       state.active_marker_count = arma::accu(comp_t > 0);
       state.maximum_absolute_effect = arma::abs(b_t).max();
       state.effect_norm = arma::norm(b_t, 2);
       state.residual_norm = arma::norm(marker_residual, 2);
       state.rebuilt_residual_norm = rebuilt_norm;
       state.maximum_residual_drift = maximum_residual_drift;
       state.relative_residual_drift = maximum_residual_drift /
        std::max(1.0, rebuilt_norm);
       state.effects_finite = b_t.is_finite();
       state.residual_finite = marker_residual.is_finite();
       state.rebuilt_residual_finite = marker_rebuilt_residual.is_finite();
       state.operator_name = op.diagnostic_name();
       state.effects = b_t;
       state.residual = marker_residual;
       state.rebuilt_residual = marker_rebuilt_residual;
       state.score = wy_t;
       state.diagonal = ww_t;
       state.component = comp_t;
       const int trajectory_start = std::max(0, it - 20);
       for (int previous = trajectory_start; previous < it; ++previous) {
        const arma::uword previous_u = static_cast<arma::uword>(previous);
        state.preceding_marker_variance.push_back(vbs_t(previous_u));
        state.preceding_genetic_variance.push_back(vgs_t(previous_u));
        state.preceding_residual_variance.push_back(ves_t(previous_u));
        const double previous_total = vgs_t(previous_u) + ves_t(previous_u);
        state.preceding_heritability.push_back(
         previous_total > 0.0 ? vgs_t(previous_u) / previous_total : NA_REAL);
        state.preceding_maximum_effect.push_back(
         maximum_effect_history[static_cast<std::size_t>(previous)]);
       }
      }
      throw std::runtime_error(
       "sampleE_ST_operator: invalid projected residual scale at iteration " +
       std::to_string(it + 1) + " (internal " + std::to_string(it) +
       "); maintained_scale=" + std::to_string(residual_scale));
     }
     sampleE_ST_operator_from_scale(nue, ve_t, residual_scale, n[t], gen_t);
    }

    if (block_residual_control.uses_block_variance()) {
     const bool retain_block = allocation_rep + 1 == allocation_repetitions &&
      it >= nburn && ((it - nburn) % nthin == 0);
     update_block_residual_variance(
      op, t, it, retain_block, static_cast<double>(n[t]),
      phenotype_variance[static_cast<std::size_t>(t)], nue, b_t, r_t,
      block_residual_control, gen_t, block_ve_state);
     ve_t = mean_block_residual_variance(block_ve_state);
    }

    vg_t = computeG_ST_operator(op, t, b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error("vg became NaN/Inf after computeG. iter=" + std::to_string(it));
    }

    if (!std::isfinite(vle_t)) {
     throw std::runtime_error("vle became NaN/Inf after computeLE. iter=" + std::to_string(it));
    }

    if (!std::isfinite(vld_t)) {
     throw std::runtime_error("vld became NaN/Inf after computeLE. iter=" + std::to_string(it));
    }

    vei_t = block_residual_control.uses_block_variance() ? ve_t :
     ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error(
       "adjusted residual variance vei became invalid. iter=" +
        std::to_string(it) +
        ", vei=" + std::to_string(vei_t)
     );
    }
    }

    if (!legacy_schedule && updateAlpha &&
        ((it + 1) % alpha_update_every == 0)) {
     for (int annotation_rep = 0; annotation_rep < annotation_updates;
          ++annotation_rep) {
      update_annotation_hierarchy();
     }
    }

    double pi_active = 0.0;
    for (int i = 0; i < m; ++i) {
     double row_active = 0.0;
     for (int k = 1; k < Kgamma; ++k) {
      row_active += snpPi_t(static_cast<arma::uword>(i), static_cast<arma::uword>(k));
     }
     pi_active += row_active;
    }
    pi_active /= static_cast<double>(m);

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_active;
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;
    summary_heritability_task(task_u, static_cast<arma::uword>(it)) =
     block_residual_control.uses_block_variance() ?
      vg_t / phenotype_variance[static_cast<std::size_t>(t)] :
      vg_t / (vg_t + ve_t);
    if (estimate_maf_effect_s) {
     maf_effect_s_task(task_u, static_cast<arma::uword>(it)) = maf_effect_s_current;
    }

    if (it>=nburn) {
     const arma::uword draw=static_cast<arma::uword>(it-nburn);
     if (context.convergence_annotations) {
      arma::uword q=0;
      for (int stick=0; stick<nstep; ++stick)
       for (int annotation=0; annotation<nAnno; ++annotation)
        convergence_alpha_task[static_cast<std::size_t>(task)](draw,q++)=
         alpha_t(static_cast<arma::uword>(annotation),static_cast<arma::uword>(stick));
     convergence_sigma_task[static_cast<std::size_t>(task)].row(draw)=sigmaSqAlpha_t.t();
     if (selection_config.enabled) {
      selection_delta_trace_task[static_cast<std::size_t>(task)].row(draw) =
       selection_state.delta.t();
      selection_pi_a_trace_task[static_cast<std::size_t>(task)](draw) =
       selection_state.pi_a;
      selection_included_trace_task[static_cast<std::size_t>(task)](draw) =
       arma::accu(selection_state.delta);
     }
     }
     for (std::size_t q=0; q<convergence_markers.size(); ++q) {
      const arma::uword marker=static_cast<arma::uword>(convergence_markers[q]);
      const int component=comp_t(marker);
      if (context.convergence_b)
       convergence_b_task[static_cast<std::size_t>(task)](draw,q)=b_t(marker);
      if (context.convergence_d)
       convergence_d_task[static_cast<std::size_t>(task)](draw,q)=component>0 ? 1 : 0;
      if (context.convergence_component)
       convergence_component_task[static_cast<std::size_t>(task)](draw,q)=component;
     }
     if (context.convergence_component)
      capture_aggregate_component_trace(
       comp_t,m,Kgamma,static_cast<int>(draw),
       convergence_aggregate_task[static_cast<std::size_t>(task)]);
    }

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      const int ci = comp_t(iu);

      if (ci < 0 || ci >= Kgamma) {
       throw std::runtime_error("component index out of range during posterior accumulation.");
      }

      bm_t(iu) += b_t(iu);
      dm_t(iu) += (ci > 0) ? 1.0 : 0.0;
      comp_prob_accum(iu, static_cast<arma::uword>(ci)) += 1.0;
      ncomp_accum(static_cast<arma::uword>(ci)) += 1.0;
     }

     alpha_accum += alpha_t;
     sigmaSqAlpha_accum += sigmaSqAlpha_t;
     if (selection_config.enabled) {
      selection_delta_accum += arma::conv_to<arma::vec>::from(
       selection_state.delta);
      selection_pi_a_accum += selection_state.pi_a;
      selection_tau2_accum += selection_state.tau2;
      selection_included_accum += arma::accu(selection_state.delta);
      for (int annotation = 0; annotation < selectable_annotation_count;
           ++annotation) {
       if (selection_state.delta(static_cast<arma::uword>(annotation)) == 1u) {
        selection_alpha_included_accum.row(static_cast<arma::uword>(annotation)) +=
         selection_state.alpha.row(static_cast<arma::uword>(annotation + 1));
        selection_alpha_included_count(static_cast<arma::uword>(annotation)) += 1.0;
       }
      }
     }
    }
   }

   if (op.uses_retained_low_rank()) {
    const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
    low_rank_residual_max_abs_drift_t = std::max(
     low_rank_residual_max_abs_drift_t, drift);
    low_rank_residual_rebuild_count_t += 1.0;
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;

   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   alpha_accum /= nsamples_t;
   sigmaSqAlpha_accum /= nsamples_t;
   comp_prob_accum /= nsamples_t;
   ncomp_accum /= nsamples_t;

   if (!bm_t.is_finite()) {
    throw std::runtime_error("posterior mean bm contains NaN/Inf.");
   }

   if (!dm_t.is_finite()) {
    throw std::runtime_error("posterior mean dm contains NaN/Inf.");
   }

   bm_task.row(task_u) = bm_t;
   dm_task.row(task_u) = dm_t;
   b_task.row(task_u)  = b_t;
   arma::rowvec marker_residual;
   op.materialize_residual(t, r_t, marker_residual);
   r_task.row(task_u)  = marker_residual;
   for (int i = 0; i < m; ++i) {
    comp_task_double(task_u, static_cast<arma::uword>(i)) =
     static_cast<double>(comp_t(static_cast<arma::uword>(i)));
   }

   vbs_task.row(task_u) = vbs_t;
   vgs_task.row(task_u) = vgs_t;
   ves_task.row(task_u) = ves_t;
   pis_task.row(task_u) = pis_t;
   vles_task.row(task_u) = vles_t;
   vlds_task.row(task_u) = vlds_t;

   final_vb_task(task_u) = vb_t;
   final_ve_task(task_u) = ve_t;
   final_vg_task(task_u) = vg_t;
   final_pi_active_task(task_u) = pis_t(static_cast<arma::uword>(nit + nburn - 1));
   for (int k = 0; k < Kgamma; ++k) {
    final_pi_component_task(task_u, static_cast<arma::uword>(k)) =
     arma::mean(snpPi_t.col(static_cast<arma::uword>(k)));
   }
   final_vle_task(task_u) = vle_t;
   final_vld_task(task_u) = vld_t;
   nsamples_task(task_u) = nsamples_t;
   ld_swap_attempted_task(task_u) = ld_swap_attempted_t;
   ld_swap_accepted_task(task_u) = ld_swap_accepted_t;
   low_rank_residual_rebuild_count_task(task_u) =
    low_rank_residual_rebuild_count_t;
   low_rank_residual_max_abs_drift_task(task_u) =
    low_rank_residual_max_abs_drift_t;
   maf_effect_s_attempted_task(task_u) = maf_effect_s_attempted_t;
   maf_effect_s_accepted_task(task_u) = maf_effect_s_accepted_t;
   if (block_residual_control.uses_block_variance()) {
    block_ve_posterior_mean_task.row(task_u) =
     posterior_mean_block_residual_variance(block_ve_state).t();
    block_ve_final_task.row(task_u) = block_ve_state.value.t();
    block_ve_resampled_task.row(task_u) =
     arma::conv_to<arma::rowvec>::from(block_ve_state.resampled.t());
    block_ve_reset_task.row(task_u) =
     arma::conv_to<arma::rowvec>::from(block_ve_state.reset_to_phenotype.t());
    if (block_residual_control.keep_history)
     block_ve_history_task[static_cast<std::size_t>(task)] =
      block_ve_state.history;
   }

   alpha_mean_task[static_cast<std::size_t>(task)] = alpha_accum;
   sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)] = sigmaSqAlpha_accum;
   comp_prob_mean_task[static_cast<std::size_t>(task)] = comp_prob_accum;
   ncomp_mean_task[static_cast<std::size_t>(task)] = ncomp_accum;
   if (selection_config.enabled) {
    selection_pip_task[static_cast<std::size_t>(task)] =
     selection_delta_accum / nsamples_t;
    selection_pi_a_mean_task[static_cast<std::size_t>(task)](0u) =
     selection_pi_a_accum / nsamples_t;
    selection_tau2_mean_task[static_cast<std::size_t>(task)] =
     selection_tau2_accum / nsamples_t;
    selection_included_mean_task[static_cast<std::size_t>(task)](0u) =
     selection_included_accum / nsamples_t;
    selection_switches_task[static_cast<std::size_t>(task)] = selection_switches;
    for (int annotation = 0; annotation < selectable_annotation_count;
         ++annotation) {
     const double count = selection_alpha_included_count(
      static_cast<arma::uword>(annotation));
     if (count > 0.0)
      selection_alpha_conditional_mean_task[static_cast<std::size_t>(task)].row(
       static_cast<arma::uword>(annotation)) =
        selection_alpha_included_accum.row(static_cast<arma::uword>(annotation)) /
        count;
    }
    selection_delta_final_task[static_cast<std::size_t>(task)] =
     selection_state.delta;
    selection_empty_stick_diagnostics_task[static_cast<std::size_t>(task)] =
     selection_empty_diagnostics;
   }

#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif

  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = e.what();
#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif
  } catch (...) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = "unknown error";
#ifdef _OPENMP
   task_seconds[static_cast<std::size_t>(task)] = omp_get_wtime() - wall_start;
#endif
  }
 }

 for (int task = 0; task < ntasks; ++task) {
  if (failed[static_cast<std::size_t>(task)]) {
   if (capture_failure_state &&
       failure_states[static_cast<std::size_t>(task)].captured)
    write_csr_sbayesrc_failure_state(
     failure_states[static_cast<std::size_t>(task)], failure_state_path);
   throw std::runtime_error(
     "stblr_cpg_omp_csr_sbayesrc failed for trait " +
      std::to_string(stblr_task_trait(task, nchains)) +
      ", chain " +
      std::to_string(stblr_task_chain(task, nchains)) +
      ": " +
      errors[static_cast<std::size_t>(task)]
   );
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);
 arma::mat ld_swap_diagnostics(nt, 3, arma::fill::zeros);
 arma::mat ld_swap_chain_diagnostics(ntasks, 5, arma::fill::zeros);
 arma::mat low_rank_residual_diagnostics(nt, 3, arma::fill::zeros);
 arma::mat low_rank_residual_chain_diagnostics(ntasks, 5, arma::fill::zeros);

 for (int task = 0; task < ntasks; ++task) {
  const arma::uword task_u = static_cast<arma::uword>(task);
  const int task_trait = stblr_task_trait(task, nchains);
  const int task_chain = stblr_task_chain(task, nchains);
  const double attempted = ld_swap_attempted_task(task_u);
  const double accepted = ld_swap_accepted_task(task_u);
  const double rate = attempted > 0.0 ? accepted / attempted : 0.0;

  ld_swap_chain_diagnostics(task_u, 0) = task_trait;
  ld_swap_chain_diagnostics(task_u, 1) = task_chain;
  ld_swap_chain_diagnostics(task_u, 2) = attempted;
  ld_swap_chain_diagnostics(task_u, 3) = accepted;
  ld_swap_chain_diagnostics(task_u, 4) = rate;

  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 0) += attempted;
  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 1) += accepted;
  low_rank_residual_chain_diagnostics(task_u, 0) = task_trait;
  low_rank_residual_chain_diagnostics(task_u, 1) = task_chain;
  low_rank_residual_chain_diagnostics(task_u, 2) =
   low_rank_residual_rebuild_every;
  low_rank_residual_chain_diagnostics(task_u, 3) =
   low_rank_residual_rebuild_count_task(task_u);
  low_rank_residual_chain_diagnostics(task_u, 4) =
   low_rank_residual_max_abs_drift_task(task_u);
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 0) =
   low_rank_residual_rebuild_every;
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 1) +=
   low_rank_residual_rebuild_count_task(task_u);
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 2) =
   std::max(
    low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 2),
    low_rank_residual_max_abs_drift_task(task_u));
 }

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  ld_swap_diagnostics(tu, 2) =
   ld_swap_diagnostics(tu, 0) > 0.0
   ? ld_swap_diagnostics(tu, 1) / ld_swap_diagnostics(tu, 0)
   : 0.0;
 }

 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);

  bm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  dm_min_mat.row(tu).fill(std::numeric_limits<double>::infinity());
  bm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());
  dm_max_mat.row(tu).fill(-std::numeric_limits<double>::infinity());

  for (int chain = 0; chain < nchains; ++chain) {
   const int task = t * nchains + chain;
   const arma::uword task_u = static_cast<arma::uword>(task);

   bm_mat.row(tu) += bm_task.row(task_u);
   dm_mat.row(tu) += dm_task.row(task_u);
   b_mat.row(tu) += b_task.row(task_u);
   r_mat.row(tu) += r_task.row(task_u);
   for (int i = 0; i < m; ++i) {
    comp_mat(tu, static_cast<arma::uword>(i)) +=
     comp_task_double(task_u, static_cast<arma::uword>(i));
    bm_min_mat(tu, static_cast<arma::uword>(i)) =
     std::min(bm_min_mat(tu, static_cast<arma::uword>(i)), bm_task(task_u, static_cast<arma::uword>(i)));
    dm_min_mat(tu, static_cast<arma::uword>(i)) =
     std::min(dm_min_mat(tu, static_cast<arma::uword>(i)), dm_task(task_u, static_cast<arma::uword>(i)));
    bm_max_mat(tu, static_cast<arma::uword>(i)) =
     std::max(bm_max_mat(tu, static_cast<arma::uword>(i)), bm_task(task_u, static_cast<arma::uword>(i)));
    dm_max_mat(tu, static_cast<arma::uword>(i)) =
     std::max(dm_max_mat(tu, static_cast<arma::uword>(i)), dm_task(task_u, static_cast<arma::uword>(i)));
   }

   vbs_mat.row(tu) += vbs_task.row(task_u);
   vgs_mat.row(tu) += vgs_task.row(task_u);
   ves_mat.row(tu) += ves_task.row(task_u);
   pis_mat.row(tu) += pis_task.row(task_u);
   vles_mat.row(tu) += vles_task.row(task_u);
   vlds_mat.row(tu) += vlds_task.row(task_u);
   if (estimate_maf_effect_s) {
    maf_effect_s_mat.row(tu) += maf_effect_s_task.row(task_u);
    maf_effect_s_attempted_vec(tu) += maf_effect_s_attempted_task(task_u);
    maf_effect_s_accepted_vec(tu) += maf_effect_s_accepted_task(task_u);
   }

   final_vb(tu) += final_vb_task(task_u);
   final_ve(tu) += final_ve_task(task_u);
   final_vg(tu) += final_vg_task(task_u);
   final_pi_active(tu) += final_pi_active_task(task_u);
   final_pi_component.row(tu) += final_pi_component_task.row(task_u);
   final_vle(tu) += final_vle_task(task_u);
   final_vld(tu) += final_vld_task(task_u);
   nsamples_vec(tu) += nsamples_task(task_u);

   alpha_mean[static_cast<std::size_t>(t)] += alpha_mean_task[static_cast<std::size_t>(task)];
   sigmaSqAlpha_mean[static_cast<std::size_t>(t)] += sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)];
   comp_prob_mean[static_cast<std::size_t>(t)] += comp_prob_mean_task[static_cast<std::size_t>(task)];
   ncomp_mean[static_cast<std::size_t>(t)] += ncomp_mean_task[static_cast<std::size_t>(task)];
  }

  bm_mat.row(tu) *= inv_chains;
  dm_mat.row(tu) *= inv_chains;
  b_mat.row(tu) *= inv_chains;
  r_mat.row(tu) *= inv_chains;
  comp_mat.row(tu) *= inv_chains;
  vbs_mat.row(tu) *= inv_chains;
  vgs_mat.row(tu) *= inv_chains;
  ves_mat.row(tu) *= inv_chains;
  pis_mat.row(tu) *= inv_chains;
  vles_mat.row(tu) *= inv_chains;
  vlds_mat.row(tu) *= inv_chains;
  if (estimate_maf_effect_s) maf_effect_s_mat.row(tu) *= inv_chains;
  final_vb(tu) *= inv_chains;
  final_ve(tu) *= inv_chains;
  final_vg(tu) *= inv_chains;
  final_pi_active(tu) *= inv_chains;
  final_pi_component.row(tu) *= inv_chains;
  final_vle(tu) *= inv_chains;
  final_vld(tu) *= inv_chains;
  nsamples_vec(tu) *= inv_chains;
  alpha_mean[static_cast<std::size_t>(t)] *= inv_chains;
  sigmaSqAlpha_mean[static_cast<std::size_t>(t)] *= inv_chains;
  comp_prob_mean[static_cast<std::size_t>(t)] *= inv_chains;
  ncomp_mean[static_cast<std::size_t>(t)] *= inv_chains;

  if (nchains > 1) {
   for (int chain = 0; chain < nchains; ++chain) {
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    arma::rowvec bm_diff = bm_task.row(task_u) - bm_mat.row(tu);
    arma::rowvec dm_diff = dm_task.row(task_u) - dm_mat.row(tu);
    bm_sd_mat.row(tu) += bm_diff % bm_diff;
    dm_sd_mat.row(tu) += dm_diff % dm_diff;
   }
   bm_sd_mat.row(tu) = arma::sqrt(bm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
   dm_sd_mat.row(tu) = arma::sqrt(dm_sd_mat.row(tu) / static_cast<double>(nchains - 1));
  }
 }

 CsrSBayesRCExecutionResult result;
 result.bm_task = std::move(bm_task);
 result.dm_task = std::move(dm_task);
 result.b_task = std::move(b_task);
 result.r_task = std::move(r_task);
 result.comp_task_double = std::move(comp_task_double);
 result.vbs_task = std::move(vbs_task);
 result.vgs_task = std::move(vgs_task);
 result.ves_task = std::move(ves_task);
 result.pis_task = std::move(pis_task);
 result.vles_task = std::move(vles_task);
 result.vlds_task = std::move(vlds_task);
 result.maf_effect_s_task = std::move(maf_effect_s_task);
 result.final_vb_task = std::move(final_vb_task);
 result.final_vg_task = std::move(final_vg_task);
 result.final_ve_task = std::move(final_ve_task);
 result.final_pi_active_task = std::move(final_pi_active_task);
 result.final_pi_component_task = std::move(final_pi_component_task);
 result.final_vle_task = std::move(final_vle_task);
 result.final_vld_task = std::move(final_vld_task);
 result.nsamples_task = std::move(nsamples_task);
 result.ld_swap_attempted_task = std::move(ld_swap_attempted_task);
 result.ld_swap_accepted_task = std::move(ld_swap_accepted_task);
 result.low_rank_residual_rebuild_count_task =
  std::move(low_rank_residual_rebuild_count_task);
 result.low_rank_residual_max_abs_drift_task =
  std::move(low_rank_residual_max_abs_drift_task);
 result.maf_effect_s_attempted_task = std::move(maf_effect_s_attempted_task);
 result.maf_effect_s_accepted_task = std::move(maf_effect_s_accepted_task);
 result.alpha_mean_task = std::move(alpha_mean_task);
 result.sigmaSqAlpha_mean_task = std::move(sigmaSqAlpha_mean_task);
 result.comp_prob_mean_task = std::move(comp_prob_mean_task);
 result.ncomp_mean_task = std::move(ncomp_mean_task);
 result.convergence_alpha_task=std::move(convergence_alpha_task);
 result.convergence_sigma_task=std::move(convergence_sigma_task);
 result.convergence_b_task=std::move(convergence_b_task);
 result.convergence_d_task=std::move(convergence_d_task);
 result.convergence_component_task=std::move(convergence_component_task);
 result.convergence_aggregate_task=std::move(convergence_aggregate_task);
 result.selection_pip_task = std::move(selection_pip_task);
 result.selection_pi_a_mean_task = std::move(selection_pi_a_mean_task);
 result.selection_tau2_mean_task = std::move(selection_tau2_mean_task);
 result.selection_included_mean_task = std::move(selection_included_mean_task);
 result.selection_switches_task = std::move(selection_switches_task);
 result.selection_alpha_conditional_mean_task =
  std::move(selection_alpha_conditional_mean_task);
 result.selection_delta_final_task = std::move(selection_delta_final_task);
 result.selection_empty_stick_diagnostics_task =
  std::move(selection_empty_stick_diagnostics_task);
 result.selection_delta_trace_task = std::move(selection_delta_trace_task);
 result.selection_pi_a_trace_task = std::move(selection_pi_a_trace_task);
 result.selection_included_trace_task = std::move(selection_included_trace_task);
 result.bm_mat = std::move(bm_mat);
 result.dm_mat = std::move(dm_mat);
 result.bm_sd_mat = std::move(bm_sd_mat);
 result.dm_sd_mat = std::move(dm_sd_mat);
 result.bm_min_mat = std::move(bm_min_mat);
 result.dm_min_mat = std::move(dm_min_mat);
 result.bm_max_mat = std::move(bm_max_mat);
 result.dm_max_mat = std::move(dm_max_mat);
 result.b_mat = std::move(b_mat);
 result.r_mat = std::move(r_mat);
 result.comp_mat = std::move(comp_mat);
 result.vbs_mat = std::move(vbs_mat);
 result.vgs_mat = std::move(vgs_mat);
 result.ves_mat = std::move(ves_mat);
 result.pis_mat = std::move(pis_mat);
 result.vles_mat = std::move(vles_mat);
 result.vlds_mat = std::move(vlds_mat);
 result.maf_effect_s_mat = std::move(maf_effect_s_mat);
 result.final_vb = std::move(final_vb);
 result.final_vg = std::move(final_vg);
 result.final_ve = std::move(final_ve);
 result.final_pi_active = std::move(final_pi_active);
 result.final_pi_component = std::move(final_pi_component);
 result.final_vle = std::move(final_vle);
 result.final_vld = std::move(final_vld);
 result.nsamples_vec = std::move(nsamples_vec);
 result.maf_effect_s_attempted_vec = std::move(maf_effect_s_attempted_vec);
 result.maf_effect_s_accepted_vec = std::move(maf_effect_s_accepted_vec);
 result.alpha_mean = std::move(alpha_mean);
 result.sigmaSqAlpha_mean = std::move(sigmaSqAlpha_mean);
 result.comp_prob_mean = std::move(comp_prob_mean);
 result.ncomp_mean = std::move(ncomp_mean);
 result.failed = std::move(failed);
 result.errors = std::move(errors);
 result.thread_used = std::move(thread_used);
 result.task_seconds = std::move(task_seconds);
 result.ld_swap_diagnostics = std::move(ld_swap_diagnostics);
 result.ld_swap_chain_diagnostics = std::move(ld_swap_chain_diagnostics);
 result.low_rank_residual_diagnostics =
  std::move(low_rank_residual_diagnostics);
 result.low_rank_residual_chain_diagnostics =
  std::move(low_rank_residual_chain_diagnostics);
 result.block_ve_posterior_mean_task =
  std::move(block_ve_posterior_mean_task);
 result.block_ve_final_task = std::move(block_ve_final_task);
 result.block_ve_resampled_task = std::move(block_ve_resampled_task);
 result.block_ve_reset_task = std::move(block_ve_reset_task);
 result.summary_heritability_task = std::move(summary_heritability_task);
 result.block_ve_history_task = std::move(block_ve_history_task);
 return result;
}

} }

#endif
