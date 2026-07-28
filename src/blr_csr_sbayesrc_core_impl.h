#ifndef SBLR_BLR_CSR_SBAYESRC_CORE_IMPL_H
#define SBLR_BLR_CSR_SBAYESRC_CORE_IMPL_H

// Implementation detail only. This header is included by
// st_sbayesrc_omp_csr.cpp after the package's Armadillo configuration and the
// concrete operator types have been established. All borrowed objects in the
// context must outlive run_csr_sbayesrc(); this is the native realization of
// the canonical storage_outlives_execution ownership contract.

namespace sblr { namespace core {

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
 const arma::rowvec& selection_s_log_h;

 // Contract views describe the same borrowed prepared inputs. The
 // numerical body continues to use the native Armadillo/operator references.
 const SBayesRCAnnotationDesignView& annotation_contract;
 const SBayesRCComponentSpec& component_contract;
 const SBayesRCAlphaSpec& alpha_contract;
 const SBayesRCProbabilityPolicy& probability_contract;
 const SBayesRCPriors& prior_contract;
 const SBayesRCControls& control_contract;
 const SBayesRCOutputControl& output_contract;

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
 bool intercept_flat;
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
 bool estimate_selection_s;
 double selection_s_initial;
 double selection_s_prior_lower;
 double selection_s_prior_upper;
 double selection_s_proposal_sd;
 const std::vector<int>& convergence_markers;
 bool convergence_annotations;
 bool convergence_b;
 bool convergence_d;
 bool convergence_component;
};

struct CsrSBayesRCExecutionResult {
 arma::mat bm_task, dm_task, b_task, r_task, comp_task_double;
 arma::mat vbs_task, vgs_task, ves_task, pis_task, vles_task, vlds_task;
 arma::mat selection_s_task;
 arma::vec final_vb_task, final_vg_task, final_ve_task;
 arma::vec final_pi_active_task, final_vle_task, final_vld_task;
 arma::mat final_pi_component_task;
 arma::vec nsamples_task, ld_swap_attempted_task, ld_swap_accepted_task;
 arma::vec selection_s_attempted_task, selection_s_accepted_task;
 std::vector<arma::mat> alpha_mean_task, comp_prob_mean_task;
 std::vector<arma::vec> sigmaSqAlpha_mean_task, ncomp_mean_task;
 std::vector<arma::mat> convergence_alpha_task, convergence_sigma_task;
 std::vector<arma::mat> convergence_b_task;
 std::vector<arma::imat> convergence_d_task, convergence_component_task;

 arma::mat bm_mat, dm_mat, bm_sd_mat, dm_sd_mat;
 arma::mat bm_min_mat, dm_min_mat, bm_max_mat, dm_max_mat;
 arma::mat b_mat, r_mat, comp_mat;
 arma::mat vbs_mat, vgs_mat, ves_mat, pis_mat, vles_mat, vlds_mat;
 arma::mat selection_s_mat;
 arma::vec final_vb, final_vg, final_ve, final_pi_active;
 arma::mat final_pi_component;
 arma::vec final_vle, final_vld, nsamples_vec;
 arma::vec selection_s_attempted_vec, selection_s_accepted_vec;
 std::vector<arma::mat> alpha_mean, comp_prob_mean;
 std::vector<arma::vec> sigmaSqAlpha_mean, ncomp_mean;
 std::vector<int> failed, thread_used;
 std::vector<std::string> errors;
 std::vector<double> task_seconds;
 arma::mat ld_swap_diagnostics, ld_swap_chain_diagnostics;
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
 const arma::rowvec& prior_scale = context.prior_scale;
 const arma::rowvec& selection_s_log_h_row = context.selection_s_log_h;
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
 const bool intercept_flat = context.intercept_flat;
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
 const bool use_selection_s_prior_scale = context.use_selection_prior_scale;
 const bool estimate_selection_s = context.estimate_selection_s;
 const double selection_s_init = context.selection_s_initial;
 const double selection_s_prior_lower = context.selection_s_prior_lower;
 const double selection_s_prior_upper = context.selection_s_prior_upper;
 const double selection_s_proposal_sd = context.selection_s_proposal_sd;
 const std::vector<int>& convergence_markers=context.convergence_markers;
 const int ntasks = stblr_num_chain_tasks(nt, nchains);

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
 arma::mat selection_s_task(ntasks, nit + nburn, arma::fill::zeros);

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
 arma::vec selection_s_attempted_task(ntasks, arma::fill::zeros);
 arma::vec selection_s_accepted_task(ntasks, arma::fill::zeros);

 std::vector<arma::mat> alpha_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> sigmaSqAlpha_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> comp_prob_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::vec> ncomp_mean_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_alpha_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_sigma_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_b_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_d_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_component_task(static_cast<std::size_t>(ntasks));

 for (int marker: convergence_markers) if (marker<0 || marker>=m)
  throw std::invalid_argument("SBayesRC convergence marker index is out of range.");

 for (int task = 0; task < ntasks; ++task) {
  alpha_mean_task[static_cast<std::size_t>(task)] = arma::mat(nAnno, nstep, arma::fill::zeros);
  sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)] = arma::vec(nstep, arma::fill::zeros);
  comp_prob_mean_task[static_cast<std::size_t>(task)] = arma::mat(m, Kgamma, arma::fill::zeros);
  ncomp_mean_task[static_cast<std::size_t>(task)] = arma::vec(Kgamma, arma::fill::zeros);
  if (context.convergence_annotations) {
   convergence_alpha_task[static_cast<std::size_t>(task)]=
    arma::mat(nit,nAnno*nstep,arma::fill::zeros);
   convergence_sigma_task[static_cast<std::size_t>(task)]=
    arma::mat(nit,nstep,arma::fill::zeros);
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
 arma::mat selection_s_mat(nt, nit + nburn, arma::fill::zeros);

 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi_active(nt, arma::fill::zeros);
 arma::mat final_pi_component(nt, Kgamma, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec selection_s_attempted_vec(nt, arma::fill::zeros);
 arma::vec selection_s_accepted_vec(nt, arma::fill::zeros);

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

   arma::rowvec r_t(m, arma::fill::zeros);
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
    op.rebuild(wy_t, b_t, r_t);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_SBayesRC_ST_csr(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   arma::mat alpha_t = alpha_init;
   arma::vec sigmaSqAlpha_t = sigmaSqAlpha_init;

   for (int j = 0; j < nstep; ++j) {
    if (!std::isfinite(sigmaSqAlpha_t(static_cast<arma::uword>(j))) ||
        sigmaSqAlpha_t(static_cast<arma::uword>(j)) <= 0.0) {
     throw std::runtime_error("sigmaSqAlpha_init contains invalid value.");
    }
   }

   arma::mat snpPi_t = st_bayesrc_compute_snp_pi(A, alpha_t, pi_floor);

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

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;
   double selection_s_current = selection_s_init;
   double selection_s_attempted_t = 0.0;
   double selection_s_accepted_t = 0.0;
   arma::rowvec dynamic_prior_scale;

   for (int it = 0; it < nit + nburn; ++it) {
    const bool use_scaled_prior = use_selection_s_prior_scale || estimate_selection_s;
    if (estimate_selection_s) {
     fill_selection_s_prior_scale_sbayesrc(
      m,
      selection_s_current,
      selection_s_log_h_row,
      dynamic_prior_scale
     );
    }
    const arma::rowvec& current_prior_scale =
     estimate_selection_s ? dynamic_prior_scale : prior_scale;

    if (use_scaled_prior) {
     for (int isort = 0; isort < m; ++isort) {
      const int i = order[static_cast<std::size_t>(isort)];

      sampleBeta_SBayesRC_ST_csr_scaled(
       i,
       snpPi_t.row(static_cast<arma::uword>(i)),
       gamma,
       vb_t,
       current_prior_scale,
       vei_t,
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
       vei_t,
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

    if (updateAlpha && ((it + 1) % alpha_update_every == 0)) {
     st_bayesrc_update_annotation_prior(
      A,
      comp_t,
      alpha_t,
      sigmaSqAlpha_t,
      intercept_flat,
      sigmaSqAlpha_a,
      sigmaSqAlpha_b,
      gen_t
     );

     snpPi_t = st_bayesrc_compute_snp_pi(A, alpha_t, pi_floor);
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

    if (estimate_selection_s) {
     selection_s_attempted_t += 1.0;
     const bool selection_s_accepted = update_selection_s_sbayesrc(
      selection_s_current,
      b_t,
      comp_t,
      vb_t,
      gamma,
      selection_s_log_h_row,
      selection_s_prior_lower,
      selection_s_prior_upper,
      selection_s_proposal_sd,
      gen_t
     );
     if (selection_s_accepted) selection_s_accepted_t += 1.0;
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      op.rebuild(wy_t, b_t, r_t);
     }

     sampleE_ST_csr(
      m,
      nue,
      ve_t,
      b_t,
      wy_t,
      r_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t)),
      n[t],
       gen_t
     );
    }

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
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

    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error(
       "adjusted residual variance vei became invalid. iter=" +
        std::to_string(it) +
        ", vei=" + std::to_string(vei_t)
     );
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
    if (estimate_selection_s) {
     selection_s_task(task_u, static_cast<arma::uword>(it)) = selection_s_current;
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
    }
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
   r_task.row(task_u)  = r_t;
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
   selection_s_attempted_task(task_u) = selection_s_attempted_t;
   selection_s_accepted_task(task_u) = selection_s_accepted_t;

   alpha_mean_task[static_cast<std::size_t>(task)] = alpha_accum;
   sigmaSqAlpha_mean_task[static_cast<std::size_t>(task)] = sigmaSqAlpha_accum;
   comp_prob_mean_task[static_cast<std::size_t>(task)] = comp_prob_accum;
   ncomp_mean_task[static_cast<std::size_t>(task)] = ncomp_accum;

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
   if (estimate_selection_s) {
    selection_s_mat.row(tu) += selection_s_task.row(task_u);
    selection_s_attempted_vec(tu) += selection_s_attempted_task(task_u);
    selection_s_accepted_vec(tu) += selection_s_accepted_task(task_u);
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
  if (estimate_selection_s) selection_s_mat.row(tu) *= inv_chains;
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
 result.selection_s_task = std::move(selection_s_task);
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
 result.selection_s_attempted_task = std::move(selection_s_attempted_task);
 result.selection_s_accepted_task = std::move(selection_s_accepted_task);
 result.alpha_mean_task = std::move(alpha_mean_task);
 result.sigmaSqAlpha_mean_task = std::move(sigmaSqAlpha_mean_task);
 result.comp_prob_mean_task = std::move(comp_prob_mean_task);
 result.ncomp_mean_task = std::move(ncomp_mean_task);
 result.convergence_alpha_task=std::move(convergence_alpha_task);
 result.convergence_sigma_task=std::move(convergence_sigma_task);
 result.convergence_b_task=std::move(convergence_b_task);
 result.convergence_d_task=std::move(convergence_d_task);
 result.convergence_component_task=std::move(convergence_component_task);
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
 result.selection_s_mat = std::move(selection_s_mat);
 result.final_vb = std::move(final_vb);
 result.final_vg = std::move(final_vg);
 result.final_ve = std::move(final_ve);
 result.final_pi_active = std::move(final_pi_active);
 result.final_pi_component = std::move(final_pi_component);
 result.final_vle = std::move(final_vle);
 result.final_vld = std::move(final_vld);
 result.nsamples_vec = std::move(nsamples_vec);
 result.selection_s_attempted_vec = std::move(selection_s_attempted_vec);
 result.selection_s_accepted_vec = std::move(selection_s_accepted_vec);
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
 return result;
}

} }

#endif
