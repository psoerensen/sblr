#include "blr_aggregate_component_trace.h"
#include "blr_csr_bayesr_policy.h"

template <class Operator>
struct CsrBayesRExecutionContext {
 Operator* op = nullptr;
 const BayesRLDLDFriends* ld_swap_friends = nullptr;
 const std::vector<double>* mixture_var = nullptr;
 const std::vector<double>* pi = nullptr;
 const std::vector<double>* alpha = nullptr;
 const std::vector<int>* n = nullptr;
 const std::vector<int>* chain_seeds = nullptr;
 const std::vector<std::vector<double>>* comp_init = nullptr;
 const std::vector<std::vector<double>>* r_init = nullptr;
 const arma::mat* wy_mat = nullptr; const arma::mat* b_init_mat = nullptr;
 const arma::mat* B = nullptr; const arma::mat* E = nullptr;
 const arma::mat* ssb_prior_mat = nullptr; const arma::mat* sse_prior_mat = nullptr;
 const arma::vec* yy_vec = nullptr;
 const arma::rowvec* prior_scale = nullptr; const arma::rowvec* maf_effect_s_log_h_row = nullptr;
 const std::vector<int>* convergence_markers = nullptr;
 const std::vector<double>* phenotype_variance = nullptr;
 sblr::core::BlockResidualControl block_residual_control;
 int K=0,m=0,nt=0,nchains=1,ncores=1,nit=0,nburn=0,nthin=1,seed=1;
 int updateE_start=0,updateE_every=1,ld_swap_moves=1;
 int low_rank_residual_rebuild_every=0;
 bool use_comp_init=false,use_r_init=false,estimate_maf_effect_s=false;
 bool use_maf_effect_s_prior_scale=false,updateLDswap=false,updateB=true,updateE=true,updatePi=true;
 bool convergence_probability=false,convergence_b=false,convergence_d=false;
 bool convergence_component=false;
 double ld_swap_prob=0.0,maf_effect_s_init=0.0,maf_effect_s_prior_lower=-3.0;
 double maf_effect_s_prior_upper=2.0,maf_effect_s_proposal_sd=0.35,nub=0.0,nue=0.0,adjE=0.0;
};

struct CsrBayesRExecutionResult {
 arma::vec mixture_var_vec;
 arma::mat bm_task,dm_task,component_mean_task,b_task,r_task,component_task;
 arma::mat vbs_task,vgs_task,ves_task,vles_task,vlds_task,pis_task,maf_effect_s_task;
 arma::mat final_pi_task,mean_pi_task;
 arma::vec final_vb_task,final_vg_task,final_ve_task;
 arma::vec maf_effect_s_attempted_task,maf_effect_s_accepted_task;
 std::vector<arma::mat> comp_prob_task;
 std::vector<arma::mat> convergence_pi_task,convergence_b_task;
 std::vector<arma::imat> convergence_d_task,convergence_component_task;
 std::vector<sblr::core::AggregateComponentTrace> convergence_aggregate_task;
 std::vector<arma::vec> ncomp_task;
 arma::mat bm,dm,bm_sd,dm_sd,bm_min,dm_min,bm_max,dm_max,component_mean;
 arma::mat b_out,r_out,component_out,vbs,vgs,ves,vle,vld,pis,maf_effect_s;
 arma::mat final_pi,mean_pi,updateE_diagnostics,ld_swap_diagnostics,ld_swap_chain_diagnostics;
 arma::mat low_rank_residual_diagnostics,low_rank_residual_chain_diagnostics;
 arma::mat block_ve_posterior_mean_task,block_ve_final_task;
 arma::mat block_ve_resampled_task,block_ve_reset_task,summary_heritability_task;
 std::vector<arma::mat> block_ve_history_task;
 arma::vec final_vb,final_vg,final_ve,maf_effect_s_attempted,maf_effect_s_accepted,nsamples;
 std::vector<arma::mat> comp_prob;
 arma::mat ncomp,covb,covg,cove,vb,vg,ve;
};

// The policy surface is deliberately limited to a marker-scale provider and
// one post-vb hook. The ordinary policy owns no state and consumes no RNG.
struct CsrBayesRNoOpPolicy {
 bool provides_prior_scale() const noexcept { return false; }
 const arma::rowvec& prior_scale() const {
  throw std::logic_error("ordinary BayesR has no policy-owned prior scale");
 }
 void after_vb_update(
   const arma::rowvec&, const arma::Row<int>&, double,
   const arma::vec&, std::mt19937&, int) noexcept {}
 void capture(int) noexcept {}
 void retain(int) noexcept {}
 void finish() noexcept {}
};

struct CsrBayesRNoOpPolicyFactory {
 CsrBayesRNoOpPolicy make(int, int, int, int) const noexcept {
  return CsrBayesRNoOpPolicy{};
 }
};

// Canonical binding-neutral, operator-aware implementation.
template <class Operator, class PolicyFactory>
CsrBayesRExecutionResult run_csr_bayesr_engine(
  CsrBayesRExecutionContext<Operator>& context,
  PolicyFactory& policy_factory) {
 Operator& op=*context.op; const BayesRLDLDFriends& ld_swap_friends=*context.ld_swap_friends;
 const std::vector<double>& mixture_var=*context.mixture_var; const std::vector<double>& pi=*context.pi;
 const std::vector<double>& alpha=*context.alpha; const std::vector<int>& n=*context.n;
 const std::vector<int>& chain_seeds=*context.chain_seeds;
 const std::vector<std::vector<double>>& comp_init=*context.comp_init;
 const std::vector<std::vector<double>>& r_init=*context.r_init;
 const arma::mat& wy_mat=*context.wy_mat; const arma::mat& b_init_mat=*context.b_init_mat;
 const arma::mat& B=*context.B; const arma::mat& E=*context.E;
 const arma::mat& ssb_prior_mat=*context.ssb_prior_mat; const arma::mat& sse_prior_mat=*context.sse_prior_mat;
 const arma::vec& yy_vec=*context.yy_vec; const arma::rowvec& prior_scale=*context.prior_scale;
 const arma::rowvec& maf_effect_s_log_h_row=*context.maf_effect_s_log_h_row;
 const std::vector<int> empty_convergence_markers;
 const std::vector<int>& convergence_markers = context.convergence_markers ?
  *context.convergence_markers : empty_convergence_markers;
 const std::vector<double> empty_phenotype_variance;
 const std::vector<double>& phenotype_variance = context.phenotype_variance ?
  *context.phenotype_variance : empty_phenotype_variance;
 const sblr::core::BlockResidualControl& block_residual_control =
  context.block_residual_control;
 const int K=context.K,m=context.m,nt=context.nt,nchains=context.nchains,ncores=context.ncores;
 const int nit=context.nit,nburn=context.nburn,nthin=context.nthin,seed=context.seed;
 const int updateE_start=context.updateE_start,updateE_every=context.updateE_every,ld_swap_moves=context.ld_swap_moves;
 const int low_rank_residual_rebuild_every=context.low_rank_residual_rebuild_every;
 const bool use_comp_init=context.use_comp_init,use_r_init=context.use_r_init;
 const bool estimate_maf_effect_s=context.estimate_maf_effect_s,use_maf_effect_s_prior_scale=context.use_maf_effect_s_prior_scale;
 const bool updateLDswap=context.updateLDswap,updateB=context.updateB,updateE=context.updateE,updatePi=context.updatePi;
 const bool convergence_probability=context.convergence_probability;
 const bool convergence_b=context.convergence_b,convergence_d=context.convergence_d;
 const bool convergence_component=context.convergence_component;
 const double ld_swap_prob=context.ld_swap_prob,maf_effect_s_init=context.maf_effect_s_init;
 const double maf_effect_s_proposal_sd=context.maf_effect_s_proposal_sd,nub=context.nub,nue=context.nue,adjE=context.adjE;
 const double maf_effect_s_prior[2]={context.maf_effect_s_prior_lower,context.maf_effect_s_prior_upper};
 arma::vec mixture_var_vec(K, arma::fill::zeros);
 for (int k = 0; k < K; ++k) mixture_var_vec(static_cast<arma::uword>(k)) = mixture_var[k];

 std::vector<int> order(static_cast<std::size_t>(m));
 std::iota(order.begin(), order.end(), 0);

 const int ntasks = stblr_num_chain_tasks(nt, nchains);
 const int nthreads = stblr_num_threads_for_tasks(ncores, ntasks);
 const int trace_len = nit + nburn;
 const int block_count = op.block_count();
 if (block_residual_control.uses_block_variance() &&
     (block_count <= 0 || static_cast<int>(phenotype_variance.size()) != nt))
  throw std::runtime_error(
   "block residual policy requires retained blocks and one phenotype variance per trait.");

 arma::mat bm_task(ntasks, m, arma::fill::zeros);
 arma::mat dm_task(ntasks, m, arma::fill::zeros);
 arma::mat component_mean_task(ntasks, m, arma::fill::zeros);
 arma::mat b_task(ntasks, m, arma::fill::zeros);
 arma::mat r_task(ntasks, m, arma::fill::zeros);
 arma::mat component_task(ntasks, m, arma::fill::zeros);
 arma::mat vbs_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vgs_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat ves_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vles_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat vlds_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat pis_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat maf_effect_s_task(ntasks, trace_len, arma::fill::zeros);
 arma::mat final_pi_task(ntasks, K, arma::fill::zeros);
 arma::mat mean_pi_task(ntasks, K, arma::fill::zeros);
 arma::vec final_vb_task(ntasks, arma::fill::zeros);
 arma::vec final_vg_task(ntasks, arma::fill::zeros);
 arma::vec final_ve_task(ntasks, arma::fill::zeros);
 arma::vec min_sse_task(ntasks, arma::fill::zeros);
 arma::ivec min_sse_iter_task(ntasks, arma::fill::value(-1));
 arma::vec min_residual_scale_task(ntasks, arma::fill::zeros);
 arma::ivec max_nonzero_components_task(ntasks, arma::fill::zeros);
 arma::vec max_abs_effect_task(ntasks, arma::fill::zeros);
 arma::vec max_fitted_quadratic_task(ntasks, arma::fill::zeros);
 arma::ivec n_updateE_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_attempted_task(ntasks, arma::fill::zeros);
 arma::vec ld_swap_accepted_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_attempted_task(ntasks, arma::fill::zeros);
 arma::vec maf_effect_s_accepted_task(ntasks, arma::fill::zeros);
 arma::ivec low_rank_residual_rebuild_count_task(ntasks, arma::fill::zeros);
 arma::vec low_rank_residual_max_abs_drift_task(ntasks, arma::fill::zeros);
 arma::vec nsamples_task(ntasks, arma::fill::zeros);
 arma::mat block_ve_posterior_mean_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_final_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_resampled_task(ntasks, block_count, arma::fill::zeros);
 arma::mat block_ve_reset_task(ntasks, block_count, arma::fill::zeros);
 arma::mat summary_heritability_task(ntasks, trace_len, arma::fill::zeros);
 std::vector<arma::mat> block_ve_history_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> comp_prob_task(static_cast<std::size_t>(ntasks));
 const int convergence_pi_count = convergence_probability ? (K == 2 ? 1 : K) : 0;
 const int convergence_marker_count = static_cast<int>(convergence_markers.size());
 std::vector<arma::mat> convergence_pi_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::mat> convergence_b_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_d_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_component_task(static_cast<std::size_t>(ntasks));
 std::vector<sblr::core::AggregateComponentTrace> convergence_aggregate_task(
  static_cast<std::size_t>(ntasks));
 for (int task = 0; task < ntasks; ++task) {
  if (convergence_pi_count > 0) convergence_pi_task[static_cast<std::size_t>(task)].zeros(nit, convergence_pi_count);
  if (convergence_b && convergence_marker_count > 0) convergence_b_task[static_cast<std::size_t>(task)].zeros(nit, convergence_marker_count);
  if (convergence_d && convergence_marker_count > 0) convergence_d_task[static_cast<std::size_t>(task)].zeros(nit, convergence_marker_count);
  if (convergence_component && convergence_marker_count > 0) convergence_component_task[static_cast<std::size_t>(task)].zeros(nit, convergence_marker_count);
  if (convergence_component) sblr::core::allocate_aggregate_component_trace(
   convergence_aggregate_task[static_cast<std::size_t>(task)],nit,K);
 }
 std::vector<arma::vec> ncomp_task(static_cast<std::size_t>(ntasks));
 std::vector<int> failed(static_cast<std::size_t>(ntasks), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(ntasks));

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int task = 0; task < ntasks; ++task) {
  const int t = stblr_task_trait(task, nchains);
  const int chain = stblr_task_chain(task, nchains);
  const arma::uword task_u = static_cast<arma::uword>(task);

  try {
   unsigned int task_seed = 0u;
   if (!chain_seeds.empty()) {
    task_seed = stblr_seed_with_chain_base(chain_seeds[static_cast<std::size_t>(chain)], t);
   } else if (nchains == 1) {
    task_seed = stblr_trait_seed(seed, t);
   } else {
    task_seed = stblr_chain_seed(seed, t, chain);
   }

   std::mt19937 gen_t(task_seed);
   std::vector<int> order_t = order;
   std::shuffle(order_t.begin(), order_t.end(), gen_t);
   auto policy = policy_factory.make(task, t, chain, m);

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   const arma::rowvec& ww_t = op.diag();
   arma::rowvec b_t = b_init_mat.row(static_cast<arma::uword>(t));
   arma::rowvec r_t(op.residual_size(), arma::fill::zeros);
   arma::Row<int> comp_t(m, arma::fill::zeros);

   if (use_comp_init) {
    for (int i = 0; i < m; ++i) {
     const double ci = comp_init[t][i];
     if (!std::isfinite(ci) || ci != std::floor(ci) || ci < 0.0 || ci >= K) {
      throw std::runtime_error("comp_init contains invalid component index.");
     }
     comp_t(static_cast<arma::uword>(i)) = static_cast<int>(ci);
    }
   } else {
    for (int i = 0; i < m; ++i) {
     comp_t(static_cast<arma::uword>(i)) = (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = r_init[t][i];
    }
    if (!r_t.is_finite()) throw std::runtime_error("r_init contains NaN/Inf.");
   } else {
    op.rebuild(t, wy_t, b_t, r_t);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   if (!std::isfinite(vb_t) || vb_t <= 0.0 || !std::isfinite(ve_t) || ve_t <= 0.0) {
    throw std::runtime_error("initial B/E diagonal values must be finite and positive.");
   }

   std::vector<double> pi_t = pi;
   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec component_mean_t(m, arma::fill::zeros);
   arma::mat comp_prob_t(m, K, arma::fill::zeros);
   arma::vec pi_mean_t(K, arma::fill::zeros);
   double nsamples_t = 0.0;
   double min_sse_t = std::numeric_limits<double>::infinity();
   int min_sse_iter_t = -1;
   double min_residual_scale_t = std::numeric_limits<double>::infinity();
   int max_nonzero_components_t = 0;
   double max_abs_effect_t = 0.0;
   double max_fitted_quadratic_t = 0.0;
   int n_updateE_t = 0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;
   double maf_effect_s_current = maf_effect_s_init;
   double maf_effect_s_attempted_t = 0.0;
   double maf_effect_s_accepted_t = 0.0;
   int low_rank_residual_rebuild_count_t = 0;
   double low_rank_residual_max_abs_drift_t = 0.0;
   arma::rowvec dynamic_prior_scale;
   sblr::core::BlockResidualChainState block_ve_state;
   if (block_residual_control.uses_block_variance())
    block_ve_state = sblr::core::make_block_residual_chain_state(
     block_count, trace_len, phenotype_variance[static_cast<std::size_t>(t)],
     block_residual_control.keep_history);

   double vg_t = computeG_ST_operator(op, t, b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_bayesr_ST_csr(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   for (int it = 0; it < trace_len; ++it) {
    if (estimate_maf_effect_s) {
     fill_maf_effect_s_prior_scale_bayesr(
      m,
      maf_effect_s_current,
      maf_effect_s_log_h_row,
      dynamic_prior_scale
     );
    }
    if (estimate_maf_effect_s || use_maf_effect_s_prior_scale ||
        policy.provides_prior_scale()) {
     const arma::rowvec& current_prior_scale =
     estimate_maf_effect_s ? dynamic_prior_scale :
      (use_maf_effect_s_prior_scale ? prior_scale : policy.prior_scale());
     for (int isort = 0; isort < m; ++isort) {
      const int marker = order_t[static_cast<std::size_t>(isort)];
      sampleBetaR_ST_csr(
       marker,
       pi_t,
       mixture_var_vec,
       vb_t,
       current_prior_scale,
       sblr::core::block_residual_variance_for_marker(
        op, marker, vei_t, block_residual_control, block_ve_state),
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
      const int marker = order_t[static_cast<std::size_t>(isort)];
      sampleBetaR_ST_csr_unscaled(
       marker,
       pi_t,
       mixture_var_vec,
       vb_t,
       sblr::core::block_residual_variance_for_marker(
        op, marker, vei_t, block_residual_control, block_ve_state),
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
       const bool accepted = (estimate_maf_effect_s || use_maf_effect_s_prior_scale ||
         policy.provides_prior_scale()) ?
        attempt_ld_swap_bayesr_ST_csr(
         m,
         t,
         chain,
         it,
         vei_t,
         vb_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
         mixture_var_vec,
         estimate_maf_effect_s ? dynamic_prior_scale :
          (use_maf_effect_s_prior_scale ? prior_scale : policy.prior_scale()),
         r_t,
         b_t,
         comp_t,
         op,
         ld_swap_friends,
         gen_t,
         attempted
        ) :
        attempt_ld_swap_bayesr_ST_csr_unscaled(
         m,
         t,
         chain,
         it,
         vei_t,
         yy_vec(static_cast<arma::uword>(t)),
         ww_t,
         wy_t,
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

    if (updateB) {
     if (estimate_maf_effect_s || use_maf_effect_s_prior_scale ||
         policy.provides_prior_scale()) {
      sampleB_bayesr_ST_csr(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       mixture_var_vec,
       estimate_maf_effect_s ? dynamic_prior_scale :
        (use_maf_effect_s_prior_scale ? prior_scale : policy.prior_scale()),
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     } else {
      sampleB_bayesr_ST_csr_unscaled(
       m,
       nub,
       vb_t,
       b_t,
       comp_t,
       mixture_var_vec,
       ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
       gen_t
      );
     }
    }

    if (estimate_maf_effect_s) {
     maf_effect_s_attempted_t += 1.0;
     const bool accepted_s = update_maf_effect_s_bayesr(
      maf_effect_s_current,
      b_t,
      comp_t,
      vb_t,
      mixture_var_vec,
      maf_effect_s_log_h_row,
      maf_effect_s_prior[0],
      maf_effect_s_prior[1],
      maf_effect_s_proposal_sd,
      gen_t
     );
     if (accepted_s) maf_effect_s_accepted_t += 1.0;
    }

    policy.after_vb_update(
     b_t, comp_t, vb_t, mixture_var_vec, gen_t, it);

    const bool do_updateE =
     !block_residual_control.uses_block_variance() &&
     updateE &&
     it >= updateE_start &&
     ((it - updateE_start) % updateE_every == 0);

    if (op.uses_retained_low_rank() &&
        low_rank_residual_rebuild_every > 0 &&
        ((it + 1) % low_rank_residual_rebuild_every == 0)) {
     const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
     low_rank_residual_max_abs_drift_t = std::max(
      low_rank_residual_max_abs_drift_t, drift
     );
     ++low_rank_residual_rebuild_count_t;
    }

    double fitted_quadratic_t = 0.0;
    if (do_updateE) {
     ensure_null_effects_bayesr_ST_csr(m, t, chain, it, b_t, comp_t);
     if (!op.uses_retained_low_rank()) op.rebuild(t, wy_t, b_t, r_t);
     const BayesRUpdateEDiagnostics diag = residual_diagnostics_bayesr_ST_csr(
      op,
      t,
      m,
      nue,
      b_t,
      r_t,
      comp_t,
      wy_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t))
     );
     if (diag.sse < min_sse_t) {
      min_sse_t = diag.sse;
      min_sse_iter_t = it;
     }
     min_residual_scale_t = std::min(min_residual_scale_t, diag.residual_scale);
     max_nonzero_components_t = std::max(max_nonzero_components_t, diag.nonzero_components);
     max_abs_effect_t = std::max(max_abs_effect_t, diag.max_abs_b);
     max_fitted_quadratic_t = std::max(
      max_fitted_quadratic_t, std::abs(diag.fitted_quadratic)
     );
     ++n_updateE_t;
     check_residual_scale_bayesr_ST_csr(
      op,
      diag,
      t,
      chain,
      it,
      ve_t,
      vb_t,
      mixture_var_vec,
      b_t,
      r_t,
      sse_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      yy_vec(static_cast<arma::uword>(t)),
      n[t],
      adjE
     );
     sampleE_ST_operator_from_scale(
      nue,
      ve_t,
      diag.residual_scale,
      n[t],
      gen_t
     );
     fitted_quadratic_t = diag.fitted_quadratic;
    } else {
     fitted_quadratic_t = op.fitted_quadratic(t, b_t, wy_t, r_t);
    }

    if (block_residual_control.uses_block_variance()) {
     const bool retain_block =
      it >= nburn && ((it - nburn) % nthin == 0);
     sblr::core::update_block_residual_variance(
      op, t, it, retain_block, static_cast<double>(n[t]),
      phenotype_variance[static_cast<std::size_t>(t)], nue, b_t, r_t,
      block_residual_control, gen_t, block_ve_state);
     ve_t = sblr::core::mean_block_residual_variance(block_ve_state);
    }

    if (updatePi) samplePi_bayesr_ST_csr(comp_t, pi_t, alpha, gen_t);

    vg_t = fitted_quadratic_t / static_cast<double>(n[t]);
    vle_t = computeLE_bayesr_ST_csr(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;
    vei_t = block_residual_control.uses_block_variance() ? ve_t :
     ve_t + adjE * vg_t;
    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error("adjusted residual variance became invalid.");
    }

    vbs_task(task_u, static_cast<arma::uword>(it)) = vb_t;
    vgs_task(task_u, static_cast<arma::uword>(it)) = vg_t;
    ves_task(task_u, static_cast<arma::uword>(it)) = ve_t;
    summary_heritability_task(task_u, static_cast<arma::uword>(it)) =
     block_residual_control.uses_block_variance() ?
      vg_t / phenotype_variance[static_cast<std::size_t>(t)] :
      vg_t / (vg_t + ve_t);
    vles_task(task_u, static_cast<arma::uword>(it)) = vle_t;
    vlds_task(task_u, static_cast<arma::uword>(it)) = vld_t;
    pis_task(task_u, static_cast<arma::uword>(it)) = 1.0 - pi_t[0];
    maf_effect_s_task(task_u, static_cast<arma::uword>(it)) = maf_effect_s_current;
    policy.capture(it);

    if (it >= nburn) {
     const arma::uword draw = static_cast<arma::uword>(it - nburn);
     if (convergence_pi_count > 0) {
      if (K == 2) {
       convergence_pi_task[static_cast<std::size_t>(task)](draw, 0) = pi_t[1];
      } else for (int k = 0; k < K; ++k) {
       convergence_pi_task[static_cast<std::size_t>(task)](draw, static_cast<arma::uword>(k)) = pi_t[static_cast<std::size_t>(k)];
      }
     }
     for (int s = 0; s < convergence_marker_count; ++s) {
      const arma::uword marker = static_cast<arma::uword>(convergence_markers[static_cast<std::size_t>(s)]);
      const int component = comp_t(marker);
      if (convergence_b) convergence_b_task[static_cast<std::size_t>(task)](draw, static_cast<arma::uword>(s)) = b_t(marker);
      if (convergence_d) convergence_d_task[static_cast<std::size_t>(task)](draw, static_cast<arma::uword>(s)) = component > 0 ? 1 : 0;
      if (convergence_component) convergence_component_task[static_cast<std::size_t>(task)](draw, static_cast<arma::uword>(s)) = component;
     }
     if (convergence_component) sblr::core::capture_aggregate_component_trace(
      comp_t,m,K,static_cast<int>(draw),
      convergence_aggregate_task[static_cast<std::size_t>(task)]);
    }

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     policy.retain(it);
     nsamples_t += 1.0;
     for (int k = 0; k < K; ++k) {
      pi_mean_t(static_cast<arma::uword>(k)) += pi_t[static_cast<std::size_t>(k)];
     }
     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      const int ci = comp_t(iu);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += ci > 0 ? 1.0 : 0.0;
      component_mean_t(iu) += static_cast<double>(ci);
      comp_prob_t(iu, static_cast<arma::uword>(ci)) += 1.0;
     }
    }
   }

   policy.finish();

   if (op.uses_retained_low_rank()) {
    const double drift = op.rebuild_and_measure_drift(t, wy_t, b_t, r_t);
    low_rank_residual_max_abs_drift_t = std::max(
     low_rank_residual_max_abs_drift_t, drift
    );
    ++low_rank_residual_rebuild_count_t;
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;
   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   component_mean_t /= nsamples_t;
   comp_prob_t /= nsamples_t;
   pi_mean_t /= nsamples_t;

   bm_task.row(task_u) = bm_t;
   dm_task.row(task_u) = dm_t;
   component_mean_task.row(task_u) = component_mean_t;
   b_task.row(task_u) = b_t;
   arma::rowvec marker_residual;
   op.materialize_residual(t, r_t, marker_residual);
   r_task.row(task_u) = marker_residual;
   for (int i = 0; i < m; ++i) {
    component_task(task_u, static_cast<arma::uword>(i)) =
     static_cast<double>(comp_t(static_cast<arma::uword>(i)));
   }
   for (int k = 0; k < K; ++k) {
    final_pi_task(task_u, static_cast<arma::uword>(k)) = pi_t[static_cast<std::size_t>(k)];
    mean_pi_task(task_u, static_cast<arma::uword>(k)) = pi_mean_t(static_cast<arma::uword>(k));
   }
   final_vb_task(task_u) = vb_t;
   final_vg_task(task_u) = vg_t;
   final_ve_task(task_u) = ve_t;
   min_sse_task(task_u) = std::isfinite(min_sse_t) ? min_sse_t : NA_REAL;
   min_sse_iter_task(task_u) = min_sse_iter_t;
   min_residual_scale_task(task_u) = std::isfinite(min_residual_scale_t) ? min_residual_scale_t : NA_REAL;
   max_nonzero_components_task(task_u) = max_nonzero_components_t;
   max_abs_effect_task(task_u) = max_abs_effect_t;
   max_fitted_quadratic_task(task_u) = max_fitted_quadratic_t;
   n_updateE_task(task_u) = n_updateE_t;
   ld_swap_attempted_task(task_u) = ld_swap_attempted_t;
   ld_swap_accepted_task(task_u) = ld_swap_accepted_t;
   maf_effect_s_attempted_task(task_u) = maf_effect_s_attempted_t;
   maf_effect_s_accepted_task(task_u) = maf_effect_s_accepted_t;
   low_rank_residual_rebuild_count_task(task_u) = low_rank_residual_rebuild_count_t;
   low_rank_residual_max_abs_drift_task(task_u) = low_rank_residual_max_abs_drift_t;
   nsamples_task(task_u) = nsamples_t;
   if (block_residual_control.uses_block_variance()) {
    block_ve_posterior_mean_task.row(task_u) =
     sblr::core::posterior_mean_block_residual_variance(block_ve_state).t();
    block_ve_final_task.row(task_u) = block_ve_state.value.t();
    block_ve_resampled_task.row(task_u) =
     arma::conv_to<arma::rowvec>::from(block_ve_state.resampled.t());
    block_ve_reset_task.row(task_u) =
     arma::conv_to<arma::rowvec>::from(block_ve_state.reset_to_phenotype.t());
    if (block_residual_control.keep_history)
     block_ve_history_task[static_cast<std::size_t>(task)] =
      block_ve_state.history;
   }
   comp_prob_task[static_cast<std::size_t>(task)] = comp_prob_t;
   ncomp_task[static_cast<std::size_t>(task)] = arma::sum(comp_prob_t, 0).t();
  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = e.what();
  } catch (...) {
   failed[static_cast<std::size_t>(task)] = 1;
   errors[static_cast<std::size_t>(task)] = "unknown error";
  }
 }

 for (int task = 0; task < ntasks; ++task) {
  if (failed[static_cast<std::size_t>(task)]) {
   throw std::runtime_error(
    "stblr_cpg_omp_csr_bayesr failed for trait " +
    std::to_string(stblr_task_trait(task, nchains)) +
    ", chain " + std::to_string(stblr_task_chain(task, nchains)) +
    ": " + errors[static_cast<std::size_t>(task)]
   );
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);
 arma::mat bm(nt, m, arma::fill::zeros);
 arma::mat dm(nt, m, arma::fill::zeros);
 arma::mat bm_sd(nt, m, arma::fill::zeros);
 arma::mat dm_sd(nt, m, arma::fill::zeros);
 arma::mat bm_min(nt, m, arma::fill::zeros);
 arma::mat dm_min(nt, m, arma::fill::zeros);
 arma::mat bm_max(nt, m, arma::fill::zeros);
 arma::mat dm_max(nt, m, arma::fill::zeros);
 arma::mat component_mean(nt, m, arma::fill::zeros);
 arma::mat b_out(nt, m, arma::fill::zeros);
 arma::mat r_out(nt, m, arma::fill::zeros);
 arma::mat component_out(nt, m, arma::fill::zeros);
 arma::mat vbs(nt, trace_len, arma::fill::zeros);
 arma::mat vgs(nt, trace_len, arma::fill::zeros);
 arma::mat ves(nt, trace_len, arma::fill::zeros);
 arma::mat vle(nt, trace_len, arma::fill::zeros);
 arma::mat vld(nt, trace_len, arma::fill::zeros);
 arma::mat pis(nt, trace_len, arma::fill::zeros);
 arma::mat maf_effect_s(nt, trace_len, arma::fill::zeros);
 arma::mat final_pi(nt, K, arma::fill::zeros);
 arma::mat mean_pi(nt, K, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::mat updateE_diagnostics(ntasks, 9, arma::fill::zeros);
 arma::mat ld_swap_diagnostics(nt, 3, arma::fill::zeros);
 arma::mat ld_swap_chain_diagnostics(ntasks, 5, arma::fill::zeros);
 arma::mat low_rank_residual_diagnostics(nt, 3, arma::fill::zeros);
 arma::mat low_rank_residual_chain_diagnostics(ntasks, 5, arma::fill::zeros);
 arma::vec maf_effect_s_attempted(nt, arma::fill::zeros);
 arma::vec maf_effect_s_accepted(nt, arma::fill::zeros);
 arma::vec nsamples(nt, arma::fill::zeros);
 std::vector<arma::mat> comp_prob(static_cast<std::size_t>(nt));
 arma::mat ncomp(nt, K, arma::fill::zeros);

 for (int task = 0; task < ntasks; ++task) {
  const arma::uword task_u = static_cast<arma::uword>(task);
  const int task_trait = stblr_task_trait(task, nchains);
  const int task_chain = stblr_task_chain(task, nchains);
  updateE_diagnostics(task_u, 0) = stblr_task_trait(task, nchains);
  updateE_diagnostics(task_u, 1) = stblr_task_chain(task, nchains);
  updateE_diagnostics(task_u, 2) = n_updateE_task(task_u);
  updateE_diagnostics(task_u, 3) = min_sse_task(task_u);
  updateE_diagnostics(task_u, 4) = min_sse_iter_task(task_u);
  updateE_diagnostics(task_u, 5) = min_residual_scale_task(task_u);
  updateE_diagnostics(task_u, 6) = max_nonzero_components_task(task_u);
  updateE_diagnostics(task_u, 7) = max_abs_effect_task(task_u);
  updateE_diagnostics(task_u, 8) = max_fitted_quadratic_task(task_u);

  const double attempted = ld_swap_attempted_task(task_u);
  const double accepted = ld_swap_accepted_task(task_u);
  ld_swap_chain_diagnostics(task_u, 0) = task_trait;
  ld_swap_chain_diagnostics(task_u, 1) = task_chain;
  ld_swap_chain_diagnostics(task_u, 2) = attempted;
  ld_swap_chain_diagnostics(task_u, 3) = accepted;
  ld_swap_chain_diagnostics(task_u, 4) =
   attempted > 0.0 ? accepted / attempted : 0.0;
  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 0) += attempted;
  ld_swap_diagnostics(static_cast<arma::uword>(task_trait), 1) += accepted;
  low_rank_residual_chain_diagnostics(task_u, 0) = task_trait;
  low_rank_residual_chain_diagnostics(task_u, 1) = task_chain;
  low_rank_residual_chain_diagnostics(task_u, 2) = low_rank_residual_rebuild_every;
  low_rank_residual_chain_diagnostics(task_u, 3) =
   low_rank_residual_rebuild_count_task(task_u);
  low_rank_residual_chain_diagnostics(task_u, 4) =
   low_rank_residual_max_abs_drift_task(task_u);
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 0) =
   low_rank_residual_rebuild_every;
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 1) =
   std::max(
    low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 1),
    static_cast<double>(low_rank_residual_rebuild_count_task(task_u))
   );
  low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 2) =
   std::max(
    low_rank_residual_diagnostics(static_cast<arma::uword>(task_trait), 2),
    low_rank_residual_max_abs_drift_task(task_u)
   );
  maf_effect_s_attempted(static_cast<arma::uword>(task_trait)) +=
   maf_effect_s_attempted_task(task_u);
  maf_effect_s_accepted(static_cast<arma::uword>(task_trait)) +=
   maf_effect_s_accepted_task(task_u);
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
  bm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  dm_min.row(tu).fill(std::numeric_limits<double>::infinity());
  bm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  dm_max.row(tu).fill(-std::numeric_limits<double>::infinity());
  comp_prob[static_cast<std::size_t>(t)] = arma::mat(m, K, arma::fill::zeros);

  for (int chain = 0; chain < nchains; ++chain) {
   const int task = t * nchains + chain;
   const arma::uword task_u = static_cast<arma::uword>(task);
   bm.row(tu) += bm_task.row(task_u);
   dm.row(tu) += dm_task.row(task_u);
   component_mean.row(tu) += component_mean_task.row(task_u);
   b_out.row(tu) += b_task.row(task_u);
   r_out.row(tu) += r_task.row(task_u);
   component_out.row(tu) += component_task.row(task_u);
   vbs.row(tu) += vbs_task.row(task_u);
   vgs.row(tu) += vgs_task.row(task_u);
   ves.row(tu) += ves_task.row(task_u);
   vle.row(tu) += vles_task.row(task_u);
   vld.row(tu) += vlds_task.row(task_u);
   pis.row(tu) += pis_task.row(task_u);
   maf_effect_s.row(tu) += maf_effect_s_task.row(task_u);
   final_pi.row(tu) += final_pi_task.row(task_u);
   mean_pi.row(tu) += mean_pi_task.row(task_u);
   final_vb(tu) += final_vb_task(task_u);
   final_vg(tu) += final_vg_task(task_u);
   final_ve(tu) += final_ve_task(task_u);
   nsamples(tu) += nsamples_task(task_u);
   comp_prob[static_cast<std::size_t>(t)] += comp_prob_task[static_cast<std::size_t>(task)];
   ncomp.row(tu) += ncomp_task[static_cast<std::size_t>(task)].t();

   for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    bm_min(tu, iu) = std::min(bm_min(tu, iu), bm_task(task_u, iu));
    dm_min(tu, iu) = std::min(dm_min(tu, iu), dm_task(task_u, iu));
    bm_max(tu, iu) = std::max(bm_max(tu, iu), bm_task(task_u, iu));
    dm_max(tu, iu) = std::max(dm_max(tu, iu), dm_task(task_u, iu));
   }
  }

  bm.row(tu) *= inv_chains;
  dm.row(tu) *= inv_chains;
  component_mean.row(tu) *= inv_chains;
  b_out.row(tu) *= inv_chains;
  r_out.row(tu) *= inv_chains;
  component_out.row(tu) *= inv_chains;
  vbs.row(tu) *= inv_chains;
  vgs.row(tu) *= inv_chains;
  ves.row(tu) *= inv_chains;
  vle.row(tu) *= inv_chains;
  vld.row(tu) *= inv_chains;
  pis.row(tu) *= inv_chains;
  maf_effect_s.row(tu) *= inv_chains;
  final_pi.row(tu) *= inv_chains;
  mean_pi.row(tu) *= inv_chains;
  final_vb(tu) *= inv_chains;
  final_vg(tu) *= inv_chains;
  final_ve(tu) *= inv_chains;
  nsamples(tu) *= inv_chains;
  comp_prob[static_cast<std::size_t>(t)] *= inv_chains;
  ncomp.row(tu) *= inv_chains;

  if (nchains > 1) {
   for (int chain = 0; chain < nchains; ++chain) {
    const int task = t * nchains + chain;
    const arma::uword task_u = static_cast<arma::uword>(task);
    arma::rowvec bm_diff = bm_task.row(task_u) - bm.row(tu);
    arma::rowvec dm_diff = dm_task.row(task_u) - dm.row(tu);
    bm_sd.row(tu) += bm_diff % bm_diff;
    dm_sd.row(tu) += dm_diff % dm_diff;
   }
   bm_sd.row(tu) = arma::sqrt(bm_sd.row(tu) / static_cast<double>(nchains - 1));
   dm_sd.row(tu) = arma::sqrt(dm_sd.row(tu) / static_cast<double>(nchains - 1));
  }
 }

 arma::mat covb(nt, nt, arma::fill::zeros);
 arma::mat covg(nt, nt, arma::fill::zeros);
 arma::mat cove(nt, nt, arma::fill::zeros);
 arma::mat vb(nt, nt, arma::fill::zeros);
 arma::mat vg(nt, nt, arma::fill::zeros);
 arma::mat ve(nt, nt, arma::fill::zeros);
 for (int t = 0; t < nt; ++t) {
  const arma::uword tu = static_cast<arma::uword>(t);
  covb(tu, tu) = final_vb(tu);
  covg(tu, tu) = final_vg(tu);
  cove(tu, tu) = final_ve(tu);
  vb(tu, tu) = final_vb(tu);
  vg(tu, tu) = final_vg(tu);
  ve(tu, tu) = final_ve(tu);
 }

 // --------------------------------------------------------------------------
 // Build named raw schema v1
 // --------------------------------------------------------------------------
 return CsrBayesRExecutionResult{
  std::move(mixture_var_vec),std::move(bm_task),std::move(dm_task),std::move(component_mean_task),
  std::move(b_task),std::move(r_task),std::move(component_task),std::move(vbs_task),std::move(vgs_task),
  std::move(ves_task),std::move(vles_task),std::move(vlds_task),std::move(pis_task),std::move(maf_effect_s_task),
  std::move(final_pi_task),std::move(mean_pi_task),std::move(final_vb_task),std::move(final_vg_task),std::move(final_ve_task),
  std::move(maf_effect_s_attempted_task),std::move(maf_effect_s_accepted_task),std::move(comp_prob_task),
  std::move(convergence_pi_task),std::move(convergence_b_task),std::move(convergence_d_task),std::move(convergence_component_task),
  std::move(convergence_aggregate_task),
  std::move(ncomp_task),
  std::move(bm),std::move(dm),std::move(bm_sd),std::move(dm_sd),std::move(bm_min),std::move(dm_min),std::move(bm_max),std::move(dm_max),
  std::move(component_mean),std::move(b_out),std::move(r_out),std::move(component_out),std::move(vbs),std::move(vgs),std::move(ves),
  std::move(vle),std::move(vld),std::move(pis),std::move(maf_effect_s),std::move(final_pi),std::move(mean_pi),
  std::move(updateE_diagnostics),std::move(ld_swap_diagnostics),std::move(ld_swap_chain_diagnostics),
  std::move(low_rank_residual_diagnostics),std::move(low_rank_residual_chain_diagnostics),
  std::move(block_ve_posterior_mean_task),std::move(block_ve_final_task),
  std::move(block_ve_resampled_task),std::move(block_ve_reset_task),
  std::move(summary_heritability_task),std::move(block_ve_history_task),
  std::move(final_vb),std::move(final_vg),std::move(final_ve),std::move(maf_effect_s_attempted),std::move(maf_effect_s_accepted),
  std::move(nsamples),std::move(comp_prob),std::move(ncomp),std::move(covb),std::move(covg),std::move(cove),std::move(vb),std::move(vg),std::move(ve)};
}

#ifdef SBLR_CSR_BAYESR_DEFINE_ORDINARY_RUNNER
template <class Operator>
CsrBayesRExecutionResult run_csr_bayesr(
  CsrBayesRExecutionContext<Operator>& context) {
 CsrBayesRNoOpPolicyFactory policy_factory;
 return run_csr_bayesr_engine(context, policy_factory);
}
#endif
