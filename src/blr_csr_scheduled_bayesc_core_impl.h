#ifndef SBLR_BLR_CSR_SCHEDULED_BAYESC_CORE_IMPL_H
#define SBLR_BLR_CSR_SCHEDULED_BAYESC_CORE_IMPL_H

#include "blr_csr_scheduled_bayesc_types.h"

#include <utility>

// Implementation detail: included only by st_cpg_omp_csr_scheduled.cpp.
namespace sblr { namespace core {

template <class Operator>
CsrScheduledBayesCExecutionResult run_csr_scheduled_bayesc(
 const CsrScheduledBayesCExecutionContext<Operator>& context
) {
 validate_csr_scheduled_bayesc_context(context);

 const Operator& ld=context.ld;
 const arma::mat& wy_mat=context.wy;
 const arma::mat& ww_mat=context.ww;
 const arma::mat& b_mat=context.initial_b;
 const arma::vec& yy_vec=context.yy;
 const arma::mat& ssb_prior_mat=context.ssb_prior;
 const arma::mat& sse_prior_mat=context.sse_prior;
 const arma::mat& B=context.initial_B;
 const arma::mat& E=context.initial_E;
 const std::vector<double>& pi=context.initial_pi;
 const std::vector<int>& n=context.sample_sizes;
 const std::vector<int>& order=context.marker_order;
 const std::vector<std::vector<double>>& d_init=context.initial_d;
 const std::vector<std::vector<double>>& r_init=context.initial_r;
 const ScheduledExecutionControl& scheduled_control=context.scheduled;
 const std::vector<int>& chain_seeds=scheduled_control.chain_seeds;
 const int m=context.marker_count;
 const int nt=context.trait_count;
 const int nit=scheduled_control.iterations;
 const int nburn=scheduled_control.burnin;
 const int nthin=scheduled_control.thinning;
 const int nchains=scheduled_control.chains;
 const int ncores=scheduled_control.cores;
 const int seed=scheduled_control.seed;
 const int full_sweep_every=scheduled_control.sweep.full_sweep_every;
 const int null_skip_base=scheduled_control.skip.base_interval;
 const int null_skip_max=scheduled_control.skip.maximum_interval;
 const bool skip_nulls_burnin_only=scheduled_control.skip.burnin_only;
 const double candidate_threshold=scheduled_control.candidate.probability_threshold;
 const int candidate_lifetime=scheduled_control.candidate.lifetime;
 const bool wakeup_ld_neighbors=scheduled_control.neighbor.enabled;
 const double wakeup_diff_threshold=
  scheduled_control.neighbor.effect_difference_threshold;
 const int wakeup_max_neighbors=scheduled_control.neighbor.maximum_neighbors;
 const double nub=context.nub;
 const double nue=context.nue;
 const double adjE=context.adjE;
 const double pi_prior_a=context.pi_prior_a;
 const double pi_prior_b=context.pi_prior_b;
 const bool use_d_init=context.use_d_init;
 const bool use_r_init=context.use_r_init;
 const bool rebuild_r_before_updateE=context.rebuild_r_before_updateE;
 const bool updateB=context.updateB;
 const bool updateE=context.updateE;
 const bool updatePi=context.updatePi;

 const int ntasks = stblr_num_chain_tasks(nt, nchains);

 arma::mat bm_task(ntasks, m, arma::fill::zeros);
 arma::mat dm_task(ntasks, m, arma::fill::zeros);
 arma::mat b_task(ntasks, m, arma::fill::zeros);
 arma::mat r_task(ntasks, m, arma::fill::zeros);
 arma::mat d_task_double(ntasks, m, arma::fill::zeros);
 arma::mat vbs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vgs_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat ves_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat pis_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vles_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::mat vlds_task(ntasks, nit + nburn, arma::fill::zeros);
 arma::vec final_vb_task(ntasks, arma::fill::zeros);
 arma::vec final_vg_task(ntasks, arma::fill::zeros);
 arma::vec final_ve_task(ntasks, arma::fill::zeros);
 arma::vec final_pi_task(ntasks, arma::fill::zeros);
 arma::vec final_vle_task(ntasks, arma::fill::zeros);
 arma::vec final_vld_task(ntasks, arma::fill::zeros);
 arma::vec nsamples_task(ntasks, arma::fill::zeros);
 const arma::uword convergence_marker_count=static_cast<arma::uword>(context.convergence_markers.size());
 std::vector<arma::mat> convergence_b_task(static_cast<std::size_t>(ntasks));
 std::vector<arma::imat> convergence_d_task(static_cast<std::size_t>(ntasks));
 for (int task=0;task<ntasks;++task) {
  if (context.convergence_b && convergence_marker_count>0) convergence_b_task[static_cast<std::size_t>(task)].zeros(nit,convergence_marker_count);
  if (context.convergence_d && convergence_marker_count>0) convergence_d_task[static_cast<std::size_t>(task)].zeros(nit,convergence_marker_count);
 }

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
 arma::mat bm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat dm_sd_mat(nt, m, arma::fill::zeros);
 arma::mat bm_min_mat(nt, m, arma::fill::zeros);
 arma::mat dm_min_mat(nt, m, arma::fill::zeros);
 arma::mat bm_max_mat(nt, m, arma::fill::zeros);
 arma::mat dm_max_mat(nt, m, arma::fill::zeros);
 arma::mat b_out_mat(nt, m, arma::fill::zeros);
 arma::mat r_out_mat(nt, m, arma::fill::zeros);
 arma::mat d_out_mat(nt, m, arma::fill::zeros);
 arma::mat vbs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vgs_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat ves_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat pis_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vles_mat(nt, nit + nburn, arma::fill::zeros);
 arma::mat vlds_mat(nt, nit + nburn, arma::fill::zeros);
 arma::vec final_vb(nt, arma::fill::zeros);
 arma::vec final_vg(nt, arma::fill::zeros);
 arma::vec final_ve(nt, arma::fill::zeros);
 arma::vec final_pi(nt, arma::fill::zeros);
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);

 // --------------------------------------------------------------------------
 // Parallel over trait-chain tasks.
 // --------------------------------------------------------------------------

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
   unsigned int task_seed;
   if (!chain_seeds.empty()) {
    task_seed = stblr_seed_with_chain_base(
     chain_seeds[static_cast<std::size_t>(chain)],
     t
    );
   } else if (nchains == 1) {
    task_seed = stblr_trait_seed(seed, t);
   } else {
    task_seed = stblr_chain_seed(seed, t, chain);
   }
   sblr::core::ScheduledChainRng chain_rng(task_seed);
   std::mt19937& gen_t=chain_rng.engine;
   std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) {
    b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));
   }

   arma::rowvec r_t(m, arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   if (use_d_init) {
    for (int i = 0; i < m; ++i) {
     d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
    }
   } else {
    for (int i = 0; i < m; ++i) {
     d_t(static_cast<arma::uword>(i)) = b_t(static_cast<arma::uword>(i)) != 0.0 ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("r_init contains NaN/Inf.");
    }
   } else {
    //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_scheduled(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   std::vector<double> pi_t = pi;
   if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
    throw std::runtime_error("invalid initial pi.");
   }
   {
    const double psum = pi_t[0] + pi_t[1];
    if (!std::isfinite(psum) || psum <= 0.0) throw std::runtime_error("invalid initial pi sum.");
    pi_t[0] /= psum;
    pi_t[1] /= psum;
   }

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);

   const int total_it = nit + nburn;
   const int bucket_count = total_it + std::max(2, null_skip_max > 0 ? null_skip_max : 4 * null_skip_base) + null_skip_base + 10;

   std::vector<std::vector<int>> scheduled(static_cast<std::size_t>(bucket_count));
   std::vector<int> scheduled_at(static_cast<std::size_t>(m), -1);
   std::vector<int> last_updated(static_cast<std::size_t>(m), -1);
   std::vector<unsigned char> is_candidate(static_cast<std::size_t>(m), 0u);
   std::vector<int> candidate_list;
   std::vector<unsigned char> in_candidate_list(static_cast<std::size_t>(m), 0u);
   std::vector<int> active_list;
   std::vector<unsigned char> in_active_list(static_cast<std::size_t>(m), 0u);
   std::vector<int> last_interesting(static_cast<std::size_t>(m), -1000000000);

   candidate_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));
   active_list.reserve(static_cast<std::size_t>(std::min(m, 10000)));

   auto add_candidate = [&](int marker, int it) {
    if (marker < 0 || marker >= m) return;
    is_candidate[static_cast<std::size_t>(marker)] = 1u;
    last_interesting[static_cast<std::size_t>(marker)] = it;
    if (in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
     candidate_list.push_back(marker);
     in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
    }
   };

   auto add_active = [&](int marker) {
    if (marker < 0 || marker >= m) return;
    if (in_active_list[static_cast<std::size_t>(marker)] == 0u) {
     active_list.push_back(marker);
     in_active_list[static_cast<std::size_t>(marker)] = 1u;
    }
   };

   auto schedule_marker = [&](int marker, int target_it) {
    if (marker < 0 || marker >= m) return;
    if (target_it >= bucket_count) target_it = bucket_count - 1;
    if (target_it < 0) target_it = 0;
    scheduled_at[static_cast<std::size_t>(marker)] = target_it;
    scheduled[static_cast<std::size_t>(target_it)].push_back(marker);
   };

   auto wakeup_neighbors = [&](int marker, int it) {
    if (!wakeup_ld_neighbors) return;
    const uint64_t start = ld.ptr[static_cast<std::size_t>(marker)];
    const uint64_t end   = ld.ptr[static_cast<std::size_t>(marker + 1)];
    int n_wake = 0;

    for (uint64_t p = start; p < end; ++p) {
     const int j = ld.idx[static_cast<std::size_t>(p)];
     add_candidate(j, it);
     scheduled_at[static_cast<std::size_t>(j)] = -1;
     ++n_wake;

     if (wakeup_max_neighbors > 0 && n_wake >= wakeup_max_neighbors) break;
    }
   };

   for (int i = 0; i < m; ++i) {
    if (d_t(static_cast<arma::uword>(i)) > 0) {
     add_active(i);
     add_candidate(i, 0);
    } else {
     const int skip = null_skip_base + (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
     schedule_marker(i, skip);
    }
   }

   auto update_one_marker = [&](int marker, int it) {
    if (marker < 0 || marker >= m) return;
    if (last_updated[static_cast<std::size_t>(marker)] == it) return;
    last_updated[static_cast<std::size_t>(marker)] = it;

    const STMarkerUpdateResult res = sampleBetaC_ST_csr_scheduled_one(
     marker,
     pi_t,
     vb_t,
     vei_t,
     ww_t,
     r_t,
     b_t,
     d_t,
     ld,
     chain_rng
    );

    if (res.d_new > 0) {
     add_active(marker);
     add_candidate(marker, it);
     scheduled_at[static_cast<std::size_t>(marker)] = -1;
    } else if (res.p1 >= candidate_threshold) {
     add_candidate(marker, it);
     scheduled_at[static_cast<std::size_t>(marker)] = -1;
    } else {
     if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
         it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
      is_candidate[static_cast<std::size_t>(marker)] = 0u;
     }

     if (is_candidate[static_cast<std::size_t>(marker)] == 0u) {
      const int skip = adaptive_skip_length_csr_scheduled(res.p1, null_skip_base, null_skip_max) +
       (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
      schedule_marker(marker, it + skip);
     }
    }

    if (std::abs(res.diff) > wakeup_diff_threshold) {
     wakeup_neighbors(marker, it);
    }
   };

   double nsamples_t = 0.0;

   for (int it = 0; it < total_it; ++it) {
    const bool skipping_allowed =
     null_skip_base > 1 &&
     (!skip_nulls_burnin_only || it < nburn);

    const bool full_sweep =
     !skipping_allowed ||
     full_sweep_every <= 0 ||
     ((it % full_sweep_every) == 0);

    if (full_sweep) {
     for (int isort = 0; isort < m; ++isort) {
      update_one_marker(order[static_cast<std::size_t>(isort)], it);
     }
    } else {
     for (int marker : active_list) {
      if (d_t(static_cast<arma::uword>(marker)) > 0) update_one_marker(marker, it);
     }

     for (int marker : candidate_list) {
      if (is_candidate[static_cast<std::size_t>(marker)] != 0u) update_one_marker(marker, it);
     }

     if (it < bucket_count) {
      const std::vector<int>& due = scheduled[static_cast<std::size_t>(it)];
      for (int marker : due) {
       if (scheduled_at[static_cast<std::size_t>(marker)] == it &&
           d_t(static_cast<arma::uword>(marker)) == 0 &&
           is_candidate[static_cast<std::size_t>(marker)] == 0u) {
        update_one_marker(marker, it);
       }
      }
     }
    }

    // Periodic compaction of stale active/candidate lists.
    if ((it + 1) % 50 == 0) {
     std::vector<int> active_new;
     active_new.reserve(active_list.size());
     std::fill(in_active_list.begin(), in_active_list.end(), 0u);

     for (int marker : active_list) {
      if (d_t(static_cast<arma::uword>(marker)) > 0 &&
          in_active_list[static_cast<std::size_t>(marker)] == 0u) {
       active_new.push_back(marker);
       in_active_list[static_cast<std::size_t>(marker)] = 1u;
      }
     }
     active_list.swap(active_new);

     std::vector<int> cand_new;
     cand_new.reserve(candidate_list.size());
     std::fill(in_candidate_list.begin(), in_candidate_list.end(), 0u);

     for (int marker : candidate_list) {
      if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
          in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
       cand_new.push_back(marker);
       in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
      }
     }
     candidate_list.swap(cand_new);
    }

    if (updateB) {
     sampleB_ST_csr_scheduled(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error("vb became invalid after sampleB. iter=" + std::to_string(it));
     }
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      //rebuild_residual_st_csr_scheduled(m, wy_t, ww_t, b_t, r_t, ld);
      rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
     }

     sampleE_ST_csr_scheduled(
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

    if (updatePi && full_sweep) {
     samplePi_ST_scheduled(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);

     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
      throw std::runtime_error("pi became invalid after samplePi. iter=" + std::to_string(it));
     }
    }

    vg_t = computeG_ST_csr_scheduled(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_scheduled(m, b_t, ww_t, n[t]);
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
     throw std::runtime_error("adjusted residual variance vei became invalid. iter=" + std::to_string(it));
    }

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_t[1];
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;
    if (it>=nburn) {
     const arma::uword draw=static_cast<arma::uword>(it-nburn);
     for (arma::uword s=0;s<convergence_marker_count;++s) {
      const arma::uword marker=static_cast<arma::uword>(context.convergence_markers[static_cast<std::size_t>(s)]);
      if (context.convergence_b) convergence_b_task[static_cast<std::size_t>(task)](draw,s)=b_t(marker);
      if (context.convergence_d) convergence_d_task[static_cast<std::size_t>(task)](draw,s)=d_t(marker);
     }
    }

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;
   bm_t /= nsamples_t;
   dm_t /= nsamples_t;

   if (!bm_t.is_finite()) throw std::runtime_error("posterior mean bm contains NaN/Inf.");
   if (!dm_t.is_finite()) throw std::runtime_error("posterior mean dm contains NaN/Inf.");

   bm_task.row(task_u) = bm_t;
   dm_task.row(task_u) = dm_t;
   b_task.row(task_u)  = b_t;
   r_task.row(task_u)  = r_t;
   for (int i = 0; i < m; ++i) {
    d_task_double(task_u, static_cast<arma::uword>(i)) =
     static_cast<double>(d_t(static_cast<arma::uword>(i)));
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
   final_vle_task(task_u) = vle_t;
   final_vld_task(task_u) = vld_t;
   final_pi_task(task_u) = pi_t[1];
   nsamples_task(task_u) = nsamples_t;

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
   const int t = stblr_task_trait(task, nchains);
   const int chain = stblr_task_chain(task, nchains);
   throw std::runtime_error(
     "stblr_cpg_omp_csr_scheduled failed for trait " +
      std::to_string(t) +
      ", chain " +
      std::to_string(chain) +
      ": " +
      errors[static_cast<std::size_t>(task)]
   );
  }
 }

 const double inv_chains = 1.0 / static_cast<double>(nchains);

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
   b_out_mat.row(tu) += b_task.row(task_u);
   r_out_mat.row(tu) += r_task.row(task_u);
   d_out_mat.row(tu) += d_task_double.row(task_u);
   vbs_mat.row(tu) += vbs_task.row(task_u);
   vgs_mat.row(tu) += vgs_task.row(task_u);
   ves_mat.row(tu) += ves_task.row(task_u);
   pis_mat.row(tu) += pis_task.row(task_u);
   vles_mat.row(tu) += vles_task.row(task_u);
   vlds_mat.row(tu) += vlds_task.row(task_u);
   final_vb(tu) += final_vb_task(task_u);
   final_vg(tu) += final_vg_task(task_u);
   final_ve(tu) += final_ve_task(task_u);
   final_vle(tu) += final_vle_task(task_u);
   final_vld(tu) += final_vld_task(task_u);
   final_pi(tu) += final_pi_task(task_u);
   nsamples_vec(tu) += nsamples_task(task_u);

   for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    bm_min_mat(tu, iu) = std::min(bm_min_mat(tu, iu), bm_task(task_u, iu));
    dm_min_mat(tu, iu) = std::min(dm_min_mat(tu, iu), dm_task(task_u, iu));
    bm_max_mat(tu, iu) = std::max(bm_max_mat(tu, iu), bm_task(task_u, iu));
    dm_max_mat(tu, iu) = std::max(dm_max_mat(tu, iu), dm_task(task_u, iu));
   }
  }

  bm_mat.row(tu) *= inv_chains;
  dm_mat.row(tu) *= inv_chains;
  b_out_mat.row(tu) *= inv_chains;
  r_out_mat.row(tu) *= inv_chains;
  d_out_mat.row(tu) *= inv_chains;
  vbs_mat.row(tu) *= inv_chains;
  vgs_mat.row(tu) *= inv_chains;
  ves_mat.row(tu) *= inv_chains;
  pis_mat.row(tu) *= inv_chains;
  vles_mat.row(tu) *= inv_chains;
  vlds_mat.row(tu) *= inv_chains;
  final_vb(tu) *= inv_chains;
  final_vg(tu) *= inv_chains;
  final_ve(tu) *= inv_chains;
  final_vle(tu) *= inv_chains;
  final_vld(tu) *= inv_chains;
  final_pi(tu) *= inv_chains;
  nsamples_vec(tu) *= inv_chains;

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

 // --------------------------------------------------------------------------

 arma::vec mean_pi(nt,arma::fill::zeros);
 for (int t=0;t<nt;++t) {
  const arma::uword tu=static_cast<arma::uword>(t);
  double sum=0.0; int count=0;
  for (int it=nburn;it<nit+nburn;++it) {
   if ((it-nburn)%nthin!=0) continue;
   sum+=pis_mat(tu,static_cast<arma::uword>(it)); ++count;
  }
  mean_pi(tu)=count>0?sum/static_cast<double>(count):final_pi(tu);
 }

 CsrScheduledBayesCExecutionResult result;
 result.bm=std::move(bm_mat);
 result.dm=std::move(dm_mat);
 result.bm_sd=std::move(bm_sd_mat);
 result.dm_sd=std::move(dm_sd_mat);
 result.bm_min=std::move(bm_min_mat);
 result.dm_min=std::move(dm_min_mat);
 result.bm_max=std::move(bm_max_mat);
 result.dm_max=std::move(dm_max_mat);
 result.b=std::move(b_out_mat);
 result.r=std::move(r_out_mat);
 result.state=std::move(d_out_mat);
 result.vbs=std::move(vbs_mat);
 result.vgs=std::move(vgs_mat);
 result.ves=std::move(ves_mat);
 result.pis=std::move(pis_mat);
 result.vle=std::move(vles_mat);
 result.vld=std::move(vlds_mat);
 result.final_vb=std::move(final_vb);
 result.final_vg=std::move(final_vg);
 result.final_ve=std::move(final_ve);
 result.final_pi=std::move(final_pi);
 result.mean_pi=std::move(mean_pi);
 result.final_vle=std::move(final_vle);
 result.final_vld=std::move(final_vld);
 result.nsamples=std::move(nsamples_vec);
 result.task_seconds=std::move(task_seconds);
 result.task_bm=std::move(bm_task);
 result.task_dm=std::move(dm_task);
 result.task_state=std::move(d_task_double);
 result.task_vbs=std::move(vbs_task);
 result.task_vgs=std::move(vgs_task);
 result.task_ves=std::move(ves_task);
 result.task_vle=std::move(vles_task);
 result.task_vld=std::move(vlds_task);
 result.task_pis=std::move(pis_task);
 result.task_final_pi=std::move(final_pi_task);
 result.task_nsamples=std::move(nsamples_task);
 result.convergence_b=std::move(convergence_b_task);
 result.convergence_d=std::move(convergence_d_task);
 result.task_mean_pi=arma::vec(ntasks,arma::fill::zeros);
 for (int task=0;task<ntasks;++task) {
  const arma::uword task_u=static_cast<arma::uword>(task);
  double sum=0.0; int count=0;
  for (int it=nburn;it<nit+nburn;++it) {
   if ((it-nburn)%nthin!=0) continue;
   sum+=result.task_pis(task_u,static_cast<arma::uword>(it)); ++count;
  }
  result.task_mean_pi(task_u)=count>0 ? sum/static_cast<double>(count) :
   result.task_final_pi(task_u);
 }
 result.marker_count=m;
 result.trait_count=nt;
 result.chain_count=nchains;
 result.task_count=ntasks;
 return result;
}

} }

#endif  // SBLR_BLR_CSR_SCHEDULED_BAYESC_CORE_IMPL_H
