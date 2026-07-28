#ifndef SBLR_BLR_BED_SCHEDULED_BAYESC_CORE_IMPL_H
#define SBLR_BLR_BED_SCHEDULED_BAYESC_CORE_IMPL_H

// Implementation detail: included only by st_cpg_omp_individual_scheduled_chains.cpp.

namespace sblr { namespace core {

template <class PackedGenotype>
BedScheduledBayesCChainExecutionResult run_bed_scheduled_bayesc_chain(
 const BedScheduledBayesCChainExecutionContext<PackedGenotype>& context
) {
 validate_bed_scheduled_bayesc_chain_context(context);
 const PackedGenotype& G=context.genotype.storage;
 const std::vector<BedScheduledBayesCMarkerMap>& marker_maps=context.marker_maps;
 const std::vector<int>& marker_order=context.marker_order;
 const arma::mat& y_mat=context.phenotype;
 const std::vector<std::vector<double>>& b_init=context.initial_effects;
 const arma::mat& B=context.initial_B;
 const arma::mat& E=context.initial_E;
 const arma::mat& ssb_prior_mat=context.ssb_prior;
 const arma::mat& sse_prior_mat=context.sse_prior;
 const std::vector<double>& pi=context.initial_pi;
 const double nub=context.nub;
 const double nue=context.nue;
 const double adjE=context.adjE;
 const double pi_prior_a=context.pi_prior_a;
 const double pi_prior_b=context.pi_prior_b;
 const int nit=context.iterations;
 const int nburn=context.burnin;
 const int nthin=context.thinning;
 const int rebuild_every=context.rebuild_every;
 const int full_sweep_every=context.sweep.full_sweep_every;
 const int null_skip_base=context.skip.base_interval;
 const int null_skip_max=context.skip.maximum_interval;
 const double candidate_threshold=context.candidate.probability_threshold;
 const int candidate_lifetime=context.candidate.lifetime;
 const bool skip_nulls_burnin_only=context.skip.burnin_only;
 const int t=context.trait_index;
 const int chain=context.chain_index;
 const int progress_every=context.progress_every;
 const bool updateB=context.updateB;
 const bool updateE=context.updateE;
 const bool updatePi=context.updatePi;
#ifdef _OPENMP
 const double wall_start = omp_get_wtime();
#else
 const double wall_start = 0.0;
#endif

 const int m = G.m;
 BedScheduledBayesCChainExecutionResult out;
 out.bm = arma::rowvec(m, arma::fill::zeros);
 out.dm = arma::rowvec(m, arma::fill::zeros);
 out.b = arma::rowvec(m, arma::fill::zeros);
 out.d_as_double = arma::rowvec(m, arma::fill::zeros);
 out.vbs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vgs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.ves = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.pis = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vles = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vlds = arma::rowvec(nit + nburn, arma::fill::zeros);
 const arma::uword convergence_marker_count=static_cast<arma::uword>(context.convergence_markers.size());
 if (context.convergence_b && convergence_marker_count>0) out.convergence_b.zeros(nit,convergence_marker_count);
 if (context.convergence_d && convergence_marker_count>0) out.convergence_d.zeros(nit,convergence_marker_count);

 try {
  BedScheduledBayesCChainRng chain_rng(context.chain_seed);
  std::mt19937& gen_t = chain_rng.engine;
  std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));

  arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
  arma::rowvec b_t(m, arma::fill::zeros);
  arma::Row<int> d_t(m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
   b_t(static_cast<arma::uword>(j)) = b_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
   d_t(static_cast<arma::uword>(j)) = b_t(static_cast<arma::uword>(j)) != 0.0 ? 1 : 0;
  }

  arma::vec xb_t = bed_xb_from_b_scheduled_chains(G, marker_maps, marker_order, b_t);
  arma::vec e_t = y_t - xb_t;

  double vb_t = B(t, t);
  double ve_t = E(t, t);
  double vg_t = computeG_sparse_scheduled_chains(y_t, e_t);
  double vei_t = ve_t + adjE * vg_t;

  std::vector<double> pi_t = pi;
  const double psum = pi_t[0] + pi_t[1];
  if (!std::isfinite(psum) || psum <= 0.0 || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
   throw std::runtime_error("invalid initial pi.");
  }
  pi_t[0] /= psum;
  pi_t[1] /= psum;

  arma::rowvec bm_t(m, arma::fill::zeros);
  arma::rowvec dm_t(m, arma::fill::zeros);
  arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
  arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);
  arma::vec log_inv_cpo_t(G.n, arma::fill::value(-std::numeric_limits<double>::infinity()));
  double log_cpo_t = NA_REAL;
  double mean_log_cpo_t = NA_REAL;

  const int total_it = nit + nburn;
  const int bucket_count = total_it +
   std::max(2, null_skip_max > 0 ? null_skip_max : 4 * null_skip_base) +
   null_skip_base + 10;

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

  auto add_candidate = [&](int marker) {
   is_candidate[static_cast<std::size_t>(marker)] = 1u;
   if (in_candidate_list[static_cast<std::size_t>(marker)] == 0u) {
    candidate_list.push_back(marker);
    in_candidate_list[static_cast<std::size_t>(marker)] = 1u;
   }
  };

  auto add_active = [&](int marker) {
   if (in_active_list[static_cast<std::size_t>(marker)] == 0u) {
    active_list.push_back(marker);
    in_active_list[static_cast<std::size_t>(marker)] = 1u;
   }
  };

  auto schedule_marker = [&](int marker, int target_it) {
   if (target_it >= bucket_count) target_it = bucket_count - 1;
   if (target_it < 0) target_it = 0;
   scheduled_at[static_cast<std::size_t>(marker)] = target_it;
   scheduled[static_cast<std::size_t>(target_it)].push_back(marker);
  };

  for (int j = 0; j < m; ++j) {
   if (d_t(static_cast<arma::uword>(j)) > 0) {
    add_active(j);
    add_candidate(j);
    last_interesting[static_cast<std::size_t>(j)] = 0;
   } else {
    const int skip = null_skip_base + (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
    schedule_marker(j, skip);
   }
  }

  auto update_one_marker = [&](int marker, int it) {
   if (marker < 0 || marker >= m) return;
   if (last_updated[static_cast<std::size_t>(marker)] == it) return;
   last_updated[static_cast<std::size_t>(marker)] = it;

   const arma::uword ju = static_cast<arma::uword>(marker);
   double bj = b_t(ju);
   int dj = d_t(ju);

   const double p1 = sample_marker_scheduled_chains(
    G,
    marker,
    marker_maps[static_cast<std::size_t>(marker)],
               pi_t,
               vb_t,
               vei_t,
               e_t,
               bj,
               dj,
               chain_rng
   );

   b_t(ju) = bj;
   d_t(ju) = dj;

   if (dj > 0) {
    add_active(marker);
    add_candidate(marker);
    last_interesting[static_cast<std::size_t>(marker)] = it;
    scheduled_at[static_cast<std::size_t>(marker)] = -1;
    return;
   }

   if (p1 >= candidate_threshold) {
    add_candidate(marker);
    last_interesting[static_cast<std::size_t>(marker)] = it;
    scheduled_at[static_cast<std::size_t>(marker)] = -1;
    return;
   }

   if (is_candidate[static_cast<std::size_t>(marker)] != 0u &&
       it - last_interesting[static_cast<std::size_t>(marker)] > candidate_lifetime) {
    is_candidate[static_cast<std::size_t>(marker)] = 0u;
   }

   if (is_candidate[static_cast<std::size_t>(marker)] == 0u) {
    const int skip = adaptive_skip_length_scheduled_chains(p1, null_skip_base, null_skip_max) +
     (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
    schedule_marker(marker, it + skip);
   }
  };

  double nsamples_t = 0.0;

  for (int it = 0; it < total_it; ++it) {
   if (progress_every > 0 &&
       (it == 0 || ((it + 1) % progress_every == 0) || it + 1 == total_it)) {
    double n_included_progress = 0.0;
    for (arma::uword jj = 0; jj < d_t.n_elem; ++jj) {
     if (d_t(jj) > 0) n_included_progress += 1.0;
    }

    const double vle_progress = computeLE_sparse_scheduled_chains(b_t, marker_maps, G.n);
    const double vld_progress = vg_t - vle_progress;

#ifdef _OPENMP
#pragma omp critical
#endif
{
 std::cout
 << "progress chain " << chain
 << ", trait " << t
 << ": iter " << (it + 1)
 << "/" << total_it
 << ", vb=" << vb_t
 << ", ve=" << ve_t
 << ", vg=" << vg_t
 << ", vle=" << vle_progress
 << ", vld=" << vld_progress
 << ", vei=" << vei_t
 << ", pi=" << pi_t[1]
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << ", n_included=" << n_included_progress
 << ", active=" << active_list.size()
 << ", candidates=" << candidate_list.size()
 << "\n";
}
   }

   const bool skipping_allowed =
    null_skip_base > 1 &&
    (!skip_nulls_burnin_only || it < nburn);

   const bool full_sweep =
    !skipping_allowed ||
    full_sweep_every <= 0 ||
    ((it % full_sweep_every) == 0);

   if (full_sweep) {
    for (int marker : marker_order) {
     update_one_marker(marker, it);
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

   if (updateE && rebuild_every > 0 && ((it + 1) % rebuild_every == 0)) {
    xb_t = bed_xb_from_b_scheduled_chains(G, marker_maps, marker_order, b_t);
    e_t = y_t - xb_t;
   }

   if (updateB) {
    sampleB_sparse_scheduled_chains(m, nub, vb_t, b_t, d_t, ssb_prior_mat(t, t), gen_t);
    if (!std::isfinite(vb_t) || vb_t <= 0.0) throw std::runtime_error("invalid vb.");
   }

   if (updateE) {
    sampleE_sparse_scheduled_chains(nue, ve_t, e_t, sse_prior_mat(t, t), gen_t);
    if (!std::isfinite(ve_t) || ve_t <= 0.0) throw std::runtime_error("invalid ve.");
   }

   // if (updatePi) {
   //  samplePi_sparse_scheduled_chains(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);
   //  if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
   //   throw std::runtime_error("invalid pi.");
   //  }
   // }
   if (updatePi && full_sweep) {
    samplePi_sparse_scheduled_chains(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);
    if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) || pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
     throw std::runtime_error("invalid pi.");
    }
   }

   vg_t = computeG_sparse_scheduled_chains(y_t, e_t);
   const double vle_t = computeLE_sparse_scheduled_chains(b_t, marker_maps, G.n);
   const double vld_t = vg_t - vle_t;
   vei_t = ve_t + adjE * vg_t;
   if (!std::isfinite(vg_t)) throw std::runtime_error("invalid vg.");
   if (!std::isfinite(vle_t)) throw std::runtime_error("invalid vle.");
   if (!std::isfinite(vld_t)) throw std::runtime_error("invalid vld.");
   if (!std::isfinite(vei_t) || vei_t <= 0.0) throw std::runtime_error("invalid adjusted residual variance.");

   vbs_t(static_cast<arma::uword>(it)) = vb_t;
   ves_t(static_cast<arma::uword>(it)) = ve_t;
   vgs_t(static_cast<arma::uword>(it)) = vg_t;
   vles_t(static_cast<arma::uword>(it)) = vle_t;
   vlds_t(static_cast<arma::uword>(it)) = vld_t;
   pis_t(static_cast<arma::uword>(it)) = pi_t[1];
   if (it>=nburn) {
    const arma::uword draw=static_cast<arma::uword>(it-nburn);
    for (arma::uword s=0;s<convergence_marker_count;++s) {
     const arma::uword marker=static_cast<arma::uword>(context.convergence_markers[static_cast<std::size_t>(s)]);
     if (context.convergence_b) out.convergence_b(draw,s)=b_t(marker);
     if (context.convergence_d) out.convergence_d(draw,s)=d_t(marker);
    }
   }

   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
    nsamples_t += 1.0;

    update_log_inv_cpo_gaussian_scheduled_chains(e_t, ve_t, log_inv_cpo_t);

    for (int j = 0; j < m; ++j) {
     const arma::uword ju = static_cast<arma::uword>(j);
     bm_t(ju) += b_t(ju);
     dm_t(ju) += static_cast<double>(d_t(ju));
    }
   }
  }

  if (nsamples_t <= 0.0) nsamples_t = 1.0;

  log_cpo_t = compute_total_log_cpo_scheduled_chains(log_inv_cpo_t, nsamples_t);
  mean_log_cpo_t = log_cpo_t / static_cast<double>(G.n);

  bm_t /= nsamples_t;
  dm_t /= nsamples_t;

  out.bm = bm_t;
  out.dm = dm_t;
  out.b = b_t;
  for (int j = 0; j < m; ++j) {
   out.d_as_double(static_cast<arma::uword>(j)) = static_cast<double>(d_t(static_cast<arma::uword>(j)));
  }
  out.vbs = vbs_t;
  out.vgs = vgs_t;
  out.ves = ves_t;
  out.pis = pis_t;
  out.vles = vles_t;
  out.vlds = vlds_t;
  out.final_vb = vb_t;
  out.final_ve = ve_t;
  out.final_vg = vg_t;
  out.final_vle = computeLE_sparse_scheduled_chains(b_t, marker_maps, G.n);
  out.final_vld = out.final_vg - out.final_vle;
  out.log_cpo = log_cpo_t;
  out.mean_log_cpo = mean_log_cpo_t;
  out.final_pi = pi_t[1];
  out.nsamples = nsamples_t;

  double mean_pi = 0.0;
  int npi = 0;
  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_t(static_cast<arma::uword>(it));
   ++npi;
  }
  out.mean_pi = npi > 0 ? mean_pi / static_cast<double>(npi) : out.final_pi;

#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#else
  out.seconds = 0.0;
#endif

 } catch (const std::exception& e) {
  out.failed = 1;
  out.error = e.what();
#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#endif
 } catch (...) {
  out.failed = 1;
  out.error = "unknown error";
#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#endif
 }

 return out;
}

} }

#endif
