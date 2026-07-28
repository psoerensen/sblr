#ifndef SBLR_BLR_BED_BAYESR_CORE_IMPL_H
#define SBLR_BLR_BED_BAYESR_CORE_IMPL_H

// Implementation detail for the public packed-BED BayesR translation unit.

namespace sblr { namespace core {

// FastPackedBedMatrix: identical definition to st_cpg_omp_individual_scheduled_chains.cpp.
// Static functions below have internal linkage so there is no ODR conflict.
// BayesR marker sampler: K-component mixture.
// Returns 1 - prob[0] (non-null probability) for the scheduling heuristic.
static inline double sample_marker_bayesr(
  const FastPackedBedMatrixBR& G,
  int marker,
  const MarkerMapBayesR& map,
  const std::vector<double>& pi,
  const std::vector<double>& c,
  double vb,
  double vei,
  arma::vec& e,
  double& b_j,
  int& d_j,
  std::mt19937& gen,
  std::uniform_real_distribution<double>& runif,
  std::normal_distribution<double>& norm01
) {
 const int K = static_cast<int>(c.size());
 const double xx = map.xx;
 const double vei_safe = std::max(vei, 1e-300);

 const double xte = br_dot_residual(G, marker, map, e.memptr());
 const double score = xte + xx * b_j;

 // K log-posteriors
 std::vector<double> logp(static_cast<std::size_t>(K));
 logp[0] = std::log(std::max(pi[0], 1e-300));
 for (int k = 1; k < K; ++k) {
  const double vbk = vb * c[static_cast<std::size_t>(k)];
  const double denom = std::max(vei_safe + xx * vbk, 1e-300);
  logp[static_cast<std::size_t>(k)] =
   std::log(std::max(pi[static_cast<std::size_t>(k)], 1e-300)) +
   0.5 * std::log(vei_safe / denom) +
   0.5 * score * score * vbk / (vei_safe * denom);
 }

 // Softmax: subtract max for numerical stability, then normalize
 double logp_max = logp[0];
 for (int k = 1; k < K; ++k)
  if (logp[static_cast<std::size_t>(k)] > logp_max)
   logp_max = logp[static_cast<std::size_t>(k)];

 std::vector<double> prob(static_cast<std::size_t>(K));
 double psum = 0.0;
 for (int k = 0; k < K; ++k) {
  prob[static_cast<std::size_t>(k)] = std::exp(logp[static_cast<std::size_t>(k)] - logp_max);
  psum += prob[static_cast<std::size_t>(k)];
 }
 for (int k = 0; k < K; ++k)
  prob[static_cast<std::size_t>(k)] /= psum;

 // Sample component via inverse CDF; default to K-1 if u falls at the tail
 int d_new = K - 1;
 const double u = runif(gen);
 double cumsum = 0.0;
 for (int k = 0; k < K - 1; ++k) {
  cumsum += prob[static_cast<std::size_t>(k)];
  if (u < cumsum) { d_new = k; break; }
 }

 double b_new = 0.0;
 if (d_new > 0) {
  const double vbk = vb * c[static_cast<std::size_t>(d_new)];
  const double lhs = xx + vei_safe / vbk;
  const double mean = score / lhs;
  const double sd = std::sqrt(vei_safe / lhs);
  b_new = mean + sd * norm01(gen);
 }

 const double diff = b_new - b_j;
 if (diff != 0.0)
  br_update_residual(G, marker, map, e.memptr(), diff);

 b_j = b_new;
 d_j = d_new;
 return 1.0 - prob[0]; // non-null probability, mirrors p1 used in scheduling
}

// BayesR vb sampler: ssb contribution for marker j is b_j^2 / c[d_j].
// BayesR pi sampler: Dirichlet(alpha[k] + n_k) via K independent Gamma samples.
static inline void samplePi_bayesr(
  const arma::Row<int>& d,
  std::vector<double>& pi,
  const std::vector<double>& alpha,
  std::mt19937& gen
) {
 const int K = static_cast<int>(pi.size());

 std::vector<double> counts(static_cast<std::size_t>(K), 0.0);
 for (arma::uword i = 0; i < d.n_elem; ++i) {
  const int di = d(i);
  if (di >= 0 && di < K)
   counts[static_cast<std::size_t>(di)] += 1.0;
 }

 double gsum = 0.0;
 std::vector<double> g(static_cast<std::size_t>(K));
 for (int k = 0; k < K; ++k) {
  const double shape = alpha[static_cast<std::size_t>(k)] + counts[static_cast<std::size_t>(k)];
  std::gamma_distribution<double> rgamma(shape, 1.0);
  g[static_cast<std::size_t>(k)] = std::max(rgamma(gen), 1e-300);
  gsum += g[static_cast<std::size_t>(k)];
 }

 for (int k = 0; k < K; ++k)
  pi[static_cast<std::size_t>(k)] = g[static_cast<std::size_t>(k)] / gsum;
}

static inline int br_adaptive_skip(double p_nonnull, int null_skip_base, int null_skip_max) {
 if (null_skip_base <= 1) return 1;

 int skip = null_skip_base;
 if (p_nonnull < 1e-6) skip = 4 * null_skip_base;
 else if (p_nonnull < 1e-5) skip = 2 * null_skip_base;
 else if (p_nonnull < 1e-4) skip = null_skip_base;
 else if (p_nonnull < 1e-3) skip = std::max(1, null_skip_base / 2);
 else skip = 1;

 if (null_skip_max > 0) skip = std::min(skip, null_skip_max);
 return std::max(1, skip);
}

template <class PackedGenotype, class MarkerMap>
BedBayesRChainExecutionResult run_bed_bayesr_chain(
 const BedBayesRChainExecutionContext<PackedGenotype,MarkerMap>& context
) {
 validate_bed_bayesr_chain_context(context);
 const PackedGenotype& G=context.genotype.storage;
 const std::vector<MarkerMap>& marker_maps=context.marker_maps;
 const std::vector<int>& marker_order=context.marker_order;
 const arma::mat& y_mat=context.phenotype;
 const std::vector<std::vector<double>>& b_init=context.initial_effects;
 const arma::mat& B=context.initial_B;
 const arma::mat& E=context.initial_E;
 const arma::mat& ssb_prior_mat=context.ssb_prior;
 const arma::mat& sse_prior_mat=context.sse_prior;
 const std::vector<double>& pi_init=context.components.initial_probabilities;
 const std::vector<double>& c=context.components.scales;
 const std::vector<double>& alpha=context.components.dirichlet_prior;
 const double nub=context.nub, nue=context.nue, adjE=context.adjE;
 const bool updateB=context.updateB, updateE=context.updateE, updatePi=context.updatePi;
 const int nit=context.iterations, nburn=context.burnin, nthin=context.thinning;
 const int rebuild_every=context.rebuild_every;
 const int full_sweep_every=context.scheduler.sweep.full_sweep_every;
 const int null_skip_base=context.scheduler.skip.base_interval;
 const int null_skip_max=context.scheduler.skip.maximum_interval;
 const double candidate_threshold=context.scheduler.candidate.probability_threshold;
 const int candidate_lifetime=context.scheduler.candidate.lifetime;
 const bool skip_nulls_burnin_only=context.scheduler.skip.burnin_only;
 const int t=context.trait_index, chain=context.chain_index;
 const int progress_every=context.progress_every;
#ifdef _OPENMP
 const double wall_start = omp_get_wtime();
#else
 const double wall_start = 0.0;
#endif

 const int m = G.m;
 const int K = static_cast<int>(c.size());
 BedBayesRChainExecutionResult out;
 out.bm = arma::rowvec(m, arma::fill::zeros);
 out.dm = arma::rowvec(m, arma::fill::zeros);
 out.component_mean = arma::rowvec(m, arma::fill::zeros);
 out.b = arma::rowvec(m, arma::fill::zeros);
 out.d_as_double = arma::rowvec(m, arma::fill::zeros);
 out.vbs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vgs = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.ves = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.pip_k = arma::mat(static_cast<arma::uword>(K), static_cast<arma::uword>(m), arma::fill::zeros);
 out.vles = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.vlds = arma::rowvec(nit + nburn, arma::fill::zeros);
 out.final_pi.assign(static_cast<std::size_t>(K), 0.0);
 out.mean_pi.assign(static_cast<std::size_t>(K), 0.0);
 const int convergence_pi_count=context.convergence_probability ? (K==2 ? 1 : K) : 0;
 const arma::uword convergence_marker_count=static_cast<arma::uword>(context.convergence_markers.size());
 if (convergence_pi_count>0) out.convergence_pi.zeros(nit,convergence_pi_count);
 if (context.convergence_b && convergence_marker_count>0) out.convergence_b.zeros(nit,convergence_marker_count);
 if (context.convergence_d && convergence_marker_count>0) out.convergence_d.zeros(nit,convergence_marker_count);
 if (context.convergence_component && convergence_marker_count>0) out.convergence_component.zeros(nit,convergence_marker_count);

 try {
  const unsigned int chain_seed=static_cast<unsigned int>(context.chain_seed);

  std::mt19937 gen_t(chain_seed);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  std::normal_distribution<double> norm01(0.0, 1.0);
  std::uniform_int_distribution<int> jitter_dist(0, std::max(0, null_skip_base - 1));

  arma::vec y_t = y_mat.col(static_cast<arma::uword>(t));
  arma::rowvec b_t(m, arma::fill::zeros);
  arma::Row<int> d_t(m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
   b_t(static_cast<arma::uword>(j)) = b_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(j)];
   d_t(static_cast<arma::uword>(j)) = b_t(static_cast<arma::uword>(j)) != 0.0 ? 1 : 0;
  }

  arma::vec xb_t = br_xb(G, marker_maps, marker_order, b_t);
  arma::vec e_t = y_t - xb_t;

  double vb_t = B(t, t);
  double ve_t = E(t, t);
  double vg_t = br_computeG(y_t, e_t);
  double vei_t = ve_t + adjE * vg_t;

  std::vector<double> pi_t = pi_init;
  {
   double psum = 0.0;
   for (int k = 0; k < K; ++k) psum += pi_t[static_cast<std::size_t>(k)];
   if (!std::isfinite(psum) || psum <= 0.0)
    throw std::runtime_error("invalid initial pi: sum <= 0.");
   for (int k = 0; k < K; ++k) {
    pi_t[static_cast<std::size_t>(k)] /= psum;
    if (!std::isfinite(pi_t[static_cast<std::size_t>(k)]) ||
        pi_t[static_cast<std::size_t>(k)] < 0.0)
     throw std::runtime_error("invalid initial pi after normalization.");
   }
  }

  arma::rowvec bm_t(m, arma::fill::zeros);
  arma::rowvec dm_t(m, arma::fill::zeros);
  arma::rowvec component_mean_t(m, arma::fill::zeros);
  arma::mat    pip_k_t(static_cast<arma::uword>(K), static_cast<arma::uword>(m), arma::fill::zeros);
  arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
  arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
  arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);
  arma::vec log_inv_cpo_t(G.n, arma::fill::value(-std::numeric_limits<double>::infinity()));

  std::vector<double> mean_pi_acc(static_cast<std::size_t>(K), 0.0);
  int npi_acc = 0;

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

   const double p_nonnull = sample_marker_bayesr(
    G,
    marker,
    marker_maps[static_cast<std::size_t>(marker)],
    pi_t,
    c,
    vb_t,
    vei_t,
    e_t,
    bj,
    dj,
    gen_t,
    runif,
    norm01
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

   if (p_nonnull >= candidate_threshold) {
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
    const int skip = br_adaptive_skip(p_nonnull, null_skip_base, null_skip_max) +
     (null_skip_base > 1 ? jitter_dist(gen_t) : 0);
    schedule_marker(marker, it + skip);
   }
  };

  double nsamples_t = 0.0;

  for (int it = 0; it < total_it; ++it) {
   if (progress_every > 0 &&
       (it == 0 || ((it + 1) % progress_every == 0) || it + 1 == total_it)) {
    double n_included_progress = 0.0;
    for (arma::uword jj = 0; jj < d_t.n_elem; ++jj)
     if (d_t(jj) > 0) n_included_progress += 1.0;

    const double vle_progress = br_computeLE(b_t, marker_maps, G.n);
    const double vld_progress = vg_t - vle_progress;
    const double pi_nonnull = 1.0 - pi_t[0];

    BedBayesRProgressEvent event;
    event.iteration=it+1; event.total_iterations=total_it;
    event.trait_index=t; event.chain_index=chain;
    event.vb=vb_t; event.ve=ve_t; event.vg=vg_t;
    event.vle=vle_progress; event.vld=vld_progress; event.vei=vei_t;
    event.pi_nonnull=pi_nonnull; event.included=n_included_progress;
    event.component_count=static_cast<std::size_t>(K);
    event.active_count=active_list.size(); event.candidate_count=candidate_list.size();
    out.progress_events.push_back(event);
   }

   const bool skipping_allowed =
    null_skip_base > 1 &&
    (!skip_nulls_burnin_only || it < nburn);

   const bool full_sweep =
    !skipping_allowed ||
    full_sweep_every <= 0 ||
    ((it % full_sweep_every) == 0);

   if (full_sweep) {
    for (int marker : marker_order)
     update_one_marker(marker, it);
   } else {
    for (int marker : active_list)
     if (d_t(static_cast<arma::uword>(marker)) > 0) update_one_marker(marker, it);

    for (int marker : candidate_list)
     if (is_candidate[static_cast<std::size_t>(marker)] != 0u) update_one_marker(marker, it);

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
    xb_t = br_xb(G, marker_maps, marker_order, b_t);
    e_t = y_t - xb_t;
   }

   if (updateB) {
    sampleB_bayesr(m, c, nub, vb_t, b_t, d_t, ssb_prior_mat(t, t), gen_t);
    if (!std::isfinite(vb_t) || vb_t <= 0.0) throw std::runtime_error("invalid vb.");
   }

   if (updateE) {
    sampleE_bayesr(nue, ve_t, e_t, sse_prior_mat(t, t), gen_t);
    if (!std::isfinite(ve_t) || ve_t <= 0.0) throw std::runtime_error("invalid ve.");
   }

   if (updatePi && full_sweep) {
    samplePi_bayesr(d_t, pi_t, alpha, gen_t);
    double psum_check = 0.0;
    for (int k = 0; k < K; ++k) psum_check += pi_t[static_cast<std::size_t>(k)];
    if (!std::isfinite(psum_check) || psum_check <= 0.0)
     throw std::runtime_error("invalid pi after sampling.");
   }

   vg_t = br_computeG(y_t, e_t);
   const double vle_t = br_computeLE(b_t, marker_maps, G.n);
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

   // Accumulate pi for mean_pi (post-burnin only, no thinning, mirrors original)
   if (it >= nburn) {
    const arma::uword draw=static_cast<arma::uword>(it-nburn);
    if (convergence_pi_count>0) {
     if (K==2) out.convergence_pi(draw,0)=pi_t[1];
     else for (int k=0;k<K;++k) out.convergence_pi(draw,static_cast<arma::uword>(k))=pi_t[static_cast<std::size_t>(k)];
    }
    for (arma::uword s=0;s<convergence_marker_count;++s) {
     const arma::uword marker=static_cast<arma::uword>(context.convergence_markers[static_cast<std::size_t>(s)]);
     const int component=d_t(marker);
     if (context.convergence_b) out.convergence_b(draw,s)=b_t(marker);
     if (context.convergence_d) out.convergence_d(draw,s)=component>0 ? 1 : 0;
     if (context.convergence_component) out.convergence_component(draw,s)=component;
    }
    for (int k = 0; k < K; ++k)
     mean_pi_acc[static_cast<std::size_t>(k)] += pi_t[static_cast<std::size_t>(k)];
    ++npi_acc;
   }

   if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
    nsamples_t += 1.0;
    br_update_log_inv_cpo(e_t, ve_t, log_inv_cpo_t);

    for (int j = 0; j < m; ++j) {
     const arma::uword ju = static_cast<arma::uword>(j);
     bm_t(ju) += b_t(ju);
     dm_t(ju) += d_t(ju) > 0 ? 1.0 : 0.0;
     component_mean_t(ju) += static_cast<double>(d_t(ju));
     pip_k_t(static_cast<arma::uword>(d_t(ju)), ju) += 1.0;
    }
   }
  }

  if (nsamples_t <= 0.0) nsamples_t = 1.0;

  const double log_cpo_t = br_compute_total_log_cpo(log_inv_cpo_t, nsamples_t);
  const double mean_log_cpo_t = log_cpo_t / static_cast<double>(G.n);

  bm_t /= nsamples_t;
  dm_t /= nsamples_t;
  component_mean_t /= nsamples_t;
  pip_k_t /= nsamples_t;

  out.bm = bm_t;
  out.dm = dm_t;
  out.component_mean = component_mean_t;
  out.b = b_t;
  for (int j = 0; j < m; ++j)
   out.d_as_double(static_cast<arma::uword>(j)) = static_cast<double>(d_t(static_cast<arma::uword>(j)));
  out.vbs = vbs_t;
  out.vgs = vgs_t;
  out.ves = ves_t;
  out.pip_k = pip_k_t;
  out.vles = vles_t;
  out.vlds = vlds_t;
  out.final_vb = vb_t;
  out.final_ve = ve_t;
  out.final_vg = vg_t;
  out.final_vle = br_computeLE(b_t, marker_maps, G.n);
  out.final_vld = out.final_vg - out.final_vle;
  out.log_cpo = log_cpo_t;
  out.mean_log_cpo = mean_log_cpo_t;
  out.nsamples = nsamples_t;

  for (int k = 0; k < K; ++k) {
   out.final_pi[static_cast<std::size_t>(k)] = pi_t[static_cast<std::size_t>(k)];
   out.mean_pi[static_cast<std::size_t>(k)] =
    npi_acc > 0
     ? mean_pi_acc[static_cast<std::size_t>(k)] / static_cast<double>(npi_acc)
     : pi_t[static_cast<std::size_t>(k)];
  }

#ifdef _OPENMP
  out.seconds = omp_get_wtime() - wall_start;
#else
  out.seconds = 0.0;
#endif

 } catch (const std::exception& ex) {
  out.failed = 1;
  out.error = ex.what();
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
