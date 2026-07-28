#ifndef SBLR_BLR_CSR_GROUP_BAYESC_CORE_IMPL_H
#define SBLR_BLR_CSR_GROUP_BAYESC_CORE_IMPL_H

// Binding-neutral group BayesC implementation detail. Include only from
// st_cpg_omp_csr_group.cpp after native helpers and concrete CSR types exist.
#ifndef SBLR_CSR_GROUP_BAYESC_CORE_IMPL_TRANSLATION_UNIT
#error "blr_csr_group_bayesc_core_impl.h may only be included by st_cpg_omp_csr_group.cpp"
#endif

#include "blr_csr_group_bayesc_types.h"
#include <iostream>

namespace sblr { namespace core {

inline CsrGroupBayesCExecutionResult run_csr_group_bayesc(
 const CsrGroupBayesCExecutionContext& context
) {
 validate_csr_group_bayesc_execution_context(context);
 const int m=context.marker_count, nt=context.trait_count, ngroup=context.group_count;
 const int nit=context.iterations, nburn=context.burnin, nthin=context.thinning;
 const int ncores=context.cores, seed=context.seed;
 const bool use_d_init=context.use_initial_inclusion;
 const bool use_r_init=context.use_initial_residual;
 const bool rebuild_r_before_updateE=context.rebuild_residual_before_update;
 const bool updateGroupVb=context.update_group_multiplier;
 const bool normalize_group_vb=context.normalize_group_multiplier;
 const bool updateB=context.update_marker_variance;
 const bool updateE=context.update_residual_variance;
 const bool updatePi=context.update_group_probability;
 const bool updateLDswap=context.update_ld_swap;
 const double nub_group=context.group_multiplier_prior_df;
 const double ssb_group_prior=context.group_multiplier_prior_scale;
 const double nub=context.marker_degrees_freedom, nue=context.residual_degrees_freedom;
 const double adjE=context.residual_adjustment;
 const double ld_swap_prob=context.ld_swap_probability;
 const int ld_swap_moves=context.ld_swap_moves;
 arma::mat& wy_mat=*context.marker_score;
 const arma::mat& ww_mat=*context.marker_diagonal;
 arma::mat& b_mat=*context.effect;
 arma::mat& r_mat=*context.residual;
 arma::Mat<int>& d_mat=*context.inclusion;
 const arma::vec& yy_vec=*context.phenotype_sum_squares;
 const arma::mat& ssb_prior_mat=*context.marker_scale_prior;
 const arma::mat& sse_prior_mat=*context.residual_scale_prior;
 const arma::mat& B=*context.marker_variance_initial;
 const arma::mat& E=*context.residual_variance_initial;
 const std::vector<double>& pi=*context.global_probability;
 (void)pi;
 const std::vector<int>& n=*context.sample_size;
 const std::vector<std::vector<double>>& d_init=*context.initial_inclusion;
 const std::vector<std::vector<double>>& r_init=*context.initial_residual;
 const std::vector<std::vector<double>>& group_pi_init=*context.initial_group_probability;
 const std::vector<std::vector<double>>& group_vb_multiplier_init=*context.initial_group_multiplier;
 const arma::rowvec& prior_a=*context.probability_prior_a;
 const arma::rowvec& prior_b=*context.probability_prior_b;
 const arma::Row<int>& group=*context.marker_group;
 const arma::rowvec& group_size=*context.group_size;
 const std::vector<int>& order=*context.marker_order;
 const STLDCSR& ld=*static_cast<const STLDCSR*>(context.ld_storage);
 const LDLDFriendsGroup& ld_swap_friends=*static_cast<const LDLDFriendsGroup*>(context.ld_friends_storage);

 arma::mat bm_mat(nt, m, arma::fill::zeros);
 arma::mat dm_mat(nt, m, arma::fill::zeros);
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
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_attempted_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_accepted_vec(nt, arma::fill::zeros);

 arma::mat group_pi_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_vb_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_nincluded_mean(nt, ngroup, arma::fill::zeros);
 arma::mat group_size_mat(nt, ngroup, arma::fill::zeros);
 std::vector<arma::mat> convergence_group_pi(static_cast<std::size_t>(nt));
 std::vector<arma::mat> convergence_group_vb(static_cast<std::size_t>(nt));
 std::vector<arma::mat> convergence_b(static_cast<std::size_t>(nt));
 std::vector<arma::imat> convergence_d(static_cast<std::size_t>(nt));
 const std::vector<int>& convergence_markers=*context.convergence_markers;
 for (int t=0;t<nt;++t) {
  if (context.convergence_annotations) {
   convergence_group_pi[static_cast<std::size_t>(t)].zeros(nit,ngroup);
   convergence_group_vb[static_cast<std::size_t>(t)].zeros(nit,ngroup);
  }
  if (context.convergence_b)
   convergence_b[static_cast<std::size_t>(t)].zeros(nit,convergence_markers.size());
  if (context.convergence_d)
   convergence_d[static_cast<std::size_t>(t)].zeros(nit,convergence_markers.size());
 }

 std::vector<int> failed(static_cast<std::size_t>(nt), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(nt));
 std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
 std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);

 int nthreads = 1;
#ifdef _OPENMP
 omp_set_dynamic(0);
 nthreads = std::max(1, std::min(ncores, nt));
 omp_set_num_threads(nthreads);
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
 for (int t = 0; t < nt; ++t) {
#ifdef _OPENMP
  const double wall_start = omp_get_wtime();
  thread_used[static_cast<std::size_t>(t)] = omp_get_thread_num();
#else
  const double wall_start = 0.0;
  thread_used[static_cast<std::size_t>(t)] = 0;
#endif

  try {
   std::mt19937 gen_t(static_cast<unsigned int>(seed + 1000003 * (t + 1)));

   arma::rowvec wy_t = wy_mat.row(static_cast<arma::uword>(t));
   arma::rowvec ww_t = ww_mat.row(static_cast<arma::uword>(t));

   arma::rowvec b_t(m, arma::fill::zeros);
   for (int i = 0; i < m; ++i) b_t(static_cast<arma::uword>(i)) = b_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(i));

   arma::rowvec r_t(m, arma::fill::zeros);
   arma::Row<int> d_t(m, arma::fill::zeros);

   if (use_d_init) {
    for (int i = 0; i < m; ++i) d_t(static_cast<arma::uword>(i)) = d_init[t][i] > 0 ? 1 : 0;
   } else {
    for (int i = 0; i < m; ++i) d_t(static_cast<arma::uword>(i)) = b_t(static_cast<arma::uword>(i)) != 0.0 ? 1 : 0;
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    if (!r_t.is_finite()) throw std::runtime_error("r_init contains NaN/Inf.");
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   arma::rowvec group_pi_t(ngroup, arma::fill::zeros);
   arma::rowvec group_vb_multiplier_t(ngroup, arma::fill::ones);
   for (int g = 0; g < ngroup; ++g) {
    const double pg = group_pi_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(g)];
    const double mg = group_vb_multiplier_init[static_cast<std::size_t>(t)][static_cast<std::size_t>(g)];
    if (!std::isfinite(pg) || pg <= 0.0 || pg >= 1.0) throw std::runtime_error("group_pi_init contains invalid value.");
    if (!std::isfinite(mg) || mg <= 0.0) throw std::runtime_error("group_vb_multiplier_init contains invalid value.");
    group_pi_t(static_cast<arma::uword>(g)) = pg;
    group_vb_multiplier_t(static_cast<arma::uword>(g)) = mg;
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_group(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);

   arma::rowvec group_pi_accum(ngroup, arma::fill::zeros);
   arma::rowvec group_vb_accum(ngroup, arma::fill::zeros);
   arma::rowvec group_nincluded_accum(ngroup, arma::fill::zeros);

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const arma::uword iu = static_cast<arma::uword>(i);
     const int g = group(iu);
     const arma::uword gu = static_cast<arma::uword>(g);

     sampleBetaC_ST_csr_group(
      i,
      group_pi_t(gu),
      vb_t,
      group_vb_multiplier_t(gu),
      vei_t,
      ww_t,
      r_t,
      b_t,
      d_t,
      ld,
      gen_t
     );
    }

    if (updateLDswap && ld_swap_moves > 0 && ld_swap_prob > 0.0) {
     std::uniform_real_distribution<double> runif(0.0, 1.0);
     if (runif(gen_t) < ld_swap_prob) {
      for (int move = 0; move < ld_swap_moves; ++move) {
       ld_swap_attempted_t += 1.0;
       if (attempt_ld_swap_st_csr_group(
            m,
            vei_t,
            yy_vec(static_cast<arma::uword>(t)),
            vb_t,
            ww_t,
            wy_t,
            group,
            group_pi_t,
            group_vb_multiplier_t,
            r_t,
            b_t,
            d_t,
            ld,
            ld_swap_friends,
            gen_t
           )) {
        ld_swap_accepted_t += 1.0;
       }
      }
     }
    }

    if (updateB) {
     sampleB_ST_csr_group(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      group,
      group_vb_multiplier_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );
    }

    if (updateGroupVb) {
     sampleGroupVbMultipliers_ST_csr_group(
      m,
      ngroup,
      nub_group,
      ssb_group_prior,
      vb_t,
      b_t,
      d_t,
      group,
      group_size,
      group_vb_multiplier_t,
      normalize_group_vb,
      gen_t
     );
    }

    if (updateE) {
     if (rebuild_r_before_updateE) rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
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

    if (updatePi) {
     samplePiGroups_ST_csr_group(d_t, group, group_pi_t, prior_a, prior_b, ngroup, gen_t);
    }

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_group(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;
    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vg_t)) throw std::runtime_error("vg became NaN/Inf.");
    if (!std::isfinite(vle_t)) throw std::runtime_error("vle became NaN/Inf.");
    if (!std::isfinite(vld_t)) throw std::runtime_error("vld became NaN/Inf.");
    if (!std::isfinite(vei_t) || vei_t <= 0.0) throw std::runtime_error("vei became invalid.");

    const double global_pi = marker_weighted_mean_group_value(group_pi_t, group_size);

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    pis_t(static_cast<arma::uword>(it)) = global_pi;
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;

    if (it>=nburn) {
     const arma::uword draw=static_cast<arma::uword>(it-nburn);
     if (context.convergence_annotations) {
      convergence_group_pi[static_cast<std::size_t>(t)].row(draw)=group_pi_t;
      convergence_group_vb[static_cast<std::size_t>(t)].row(draw)=group_vb_multiplier_t;
     }
     for (std::size_t q=0;q<convergence_markers.size();++q) {
      const arma::uword marker=static_cast<arma::uword>(convergence_markers[q]);
      if (context.convergence_b)
       convergence_b[static_cast<std::size_t>(t)](draw,q)=b_t(marker);
      if (context.convergence_d)
       convergence_d[static_cast<std::size_t>(t)](draw,q)=d_t(marker);
     }
    }

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;
     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }
     group_pi_accum += group_pi_t;
     group_vb_accum += group_vb_multiplier_t;
     group_nincluded_accum += count_group_inclusions(d_t, group, ngroup);
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;
   bm_t /= nsamples_t;
   dm_t /= nsamples_t;
   group_pi_accum /= nsamples_t;
   group_vb_accum /= nsamples_t;
   group_nincluded_accum /= nsamples_t;

   bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
   dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
   b_mat.row(static_cast<arma::uword>(t)) = b_t;
   r_mat.row(static_cast<arma::uword>(t)) = r_t;
   d_mat.row(static_cast<arma::uword>(t)) = d_t;
   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
   vles_mat.row(static_cast<arma::uword>(t)) = vles_t;
   vlds_mat.row(static_cast<arma::uword>(t)) = vlds_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_pi(static_cast<arma::uword>(t)) = marker_weighted_mean_group_value(group_pi_t, group_size);
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
   ld_swap_attempted_vec(static_cast<arma::uword>(t)) = ld_swap_attempted_t;
   ld_swap_accepted_vec(static_cast<arma::uword>(t)) = ld_swap_accepted_t;

   group_pi_mean.row(static_cast<arma::uword>(t)) = group_pi_accum;
   group_vb_mean.row(static_cast<arma::uword>(t)) = group_vb_accum;
   group_nincluded_mean.row(static_cast<arma::uword>(t)) = group_nincluded_accum;
   group_size_mat.row(static_cast<arma::uword>(t)) = group_size;

#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif

  } catch (const std::exception& e) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = e.what();
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  } catch (...) {
   failed[static_cast<std::size_t>(t)] = 1;
   errors[static_cast<std::size_t>(t)] = "unknown error";
#ifdef _OPENMP
   trait_seconds[static_cast<std::size_t>(t)] = omp_get_wtime() - wall_start;
#endif
  }
 }

 for (int t = 0; t < nt; ++t) {
  if (failed[static_cast<std::size_t>(t)]) {
   throw std::runtime_error("stblr_cpg_omp_csr_group_annot failed for trait " + std::to_string(t) + ": " + errors[static_cast<std::size_t>(t)]);
  }
 }

 std::vector<std::vector<std::vector<double>>> result(28);
 for (int k = 0; k < 28; ++k) result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int k = 0; k <= 6; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(m));
  for (int k = 7; k <= 9; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nit + nburn));
  for (int k = 10; k <= 15; ++k) result[static_cast<std::size_t>(k)][ts].resize(static_cast<std::size_t>(nt));
  result[16][ts].resize(2);
  result[17][ts].resize(2);
  result[18][ts].resize(4);
  result[19][ts].resize(2);
  result[20][ts].resize(static_cast<std::size_t>(nit + nburn));
  result[21][ts].resize(static_cast<std::size_t>(nit + nburn));
  result[22][ts].resize(static_cast<std::size_t>(ngroup));
  result[23][ts].resize(static_cast<std::size_t>(ngroup));
  result[24][ts].resize(static_cast<std::size_t>(ngroup));
  result[25][ts].resize(static_cast<std::size_t>(ngroup));
  result[26][ts].resize(4);
  result[27][ts].resize(static_cast<std::size_t>(nit + nburn));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int i = 0; i < m; ++i) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword iu = static_cast<arma::uword>(i);
   const std::size_t is = static_cast<std::size_t>(i);
   result[0][ts][is] = bm_mat(tu, iu);
   result[1][ts][is] = dm_mat(tu, iu);
   result[2][ts][is] = wy_mat(tu, iu);
   result[3][ts][is] = r_mat(tu, iu);
   result[4][ts][is] = b_mat(tu, iu);
   result[5][ts][is] = static_cast<double>(d_mat(tu, iu));
   result[6][ts][is] = static_cast<double>(i);
  }
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  for (int it = 0; it < nit + nburn; ++it) {
   const arma::uword tu = static_cast<arma::uword>(t);
   const arma::uword itu = static_cast<arma::uword>(it);
   const std::size_t its = static_cast<std::size_t>(it);
   result[7][ts][its] = vbs_mat(tu, itu);
   result[8][ts][its] = vgs_mat(tu, itu);
   result[9][ts][its] = ves_mat(tu, itu);
   result[20][ts][its] = vles_mat(tu, itu);
   result[21][ts][its] = vlds_mat(tu, itu);
   result[27][ts][its] = pis_mat(tu, itu);
  }
 }

 for (int t1 = 0; t1 < nt; ++t1) {
  const std::size_t t1s = static_cast<std::size_t>(t1);
  for (int t2 = 0; t2 < nt; ++t2) {
   const std::size_t t2s = static_cast<std::size_t>(t2);
   result[10][t1s][t2s] = 0.0;
   result[11][t1s][t2s] = 0.0;
   result[12][t1s][t2s] = 0.0;
   result[13][t1s][t2s] = 0.0;
   result[14][t1s][t2s] = 0.0;
   result[15][t1s][t2s] = 0.0;
  }
  result[10][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[11][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[12][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
  result[13][t1s][t1s] = final_vb(static_cast<arma::uword>(t1));
  result[14][t1s][t1s] = final_vg(static_cast<arma::uword>(t1));
  result[15][t1s][t1s] = final_ve(static_cast<arma::uword>(t1));
 }

 for (int t = 0; t < nt; ++t) {
  const std::size_t ts = static_cast<std::size_t>(t);
  result[16][ts][0] = 1.0 - final_pi(static_cast<arma::uword>(t));
  result[16][ts][1] = final_pi(static_cast<arma::uword>(t));

  double mean_pi = 0.0;
  int npi = 0;
  for (int it = nburn; it < nit + nburn; ++it) {
   mean_pi += pis_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(it));
   ++npi;
  }
  if (npi > 0) mean_pi /= static_cast<double>(npi);
  else mean_pi = final_pi(static_cast<arma::uword>(t));
  result[17][ts][0] = 1.0 - mean_pi;
  result[17][ts][1] = mean_pi;

  result[18][ts][0] = 0.0;
  result[18][ts][1] = 0.0;
  result[18][ts][2] = trait_seconds[static_cast<std::size_t>(t)];
  result[18][ts][3] = trait_seconds[static_cast<std::size_t>(t)];
  result[19][ts][0] = nsamples_vec(static_cast<arma::uword>(t));
  result[19][ts][1] = static_cast<double>(n[t]);

  for (int g = 0; g < ngroup; ++g) {
   const std::size_t gs = static_cast<std::size_t>(g);
   result[22][ts][gs] = group_pi_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[23][ts][gs] = group_vb_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[24][ts][gs] = group_nincluded_mean(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
   result[25][ts][gs] = group_size_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(g));
  }

  const double attempted = ld_swap_attempted_vec(static_cast<arma::uword>(t));
  const double accepted = ld_swap_accepted_vec(static_cast<arma::uword>(t));
  result[26][ts][0] = attempted;
  result[26][ts][1] = accepted;
  result[26][ts][2] = attempted > 0.0 ? accepted / attempted : 0.0;
  result[26][ts][3] = 1.0;
 }

 CsrGroupBayesCExecutionResult execution_result;
 execution_result.raw=std::move(result);
 execution_result.convergence_group_pi=std::move(convergence_group_pi);
 execution_result.convergence_group_vb=std::move(convergence_group_vb);
 execution_result.convergence_b=std::move(convergence_b);
 execution_result.convergence_d=std::move(convergence_d);
 return execution_result;

}

} }

#endif
