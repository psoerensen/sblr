#ifndef SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_H
#define SBLR_BLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_H

// Binding-neutral learned-annotation BayesC implementation detail. Include only
// from st_cpg_omp_csr_annot.cpp after native helpers and concrete CSR types exist.
#ifndef SBLR_CSR_LEARNED_ANNOTATION_BAYESC_CORE_IMPL_TRANSLATION_UNIT
#error "blr_csr_learned_annotation_bayesc_core_impl.h may only be included by st_cpg_omp_csr_annot.cpp"
#endif

#include "blr_csr_learned_annotation_bayesc_types.h"
#include <iostream>

namespace sblr { namespace core {

inline CsrLearnedAnnotationBayesCExecutionResult
run_csr_learned_annotation_bayesc(
 const CsrLearnedAnnotationBayesCExecutionContext& context
) {
 validate_csr_learned_annotation_bayesc_execution_context(context);
 const int m=context.marker_count, nt=context.trait_count, K=context.annotation_count;
 const int nit=context.iterations, nburn=context.burnin, nthin=context.thinning;
 const int ncores=context.cores, seed=context.seed;
 const bool use_d_init=context.use_initial_inclusion;
 const bool use_r_init=context.use_initial_residual;
 const bool rebuild_r_before_updateE=context.rebuild_residual_before_update;
 const bool learn_pi_annot=context.learn_probability;
 const bool learn_vb_annot=context.learn_multiplier;
 const bool updateB=context.update_marker_variance;
 const bool updateE=context.update_residual_variance;
 const bool updatePi=context.update_global_probability;
 const bool updateLDswap=context.update_ld_swap;
 const int annot_update_every=context.annotation_update_every;
 const int ld_swap_moves=context.ld_swap_moves;
 const double sigma_eta_pi=context.probability_prior_sd;
 const double sigma_eta_vb=context.multiplier_prior_sd;
 const double rw_sd_eta_pi=context.probability_proposal_sd;
 const double rw_sd_eta_vb=context.multiplier_proposal_sd;
 const double pi_min=context.probability_min, pi_max=context.probability_max;
 const double vb_multiplier_min=context.multiplier_min;
 const double vb_multiplier_max=context.multiplier_max;
 const double nub=context.marker_degrees_freedom, nue=context.residual_degrees_freedom;
 const double adjE=context.residual_adjustment;
 const double pi_prior_a=context.global_probability_prior_a;
 const double pi_prior_b=context.global_probability_prior_b;
 const double ld_swap_prob=context.ld_swap_probability;
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
 const std::vector<int>& n=*context.sample_size;
 const std::vector<std::vector<double>>& d_init=*context.initial_inclusion;
 const std::vector<std::vector<double>>& r_init=*context.initial_residual;
 const arma::mat& A=*context.annotation;
 const arma::mat& eta_pi_init=*context.initial_probability_coefficient;
 const arma::mat& eta_vb_init=*context.initial_multiplier_coefficient;
 const std::vector<int>& order=*context.marker_order;
 const STLDCSR& ld=*static_cast<const STLDCSR*>(context.ld_storage);
 const LDLDFriendsAnnot& ld_swap_friends=*static_cast<const LDLDFriendsAnnot*>(context.ld_friends_storage);

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
 arma::vec final_vle(nt, arma::fill::zeros);
 arma::vec final_vld(nt, arma::fill::zeros);
 arma::vec nsamples_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_attempted_vec(nt, arma::fill::zeros);
 arma::vec ld_swap_accepted_vec(nt, arma::fill::zeros);

 arma::mat eta_pi_mean(K, nt, arma::fill::zeros);
 arma::mat eta_vb_mean(K, nt, arma::fill::zeros);
 arma::mat eta_pi_final(K, nt, arma::fill::zeros);
 arma::mat eta_vb_final(K, nt, arma::fill::zeros);
 arma::mat eta_pi_accept(nt, 2, arma::fill::zeros);
 arma::mat eta_vb_accept(nt, 2, arma::fill::zeros);

 std::vector<int> failed(static_cast<std::size_t>(nt), 0);
 std::vector<std::string> errors(static_cast<std::size_t>(nt));
 std::vector<int> thread_used(static_cast<std::size_t>(nt), 0);
 std::vector<double> trait_seconds(static_cast<std::size_t>(nt), 0.0);

 int nthreads = 1;

#ifdef _OPENMP
 omp_set_dynamic(0);
 nthreads = std::max(1, std::min(ncores, nt));
 omp_set_num_threads(nthreads);

 std::cout
 << "STBLR annotation CSR OpenMP requested threads = "
 << nthreads
 << ", omp_get_max_threads = "
 << omp_get_max_threads()
 << ", num procs = "
 << omp_get_num_procs()
 << "\n";
#endif

 std::cout
 << "STBLR annotation CSR: K=" << K
 << ", learn_pi_annot=" << learn_pi_annot
 << ", learn_vb_annot=" << learn_vb_annot
 << ", annot_update_every=" << annot_update_every
 << ", pi bounds=[" << pi_min << "," << pi_max << "]"
 << ", vb multiplier bounds=[" << vb_multiplier_min << "," << vb_multiplier_max << "]"
 << ", pi_prior_a=" << pi_prior_a
 << ", pi_prior_b=" << pi_prior_b
 << "\n";

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
     d_t(static_cast<arma::uword>(i)) =
      (b_t(static_cast<arma::uword>(i)) != 0.0) ? 1 : 0;
    }
   }

   if (use_r_init) {
    for (int i = 0; i < m; ++i) {
     r_t(static_cast<arma::uword>(i)) = static_cast<double>(r_init[t][i]);
    }

    if (!r_t.is_finite()) {
     throw std::runtime_error("stblr_cpg_omp_csr_annot: r_init contains NaN/Inf.");
    }
   } else {
    rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
   }

   double vb_t = B(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double ve_t = E(static_cast<arma::uword>(t), static_cast<arma::uword>(t));
   double vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
   double vle_t = computeLE_ST_csr_annot(m, b_t, ww_t, n[t]);
   double vld_t = vg_t - vle_t;
   double vei_t = ve_t + adjE * vg_t;

   std::vector<double> pi_t = pi;

   if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
       pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
    throw std::runtime_error(
      "invalid initial pi: pi0=" + std::to_string(pi_t[0]) +
       ", pi1=" + std::to_string(pi_t[1])
    );
   }

   {
    const double psum = pi_t[0] + pi_t[1];

    if (!std::isfinite(psum) || psum <= 0.0) {
     throw std::runtime_error("invalid initial pi sum.");
    }

    pi_t[0] /= psum;
    pi_t[1] /= psum;
   }

   arma::vec eta_pi_t(K, arma::fill::zeros);
   arma::vec eta_vb_t(K, arma::fill::zeros);

   if (learn_pi_annot) {
    eta_pi_t = eta_pi_init.col(static_cast<arma::uword>(t));
   }

   if (learn_vb_annot) {
    eta_vb_t = eta_vb_init.col(static_cast<arma::uword>(t));
   }

   arma::rowvec pi_marker_t(m, arma::fill::zeros);
   arma::rowvec vb_multiplier_t(m, arma::fill::ones);

   if (learn_pi_annot) {
    pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
   } else {
    pi_marker_t.fill(clamp_double(pi_t[1], pi_min, pi_max));
   }

   if (learn_vb_annot) {
    vb_multiplier_t = make_vb_multiplier_from_annotation(
     A,
     eta_vb_t,
     vb_multiplier_min,
     vb_multiplier_max
    );
   }

   arma::rowvec bm_t(m, arma::fill::zeros);
   arma::rowvec dm_t(m, arma::fill::zeros);
   arma::rowvec vbs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vgs_t(nit + nburn, arma::fill::zeros);
   arma::rowvec ves_t(nit + nburn, arma::fill::zeros);
   arma::rowvec pis_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vles_t(nit + nburn, arma::fill::zeros);
   arma::rowvec vlds_t(nit + nburn, arma::fill::zeros);

   arma::vec eta_pi_accum(K, arma::fill::zeros);
   arma::vec eta_vb_accum(K, arma::fill::zeros);

   int eta_pi_accepted = 0;
   int eta_pi_proposed = 0;
   int eta_vb_accepted = 0;
   int eta_vb_proposed = 0;

   double nsamples_t = 0.0;
   double ld_swap_attempted_t = 0.0;
   double ld_swap_accepted_t = 0.0;

   for (int it = 0; it < nit + nburn; ++it) {
    for (int isort = 0; isort < m; ++isort) {
     const int i = order[static_cast<std::size_t>(isort)];
     const arma::uword iu = static_cast<arma::uword>(i);

     sampleBetaC_ST_csr_prior(
      i,
      pi_marker_t(iu),
      vb_t,
      vb_multiplier_t(iu),
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
       if (attempt_ld_swap_st_csr_annot(
            m,
            vei_t,
            yy_vec(static_cast<arma::uword>(t)),
            vb_t,
            ww_t,
            wy_t,
            pi_marker_t,
            learn_pi_annot,
            vb_multiplier_t,
            learn_vb_annot,
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

    if ((it + 1) % annot_update_every == 0) {
     if (learn_pi_annot) {
      rw_update_eta_pi(
       A,
       eta_pi_t,
       d_t,
       pi_t[1],
           pi_min,
           pi_max,
           sigma_eta_pi,
           rw_sd_eta_pi,
           gen_t,
           eta_pi_accepted,
           eta_pi_proposed
      );

      pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
     }

     if (learn_vb_annot) {
      rw_update_eta_vb(
       A,
       eta_vb_t,
       b_t,
       d_t,
       vb_t,
       vb_multiplier_min,
       vb_multiplier_max,
       sigma_eta_vb,
       rw_sd_eta_vb,
       gen_t,
       eta_vb_accepted,
       eta_vb_proposed
      );

      vb_multiplier_t = make_vb_multiplier_from_annotation(
       A,
       eta_vb_t,
       vb_multiplier_min,
       vb_multiplier_max
      );
     }
    }

    if (updateB) {
     sampleB_ST_csr_prior(
      m,
      nub,
      vb_t,
      b_t,
      d_t,
      vb_multiplier_t,
      ssb_prior_mat(static_cast<arma::uword>(t), static_cast<arma::uword>(t)),
      gen_t
     );

     if (!std::isfinite(vb_t) || vb_t <= 0.0) {
      throw std::runtime_error(
        "vb became invalid after sampleB. iter=" +
         std::to_string(it) +
         ", vb=" + std::to_string(vb_t)
      );
     }
    }

    if (updateE) {
     if (rebuild_r_before_updateE) {
      rebuild_residual_st_csr(m, wy_t, ww_t, b_t, r_t, ld);
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

    if (updatePi) {
     samplePi_ST_annot(d_t, pi_t, pi_prior_a, pi_prior_b, gen_t);

     if (!std::isfinite(pi_t[0]) || !std::isfinite(pi_t[1]) ||
         pi_t[0] <= 0.0 || pi_t[1] <= 0.0) {
      throw std::runtime_error(
        "pi became invalid after samplePi. iter=" +
         std::to_string(it) +
         ", pi0=" + std::to_string(pi_t[0]) +
         ", pi1=" + std::to_string(pi_t[1])
      );
     }

     if (learn_pi_annot) {
      pi_marker_t = make_pi_from_annotation(A, eta_pi_t, pi_t[1], pi_min, pi_max);
     } else {
      pi_marker_t.fill(clamp_double(pi_t[1], pi_min, pi_max));
     }
    }

    vg_t = computeG_ST_csr(b_t, wy_t, r_t, n[t]);
    vle_t = computeLE_ST_csr_annot(m, b_t, ww_t, n[t]);
    vld_t = vg_t - vle_t;

    if (!std::isfinite(vg_t)) {
     throw std::runtime_error(
       "vg became NaN/Inf after computeG. iter=" +
        std::to_string(it)
     );
    }

    if (!std::isfinite(vle_t)) {
     throw std::runtime_error(
       "vle became NaN/Inf after computeLE. iter=" +
        std::to_string(it)
     );
    }

    if (!std::isfinite(vld_t)) {
     throw std::runtime_error(
       "vld became NaN/Inf after computeLE. iter=" +
        std::to_string(it)
     );
    }

    vei_t = ve_t + adjE * vg_t;

    if (!std::isfinite(vei_t) || vei_t <= 0.0) {
     throw std::runtime_error(
       "adjusted residual variance vei became invalid. iter=" +
        std::to_string(it) +
        ", vei=" + std::to_string(vei_t)
     );
    }

    vbs_t(static_cast<arma::uword>(it)) = vb_t;
    ves_t(static_cast<arma::uword>(it)) = ve_t;
    vgs_t(static_cast<arma::uword>(it)) = vg_t;
    pis_t(static_cast<arma::uword>(it)) = pi_t[1];
    vles_t(static_cast<arma::uword>(it)) = vle_t;
    vlds_t(static_cast<arma::uword>(it)) = vld_t;

    if ((it >= nburn) && ((it - nburn) % nthin == 0)) {
     nsamples_t += 1.0;

     for (int i = 0; i < m; ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      bm_t(iu) += b_t(iu);
      dm_t(iu) += static_cast<double>(d_t(iu));
     }

     if (K > 0) {
      eta_pi_accum += eta_pi_t;
      eta_vb_accum += eta_vb_t;
     }
    }
   }

   if (nsamples_t <= 0.0) nsamples_t = 1.0;

   bm_t /= nsamples_t;
   dm_t /= nsamples_t;

   if (K > 0) {
    eta_pi_accum /= nsamples_t;
    eta_vb_accum /= nsamples_t;
   }

   if (!bm_t.is_finite()) {
    throw std::runtime_error("posterior mean bm contains NaN/Inf.");
   }

   if (!dm_t.is_finite()) {
    throw std::runtime_error("posterior mean dm contains NaN/Inf.");
   }

   bm_mat.row(static_cast<arma::uword>(t)) = bm_t;
   dm_mat.row(static_cast<arma::uword>(t)) = dm_t;
   b_mat.row(static_cast<arma::uword>(t))  = b_t;
   r_mat.row(static_cast<arma::uword>(t))  = r_t;
   d_mat.row(static_cast<arma::uword>(t))  = d_t;

   vbs_mat.row(static_cast<arma::uword>(t)) = vbs_t;
   vgs_mat.row(static_cast<arma::uword>(t)) = vgs_t;
   ves_mat.row(static_cast<arma::uword>(t)) = ves_t;
   pis_mat.row(static_cast<arma::uword>(t)) = pis_t;
   vles_mat.row(static_cast<arma::uword>(t)) = vles_t;
   vlds_mat.row(static_cast<arma::uword>(t)) = vlds_t;

   final_vb(static_cast<arma::uword>(t)) = vb_t;
   final_ve(static_cast<arma::uword>(t)) = ve_t;
   final_vg(static_cast<arma::uword>(t)) = vg_t;
   final_pi(static_cast<arma::uword>(t)) = pi_t[1];
   final_vle(static_cast<arma::uword>(t)) = vle_t;
   final_vld(static_cast<arma::uword>(t)) = vld_t;
   nsamples_vec(static_cast<arma::uword>(t)) = nsamples_t;
   ld_swap_attempted_vec(static_cast<arma::uword>(t)) = ld_swap_attempted_t;
   ld_swap_accepted_vec(static_cast<arma::uword>(t)) = ld_swap_accepted_t;

   if (K > 0) {
    eta_pi_mean.col(static_cast<arma::uword>(t)) = eta_pi_accum;
    eta_vb_mean.col(static_cast<arma::uword>(t)) = eta_vb_accum;
    eta_pi_final.col(static_cast<arma::uword>(t)) = eta_pi_t;
    eta_vb_final.col(static_cast<arma::uword>(t)) = eta_vb_t;
   }

   eta_pi_accept(static_cast<arma::uword>(t), 0) = eta_pi_accepted;
   eta_pi_accept(static_cast<arma::uword>(t), 1) = eta_pi_proposed;
   eta_vb_accept(static_cast<arma::uword>(t), 0) = eta_vb_accepted;
   eta_vb_accept(static_cast<arma::uword>(t), 1) = eta_vb_proposed;

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

#ifdef _OPENMP
 for (int t = 0; t < nt; ++t) {
  std::cout
  << "trait " << t
  << " used thread " << thread_used[static_cast<std::size_t>(t)]
  << ", seconds = " << trait_seconds[static_cast<std::size_t>(t)]
  << "\n";
 }
#endif

 for (int t = 0; t < nt; ++t) {
  if (failed[static_cast<std::size_t>(t)]) {
   throw std::runtime_error(
     "stblr_cpg_omp_csr_annot failed for trait " +
      std::to_string(t) +
      ": " +
      errors[static_cast<std::size_t>(t)]
   );
  }
 }

 std::vector<std::vector<std::vector<double>>> result(24);

 for (int k = 0; k < 24; ++k) {
  result[static_cast<std::size_t>(k)].resize(static_cast<std::size_t>(nt));
 }

 for (int t = 0; t < nt; ++t) {
  result[0][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[1][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[2][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[3][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[4][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[5][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));
  result[6][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(m));

  result[7][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[8][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[9][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));

  result[10][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[11][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[12][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[13][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[14][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));
  result[15][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nt));

  result[16][static_cast<std::size_t>(t)].resize(2);
  result[17][static_cast<std::size_t>(t)].resize(2);

  result[18][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(std::max(K, 1)));
  result[19][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(std::max(K, 1)));
  result[20][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[21][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
  result[22][static_cast<std::size_t>(t)].resize(4);
  result[23][static_cast<std::size_t>(t)].resize(static_cast<std::size_t>(nit + nburn));
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
   result[23][ts][its] = pis_mat(tu, itu);
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

  if (K > 0) {
   for (int k = 0; k < K; ++k) {
    const std::size_t ks = static_cast<std::size_t>(k);
    result[18][ts][ks] = eta_pi_mean(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
    result[19][ts][ks] = eta_vb_mean(static_cast<arma::uword>(k), static_cast<arma::uword>(t));
   }
  } else {
   result[18][ts][0] = 0.0;
   result[19][ts][0] = 0.0;
  }

  const double attempted = ld_swap_attempted_vec(static_cast<arma::uword>(t));
  const double accepted = ld_swap_accepted_vec(static_cast<arma::uword>(t));
  result[22][ts][0] = attempted;
  result[22][ts][1] = accepted;
  result[22][ts][2] = attempted > 0.0 ? accepted / attempted : 0.0;
  result[22][ts][3] = 1.0;
 }

 std::cout << "Annotation eta_pi MH accepted/proposed by trait:\n";
 for (int t = 0; t < nt; ++t) {
  std::cout
  << "trait " << t
  << ": " << eta_pi_accept(static_cast<arma::uword>(t), 0)
  << "/" << eta_pi_accept(static_cast<arma::uword>(t), 1)
  << "\n";
 }

 std::cout << "Annotation eta_vb MH accepted/proposed by trait:\n";
 for (int t = 0; t < nt; ++t) {
  std::cout
  << "trait " << t
  << ": " << eta_vb_accept(static_cast<arma::uword>(t), 0)
  << "/" << eta_vb_accept(static_cast<arma::uword>(t), 1)
  << "\n";
 }

 CsrLearnedAnnotationBayesCExecutionResult execution_result;
 execution_result.raw=std::move(result);
 return execution_result;

}

} }

#endif
