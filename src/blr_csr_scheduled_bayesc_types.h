#ifndef SBLR_BLR_CSR_SCHEDULED_BAYESC_TYPES_H
#define SBLR_BLR_CSR_SCHEDULED_BAYESC_TYPES_H

#include "blr_scheduled_execution_types.h"

#include <armadillo>
#include <stdexcept>
#include <vector>

namespace sblr { namespace core {

// All referenced storage is borrowed, immutable, and must outlive
// run_csr_scheduled_bayesc(). Mutable sampler, scheduler, and RNG state is
// allocated by the callable core for each logical trait-chain task.
template <class Operator>
struct CsrScheduledBayesCExecutionContext {
 const Operator& ld;
 const arma::mat& wy;
 const arma::mat& ww;
 const arma::mat& initial_b;
 const arma::vec& yy;
 const arma::mat& ssb_prior;
 const arma::mat& sse_prior;
 const arma::mat& initial_B;
 const arma::mat& initial_E;
 const std::vector<double>& initial_pi;
 const std::vector<int>& sample_sizes;
 const std::vector<int>& marker_order;
 const std::vector<std::vector<double>>& initial_d;
 const std::vector<std::vector<double>>& initial_r;
 const ScheduledExecutionControl& scheduled;

 int marker_count;
 int trait_count;
 double nub;
 double nue;
 double adjE;
 double pi_prior_a;
 double pi_prior_b;
 bool use_d_init;
 bool use_r_init;
 bool rebuild_r_before_updateE;
 bool updateB;
 bool updateE;
 bool updatePi;
};

struct CsrScheduledBayesCExecutionResult {
 arma::mat bm;
 arma::mat dm;
 arma::mat bm_sd;
 arma::mat dm_sd;
 arma::mat bm_min;
 arma::mat dm_min;
 arma::mat bm_max;
 arma::mat dm_max;
 arma::mat b;
 arma::mat r;
 arma::mat state;
 arma::mat vbs;
 arma::mat vgs;
 arma::mat ves;
 arma::mat pis;
 arma::mat vle;
 arma::mat vld;
 arma::vec final_vb;
 arma::vec final_vg;
 arma::vec final_ve;
 arma::vec final_pi;
 arma::vec mean_pi;
 arma::vec final_vle;
 arma::vec final_vld;
 arma::vec nsamples;
 std::vector<double> task_seconds;
 arma::mat task_bm;
 arma::mat task_dm;
 arma::mat task_state;
 arma::mat task_vbs;
 arma::mat task_vgs;
 arma::mat task_ves;
 arma::mat task_vle;
 arma::mat task_vld;
 arma::mat task_pis;
 arma::vec task_final_pi;
 arma::vec task_mean_pi;
 arma::vec task_nsamples;
 int marker_count=0;
 int trait_count=0;
 int chain_count=0;
 int task_count=0;
};

template <class Operator>
inline void validate_csr_scheduled_bayesc_context(
 const CsrScheduledBayesCExecutionContext<Operator>& x
) {
 const ScheduledExecutionControl& s=x.scheduled;
 if (x.marker_count<=0) throw std::invalid_argument("scheduled marker_count must be positive");
 if (x.trait_count<=0) throw std::invalid_argument("scheduled trait_count must be positive");
 if (s.marker_count!=static_cast<std::size_t>(x.marker_count))
  throw std::invalid_argument("scheduled marker_count does not match prepared data");
 if (s.trait_count!=static_cast<std::size_t>(x.trait_count))
  throw std::invalid_argument("scheduled trait_count does not match prepared data");
 if (s.iterations<=0) throw std::invalid_argument("scheduled iterations must be positive");
 if (s.burnin<0) throw std::invalid_argument("scheduled burnin must be non-negative");
 if (s.thinning<=0) throw std::invalid_argument("scheduled thinning must be positive");
 if (s.chains<=0) throw std::invalid_argument("scheduled chains must be positive");
 if (s.cores<=0) throw std::invalid_argument("scheduled cores must be positive");
 if (!s.chain_seeds.empty() && static_cast<int>(s.chain_seeds.size())!=s.chains)
  throw std::invalid_argument("scheduled chain_seeds length must equal chains");
 if (s.sweep.full_sweep_every<0)
  throw std::invalid_argument("scheduled full_sweep_every must be non-negative");
 if (s.skip.base_interval<=0)
  throw std::invalid_argument("scheduled null_skip_base must be positive");
 if (s.skip.maximum_interval<0)
  throw std::invalid_argument("scheduled null_skip_max must be non-negative");
 if (!(s.candidate.probability_threshold>=0.0 &&
       s.candidate.probability_threshold<=1.0))
  throw std::invalid_argument("scheduled candidate_threshold must be in [0,1]");
 if (s.candidate.lifetime<0)
  throw std::invalid_argument("scheduled candidate_lifetime must be non-negative");
 if (!(s.neighbor.effect_difference_threshold>=0.0))
  throw std::invalid_argument("scheduled wakeup_diff_threshold must be non-negative");
 if (s.neighbor.maximum_neighbors<0)
  throw std::invalid_argument("scheduled wakeup_max_neighbors must be non-negative");
 if (!s.neighbor.shared_read_only || !s.neighbor.storage_outlives_execution)
  throw std::invalid_argument("scheduled neighbor data must be borrowed immutable storage");
 if (s.neighbor.friend_marker_count!=static_cast<std::size_t>(x.marker_count))
  throw std::invalid_argument("scheduled friend-list marker dimension mismatch");
 if (x.wy.n_rows!=static_cast<arma::uword>(x.trait_count) ||
     x.wy.n_cols!=static_cast<arma::uword>(x.marker_count) ||
     x.ww.n_rows!=x.wy.n_rows || x.ww.n_cols!=x.wy.n_cols ||
     x.initial_b.n_rows!=x.wy.n_rows || x.initial_b.n_cols!=x.wy.n_cols)
  throw std::invalid_argument("scheduled summary-statistic dimensions are inconsistent");
 if (x.yy.n_elem!=static_cast<arma::uword>(x.trait_count) ||
     x.sample_sizes.size()!=static_cast<std::size_t>(x.trait_count))
  throw std::invalid_argument("scheduled trait metadata dimensions are inconsistent");
 if (x.initial_B.n_rows!=x.wy.n_rows || x.initial_B.n_cols!=x.wy.n_rows ||
     x.initial_E.n_rows!=x.wy.n_rows || x.initial_E.n_cols!=x.wy.n_rows ||
     x.ssb_prior.n_rows!=x.wy.n_rows || x.ssb_prior.n_cols!=x.wy.n_rows ||
     x.sse_prior.n_rows!=x.wy.n_rows || x.sse_prior.n_cols!=x.wy.n_rows)
  throw std::invalid_argument("scheduled prior dimensions are inconsistent");
 if (x.initial_pi.size()!=2)
  throw std::invalid_argument("scheduled initial_pi must have length two");
 if (x.marker_order.size()!=static_cast<std::size_t>(x.marker_count))
  throw std::invalid_argument("scheduled marker_order dimension mismatch");
 if (x.use_d_init && x.initial_d.size()!=static_cast<std::size_t>(x.trait_count))
  throw std::invalid_argument("scheduled initial_d trait dimension mismatch");
 if (x.use_r_init && x.initial_r.size()!=static_cast<std::size_t>(x.trait_count))
  throw std::invalid_argument("scheduled initial_r trait dimension mismatch");
 if (s.rng_ownership.engine_owner!="chain" ||
     s.rng_ownership.distribution_owner!="chain" ||
     s.rng_ownership.lifetime!="one_chain_execution" ||
     s.rng_ownership.worker_thread_owner!="none" ||
     s.rng_ownership.fit_persistent_distribution_state)
  throw std::invalid_argument("scheduled RNG ownership must be chain-local and fit-bounded");
}

} }

#endif
