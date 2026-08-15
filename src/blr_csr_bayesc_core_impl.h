#ifndef SBLR_CORE_BLR_CSR_BAYESC_CORE_IMPL_H
#define SBLR_CORE_BLR_CSR_BAYESC_CORE_IMPL_H

// Binding-neutral implementation header. Translation units must select the
// package's established Armadillo configuration before including it. The
// ordinary binding defines SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT to emit
// the legacy run_csr_bayesc() entry point; scientific adapters instantiate the
// policy engine directly.

#include "blr_csr_bayesc_types.h"
#include "blr_scalar_execution.h"

#include "st_chain_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sblr {
namespace core {
namespace {

inline void sample_marker_scaled(
  int i, const std::vector<double>& pi, double vb,
  const arma::rowvec& prior_scale, double vei,
  const arma::rowvec& ww, arma::rowvec& r, arma::rowvec& b,
  arma::Row<int>& d, const CsrBayesCDataView& op, std::mt19937& gen
) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double wi = ww(iu);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double pi0 = std::max(pi[0], 1e-300);
  const double pi1 = std::max(pi[1], 1e-300);
  const double vei_safe = std::max(vei, 1e-300);
  const double vbi = std::max(vb * prior_scale(iu), 1e-300);
  const double score = r(iu) + wi * b(iu);
  const double denom = std::max(vei_safe + wi * vbi, 1e-300);
  const double log_bf =
    0.5 * std::log(vei_safe / denom) +
    0.5 * score * score * vbi / (vei_safe * denom);
  const double delta_log = std::log(pi0) - (std::log(pi1) + log_bf);
  double p1 = 0.0;
  if (delta_log > 35.0) p1 = 0.0;
  else if (delta_log < -35.0) p1 = 1.0;
  else p1 = 1.0 / (1.0 + std::exp(delta_log));
  const int di = (runif(gen) < p1) ? 1 : 0;
  double b_new = 0.0;
  if (di == 1) {
    const double lhs = wi + vei_safe / vbi;
    const double mean = score / lhs;
    const double sd = std::sqrt(vei_safe / lhs);
    b_new = mean + sd * norm01(gen);
  }
  const double diff = b_new - b(iu);
  if (diff != 0.0) {
    r(iu) -= wi * diff;
    op.ld.apply_offdiag(i, diff, r);
  }
  b(iu) = b_new;
  d(iu) = di;
}

inline void sample_marker_unscaled(
  int i, const std::vector<double>& pi, double vb, double vei,
  const arma::rowvec& ww, arma::rowvec& r, arma::rowvec& b,
  arma::Row<int>& d, const CsrBayesCDataView& op, std::mt19937& gen
) {
  const arma::uword iu = static_cast<arma::uword>(i);
  const double wi = ww(iu);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  std::normal_distribution<double> norm01(0.0, 1.0);
  const double pi0 = std::max(pi[0], 1e-300);
  const double pi1 = std::max(pi[1], 1e-300);
  const double vei_safe = std::max(vei, 1e-300);
  const double vbi = std::max(vb, 1e-300);
  const double score = r(iu) + wi * b(iu);
  const double denom = std::max(vei_safe + wi * vbi, 1e-300);
  const double log_bf =
    0.5 * std::log(vei_safe / denom) +
    0.5 * score * score * vbi / (vei_safe * denom);
  const double delta_log = std::log(pi0) - (std::log(pi1) + log_bf);
  double p1 = 0.0;
  if (delta_log > 35.0) p1 = 0.0;
  else if (delta_log < -35.0) p1 = 1.0;
  else p1 = 1.0 / (1.0 + std::exp(delta_log));
  const int di = (runif(gen) < p1) ? 1 : 0;
  double b_new = 0.0;
  if (di == 1) {
    const double lhs = wi + vei_safe / vbi;
    const double mean = score / lhs;
    const double sd = std::sqrt(vei_safe / lhs);
    b_new = mean + sd * norm01(gen);
  }
  const double diff = b_new - b(iu);
  if (diff != 0.0) {
    r(iu) -= wi * diff;
    op.ld.apply_offdiag(i, diff, r);
  }
  b(iu) = b_new;
  d(iu) = di;
}

inline void fill_maf_effect_scale(int m, double s, const arma::rowvec& log_h,
                          arma::rowvec& scale) {
  scale.set_size(static_cast<arma::uword>(m));
  for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    const double value = std::exp((s + 1.0) * log_h(iu));
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::runtime_error("dynamic maf_effect_s prior scale became invalid.");
    }
    scale(iu) = value;
  }
}

inline void sample_marker_variance_scaled(
  int m, double degrees_freedom, double& variance,
  const arma::rowvec& effects, const arma::Row<int>& state,
  const arma::rowvec& scale, double prior, std::mt19937& gen
) {
  double sum_squares = 0.0;
  double active = 0.0;
  for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    if (state(iu) > 0) {
      sum_squares += effects(iu) * effects(iu) / scale(iu);
      active += 1.0;
    }
  }
  const double posterior_scale = sum_squares + degrees_freedom * prior;
  std::chi_squared_distribution<double> rchisq(active + degrees_freedom);
  variance = std::max(
    posterior_scale / std::max(rchisq(gen), 1e-300), 1e-12
  );
}

inline void sample_marker_variance_unscaled(
  int m, double degrees_freedom, double& variance,
  const arma::rowvec& effects, const arma::Row<int>& state,
  double prior, std::mt19937& gen
) {
  double sum_squares = 0.0;
  double active = 0.0;
  for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    if (state(iu) > 0) {
      sum_squares += effects(iu) * effects(iu);
      active += 1.0;
    }
  }
  const double posterior_scale = sum_squares + degrees_freedom * prior;
  std::chi_squared_distribution<double> rchisq(active + degrees_freedom);
  variance = std::max(
    posterior_scale / std::max(rchisq(gen), 1e-300), 1e-12
  );
}

inline double selection_log_posterior(
  double s, const arma::rowvec& effects, const arma::Row<int>& state,
  double variance, const arma::rowvec& log_h, double lower, double upper
) {
  if (!std::isfinite(s) || s < lower || s > upper) {
    return -std::numeric_limits<double>::infinity();
  }
  const double variance_safe = std::max(variance, 1e-300);
  double value = 0.0;
  for (arma::uword marker = 0; marker < effects.n_elem; ++marker) {
    if (state(marker) <= 0) continue;
    const double log_scale = (s + 1.0) * log_h(marker);
    const double scale = std::exp(log_scale);
    if (!std::isfinite(scale) || scale <= 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    value += -0.5 * (
      log_scale + effects(marker) * effects(marker) /
      (variance_safe * scale)
    );
  }
  return value;
}

inline bool update_maf_effect_s(
  double& current, const arma::rowvec& effects,
  const arma::Row<int>& state, double variance,
  const arma::rowvec& log_h, double lower, double upper,
  double proposal_sd, std::mt19937& gen
) {
  std::normal_distribution<double> proposal(0.0, proposal_sd);
  const double candidate = current + proposal(gen);
  if (candidate < lower || candidate > upper || !std::isfinite(candidate)) {
    return false;
  }
  const double log_alpha =
    selection_log_posterior(candidate, effects, state, variance, log_h,
                            lower, upper) -
    selection_log_posterior(current, effects, state, variance, log_h,
                            lower, upper);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  if (std::log(std::max(runif(gen), 1e-300)) < log_alpha) {
    current = candidate;
    return true;
  }
  return false;
}

inline void sample_inclusion_probability(
  const arma::Row<int>& state, std::vector<double>& pi,
  double prior_active, double prior_null, std::mt19937& gen
) {
  double active = prior_active;
  double null = prior_null;
  for (arma::uword marker = 0; marker < state.n_elem; ++marker) {
    if (state(marker) > 0) active += 1.0;
    else null += 1.0;
  }
  std::gamma_distribution<double> draw_null(null, 1.0);
  std::gamma_distribution<double> draw_active(active, 1.0);
  const double value_null = std::max(draw_null(gen), 1e-300);
  const double value_active = std::max(draw_active(gen), 1e-300);
  const double total = value_null + value_active;
  pi[0] = value_null / total;
  pi[1] = value_active / total;
}

inline double compute_le(int m, const arma::rowvec& effects,
                  const arma::rowvec& diagonal, int sample_size) {
  double value = 0.0;
  for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    if (effects(iu) != 0.0) {
      value += diagonal(iu) * effects(iu) * effects(iu);
    }
  }
  return value / static_cast<double>(sample_size);
}

inline double residual_sse(int m, const arma::rowvec& effects,
                    const arma::rowvec& wy, const arma::rowvec& residual,
                    double yy) {
  double product = 0.0;
  for (int i = 0; i < m; ++i) {
    const arma::uword iu = static_cast<arma::uword>(i);
    product += effects(iu) * (residual(iu) + wy(iu));
  }
  return yy - product;
}

inline double compute_genetic_variance(const arma::rowvec& effects,
                                const arma::rowvec& wy,
                                const arma::rowvec& residual,
                                int sample_size) {
  double value = 0.0;
  for (arma::uword marker = 0; marker < effects.n_elem; ++marker) {
    value += effects(marker) * (wy(marker) - residual(marker));
  }
  return value / static_cast<double>(sample_size);
}

inline void sample_residual_variance(
  int m, double degrees_freedom, double& variance,
  const arma::rowvec& effects, const arma::rowvec& wy,
  const arma::rowvec& residual, double prior, double yy,
  int sample_size, std::mt19937& gen
) {
  const double sse = residual_sse(m, effects, wy, residual, yy);
  const double scale = sse + degrees_freedom * prior;
  if (!std::isfinite(scale) || scale <= 0.0) {
    throw std::runtime_error("sampleE_ST_csr: invalid residual scale.");
  }
  std::chi_squared_distribution<double> rchisq(sample_size + degrees_freedom);
  const double sampled = scale / std::max(rchisq(gen), 1e-300);
  if (!std::isfinite(sampled) || sampled <= 0.0) {
    throw std::runtime_error("sampleE_ST_csr: sampled ve is invalid.");
  }
  variance = std::max(sampled, 1e-12);
}

inline int count_excluded_friends(int marker, const arma::Row<int>& state,
                           const CsrBayesCLdFriendsView& friends) {
  int count = 0;
  const std::uint64_t start = friends.row_ptr[static_cast<std::size_t>(marker)];
  const std::uint64_t end = friends.row_ptr[static_cast<std::size_t>(marker + 1)];
  for (std::uint64_t position = start; position < end; ++position) {
    const int neighbour = friends.index[static_cast<std::size_t>(position)];
    if (state(static_cast<arma::uword>(neighbour)) == 0) ++count;
  }
  return count;
}

inline int collect_swap_candidates(
  int m, const arma::Row<int>& state,
  const CsrBayesCLdFriendsView& friends,
  std::vector<int>& candidates, std::vector<int>& excluded_counts
) {
  candidates.clear();
  excluded_counts.clear();
  for (int marker = 0; marker < m; ++marker) {
    if (state(static_cast<arma::uword>(marker)) <= 0) continue;
    const int count = count_excluded_friends(marker, state, friends);
    if (count > 0) {
      candidates.push_back(marker);
      excluded_counts.push_back(count);
    }
  }
  return static_cast<int>(candidates.size());
}

inline void set_marker_effect(
  int marker, double new_effect, int new_state,
  const arma::rowvec& diagonal, arma::rowvec& residual,
  arma::rowvec& effects, arma::Row<int>& state,
  const CsrBayesCDataView& op
) {
  const arma::uword marker_u = static_cast<arma::uword>(marker);
  const double difference = new_effect - effects(marker_u);
  if (difference != 0.0) {
    residual(marker_u) -= diagonal(marker_u) * difference;
    op.ld.apply_offdiag(marker, difference, residual);
  }
  effects(marker_u) = new_effect;
  state(marker_u) = new_state;
}

inline bool attempt_ld_swap(
  int m, double vei, double marker_variance, double yy,
  const arma::rowvec& diagonal, const arma::rowvec& wy,
  const arma::rowvec* maf_effect_scale, arma::rowvec& residual,
  arma::rowvec& effects, arma::Row<int>& state,
  const CsrBayesCDataView& op, const CsrBayesCLdFriendsView& friends,
  std::mt19937& gen
) {
  if (!std::isfinite(vei) || vei <= 0.0) return false;
  std::vector<int> candidates;
  std::vector<int> excluded_counts;
  const int candidate_count = collect_swap_candidates(
    m, state, friends, candidates, excluded_counts
  );
  if (candidate_count <= 0) return false;
  std::uniform_int_distribution<int> pick_candidate(0, candidate_count - 1);
  const int candidate_position = pick_candidate(gen);
  const int old_marker = candidates[static_cast<std::size_t>(candidate_position)];
  const int forward_friend_count =
    excluded_counts[static_cast<std::size_t>(candidate_position)];
  if (forward_friend_count <= 0) return false;
  std::uniform_int_distribution<int> pick_friend(0, forward_friend_count - 1);
  const int friend_position = pick_friend(gen);
  int new_marker = -1;
  int seen = 0;
  const std::uint64_t start =
    friends.row_ptr[static_cast<std::size_t>(old_marker)];
  const std::uint64_t end =
    friends.row_ptr[static_cast<std::size_t>(old_marker + 1)];
  for (std::uint64_t position = start; position < end; ++position) {
    const int neighbour = friends.index[static_cast<std::size_t>(position)];
    if (state(static_cast<arma::uword>(neighbour)) != 0) continue;
    if (seen == friend_position) {
      new_marker = neighbour;
      break;
    }
    ++seen;
  }
  if (new_marker < 0) return false;
  const arma::uword old_u = static_cast<arma::uword>(old_marker);
  const arma::uword new_u = static_cast<arma::uword>(new_marker);
  const double old_effect = effects(old_u);
  const double new_effect = effects(new_u);
  const int old_state = state(old_u);
  const int new_state = state(new_u);
  if (old_state <= 0 || new_state != 0 || old_effect == 0.0) return false;
  const double sse_old = residual_sse(m, effects, wy, residual, yy);
  if (!std::isfinite(sse_old)) return false;
  const arma::rowvec residual_old = residual;
  set_marker_effect(old_marker, 0.0, 0, diagonal, residual, effects, state, op);
  set_marker_effect(new_marker, old_effect, 1, diagonal, residual, effects,
                    state, op);
  const double sse_new = residual_sse(m, effects, wy, residual, yy);
  bool accept = false;
  if (std::isfinite(sse_new)) {
    std::vector<int> reverse_candidates;
    std::vector<int> reverse_counts;
    const int reverse_candidate_count = collect_swap_candidates(
      m, state, friends, reverse_candidates, reverse_counts
    );
    int reverse_friend_count = 0;
    for (std::size_t position = 0; position < reverse_candidates.size(); ++position) {
      if (reverse_candidates[position] == new_marker) {
        reverse_friend_count = reverse_counts[position];
        break;
      }
    }
    if (reverse_candidate_count > 0 && reverse_friend_count > 0) {
      const double log_forward =
        -std::log(static_cast<double>(candidate_count)) -
        std::log(static_cast<double>(forward_friend_count));
      const double log_reverse =
        -std::log(static_cast<double>(reverse_candidate_count)) -
        std::log(static_cast<double>(reverse_friend_count));
      double log_prior_ratio = 0.0;
      if (maf_effect_scale != nullptr) {
        const double variance_old = std::max(
          marker_variance * (*maf_effect_scale)(old_u), 1e-300
        );
        const double variance_new = std::max(
          marker_variance * (*maf_effect_scale)(new_u), 1e-300
        );
        log_prior_ratio = -0.5 * (
          std::log(variance_new / variance_old) +
          old_effect * old_effect *
          (1.0 / variance_new - 1.0 / variance_old)
        );
      }
      const double log_alpha =
        -0.5 * (sse_new - sse_old) / vei + log_prior_ratio +
        log_reverse - log_forward;
      std::uniform_real_distribution<double> runif(0.0, 1.0);
      accept = std::log(std::max(runif(gen), 1e-300)) < log_alpha;
    }
  }
  if (!accept) {
    residual = residual_old;
    effects(old_u) = old_effect;
    effects(new_u) = new_effect;
    state(old_u) = old_state;
    state(new_u) = new_state;
  }
  return accept;
}

}  // namespace

inline void validate_csr_bayesc_execution_input(
  const CsrBayesCExecutionInput& input
) {
  validate_resolved_spec(input.specification);
  const std::size_t markers = input.data.marker_count;
  const std::size_t traits = input.data.trait_count;
  if (markers == 0 || traits == 0 || input.data.wy == nullptr ||
      input.data.yy == nullptr || input.data.sample_size == nullptr) {
    throw std::invalid_argument("csr_bayesc$data view is incomplete");
  }
  if (input.data.ld.marker_count != markers) {
    throw std::invalid_argument("csr_bayesc$data LD marker count is inconsistent");
  }
  validate_sparse_ld_csr_view(input.data.ld);
  if (input.data.wy->n_rows != traits || input.data.wy->n_cols != markers ||
      input.data.yy->n_elem != traits ||
      input.data.sample_size->size() != traits) {
    throw std::invalid_argument("csr_bayesc$data dimensions are inconsistent");
  }
  if (input.priors.marker_variance == nullptr ||
      input.priors.residual_variance == nullptr ||
      input.priors.marker_scale_prior == nullptr ||
      input.priors.residual_scale_prior == nullptr ||
      input.initial.effects == nullptr || input.marker_order == nullptr) {
    throw std::invalid_argument("csr_bayesc typed input is incomplete");
  }
  if (input.marker_order->size() != markers) {
    throw std::invalid_argument("csr_bayesc marker order must have length m");
  }
  const CsrBayesCControls& controls = input.controls;
  if (controls.nit <= 0 || controls.nburn < 0 || controls.nthin <= 0 ||
      controls.nchains <= 0 || controls.ncores <= 0) {
    throw std::invalid_argument("csr_bayesc MCMC controls are invalid");
  }
  if (!controls.chain_seeds.empty() &&
      controls.chain_seeds.size() != static_cast<std::size_t>(controls.nchains)) {
    throw std::invalid_argument("csr_bayesc chain seeds must match nchains");
  }
  const std::size_t task_count = traits *
    static_cast<std::size_t>(controls.nchains);
  if (controls.seed_contract_version == 1 &&
      controls.task_seeds.size() != task_count) {
    throw std::invalid_argument(
      "csr_bayesc Phase 3 task seeds must match the logical task count"
    );
  }
  if (controls.retention_contract_version == 1) {
    int previous = 0;
    for (int index : controls.retained_transition_indices) {
      if (index <= previous || index < 1 || index > controls.nit) {
        throw std::invalid_argument(
          "csr_bayesc Phase 3 retained indices must be ordered post-burn transitions"
        );
      }
      previous = index;
    }
  }
  if (controls.use_fixed_maf_effect_scale &&
      (controls.fixed_maf_effect_scale == nullptr ||
       controls.fixed_maf_effect_scale->n_elem != markers)) {
    throw std::invalid_argument("csr_bayesc fixed selection scale is invalid");
  }
  if (controls.estimate_maf_effect_s &&
      (controls.maf_effect_s_log_h == nullptr ||
       controls.maf_effect_s_log_h->n_elem != markers)) {
    throw std::invalid_argument("csr_bayesc maf_effect_s log-h view is invalid");
  }
  if (controls.update_ld_swap &&
      (input.ld_friends.row_ptr == nullptr ||
       input.ld_friends.row_ptr_size != markers + 1 ||
       (input.ld_friends.friend_count > 0 && input.ld_friends.index == nullptr))) {
    throw std::invalid_argument("csr_bayesc LD-swap friend view is invalid");
  }
  for (int marker : controls.convergence_markers) if (marker < 0 ||
      static_cast<std::size_t>(marker) >= markers) {
    throw std::invalid_argument("csr_bayesc convergence marker index is out of range");
  }
}

// The policy surface is intentionally narrow. It supplies an optional dynamic
// marker scale and one post-vb update point. Ordinary BayesC compiles against
// the no-op policy below, which performs no draws and owns no mutable state.
struct CsrBayesCNoOpPolicy {
  bool provides_prior_scale() const noexcept { return false; }
  const arma::rowvec& prior_scale() const {
    throw std::logic_error("ordinary BayesC has no policy-owned prior scale");
  }
  void after_vb_update(
    const arma::rowvec&, const arma::Row<int>&, double, std::mt19937&, int
  ) noexcept {}
  void capture(int) noexcept {}
  void retain(int) noexcept {}
  void finish() noexcept {}
};

struct CsrBayesCNoOpPolicyFactory {
  CsrBayesCNoOpPolicy make(int, int, int, int) const noexcept {
    return CsrBayesCNoOpPolicy{};
  }
};

template <class PolicyFactory>
CsrBayesCResult run_csr_bayesc_engine(
  const CsrBayesCExecutionInput& input,
  PolicyFactory& policy_factory
) {
  validate_csr_bayesc_execution_input(input);
  const int m = static_cast<int>(input.data.marker_count);
  const int nt = static_cast<int>(input.data.trait_count);
  const int n_trace = input.controls.nit + input.controls.nburn;
  const std::vector<ScalarChainTask> tasks = make_scalar_chain_tasks(
    input.data.trait_count,
    static_cast<std::size_t>(input.controls.nchains)
  );
  const int ntasks = static_cast<int>(tasks.size());
  std::vector<CsrBayesCChainResult> chains(static_cast<std::size_t>(ntasks));
  std::vector<ScalarChainExecutionStatus> statuses(
    static_cast<std::size_t>(ntasks)
  );
  int nthreads = 1;
#ifdef _OPENMP
  omp_set_dynamic(0);
  nthreads = stblr_num_threads_for_tasks(input.controls.ncores, ntasks);
  omp_set_num_threads(nthreads);
#endif

#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
  for (int task = 0; task < ntasks; ++task) {
    const ScalarChainTask& identity = tasks[static_cast<std::size_t>(task)];
    const int trait = static_cast<int>(identity.trait_index);
    const int chain = static_cast<int>(identity.chain_index);
    ScalarChainExecutionStatus& status =
      statuses[static_cast<std::size_t>(task)];
    status.task = identity;
#ifdef _OPENMP
    const double wall_start = omp_get_wtime();
#else
    const double wall_start = 0.0;
#endif
    CsrBayesCChainResult& output = chains[static_cast<std::size_t>(task)];
#ifdef _OPENMP
    output.thread_used = omp_get_thread_num();
    output.team_size_used = omp_get_num_threads();
#endif
    try {
      const unsigned int task_seed =
        input.controls.seed_contract_version == 1
          ? input.controls.task_seeds[static_cast<std::size_t>(task)]
          : resolve_scalar_chain_seed(
              input.controls.seed,
              static_cast<std::size_t>(input.controls.nchains),
              input.controls.chain_seeds,
              identity
            );
      status.seed = task_seed;
      std::mt19937 gen(task_seed);
      auto policy = policy_factory.make(task, trait, chain, m);
      arma::rowvec wy = input.data.wy->row(static_cast<arma::uword>(trait));
      const arma::rowvec& diagonal = input.data.ld.diag();
      arma::rowvec effects(m, arma::fill::zeros);
      for (int marker = 0; marker < m; ++marker) {
        effects(static_cast<arma::uword>(marker)) =
          (*input.initial.effects)(static_cast<arma::uword>(trait),
                                   static_cast<arma::uword>(marker));
      }
      arma::rowvec residual(m, arma::fill::zeros);
      arma::Row<int> state(m, arma::fill::zeros);
      if (input.initial.use_inclusion) {
        for (int marker = 0; marker < m; ++marker) {
          state(static_cast<arma::uword>(marker)) =
            (*input.initial.inclusion)[static_cast<std::size_t>(trait)]
              [static_cast<std::size_t>(marker)] > 0 ? 1 : 0;
        }
      } else {
        for (int marker = 0; marker < m; ++marker) {
          state(static_cast<arma::uword>(marker)) =
            effects(static_cast<arma::uword>(marker)) != 0.0 ? 1 : 0;
        }
      }
      if (input.initial.use_residual) {
        for (int marker = 0; marker < m; ++marker) {
          residual(static_cast<arma::uword>(marker)) =
            (*input.initial.residual)[static_cast<std::size_t>(trait)]
              [static_cast<std::size_t>(marker)];
        }
        if (!residual.is_finite()) {
          throw std::runtime_error(
            "stblr_cpg_omp_csr_state: r_init contains NaN/Inf."
          );
        }
      } else {
        input.data.ld.rebuild(wy, effects, residual);
      }
      double marker_variance = (*input.priors.marker_variance)(trait, trait);
      double residual_variance = (*input.priors.residual_variance)(trait, trait);
      double genetic_variance = 0.0;
      double le_variance = compute_le(
        m, effects, diagonal,
        (*input.data.sample_size)[static_cast<std::size_t>(trait)]
      );
      double ld_variance = genetic_variance - le_variance;
      double adjusted_residual = residual_variance +
        input.controls.residual_adjustment * genetic_variance;
      std::vector<double> pi = input.priors.inclusion_probability;
      if (!std::isfinite(pi[0]) || !std::isfinite(pi[1]) ||
          pi[0] <= 0.0 || pi[1] <= 0.0) {
        throw std::runtime_error(
          "invalid initial pi: pi0=" + std::to_string(pi[0]) +
          ", pi1=" + std::to_string(pi[1])
        );
      }
      const double pi_sum = pi[0] + pi[1];
      if (!std::isfinite(pi_sum) || pi_sum <= 0.0) {
        throw std::runtime_error("invalid initial pi sum.");
      }
      pi[0] /= pi_sum;
      pi[1] /= pi_sum;

      output.marker_mean.zeros(m);
      output.marker_pip.zeros(m);
      output.marker_variance_trace.zeros(n_trace);
      output.genetic_variance_trace.zeros(n_trace);
      output.residual_variance_trace.zeros(n_trace);
      output.inclusion_trace.zeros(n_trace);
      output.le_variance_trace.zeros(n_trace);
      output.ld_variance_trace.zeros(n_trace);
      output.maf_effect_s_trace.zeros(n_trace);
      const arma::uword selected_count=static_cast<arma::uword>(input.controls.convergence_markers.size());
      if (input.controls.convergence_b && selected_count>0) output.convergence_b.zeros(input.controls.nit,selected_count);
      if (input.controls.convergence_d && selected_count>0) output.convergence_d.zeros(input.controls.nit,selected_count);
      arma::rowvec dynamic_maf_effect_scale;
      double maf_effect_s = input.controls.maf_effect_s_initial;

      for (int iteration = 0; iteration < n_trace; ++iteration) {
        if (input.controls.estimate_maf_effect_s) {
          fill_maf_effect_scale(
            m, maf_effect_s, *input.controls.maf_effect_s_log_h,
            dynamic_maf_effect_scale
          );
        }
        if (input.controls.estimate_maf_effect_s ||
            input.controls.use_fixed_maf_effect_scale ||
            policy.provides_prior_scale()) {
          const arma::rowvec& active_scale = input.controls.estimate_maf_effect_s
            ? dynamic_maf_effect_scale
            : (input.controls.use_fixed_maf_effect_scale
                ? *input.controls.fixed_maf_effect_scale
                : policy.prior_scale());
          // Marker sweep begin.
          for (int sorted = 0; sorted < m; ++sorted) {
            const int marker =
              (*input.marker_order)[static_cast<std::size_t>(sorted)];
            sample_marker_scaled(
              marker, pi, marker_variance, active_scale, adjusted_residual,
              diagonal, residual, effects, state, input.data, gen
            );
          }
          // Marker sweep end.
        } else {
          // Marker sweep begin.
          for (int sorted = 0; sorted < m; ++sorted) {
            const int marker =
              (*input.marker_order)[static_cast<std::size_t>(sorted)];
            sample_marker_unscaled(
              marker, pi, marker_variance, adjusted_residual, diagonal,
              residual, effects, state, input.data, gen
            );
          }
          // Marker sweep end.
        }
        if (input.controls.update_ld_swap && input.controls.ld_swap_moves > 0 &&
            input.controls.ld_swap_probability > 0.0) {
          std::uniform_real_distribution<double> runif(0.0, 1.0);
          if (runif(gen) < input.controls.ld_swap_probability) {
            for (int move = 0; move < input.controls.ld_swap_moves; ++move) {
              output.ld_swap_attempted += 1.0;
              const arma::rowvec* active_scale = nullptr;
              if (input.controls.estimate_maf_effect_s) {
                active_scale = &dynamic_maf_effect_scale;
              } else if (input.controls.use_fixed_maf_effect_scale) {
                active_scale = input.controls.fixed_maf_effect_scale;
              } else if (policy.provides_prior_scale()) {
                active_scale = &policy.prior_scale();
              }
              if (attempt_ld_swap(
                    m, adjusted_residual, marker_variance,
                    (*input.data.yy)(static_cast<arma::uword>(trait)),
                    diagonal, wy, active_scale, residual, effects, state,
                    input.data, input.ld_friends, gen)) {
                output.ld_swap_accepted += 1.0;
              }
            }
          }
        }
        if (input.controls.update_marker_variance) {
          if (input.controls.estimate_maf_effect_s ||
              input.controls.use_fixed_maf_effect_scale ||
              policy.provides_prior_scale()) {
            const arma::rowvec& active_scale = input.controls.estimate_maf_effect_s
              ? dynamic_maf_effect_scale
              : (input.controls.use_fixed_maf_effect_scale
                  ? *input.controls.fixed_maf_effect_scale
                  : policy.prior_scale());
            sample_marker_variance_scaled(
              m, input.priors.marker_degrees_freedom, marker_variance,
              effects, state, active_scale,
              (*input.priors.marker_scale_prior)(trait, trait), gen
            );
          } else {
            sample_marker_variance_unscaled(
              m, input.priors.marker_degrees_freedom, marker_variance,
              effects, state,
              (*input.priors.marker_scale_prior)(trait, trait), gen
            );
          }
          if (!std::isfinite(marker_variance) || marker_variance <= 0.0) {
            throw std::runtime_error(
              "vb became invalid after sampleB. iter=" +
              std::to_string(iteration) + ", vb=" +
              std::to_string(marker_variance)
            );
          }
        }
        if (input.controls.estimate_maf_effect_s) {
          output.maf_effect_s_attempted += 1.0;
          if (update_maf_effect_s(
                maf_effect_s, effects, state, marker_variance,
                *input.controls.maf_effect_s_log_h,
                input.controls.maf_effect_s_prior_lower,
                input.controls.maf_effect_s_prior_upper,
                input.controls.maf_effect_s_proposal_sd, gen)) {
            output.maf_effect_s_accepted += 1.0;
          }
        }
        policy.after_vb_update(
          effects, state, marker_variance, gen, iteration
        );
        if (input.controls.update_residual_variance) {
          if (input.controls.rebuild_residual_before_update) {
            input.data.ld.rebuild(wy, effects, residual);
          }
          sample_residual_variance(
            m, input.priors.residual_degrees_freedom, residual_variance,
            effects, wy, residual,
            (*input.priors.residual_scale_prior)(trait, trait),
            (*input.data.yy)(static_cast<arma::uword>(trait)),
            (*input.data.sample_size)[static_cast<std::size_t>(trait)], gen
          );
        }
        if (input.controls.update_inclusion_probability) {
          sample_inclusion_probability(
            state, pi, input.priors.inclusion_prior_active,
            input.priors.inclusion_prior_null, gen
          );
          if (!std::isfinite(pi[0]) || !std::isfinite(pi[1]) ||
              pi[0] <= 0.0 || pi[1] <= 0.0) {
            throw std::runtime_error(
              "pi became invalid after samplePi. iter=" +
              std::to_string(iteration) + ", pi0=" +
              std::to_string(pi[0]) + ", pi1=" + std::to_string(pi[1])
            );
          }
        }
        genetic_variance = compute_genetic_variance(
          effects, wy, residual,
          (*input.data.sample_size)[static_cast<std::size_t>(trait)]
        );
        le_variance = compute_le(
          m, effects, diagonal,
          (*input.data.sample_size)[static_cast<std::size_t>(trait)]
        );
        ld_variance = genetic_variance - le_variance;
        if (!std::isfinite(genetic_variance)) {
          throw std::runtime_error(
            "vg became NaN/Inf after computeG. iter=" +
            std::to_string(iteration)
          );
        }
        if (!std::isfinite(le_variance)) {
          throw std::runtime_error(
            "vle became NaN/Inf after computeLE. iter=" +
            std::to_string(iteration)
          );
        }
        if (!std::isfinite(ld_variance)) {
          throw std::runtime_error(
            "vld became NaN/Inf after computeLE. iter=" +
            std::to_string(iteration)
          );
        }
        adjusted_residual = residual_variance +
          input.controls.residual_adjustment * genetic_variance;
        if (!std::isfinite(adjusted_residual) || adjusted_residual <= 0.0) {
          throw std::runtime_error(
            "adjusted residual variance vei became invalid. iter=" +
            std::to_string(iteration) + ", vei=" +
            std::to_string(adjusted_residual)
          );
        }
        const arma::uword iteration_u = static_cast<arma::uword>(iteration);
        output.marker_variance_trace(iteration_u) = marker_variance;
        output.residual_variance_trace(iteration_u) = residual_variance;
        output.genetic_variance_trace(iteration_u) = genetic_variance;
        output.inclusion_trace(iteration_u) = pi[1];
        output.le_variance_trace(iteration_u) = le_variance;
        output.ld_variance_trace(iteration_u) = ld_variance;
        output.maf_effect_s_trace(iteration_u) = maf_effect_s;
        policy.capture(iteration);
        if (iteration >= input.controls.nburn) {
          const arma::uword draw=static_cast<arma::uword>(iteration-input.controls.nburn);
          for (arma::uword s=0;s<selected_count;++s) {
            const arma::uword marker=static_cast<arma::uword>(input.controls.convergence_markers[static_cast<std::size_t>(s)]);
            if (input.controls.convergence_b) output.convergence_b(draw,s)=effects(marker);
            if (input.controls.convergence_d) output.convergence_d(draw,s)=state(marker);
          }
        }
        const int post_burn = iteration - input.controls.nburn + 1;
        const bool retained = input.controls.retention_contract_version == 1
          ? std::binary_search(
              input.controls.retained_transition_indices.begin(),
              input.controls.retained_transition_indices.end(),
              post_burn
            )
          : scalar_iteration_is_retained(
              iteration, input.controls.nburn, input.controls.nthin
            );
        if (retained) {
          policy.retain(iteration);
          output.retained_samples += 1.0;
          for (int marker = 0; marker < m; ++marker) {
            const arma::uword marker_u = static_cast<arma::uword>(marker);
            output.marker_mean(marker_u) += effects(marker_u);
            output.marker_pip(marker_u) += static_cast<double>(state(marker_u));
          }
        }
      }
      if (output.retained_samples <= 0.0) output.retained_samples = 1.0;
      output.marker_mean /= output.retained_samples;
      output.marker_pip /= output.retained_samples;
      if (!output.marker_mean.is_finite()) {
        throw std::runtime_error("posterior mean bm contains NaN/Inf.");
      }
      if (!output.marker_pip.is_finite()) {
        throw std::runtime_error("posterior mean dm contains NaN/Inf.");
      }
      output.final_effect = effects;
      output.final_residual = residual;
      output.final_state = state;
      output.final_marker_variance = marker_variance;
      output.final_residual_variance = residual_variance;
      output.final_genetic_variance = genetic_variance;
      output.final_le_variance = le_variance;
      output.final_ld_variance = ld_variance;
      output.final_inclusion_probability = pi[1];
      policy.finish();
      status.retained_samples = output.retained_samples;
#ifdef _OPENMP
      output.seconds = omp_get_wtime() - wall_start;
      status.elapsed_seconds = output.seconds;
#endif
    } catch (const std::exception& error) {
      status.failed = true;
      status.failure_message = error.what();
#ifdef _OPENMP
      output.seconds = omp_get_wtime() - wall_start;
      status.elapsed_seconds = output.seconds;
#endif
    } catch (...) {
      status.failed = true;
      status.failure_message = "unknown error";
#ifdef _OPENMP
      output.seconds = omp_get_wtime() - wall_start;
      status.elapsed_seconds = output.seconds;
#endif
    }
  }

  for (int task = 0; task < ntasks; ++task) {
    const ScalarChainExecutionStatus& status =
      statuses[static_cast<std::size_t>(task)];
    if (status.failed) {
      throw std::runtime_error(
        "stblr_cpg_omp_csr failed for trait " +
        std::to_string(status.task.trait_index) +
        ", chain " +
        std::to_string(status.task.chain_index) +
        ": " + status.failure_message
      );
    }
  }

  CsrBayesCResult result;
  result.marker_count = input.data.marker_count;
  result.trait_count = input.data.trait_count;
  result.trace_count = static_cast<std::size_t>(n_trace);
  result.chain_count = input.controls.nchains;
  result.marker_mean.zeros(nt, m);
  result.marker_pip.zeros(nt, m);
  result.marker_score = *input.data.wy;
  result.final_residual.zeros(nt, m);
  result.final_effect = *input.initial.effects;
  result.final_state.zeros(nt, m);
  result.marker_mean_sd.zeros(nt, m);
  result.marker_pip_sd.zeros(nt, m);
  result.marker_mean_min.zeros(nt, m);
  result.marker_pip_min.zeros(nt, m);
  result.marker_mean_max.zeros(nt, m);
  result.marker_pip_max.zeros(nt, m);
  result.marker_variance_trace.zeros(nt, n_trace);
  result.genetic_variance_trace.zeros(nt, n_trace);
  result.residual_variance_trace.zeros(nt, n_trace);
  result.inclusion_trace.zeros(nt, n_trace);
  result.le_variance_trace.zeros(nt, n_trace);
  result.ld_variance_trace.zeros(nt, n_trace);
  result.maf_effect_s_trace.zeros(nt, n_trace);
  result.final_marker_variance.zeros(nt);
  result.final_genetic_variance.zeros(nt);
  result.final_residual_variance.zeros(nt);
  result.final_le_variance.zeros(nt);
  result.final_ld_variance.zeros(nt);
  result.final_inclusion_probability.zeros(nt);
  result.retained_samples.zeros(nt);
  result.ld_swap_attempted.zeros(nt);
  result.ld_swap_accepted.zeros(nt);
  result.maf_effect_s_attempted.zeros(nt);
  result.maf_effect_s_accepted.zeros(nt);
  result.chains = std::move(chains);
  const double inverse_chains = 1.0 / static_cast<double>(input.controls.nchains);
  for (int trait = 0; trait < nt; ++trait) {
    const arma::uword trait_u = static_cast<arma::uword>(trait);
    result.marker_mean_min.row(trait_u).fill(
      std::numeric_limits<double>::infinity()
    );
    result.marker_pip_min.row(trait_u).fill(
      std::numeric_limits<double>::infinity()
    );
    result.marker_mean_max.row(trait_u).fill(
      -std::numeric_limits<double>::infinity()
    );
    result.marker_pip_max.row(trait_u).fill(
      -std::numeric_limits<double>::infinity()
    );
    for (int chain = 0; chain < input.controls.nchains; ++chain) {
      const int task = trait * input.controls.nchains + chain;
      const CsrBayesCChainResult& current =
        result.chains[static_cast<std::size_t>(task)];
      result.marker_mean.row(trait_u) += current.marker_mean;
      result.marker_pip.row(trait_u) += current.marker_pip;
      result.final_effect.row(trait_u) += current.final_effect;
      result.final_residual.row(trait_u) += current.final_residual;
      result.final_state.row(trait_u) += arma::conv_to<arma::rowvec>::from(
        current.final_state
      );
      result.marker_variance_trace.row(trait_u) +=
        current.marker_variance_trace;
      result.genetic_variance_trace.row(trait_u) +=
        current.genetic_variance_trace;
      result.residual_variance_trace.row(trait_u) +=
        current.residual_variance_trace;
      result.inclusion_trace.row(trait_u) += current.inclusion_trace;
      result.le_variance_trace.row(trait_u) += current.le_variance_trace;
      result.ld_variance_trace.row(trait_u) += current.ld_variance_trace;
      result.maf_effect_s_trace.row(trait_u) += current.maf_effect_s_trace;
      result.final_marker_variance(trait_u) += current.final_marker_variance;
      result.final_genetic_variance(trait_u) += current.final_genetic_variance;
      result.final_residual_variance(trait_u) += current.final_residual_variance;
      result.final_le_variance(trait_u) += current.final_le_variance;
      result.final_ld_variance(trait_u) += current.final_ld_variance;
      result.final_inclusion_probability(trait_u) +=
        current.final_inclusion_probability;
      result.retained_samples(trait_u) += current.retained_samples;
      result.ld_swap_attempted(trait_u) += current.ld_swap_attempted;
      result.ld_swap_accepted(trait_u) += current.ld_swap_accepted;
      result.maf_effect_s_attempted(trait_u) += current.maf_effect_s_attempted;
      result.maf_effect_s_accepted(trait_u) += current.maf_effect_s_accepted;
      for (int marker = 0; marker < m; ++marker) {
        const arma::uword marker_u = static_cast<arma::uword>(marker);
        result.marker_mean_min(trait_u, marker_u) = std::min(
          result.marker_mean_min(trait_u, marker_u),
          current.marker_mean(marker_u)
        );
        result.marker_pip_min(trait_u, marker_u) = std::min(
          result.marker_pip_min(trait_u, marker_u),
          current.marker_pip(marker_u)
        );
        result.marker_mean_max(trait_u, marker_u) = std::max(
          result.marker_mean_max(trait_u, marker_u),
          current.marker_mean(marker_u)
        );
        result.marker_pip_max(trait_u, marker_u) = std::max(
          result.marker_pip_max(trait_u, marker_u),
          current.marker_pip(marker_u)
        );
      }
    }
    result.marker_mean.row(trait_u) *= inverse_chains;
    result.marker_pip.row(trait_u) *= inverse_chains;
    result.final_effect.row(trait_u) *= inverse_chains;
    result.final_residual.row(trait_u) *= inverse_chains;
    result.final_state.row(trait_u) *= inverse_chains;
    result.marker_variance_trace.row(trait_u) *= inverse_chains;
    result.genetic_variance_trace.row(trait_u) *= inverse_chains;
    result.residual_variance_trace.row(trait_u) *= inverse_chains;
    result.inclusion_trace.row(trait_u) *= inverse_chains;
    result.le_variance_trace.row(trait_u) *= inverse_chains;
    result.ld_variance_trace.row(trait_u) *= inverse_chains;
    result.maf_effect_s_trace.row(trait_u) *= inverse_chains;
    result.final_marker_variance(trait_u) *= inverse_chains;
    result.final_genetic_variance(trait_u) *= inverse_chains;
    result.final_residual_variance(trait_u) *= inverse_chains;
    result.final_le_variance(trait_u) *= inverse_chains;
    result.final_ld_variance(trait_u) *= inverse_chains;
    result.final_inclusion_probability(trait_u) *= inverse_chains;
    result.retained_samples(trait_u) *= inverse_chains;
    if (input.controls.nchains > 1) {
      for (int chain = 0; chain < input.controls.nchains; ++chain) {
        const CsrBayesCChainResult& current = result.chains[
          static_cast<std::size_t>(trait * input.controls.nchains + chain)
        ];
        const arma::rowvec marker_mean_difference =
          current.marker_mean - result.marker_mean.row(trait_u);
        const arma::rowvec marker_pip_difference =
          current.marker_pip - result.marker_pip.row(trait_u);
        result.marker_mean_sd.row(trait_u) +=
          marker_mean_difference % marker_mean_difference;
        result.marker_pip_sd.row(trait_u) +=
          marker_pip_difference % marker_pip_difference;
      }
      result.marker_mean_sd.row(trait_u) = arma::sqrt(
        result.marker_mean_sd.row(trait_u) /
        static_cast<double>(input.controls.nchains - 1)
      );
      result.marker_pip_sd.row(trait_u) = arma::sqrt(
        result.marker_pip_sd.row(trait_u) /
        static_cast<double>(input.controls.nchains - 1)
      );
    }
  }
  return result;
}

#ifdef SBLR_CSR_BAYESC_CORE_IMPL_TRANSLATION_UNIT
inline CsrBayesCResult run_csr_bayesc(
  const CsrBayesCExecutionInput& input
) {
  CsrBayesCNoOpPolicyFactory policy_factory;
  return run_csr_bayesc_engine(input, policy_factory);
}
#endif

}  // namespace core
}  // namespace sblr

#endif
