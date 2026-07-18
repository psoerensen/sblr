#ifndef SBLR_BLR_SCHEDULED_EXECUTION_TYPES_H
#define SBLR_BLR_SCHEDULED_EXECUTION_TYPES_H

#include <cstddef>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace core {

// One instance belongs to one logical trait-chain task. It is constructed after
// seed resolution and destroyed when that chain execution ends; worker threads
// never own or share its engine or distributions.
struct ScheduledChainRng {
 std::mt19937 engine;
 std::normal_distribution<double> normal;
 std::uniform_real_distribution<double> uniform;

 explicit ScheduledChainRng(std::uint64_t seed)
  : engine(static_cast<std::mt19937::result_type>(seed)),
    normal(0.0,1.0), uniform(0.0,1.0) {}
};

struct ScheduledRngOwnership {
 std::string engine_owner="chain";
 std::string distribution_owner="chain";
 std::string lifetime="one_chain_execution";
 std::string worker_thread_owner="none";
 bool fit_persistent_distribution_state=false;
};

struct ScheduledSweepControl {
 // Zero is accepted and means a full sweep on every iteration. Positive
 // values request the established periodic schedule; negative values reject.
 int full_sweep_every=1;
 bool iteration_zero_is_full=true;
};

struct NullSkipControl {
 int base_interval=1;
 int maximum_interval=1;
 bool burnin_only=false;
 std::string growth_rule="probability_adaptive";
};

struct CandidateControl {
 double probability_threshold=0.0;
 int lifetime=0;
};

struct NeighborWakeupControl {
 bool enabled=false;
 double effect_difference_threshold=0.0;
 int maximum_neighbors=0;
 const void* friend_data=nullptr;
 std::size_t friend_marker_count=0;
 bool shared_read_only=true;
 bool storage_outlives_execution=true;
};

struct ScheduledExecutionControl {
 std::size_t marker_count=0;
 std::size_t trait_count=0;
 int iterations=0;
 int burnin=0;
 int thinning=1;
 int chains=1;
 int cores=1;
 int seed=0;
 std::vector<int> chain_seeds;
 bool keep_chains=false;
 ScheduledRngOwnership rng_ownership;
 ScheduledSweepControl sweep;
 NullSkipControl skip;
 CandidateControl candidate;
 NeighborWakeupControl neighbor;
};

// Current production vocabulary. All mutable vectors are chain-owned and have
// marker_count entries; the shared CSR and friend data are not represented here.
struct ScheduledMarkerState {
 std::vector<int> scheduled_at;
 std::vector<int> last_updated;
 std::vector<unsigned char> candidate;
 std::vector<unsigned char> in_candidate_list;
 std::vector<unsigned char> in_active_list;
 std::vector<int> last_interesting;
};

struct ScheduledDiagnostics {
 std::size_t attempted_updates=0;
 std::size_t completed_updates=0;
 std::size_t skipped_updates=0;
 std::size_t full_sweeps=0;
 std::size_t candidate_entries=0;
 std::size_t candidate_expiries=0;
 std::size_t neighbor_wakeups=0;
};

struct ScheduledResultVocabulary {
 std::size_t marker_count=0;
 std::size_t trait_count=0;
 std::size_t trace_length=0;
 bool has_marker_summaries=true;
 bool has_variance_traces=true;
 bool has_scheduler_diagnostics=false;
 bool has_optional_chain_payloads=false;
 ScheduledDiagnostics diagnostics;
};

inline void validate_scheduled_execution_control(
 const ScheduledExecutionControl& x
) {
 if (x.marker_count==0) throw std::invalid_argument("marker_count must be positive");
 if (x.trait_count==0) throw std::invalid_argument("trait_count must be positive");
 if (x.iterations<=0) throw std::invalid_argument("iterations must be positive");
 if (x.burnin<0) throw std::invalid_argument("burnin must be non-negative");
 if (x.thinning<=0) throw std::invalid_argument("thinning must be positive");
 if (x.chains<=0) throw std::invalid_argument("chains must be positive");
 if (x.cores<=0) throw std::invalid_argument("cores must be positive");
 if (!x.chain_seeds.empty() && static_cast<int>(x.chain_seeds.size())!=x.chains)
  throw std::invalid_argument("chain_seeds length must equal chains");
 if (x.rng_ownership.engine_owner!="chain" ||
     x.rng_ownership.distribution_owner!="chain" ||
     x.rng_ownership.lifetime!="one_chain_execution" ||
     x.rng_ownership.worker_thread_owner!="none" ||
     x.rng_ownership.fit_persistent_distribution_state)
  throw std::invalid_argument("scheduled RNG ownership must be chain-local and fit-bounded");
 if (x.sweep.full_sweep_every<0)
  throw std::invalid_argument("full_sweep_every must be non-negative");
 if (x.skip.base_interval<=0)
  throw std::invalid_argument("null_skip_base must be positive");
 if (x.skip.maximum_interval<0)
  throw std::invalid_argument("null_skip_max must be non-negative");
 if (x.skip.maximum_interval>0 && x.skip.base_interval>x.skip.maximum_interval)
  throw std::invalid_argument("null_skip_base must not exceed null_skip_max");
 if (x.skip.growth_rule!="probability_adaptive")
  throw std::invalid_argument("null skip growth_rule is unsupported");
 if (!(x.candidate.probability_threshold>=0.0 &&
       x.candidate.probability_threshold<=1.0))
  throw std::invalid_argument("candidate_threshold must be in [0,1]");
 if (x.candidate.lifetime<0)
  throw std::invalid_argument("candidate_lifetime must be non-negative");
 if (!(x.neighbor.effect_difference_threshold>=0.0))
  throw std::invalid_argument("wakeup_diff_threshold must be non-negative");
 if (x.neighbor.maximum_neighbors<0)
  throw std::invalid_argument("wakeup_max_neighbors must be non-negative");
 if (!x.neighbor.shared_read_only || !x.neighbor.storage_outlives_execution)
  throw std::invalid_argument("neighbor friend data must be borrowed immutable storage");
 if (x.neighbor.enabled && x.neighbor.friend_marker_count!=x.marker_count)
  throw std::invalid_argument("neighbor friend-list marker dimension mismatch");
}

inline void validate_scheduled_marker_state(
 const ScheduledMarkerState& state,
 std::size_t marker_count
) {
 if (state.scheduled_at.size()!=marker_count ||
     state.last_updated.size()!=marker_count ||
     state.candidate.size()!=marker_count ||
     state.in_candidate_list.size()!=marker_count ||
     state.in_active_list.size()!=marker_count ||
     state.last_interesting.size()!=marker_count)
  throw std::invalid_argument("scheduler state dimensions must equal marker_count");
}

} }

#endif
