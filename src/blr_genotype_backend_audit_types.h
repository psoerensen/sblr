#ifndef SBLR_BLR_GENOTYPE_BACKEND_AUDIT_TYPES_H
#define SBLR_BLR_GENOTYPE_BACKEND_AUDIT_TYPES_H

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace sblr { namespace audit {

enum class GenotypeAccessKind { packed_bed_marker };
enum class SchedulerKind { sparse_probability, adaptive_due, full_sweep_only };
enum class RngOwner { marker_call, chain, worker_thread_persistent };

struct PackedBedSourceView {
 std::vector<std::string> bed_paths;
 std::size_t sample_count=0;
 std::size_t marker_count=0;
 bool snp_major=true;
 bool borrowed_immutable_paths=true;
 bool fit_local_decoded_storage=true;
 std::size_t read_block_size=0;
};

struct GenotypeScalingControl {
 bool center=true;
 bool scale=true;
 bool mean_impute_missing=true;
};

struct AdaptiveBedSchedulerControl {
 SchedulerKind kind=SchedulerKind::adaptive_due;
 int full_sweep_every=1;
 int null_skip_base=1;
 int null_skip_max=1;
 double candidate_threshold=0.0;
 int candidate_lifetime=0;
 bool skip_nulls_burnin_only=false;
 bool has_neighbor_wakeup=false;
};

struct BedChainSpec {
 int chains=1;
 int cores=1;
 int base_seed=0;
 std::vector<int> explicit_chain_seeds;
 RngOwner engine_owner=RngOwner::chain;
 RngOwner distribution_owner=RngOwner::chain;
 bool fit_persistent_distribution_state=false;
 bool worker_owns_distribution_state=false;
 bool distribution_lifetime_is_one_chain=true;
};

struct BedExecutionAuditContract {
 std::string model;
 PackedBedSourceView genotype;
 GenotypeScalingControl scaling;
 AdaptiveBedSchedulerControl scheduler;
 BedChainSpec chains;
 int traits=1;
 int iterations=1;
 int burnin=0;
 int thinning=1;
 std::size_t annotation_rows=0;
 std::size_t annotation_columns=0;
};

struct BedExecutionResultVocabulary {
 bool marker_posterior=true;
 bool mixture_assignments=false;
 bool mixture_probabilities=false;
 bool annotation_coefficients=false;
 bool chain_summaries=true;
 bool optional_chain_payloads=false;
 bool scheduler_diagnostics=false;
 bool timing=true;
 bool failure_state=true;
};

inline void validate_bed_execution_audit_contract(
 const BedExecutionAuditContract& x
) {
 if (x.model!="bayesc" && x.model!="bayesr" && x.model!="bayesrc")
  throw std::invalid_argument("bed audit model is unsupported");
 if (x.genotype.bed_paths.empty()) throw std::invalid_argument("bed_paths must not be empty");
 if (x.genotype.sample_count==0) throw std::invalid_argument("sample_count must be positive");
 if (x.genotype.marker_count==0) throw std::invalid_argument("marker_count must be positive");
 if (!x.genotype.snp_major) throw std::invalid_argument("packed BED must be SNP-major");
 if (!x.genotype.borrowed_immutable_paths)
  throw std::invalid_argument("BED paths must be borrowed immutable metadata");
 if (x.traits<=0 || x.iterations<=0 || x.burnin<0 || x.thinning<=0)
  throw std::invalid_argument("invalid BED execution dimensions or MCMC controls");
 if (x.chains.chains<=0 || x.chains.cores<=0)
  throw std::invalid_argument("BED chains and cores must be positive");
 if (x.chains.engine_owner!=RngOwner::chain ||
     x.chains.distribution_owner!=RngOwner::chain ||
     x.chains.fit_persistent_distribution_state ||
     x.chains.worker_owns_distribution_state ||
     !x.chains.distribution_lifetime_is_one_chain)
  throw std::invalid_argument(
   "scheduled BED RNG engines and distributions must be one-chain owned");
 if (!x.chains.explicit_chain_seeds.empty() &&
     static_cast<int>(x.chains.explicit_chain_seeds.size())!=x.chains.chains)
  throw std::invalid_argument("BED explicit chain seeds must match chains");
 if (x.scheduler.kind==SchedulerKind::adaptive_due) {
  if (x.scheduler.full_sweep_every<=0)
   throw std::invalid_argument("full_sweep_every must be positive");
  if (x.scheduler.null_skip_base<=0 ||
      x.scheduler.null_skip_max<x.scheduler.null_skip_base)
   throw std::invalid_argument("invalid adaptive BED skip controls");
  if (x.scheduler.candidate_threshold<0.0 || x.scheduler.candidate_threshold>1.0 ||
      x.scheduler.candidate_lifetime<0)
   throw std::invalid_argument("invalid adaptive BED candidate controls");
 }
 if (x.model=="bayesrc" &&
     (x.annotation_rows!=x.genotype.marker_count || x.annotation_columns==0))
  throw std::invalid_argument(
   "BayesRC annotations must have one row per marker and positive columns");
}

} }

#endif
