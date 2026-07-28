#ifndef SBLR_BLR_MT_BED_CHAINS_EXECUTION_IMPL_H
#define SBLR_BLR_MT_BED_CHAINS_EXECUTION_IMPL_H

#include "blr_mt_bed_chains_types.h"
#include "blr_mt_bed_core_impl.h"

#include <chrono>
#include <cstdint>
#include <exception>
#include <vector>

namespace sblr {
namespace mt {

inline std::uint32_t resolve_mt_bed_chain_seed(
    int seed,
    int zero_based_chain) {
 if (zero_based_chain<0) {
  throw std::invalid_argument("chain index must be nonnegative");
 }
 const std::uint64_t base=static_cast<std::uint32_t>(seed);
 const std::uint64_t offset=9176ULL*
  static_cast<std::uint64_t>(zero_based_chain);
 return static_cast<std::uint32_t>((base+offset)&0xffffffffULL);
}

inline std::uint32_t resolve_mt_bed_explicit_chain_seed(int seed) {
 return static_cast<std::uint32_t>(seed);
}

template <class PackedGenotype>
MtBedChainExecutionResult run_mt_bed_chain_task(
    const MtBedChainTask& task,
    const MtBedDataView<PackedGenotype>& data,
    const MtBedInitialState& initial,
    const std::vector<std::vector<int>>& sets,
    const std::vector<std::vector<double>>& ssb_prior,
    const std::vector<std::vector<double>>& sse_prior,
    const std::vector<std::vector<int>>& models,
    double nub,
    double nue,
    const MtBedExecutionSpec& base_execution,
    const MtJointStateSpec* joint=nullptr,
    const std::vector<double>* marker_scale=nullptr,
    const std::vector<double>* pi_prior=nullptr,
    const MtBayesRCSpec* bayesrc=nullptr,
    const MtExtendedTraceSpec* convergence=nullptr) {
 MtBedChainExecutionResult result;
 result.chain=task.chain;
 result.seed=task.seed;
 const auto start=std::chrono::steady_clock::now();
 try {
  MtBedExecutionSpec execution=base_execution;
  execution.seed=task.seed;
  result.core=run_mt_bed_bayesc_core(
   data, initial, sets, ssb_prior, sse_prior, models, nub, nue,
   execution,joint,marker_scale,pi_prior,bayesrc,convergence);
  result.failed=false;
  result.error.clear();
 } catch (const std::exception& error) {
  result.failed=true;
  result.error=error.what();
 } catch (...) {
  result.failed=true;
  result.error="unknown C++ exception";
 }
 result.seconds=std::chrono::duration<double>(
  std::chrono::steady_clock::now()-start).count();
 return result;
}

template <class PackedGenotype>
std::vector<MtBedChainExecutionResult> dispatch_mt_bed_chain_tasks(
    const std::vector<MtBedChainTask>& tasks,
    int worker_count,
    const MtBedDataView<PackedGenotype>& data,
    const MtBedInitialState& initial,
    const std::vector<std::vector<int>>& sets,
    const std::vector<std::vector<double>>& ssb_prior,
    const std::vector<std::vector<double>>& sse_prior,
    const std::vector<std::vector<int>>& models,
    double nub,
    double nue,
    const MtBedExecutionSpec& base_execution,
    const MtJointStateSpec* joint=nullptr,
    const std::vector<double>* marker_scale=nullptr,
    const std::vector<double>* pi_prior=nullptr,
    const MtBayesRCSpec* bayesrc=nullptr,
    const MtExtendedTraceSpec* convergence=nullptr) {
 std::vector<MtBedChainExecutionResult> results(tasks.size());
 for (std::size_t chain=0; chain<tasks.size(); ++chain) {
  results[chain].chain=tasks[chain].chain;
  results[chain].seed=tasks[chain].seed;
 }
 if (worker_count<=1) {
  for (std::size_t chain=0; chain<tasks.size(); ++chain) {
   results[chain]=run_mt_bed_chain_task(
    tasks[chain], data, initial, sets, ssb_prior, sse_prior,
    models, nub, nue, base_execution,joint,marker_scale,pi_prior,bayesrc,
    convergence);
  }
  return results;
 }
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
 for (int chain=0; chain<static_cast<int>(tasks.size()); ++chain) {
  results[static_cast<std::size_t>(chain)]=run_mt_bed_chain_task(
   tasks[static_cast<std::size_t>(chain)], data, initial, sets,
   ssb_prior, sse_prior, models, nub, nue, base_execution,
   joint,marker_scale,pi_prior,bayesrc,convergence);
 }
#else
 for (std::size_t chain=0; chain<tasks.size(); ++chain) {
  results[chain]=run_mt_bed_chain_task(
   tasks[chain], data, initial, sets, ssb_prior, sse_prior,
   models, nub, nue, base_execution,joint,marker_scale,pi_prior,bayesrc,
   convergence);
 }
#endif
 return results;
}

}  // namespace mt
}  // namespace sblr

#endif
