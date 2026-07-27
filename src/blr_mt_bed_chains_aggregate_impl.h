#ifndef SBLR_BLR_MT_BED_CHAINS_AGGREGATE_IMPL_H
#define SBLR_BLR_MT_BED_CHAINS_AGGREGATE_IMPL_H

#include "blr_mt_bed_chains_types.h"
#include "blr_mt_default_finalize_impl.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace sblr {
namespace mt {

inline MtDefaultFinalResult finalize_mt_bed_chains_result(
    MtDefaultCoreResult&& pooled) {
 return finalize_mt_default_result(std::move(pooled));
}

inline void mt_bed_require_same_shape(
    const std::vector<std::vector<double>>& left,
    const std::vector<std::vector<double>>& right,
    const char* label) {
 if (left.size()!=right.size()) {
  throw std::invalid_argument(std::string("inconsistent chain ")+label+
                              " dimensions");
 }
 for (std::size_t i=0; i<left.size(); ++i) {
  if (left[i].size()!=right[i].size()) {
   throw std::invalid_argument(std::string("inconsistent chain ")+label+
                               " dimensions");
  }
 }
}

inline void mt_bed_add_nested(
    std::vector<std::vector<double>>& target,
    const std::vector<std::vector<double>>& source) {
 for (std::size_t i=0; i<target.size(); ++i) {
  for (std::size_t j=0; j<target[i].size(); ++j) {
   target[i][j]+=source[i][j];
  }
 }
}

inline void mt_bed_scale_nested(
    std::vector<std::vector<double>>& target,
    double denominator) {
 for (std::vector<double>& row : target) {
  for (double& value : row) value/=denominator;
 }
}

inline std::vector<std::vector<double>> mt_bed_mean_nested(
    const std::vector<std::vector<double>>& values,
    double count) {
 std::vector<std::vector<double>> result=values;
 if (count>0.0) mt_bed_scale_nested(result, count);
 else {
  for (std::vector<double>& row : result) {
   std::fill(row.begin(), row.end(), 0.0);
  }
 }
 return result;
}

inline std::vector<double> mt_bed_mean_vector(
    const std::vector<double>& values,
    double count) {
 std::vector<double> result=values;
 if (count>0.0) {
  for (double& value : result) value/=count;
 } else {
  std::fill(result.begin(), result.end(), 0.0);
 }
 return result;
}

inline void validate_mt_bed_chain_consistency(
    const MtDefaultCoreResult& reference,
    const MtDefaultCoreResult& candidate) {
 if (candidate.nt!=reference.nt || candidate.m!=reference.m ||
     candidate.nmodels!=reference.nmodels ||
     candidate.order!=reference.order ||
     candidate.marker_retained_count<=0.0) {
  throw std::invalid_argument("inconsistent MT BED chain result");
 }
 mt_bed_require_same_shape(reference.bm, candidate.bm, "bm");
 mt_bed_require_same_shape(reference.dm, candidate.dm, "dm");
 mt_bed_require_same_shape(reference.r, candidate.r, "r");
 mt_bed_require_same_shape(reference.b, candidate.b, "b");
 mt_bed_require_same_shape(reference.vbs, candidate.vbs, "vbs");
 mt_bed_require_same_shape(reference.vgs, candidate.vgs, "vgs");
 mt_bed_require_same_shape(reference.ves, candidate.ves, "ves");
 mt_bed_require_same_shape(reference.vle, candidate.vle, "vle");
 mt_bed_require_same_shape(reference.vld, candidate.vld, "vld");
 mt_bed_require_same_shape(reference.cvbm, candidate.cvbm, "covB");
 mt_bed_require_same_shape(reference.cvgm, candidate.cvgm, "covG");
 mt_bed_require_same_shape(reference.cvem, candidate.cvem, "covE");
 if (candidate.d.size()!=reference.d.size() ||
     candidate.B.n_rows!=reference.B.n_rows ||
     candidate.B.n_cols!=reference.B.n_cols ||
     candidate.G.n_rows!=reference.G.n_rows ||
     candidate.G.n_cols!=reference.G.n_cols ||
     candidate.E.n_rows!=reference.E.n_rows ||
     candidate.E.n_cols!=reference.E.n_cols ||
     candidate.pi.size()!=reference.pi.size() ||
     candidate.pis.size()!=reference.pis.size()) {
  throw std::invalid_argument("inconsistent MT BED chain result dimensions");
 }
 for (std::size_t trait=0; trait<reference.d.size(); ++trait) {
  if (candidate.d[trait].size()!=reference.d[trait].size()) {
   throw std::invalid_argument("inconsistent MT BED chain state dimensions");
  }
 }
}

inline MtBedChainSummary summarize_mt_bed_chain(
    const MtBedChainExecutionResult& result) {
 const MtDefaultCoreResult& core=result.core.posterior;
 MtBedChainSummary summary;
 summary.chain=result.chain;
 summary.seed=result.seed;
 summary.seconds=result.seconds;
 summary.bm=mt_bed_mean_nested(core.bm, core.marker_retained_count);
 summary.dm=mt_bed_mean_nested(core.dm, core.marker_retained_count);
 summary.b=core.b;
 summary.state=core.d;
 summary.vbs=core.vbs;
 summary.vgs=core.vgs;
 summary.ves=core.ves;
 summary.vle=core.vle;
 summary.vld=core.vld;
 summary.covb=mt_bed_mean_nested(core.cvbm, core.covb_retained_count);
 summary.covg=mt_bed_mean_nested(core.cvgm, core.covg_retained_count);
 summary.cove=mt_bed_mean_nested(core.cvem, core.cove_retained_count);
 summary.B=core.B;
 summary.G=core.G;
 summary.E=core.E;
 summary.pi_final=core.pi;
 summary.pi_mean=mt_bed_mean_vector(core.pis, core.pi_retained_count);
 summary.diagnostics=result.core.diagnostics;
 return summary;
}

inline MtBedChainsAggregateResult aggregate_mt_bed_chains(
    const std::vector<MtBedChainExecutionResult>& results,
    bool keep_chains) {
 if (results.empty()) {
  throw std::invalid_argument("at least one MT BED chain result is required");
 }
 for (std::size_t chain=0; chain<results.size(); ++chain) {
  if (results[chain].failed || results[chain].chain!=static_cast<int>(chain)) {
   throw std::invalid_argument("cannot aggregate failed or unordered MT BED chains");
  }
 }
 const MtDefaultCoreResult& reference=results[0].core.posterior;
 if (reference.marker_retained_count<=0.0) {
  throw std::invalid_argument("every MT BED chain must retain markers");
 }
 for (std::size_t chain=1; chain<results.size(); ++chain) {
  validate_mt_bed_chain_consistency(reference, results[chain].core.posterior);
 }

 MtBedChainsAggregateResult aggregate;
 aggregate.pooled=reference;
 MtDefaultCoreResult& pooled=aggregate.pooled;
 for (std::vector<double>& row : pooled.bm) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.dm) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.cvbm) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.cvgm) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.cvem) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.vbs) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.vgs) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.ves) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.vle) std::fill(row.begin(), row.end(), 0.0);
 for (std::vector<double>& row : pooled.vld) std::fill(row.begin(), row.end(), 0.0);
 std::fill(pooled.pis.begin(), pooled.pis.end(), 0.0);
 pooled.marker_retained_count=0.0;
 pooled.covb_retained_count=0.0;
 pooled.covg_retained_count=0.0;
 pooled.cove_retained_count=0.0;
 pooled.pi_retained_count=0.0;

 std::vector<MtBedChainSummary> summaries;
 summaries.reserve(results.size());
 for (const MtBedChainExecutionResult& result : results) {
  const MtDefaultCoreResult& core=result.core.posterior;
  mt_bed_add_nested(pooled.bm, core.bm);
  mt_bed_add_nested(pooled.dm, core.dm);
  mt_bed_add_nested(pooled.cvbm, core.cvbm);
  mt_bed_add_nested(pooled.cvgm, core.cvgm);
  mt_bed_add_nested(pooled.cvem, core.cvem);
  mt_bed_add_nested(pooled.vbs, core.vbs);
  mt_bed_add_nested(pooled.vgs, core.vgs);
  mt_bed_add_nested(pooled.ves, core.ves);
  mt_bed_add_nested(pooled.vle, core.vle);
  mt_bed_add_nested(pooled.vld, core.vld);
  for (std::size_t model=0; model<pooled.pis.size(); ++model) {
   pooled.pis[model]+=core.pis[model];
  }
  pooled.marker_retained_count+=core.marker_retained_count;
  pooled.covb_retained_count+=core.covb_retained_count;
  pooled.covg_retained_count+=core.covg_retained_count;
  pooled.cove_retained_count+=core.cove_retained_count;
  pooled.pi_retained_count+=core.pi_retained_count;
  summaries.push_back(summarize_mt_bed_chain(result));

  aggregate.diagnostics.marker_cholesky_jitter_attempts+=
   result.core.diagnostics.marker_cholesky_jitter_attempts;
  aggregate.diagnostics.marker_cholesky_max_increment=std::max(
   aggregate.diagnostics.marker_cholesky_max_increment,
   result.core.diagnostics.marker_cholesky_max_increment);
  aggregate.diagnostics.full_e_updates+=result.core.diagnostics.full_e_updates;
  aggregate.diagnostics.diagonal_e_updates+=
   result.core.diagnostics.diagonal_e_updates;
  aggregate.diagnostics.chain_marker_cholesky_jitter_attempts.push_back(
   static_cast<double>(result.core.diagnostics.marker_cholesky_jitter_attempts));
  aggregate.diagnostics.chain_marker_cholesky_max_increment.push_back(
   result.core.diagnostics.marker_cholesky_max_increment);
  aggregate.diagnostics.chain_full_e_updates.push_back(
   static_cast<double>(result.core.diagnostics.full_e_updates));
  aggregate.diagnostics.chain_diagonal_e_updates.push_back(
   static_cast<double>(result.core.diagnostics.diagonal_e_updates));
 }
 mt_bed_scale_nested(pooled.vbs, static_cast<double>(results.size()));
 mt_bed_scale_nested(pooled.vgs, static_cast<double>(results.size()));
 mt_bed_scale_nested(pooled.ves, static_cast<double>(results.size()));
 mt_bed_scale_nested(pooled.vle, static_cast<double>(results.size()));
 mt_bed_scale_nested(pooled.vld, static_cast<double>(results.size()));

 const std::size_t nt=reference.bm.size();
 const std::size_t m=nt==0 ? 0 : reference.bm[0].size();
 auto zeros=[&]() { return std::vector<std::vector<double>>(
  nt, std::vector<double>(m, 0.0)); };
 aggregate.bm_sd=zeros(); aggregate.bm_min=zeros(); aggregate.bm_max=zeros();
 aggregate.dm_sd=zeros(); aggregate.dm_min=zeros(); aggregate.dm_max=zeros();
 for (std::size_t trait=0; trait<nt; ++trait) {
  for (std::size_t marker=0; marker<m; ++marker) {
   double bm_sum=0.0, dm_sum=0.0;
   double bm_min=std::numeric_limits<double>::infinity();
   double bm_max=-std::numeric_limits<double>::infinity();
   double dm_min=std::numeric_limits<double>::infinity();
   double dm_max=-std::numeric_limits<double>::infinity();
   for (const MtBedChainSummary& summary : summaries) {
    const double bm=summary.bm[trait][marker];
    const double dm=summary.dm[trait][marker];
    bm_sum+=bm; dm_sum+=dm;
    bm_min=std::min(bm_min, bm); bm_max=std::max(bm_max, bm);
    dm_min=std::min(dm_min, dm); dm_max=std::max(dm_max, dm);
   }
   const double bm_mean=bm_sum/static_cast<double>(summaries.size());
   const double dm_mean=dm_sum/static_cast<double>(summaries.size());
   double bm_ss=0.0, dm_ss=0.0;
   for (const MtBedChainSummary& summary : summaries) {
    bm_ss+=(summary.bm[trait][marker]-bm_mean)*
     (summary.bm[trait][marker]-bm_mean);
    dm_ss+=(summary.dm[trait][marker]-dm_mean)*
     (summary.dm[trait][marker]-dm_mean);
   }
   aggregate.bm_sd[trait][marker]=summaries.size()==1 ? 0.0 :
    std::sqrt(bm_ss/static_cast<double>(summaries.size()-1));
   aggregate.dm_sd[trait][marker]=summaries.size()==1 ? 0.0 :
    std::sqrt(dm_ss/static_cast<double>(summaries.size()-1));
   aggregate.bm_min[trait][marker]=bm_min;
   aggregate.bm_max[trait][marker]=bm_max;
   aggregate.dm_min[trait][marker]=dm_min;
   aggregate.dm_max[trait][marker]=dm_max;
  }
 }
 if (keep_chains) aggregate.chains=std::move(summaries);
 return aggregate;
}

}  // namespace mt
}  // namespace sblr

#endif
