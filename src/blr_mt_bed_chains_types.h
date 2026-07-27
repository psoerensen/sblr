#ifndef SBLR_BLR_MT_BED_CHAINS_TYPES_H
#define SBLR_BLR_MT_BED_CHAINS_TYPES_H

#include "blr_mt_bed_types.h"

#include <cstdint>
#include <string>
#include <vector>

namespace sblr {
namespace mt {

struct MtBedChainTask {
 int chain=-1;
 std::uint32_t seed=0;
};

struct MtBedChainExecutionResult {
 int chain=-1;
 std::uint32_t seed=0;
 bool failed=true;
 std::string error="chain did not execute";
 double seconds=0.0;
 MtBedCoreResult core;
};

struct MtBedChainSummary {
 int chain=-1;
 std::uint32_t seed=0;
 double seconds=0.0;
 std::vector<std::vector<double>> bm;
 std::vector<std::vector<double>> dm;
 std::vector<std::vector<double>> b;
 std::vector<std::vector<int>> state;
 std::vector<int> component_final;
 std::vector<std::vector<double>> component_probabilities;
 std::vector<std::vector<double>> vbs;
 std::vector<std::vector<double>> vgs;
 std::vector<std::vector<double>> ves;
 std::vector<std::vector<double>> vle;
 std::vector<std::vector<double>> vld;
 std::vector<std::vector<double>> covb;
 std::vector<std::vector<double>> covg;
 std::vector<std::vector<double>> cove;
 arma::mat B;
 arma::mat G;
 arma::mat E;
 std::vector<double> pi_final;
 std::vector<double> pi_mean;
 std::vector<std::vector<double>> pi_trace;
 arma::mat annotation_alpha_final;
 arma::mat annotation_alpha_mean;
 arma::vec annotation_sigma_final;
 arma::vec annotation_sigma_mean;
 std::vector<double> pattern_pi_final;
 std::vector<double> pattern_pi_mean;
 std::vector<std::vector<double>> pattern_pi_trace;
 arma::mat prior_component_probabilities;
 int annotation_updates_attempted=0;
 int annotation_updates_completed=0;
 MtBedCoreDiagnostics diagnostics;
};

struct MtBedChainsDiagnostics {
 std::size_t marker_cholesky_jitter_attempts=0;
 double marker_cholesky_max_increment=0.0;
 std::size_t full_e_updates=0;
 std::size_t diagonal_e_updates=0;
 std::vector<double> chain_marker_cholesky_jitter_attempts;
 std::vector<double> chain_marker_cholesky_max_increment;
 std::vector<double> chain_full_e_updates;
 std::vector<double> chain_diagonal_e_updates;
};

struct MtBedChainsAggregateResult {
 MtDefaultCoreResult pooled;
 MtBedChainsDiagnostics diagnostics;
 std::vector<std::vector<double>> bm_sd;
 std::vector<std::vector<double>> bm_min;
 std::vector<std::vector<double>> bm_max;
 std::vector<std::vector<double>> dm_sd;
 std::vector<std::vector<double>> dm_min;
 std::vector<std::vector<double>> dm_max;
 std::vector<MtBedChainSummary> chains;
};

}  // namespace mt
}  // namespace sblr

#endif
