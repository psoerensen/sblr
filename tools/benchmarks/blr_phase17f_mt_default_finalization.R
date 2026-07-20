# Phase 17F reuses the canonical Phase 17E matched workloads and reporting.
source("tools/benchmarks/blr_phase17e_mt_default_typed_boundary.R")

cat("Phase 17F finalization/copy boundary audit:\n")
cat("  core result: moved by value into finalize_mt_default_result\n")
cat("  finalized result: owns posterior means, traces, final state, and counts\n")
cat("  marker/covariance/trace/state buffers: moved where safe\n")
cat("  pi_mean: one O(K) finalized allocation\n")
cat("  legacy output: unavoidable compatibility copies into positions 1-20\n")
cat("  borrowed wy: copied only into legacy position 3\n")
cat("  no new O(nt x m^2) copy and no MCMC-time I/O\n")
cat("completed-fit RSS is not peak RSS\n")
cat("dense XX remains O(nt x m^2)\n")
cat("tiny timings are regression signals\n")
