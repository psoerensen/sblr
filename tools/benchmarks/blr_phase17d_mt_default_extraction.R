# Phase 17D reuses the canonical Phase 17C workload/reporting implementation.
source("tools/benchmarks/blr_phase17c_mt_default_corrected.R")

cat("Phase 17C pre-extraction matched baseline (seconds):\n")
cat("  updated: 0.21,0.00,0.01,0.00,0.02; mean 0.048; median 0.01; RSS 120.58 MiB\n")
cat("  fixed-B: 0.00,0.01,0.02,0.00,0.02; mean 0.010; median 0.01; RSS 120.63 MiB\n")
cat("  explicit sets: 0.01,0.04,0.00,0.00,0.01; mean 0.012; median 0.01; RSS 120.63 MiB\n")
cat("  moderate dense: 0.07,0.05,0.05,0.03,0.06; mean 0.052; median 0.05; RSS 121.52 MiB\n")
cat("Phase 17D is an extraction comparison within each matched workload only.\n")
