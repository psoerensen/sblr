# Phase 10B post-correction baseline. Workloads and measurement method are
# intentionally identical to Phase 10A for direct comparison.
cat("Phase 10B chain-owned RNG scheduled CSR benchmark\n")
cat("RNG ownership: one engine, normal distribution, and uniform distribution per trait-chain task\n")
cat("Scheduler event counters remain unavailable in the unchanged public result schema\n")
source(file.path("tools","benchmarks","blr_phase10a_scheduled_csr.R"),chdir=FALSE)
