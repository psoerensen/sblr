# Phase 10C3 final migrated scheduled CSR BayesC benchmark.
# The Phase 10B driver supplies identical checkpoint workloads, warm-up, five
# repetitions, timing distributions, and completed-fit RSS. The 2,000-marker
# representative workload is retained so Phase 10B--10C3 remain comparable.
cat("Phase 10C3 migrated scheduled CSR BayesC benchmark\n")
cat("Execution: typed context -> callable core -> typed result -> named binding converter\n")
cat("RNG ownership: one engine, normal distribution, and uniform distribution per trait-chain task\n")
cat("Memory: completed-fit whole-process RSS; not peak RSS\n")
source(file.path("tools", "benchmarks", "blr_phase10b_scheduled_csr.R"), chdir = FALSE)
