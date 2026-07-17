# Phase 10C2 typed-boundary scheduled CSR benchmark.
# The Phase 10B workload driver supplies warm-up, five primary repetitions,
# dense and skipping configurations, chain/core combinations, scheduler
# controls, elapsed summaries, and completed-fit RSS. This wrapper deliberately
# uses identical workloads for direct checkpoint comparison.
cat("Phase 10C2 typed scheduled CSR BayesC benchmark\n")
cat("Execution: CsrScheduledBayesCExecutionContext -> run_csr_scheduled_bayesc -> typed result\n")
cat("RNG ownership: one engine, normal distribution, and uniform distribution per trait-chain task\n")
cat("Memory method and scheduler-counter limitations are reported by the Phase 10A driver\n")
source(file.path("tools", "benchmarks", "blr_phase10b_scheduled_csr.R"), chdir = FALSE)

