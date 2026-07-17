# Canonical scheduled ordinary-CSR BayesC baseline.
# The stable Phase 10C3/10B driver provides the 2,000-marker representative
# workloads, warm-up, five repetitions, all supported scheduler modes, timing
# distributions, and completed-fit whole-process RSS. It is intentionally
# reused to keep the canonical checkpoint directly comparable.
cat("Phase 10D canonical scheduled ordinary-CSR BayesC baseline\n")
cat("Architecture: typed context -> canonical core -> typed result -> named converter\n")
cat("RNG ownership: fit-bounded chain engine and distributions; no worker ownership\n")
cat("Memory method: completed-fit whole-process RSS, not peak RSS\n")
source(file.path("tools", "benchmarks", "blr_phase10b_scheduled_csr.R"), chdir = FALSE)
