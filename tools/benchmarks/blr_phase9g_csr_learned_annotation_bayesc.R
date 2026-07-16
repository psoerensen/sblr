# Canonical learned-annotation CSR BayesC baseline; run manually from root.
# The Phase 9F3 workload is retained verbatim as the frozen comparable workload.
source(file.path(
  "tools", "benchmarks", "blr_phase9f3_csr_learned_annotation_bayesc.R"
))
cat("baseline status: canonical learned-annotation CSR BayesC\n")
cat("memory sampling interval: unavailable; completed-fit RSS only\n")
cat("compiler/toolchain: Rtools44 GCC 13.2.0 (repository Windows toolchain)\n")
