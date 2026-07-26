root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
contract <- read("docs/dev/blr_mt_bed_multichain_contract.md")
test <- read("tests/testthat/test-blr-framework-phase17q.R")
helper <- read("tests/testthat/helper-mtblr-bed-multichain-contract.R")
production <- paste(read("R/mtblr-bed.R"), read("src/blr_mt_bed_core_impl.h"),
                    read("src/mtblr.cpp"), read("NAMESPACE"))
tokens <- c(
  "task_count = nchains", "Completion order", "9176*c", "worker identity",
  "supplied order", "packed owner is counted once", "one stationary packed owner",
  "centered Y", "sample residual R", "decoded-marker workspace", "one RNG",
  "B, E, G, pi", "no R API", "Workers emit no progress", "schedule(static)",
  "ascending chain order", "no partial aggregation", "retained counts",
  "iterationwise arithmetic", "Binary states are never averaged",
  "primary_chain = 1", "omit packed bytes", "Chain-stability summaries",
  "mtblr_raw version 1", "oversubscription", "Phase 17Q changes no formals"
)
detected <- vapply(tokens, function(x) grepl(x, contract, fixed = TRUE), logical(1))
guards <- c(
  detected,
  NO_MULTICHAIN_ROUTE = !grepl("mtblr_bed_chains_internal", production, fixed = TRUE),
  NO_MT_OPENMP = !grepl("#pragma omp", read("src/blr_mt_bed_core_impl.h"), fixed = TRUE),
  ORACLES_COVER_FAILURE = grepl("chain 1: first", test, fixed = TRUE),
  ORACLES_COVER_MEMORY = grepl("shared_bytes", helper, fixed = TRUE),
  NO_PUBLIC_CONTROLS = !grepl("nchains\\s*=|ncores\\s*=|chain_seeds\\s*=|keep_chains\\s*=",
                              read("R/mtblr-bed.R"))
)
cat("PHASE17Q_MUTATIONS_DETECTED=", sum(guards), "/", length(guards), "\n", sep = "")
if (!all(guards)) stop("Undetected Phase 17Q mutation(s): ",
                       paste(names(guards)[!guards], collapse = ", "))
