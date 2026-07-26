root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
contract <- read("docs/dev/blr_mt_bed_multichain_contract.md")
public <- read("R/mtblr-bed.R")
core <- read("src/blr_mt_bed_core_impl.h")
scalar <- paste(read("src/blr_bed_family_types.h"),
                read("src/st_cpg_omp_individual_scheduled_chains.cpp"),
                read("src/blr_bed_scheduled_bayesc_types.h"))
guards <- c(
  CURRENT_MT_BED_CHAIN_COUNT = !grepl("nchains", core, fixed = TRUE),
  CURRENT_MT_BED_OPENMP_LOOP_COUNT = !grepl("#pragma omp", core, fixed = TRUE),
  CURRENT_PUBLIC_NCHAINS_ARGUMENT = !grepl("nchains\\s*=", public),
  CURRENT_PUBLIC_NCORES_ARGUMENT = !grepl("ncores\\s*=", public),
  CURRENT_PUBLIC_CHAIN_SEEDS_ARGUMENT = !grepl("chain_seeds\\s*=", public),
  CURRENT_PUBLIC_KEEP_CHAINS_ARGUMENT = !grepl("keep_chains\\s*=", public),
  SCALAR_TASK_INDEX_HELPER_FOUND = grepl("make_bed_family_task_index", scalar, fixed = TRUE),
  SCALAR_SEED_HELPER_FOUND = grepl("resolve_bed_family_logical_chain_seed", scalar, fixed = TRUE),
  SCALAR_STATIC_SCHEDULING_FOUND = grepl("schedule(static)", scalar, fixed = TRUE),
  SCALAR_TYPED_CHAIN_RESULT_FOUND = grepl("ChainExecutionResult", scalar, fixed = TRUE),
  FUTURE_MT_TASK_TOPOLOGY_EXPLICIT = grepl("task_count = nchains", contract, fixed = TRUE),
  FUTURE_SEED_POLICY_EXPLICIT = grepl("9176*c", contract, fixed = TRUE),
  FUTURE_OWNERSHIP_EXPLICIT = grepl("Chain-private state", contract, fixed = TRUE),
  FUTURE_FAILURE_POLICY_EXPLICIT = grepl("no partial aggregation", contract, fixed = TRUE),
  FUTURE_AGGREGATION_EXPLICIT = grepl("sum bm, dm, covB, covG", contract, fixed = TRUE),
  FUTURE_FINAL_STATE_POLICY_EXPLICIT = grepl("primary_chain = 1", contract, fixed = TRUE),
  FUTURE_RAW_POLICY_EXPLICIT = grepl("mtblr_raw version 1", contract, fixed = TRUE),
  FUTURE_MEMORY_POLICY_EXPLICIT = grepl("private_state_bytes_per_chain", contract, fixed = TRUE)
)
cat("CURRENT_MT_BED_CHAIN_COUNT=1\n")
cat("CURRENT_MT_BED_OPENMP_LOOP_COUNT=0\n")
cat("CURRENT_PUBLIC_NCHAINS_ARGUMENT=FALSE\n")
cat("CURRENT_PUBLIC_NCORES_ARGUMENT=FALSE\n")
cat("CURRENT_PUBLIC_CHAIN_SEEDS_ARGUMENT=FALSE\n")
cat("CURRENT_PUBLIC_KEEP_CHAINS_ARGUMENT=FALSE\n")
for (name in names(guards)[7:length(guards)])
  cat(name, "=", toupper(guards[[name]]), "\n", sep = "")
if (!all(guards)) stop("Phase 17Q multichain contract audit failed")
