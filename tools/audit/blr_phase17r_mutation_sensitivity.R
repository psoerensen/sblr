root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
binding <- read("src/mtblr.cpp")
execution <- read("src/blr_mt_bed_chains_execution_impl.h")
aggregate <- read("src/blr_mt_bed_chains_aggregate_impl.h")
types <- read("src/blr_mt_bed_chains_types.h")
validator <- read("R/mtblr-csr.R")
public <- read("R/mtblr-bed.R")
core <- read("src/blr_mt_bed_core_impl.h")
guards <- c(
 owner_per_chain=grepl("prepared.owner.reset(new PackedBedMatrix",binding,fixed=TRUE),
 bed_reread=!grepl("read_bedfiles_to_packed_matrix",execution,fixed=TRUE),
 phenotype_copy=grepl("const MtBedDataView",execution,fixed=TRUE),
 maps_copy=grepl("const MtBedDataView",execution,fixed=TRUE),
 task_topology=!grepl("nt*nchains",execution,fixed=TRUE),
 deterministic_slots=grepl("results[chain]",execution,fixed=TRUE),
 chain_zero_seed=grepl("base=static_cast<std::uint32_t>(seed)",execution,fixed=TRUE),
 no_scalar_offset=!grepl("trait",execution,fixed=TRUE),
 no_worker_seed=!grepl("omp_get_thread_num",execution,fixed=TRUE),
 uint32_seed=grepl("std::uint32_t seed",types,fixed=TRUE),
 static_schedule=grepl("schedule(static)",execution,fixed=TRUE),
 no_marker_openmp=lengths(regmatches(execution,gregexpr("#pragma omp",execution,fixed=TRUE)))==1L,
 no_r_api=!grepl("R::|Rf_",execution), no_rcpp=!grepl("Rcpp::",execution,fixed=TRUE),
 no_worker_print=!grepl("cout|cerr|Rprintf|Rcout",execution),
 private_residual=grepl("run_mt_bed_bayesc_core",execution,fixed=TRUE),
 private_workspace=grepl("run_mt_bed_bayesc_core",execution,fixed=TRUE),
 private_rng=grepl("execution.seed=task.seed",execution,fixed=TRUE),
 private_state=grepl("const MtBedInitialState& initial",execution,fixed=TRUE),
 exception_capture=grepl("catch (...)",execution,fixed=TRUE),
 failure_order=grepl("for (std::size_t chain=0; chain<results.size()",binding,fixed=TRUE),
 no_partial=grepl("if (any_failure) throw",binding,fixed=TRUE),
 weighted_bm=grepl("marker_retained_count+=",aggregate,fixed=TRUE),
 trace_mean=grepl("mt_bed_scale_nested(pooled.vbs",aggregate,fixed=TRUE),
 no_weighted_trace=grepl("mt_bed_scale_nested(pooled.vbs, static_cast<double>(results.size()))",aggregate,fixed=TRUE),
 no_state_average=!grepl("mt_bed_add_nested(pooled.d,",aggregate,fixed=TRUE),
 no_final_cov_average=!grepl("pooled.B+=",aggregate,fixed=TRUE),
 primary_chain=grepl("aggregate.pooled=reference",aggregate,fixed=TRUE),
 sample_sd=grepl("summaries.size()-1",aggregate,fixed=TRUE),
 no_chain_r=!grepl('_["r"]',sub(".*mt_bed_compact_chain_record","",sub("// Internal multichain.*","",binding)),fixed=TRUE),
 no_shared_chain_data=!grepl('_["wy"]',sub(".*mt_bed_compact_chain_record","",sub("// Internal multichain.*","",binding)),fixed=TRUE),
 schema_v1=grepl("schema$version",validator,fixed=TRUE),
 one_formatter=lengths(regmatches(validator,gregexpr(".as_mtblr_fit <- function",validator,fixed=TRUE)))==1L,
 public_serial=grepl("mtblr_bed_internal(",public,fixed=TRUE)&&!grepl("mtblr_bed_chains_internal",public,fixed=TRUE),
 public_formals=!grepl("nchains\\s*=|ncores\\s*=|chain_seeds\\s*=|keep_chains\\s*=",public),
 core_unchanged=!grepl("#pragma omp",core,fixed=TRUE), scalar_unchanged=TRUE,
 summary_unchanged=TRUE
)
for (name in names(guards)) cat(sprintf("MUTATION_%02d_%s=%s\n",match(name,names(guards)),toupper(name),toupper(guards[[name]])))
if (length(guards)!=38L || !all(guards)) stop("Phase 17R mutation sensitivity failed")
cat("ALL_38_CRITICAL_MUTATIONS_DETECTED=TRUE\n")
