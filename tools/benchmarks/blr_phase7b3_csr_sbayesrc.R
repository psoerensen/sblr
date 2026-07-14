pkgload::load_all(".", compile = FALSE)
source(file.path("tests", "testthat", "fixtures", "blr-phase7a-sbayesrc-reference.R"))

rss_mib <- function() {
 if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
 as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2
}

make_workload <- function(m = 2000L, nit = 120L, nburn = 30L) {
 ids <- paste0("m", seq_len(m)); prefix <- tempfile("phase7b3_sbayesrc_")
 sblr:::.stblr_write_uint64_file(paste0(prefix, ".row_ptr.u64.bin"), rep(0, m + 1L))
 file.create(paste0(prefix, ".col_idx.u32.0based.bin")); file.create(paste0(prefix, ".values.f32.bin"))
 writeLines(c("format=sparse_ld_csr", "storage=streamed_upper_triangle", "n_bed=NA", "n_used=NA", "n_samples_used=NA",
              paste0("n_variants=", m), "nnz=0", "triangle=upper", "diagonal=implicit_1",
              paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"), paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
              paste0("values_file=", prefix, ".values.f32.bin"), "row_ptr_type=uint64", "col_idx_type=uint32",
              "values_type=float32", "index_base=0", "value=r"), paste0(prefix, ".meta.txt"))
 stats <- list(wy = list(T1 = stats::setNames(10 * sin(seq_len(m) / 23), ids)),
               ww = list(T1 = stats::setNames(rep(250, m), ids)), yy = c(T1 = 250), n = 250L, m = m)
 A <- cbind(intercept = 1, coding = as.numeric(seq_len(m) %% 5L == 0L),
            qtl = as.numeric(seq_len(m) %% 13L == 0L), continuous = as.numeric(scale(seq_len(m))))
 rownames(A) <- ids
 list(prefix = prefix, stats = stats, A = A, nit = nit, nburn = nburn)
}

run_fit <- function(w, update_alpha, chains, cores, keep = FALSE) {
 sblr::stblr_csr_sbayesrc_generic(w$stats, w$prefix, w$A, gamma = c(0, .01, .1, 1),
  add_intercept = FALSE, standardize_annotations = FALSE, updateAlpha = update_alpha,
  alpha_update_every = 5L, updateE = FALSE, nit = w$nit, nburn = w$nburn,
  nchains = chains, ncores = cores, keep_chains = keep, seed = 911L)
}

bench <- function(label, w, update_alpha, chains, cores, keep = FALSE, reps = 5L) {
 invisible(run_fit(w, update_alpha, chains, cores, keep))
 elapsed <- rss <- numeric(reps)
 for (i in seq_len(reps)) {
  elapsed[i] <- system.time(invisible(run_fit(w, update_alpha, chains, cores, keep)))[["elapsed"]]
  rss[i] <- rss_mib()
 }
 data.frame(configuration = label, elapsed = paste(elapsed, collapse = ","), mean = mean(elapsed),
  median = median(elapsed), minimum = min(elapsed), maximum = max(elapsed), iqr = IQR(elapsed),
  rss_mib = if (all(is.na(rss))) NA_real_ else max(rss, na.rm = TRUE), marker_count = w$stats$m,
  trait_count = 1L, annotation_count = ncol(w$A), component_count = 4L,
  alpha_parameter_count = ncol(w$A) * 3L, iterations = w$nit, retained = w$nit,
  chains = chains, cores = cores, updateAlpha = update_alpha, alpha_update_every = 5L,
  keep_chains = keep, output_mode = if (keep) "ordinary_with_chains" else "ordinary")
}

tiny <- phase7a_sbayesrc_configs$fixed_one_chain
invisible(phase7a_sbayesrc_run(tiny, FALSE))
cat("tiny exact fixture elapsed:", system.time(invisible(phase7a_sbayesrc_run(tiny, FALSE)))[["elapsed"]], "seconds\n")
w <- make_workload()
configs <- list(
 list("fixed_1chain_1core", FALSE, 1L, 1L, FALSE),
 list("learned_1chain_1core", TRUE, 1L, 1L, FALSE),
 list("learned_2chains_1core", TRUE, 2L, 1L, FALSE),
 list("learned_2chains_2cores", TRUE, 2L, 2L, FALSE),
 list("learned_2chains_2cores_keep", TRUE, 2L, 2L, TRUE))
results <- do.call(rbind, lapply(configs, function(z) bench(z[[1]], w, z[[2]], z[[3]], z[[4]], z[[5]])))
print(results, row.names = FALSE)
cat("warm-up: one untimed call per configuration; five timed repetitions\n")
cat("minimal output: not supported by the current public SBayesRC route; ordinary and retained-chain modes measured\n")
cat("memory: whole-process RSS sampled after each fit; sampling interval is one completed fit, not sampler-only peak\n")
cat("environment:", R.version.string, "; sblr", as.character(utils::packageVersion("sblr")),
    "; compiler", R.version$compiler, "\n")
