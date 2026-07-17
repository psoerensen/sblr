# Phase 11D canonical packed-BED BayesC baseline.
# Completed-fit RSS is not peak memory; repeated reads may use the OS page cache.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11a-bed-reference.R"))

rss_mb <- function() if (requireNamespace("ps", quietly = TRUE))
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2 else NA_real_

run_fit <- function(cores, chains, seed, control) {
  x <- phase11a_fixture()
  do.call(sblr::stblr_bed, c(list(y = x$y, Glist = x$Glist, method = "bayesc",
    nit = 20L, nburn = 5L, nthin = 1L, seed = seed, ncores = cores,
    nchains = chains, updateB = FALSE, updateE = FALSE, rebuild_every = 5L,
    read_block_size = 64L, pi_init = .5, pi_prior_mean = .5,
    pi_prior_strength = 4), control))
}

bench <- function(workload, cores, chains, control, markers, samples, reps = 5L) {
  invisible(run_fit(cores, chains, 1100L, control))
  before <- rss_mb()
  times <- vapply(seq_len(reps), function(i)
    system.time(run_fit(cores, chains, 1100L + i, control))[["elapsed"]], numeric(1))
  data.frame(workload, markers, samples, chains, cores,
    scheduler = paste(names(control), unlist(control), sep = "=", collapse = ","),
    times = paste(times, collapse = ","), mean = mean(times), median = median(times),
    minimum = min(times), maximum = max(times), range = diff(range(times)),
    rss_before_mb = before, rss_after_mb = rss_mb(),
    memory_method = "ps whole-process RSS after completed fits; not peak RSS")
}

dense <- list(full_sweep_every = 1L, null_skip_base = 1L,
  null_skip_max = 1L, candidate_threshold = 0, candidate_lifetime = 0L)
aggressive <- list(full_sweep_every = 20L, null_skip_base = 50L,
  null_skip_max = 200L, candidate_threshold = 1e-4, candidate_lifetime = 10L)
conservative <- list(full_sweep_every = 5L, null_skip_base = 2L,
  null_skip_max = 10L, candidate_threshold = 1e-3, candidate_lifetime = 20L)

tiny <- rbind(
  bench("dense_1x1", 1L, 1L, dense, 2L, 6L),
  bench("dense_2x1", 1L, 2L, dense, 2L, 6L),
  bench("dense_2x2", 2L, 2L, dense, 2L, 6L),
  bench("aggressive_2x2", 2L, 2L, aggressive, 2L, 6L),
  bench("conservative_2x2", 2L, 2L, conservative, 2L, 6L))

original_fixture <- phase11a_fixture
set.seed(1113)
dosage <- matrix(sample(0:2, 2000L * 200L, replace = TRUE), 2000L, 200L)
bed <- tempfile(fileext = ".bed")
phase11a_write_bed(bed, dosage)
moderate <- list(Glist = list(n = 200L, ids = paste0("id", 1:200), bedfiles = bed,
  rsids = list(paste0("rs", 1:2000)), rsidsLD = list(paste0("rs", 1:2000)),
  chr = list(rep(1L, 2000)), pos = list(seq_len(2000)),
  af = list(rowMeans(dosage) / 2)), y = matrix(scale(seq_len(200)), ncol = 1,
  dimnames = list(NULL, "D1")))
phase11a_fixture <- function() moderate
moderate_results <- rbind(
  bench("moderate_dense_1x1", 1L, 1L, dense, 2000L, 200L),
  bench("moderate_aggressive_2x2", 2L, 2L, aggressive, 2000L, 200L))
phase11a_fixture <- original_fixture

cat("Phase 11D canonical packed-BED BayesC baseline\n")
cat("BED bytes:", file.info(bed)$size,
  "; first calls warm code/data and later calls may use the OS page cache.\n")
cat("R:", R.version.string, "; sblr:", as.character(packageVersion("sblr")), "\n")
cat("Compiler/toolchain: R package build toolchain reported by sessionInfo/build logs.\n")
print(rbind(tiny, moderate_results), row.names = FALSE)


