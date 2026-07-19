root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)
setwd(root)
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")
source(paste0("tests/testthat/fixtures/blr_phase17c_mt_default_corrected/",
  "blr-phase17c-mt-default-corrected-reference.R"))

`%||%` <- function(x, y) if (is.null(x)) y else x
rss_mib <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^2
}
run_case <- function(name, config, reps = 5L) {
  config$r_seed <- config$r_seed %||% 1701L
  set.seed(config$r_seed)
  invisible(do.call(sblr, phase17b_mt_public_args(config)))
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    set.seed(config$r_seed)
    times[i] <- system.time(do.call(sblr,
      phase17b_mt_public_args(config)))[["elapsed"]]
  }
  counts <- phase17c_mt_retained_counts(config)
  data.frame(case = name, samples = paste(config$n, collapse = "/"),
    markers = length(config$Xy[[1]]), traits = length(config$Xy),
    iterations = config$nit, burnin = config$nburn, thinning = config$nthin,
    marker_retained = counts$marker_retained_count,
    covb_retained = counts$covb_retained_count,
    covg_retained = counts$covg_retained_count,
    cove_retained = counts$cove_retained_count,
    pi_retained = counts$pi_retained_count,
    updateB = config$updateB, updateE = config$updateE,
    updatePi = config$updatePi,
    times = paste(format(times, digits = 4), collapse = ","),
    mean = mean(times), median = median(times), minimum = min(times),
    maximum = max(times), range = diff(range(times)),
    completed_fit_rss_mib = rss_mib(), stringsAsFactors = FALSE)
}

set.seed(17)
m <- 80L; nt <- 3L
ids <- paste0("M", seq_len(m)); traits <- paste0("T", seq_len(nt))
A <- matrix(rnorm(m * 12), m, 12)
ld <- cov2cor(crossprod(t(A)) + diag(m) * 0.25)
dimnames(ld) <- list(ids, ids)
moderate <- list(yy = setNames(rep(150, nt), traits),
  Xy = setNames(lapply(seq_len(nt), function(t) setNames(rnorm(m), ids)), traits),
  XX = setNames(replicate(nt, ld, simplify = FALSE), traits), n = rep(120L, nt),
  nit = 12L, nburn = 5L, nthin = 2L, verbose = FALSE, method = "bayesC",
  algorithm = "default", pi = 0.05, nub = 4, nue = 4, r_seed = 1710L,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE)

result <- do.call(rbind, list(
  run_case("updated-covariance-probability", phase17c_mt_config(1L)),
  run_case("corrected-fixed-B", phase17c_mt_config(2L)),
  run_case("explicit-sets", phase17c_mt_config(3L)),
  run_case("moderate-dense", moderate)))
print(result, row.names = FALSE)
cat("Phase 17B matched tiny baselines (seconds): updated mean 0.026, median 0.020;",
  "fixed mean 0.008, median 0.010.\n")
cat("R:", R.version.string, "\n")
cat("sblr:", as.character(utils::packageVersion("sblr")), "\n")
cat("Dense XX sizes (MiB): tiny=",
  format(as.numeric(object.size(phase17c_mt_config(1L)$XX)) / 1024^2, digits = 4),
  " moderate=", format(as.numeric(object.size(moderate$XX)) / 1024^2,
    digits = 4), "\n")
cat("sampled peak RSS: not measured (opt-in tooling only)\n")
cat("completed-fit RSS is not peak RSS\n")
cat("dense XX remains O(nt x m^2)\n")
cat("tiny timings are regression signals\n")
cat("all inputs are in-memory; no MCMC-time file I/O or page-cache dependency\n")
