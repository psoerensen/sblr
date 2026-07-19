root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)
setwd(root)
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/fixtures/blr_phase17b_mt_default/blr-phase17b-mt-default-reference.R")

rss_mib <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^2
}
run_case <- function(name, config, reps = 5L) {
  config$r_seed <- config$r_seed %||% 1701L
  set.seed(config$r_seed); invisible(do.call(sblr, phase17b_mt_public_args(config)))
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    set.seed(config$r_seed)
    times[i] <- system.time(do.call(sblr, phase17b_mt_public_args(config)))[["elapsed"]]
  }
  data.frame(case = name, samples = paste(config$n, collapse = "/"),
    markers = length(config$Xy[[1]]), traits = length(config$Xy),
    iterations = config$nit, burnin = config$nburn, thinning = config$nthin,
    cores = "unsupported (single-chain default)", updateB = config$updateB,
    updateE = config$updateE, updatePi = config$updatePi,
    times = paste(format(times, digits = 4), collapse = ","),
    mean = mean(times), median = median(times), minimum = min(times),
    maximum = max(times), range = diff(range(times)), completed_fit_rss_mib = rss_mib(),
    stringsAsFactors = FALSE)
}
`%||%` <- function(x, y) if (is.null(x)) y else x

tiny <- phase17b_mt_config(1)
fixed <- phase17b_mt_config(2)
set.seed(17)
m <- 80L; nt <- 3L
ids <- paste0("M", seq_len(m)); traits <- paste0("T", seq_len(nt))
A <- matrix(rnorm(m * 12), m, 12); ld <- crossprod(t(A));
ld <- cov2cor(ld + diag(m) * 0.25); dimnames(ld) <- list(ids, ids)
moderate <- list(yy = setNames(rep(150, nt), traits),
  Xy = setNames(lapply(seq_len(nt), function(t) setNames(rnorm(m), ids)), traits),
  XX = setNames(replicate(nt, ld, simplify = FALSE), traits), n = rep(120L, nt),
  nit = 12L, nburn = 5L, nthin = 2L, verbose = FALSE, method = "bayesC",
  algorithm = "default", pi = 0.05, nub = 4, nue = 4, r_seed = 1710L,
  updateB = TRUE, updateE = TRUE, updatePi = TRUE)

result <- do.call(rbind, list(run_case("tiny-updated", tiny),
  run_case("tiny-fixed", fixed), run_case("moderate-dense", moderate)))
print(result, row.names = FALSE)
cat("R:", R.version.string, "\n")
cat("sblr:", as.character(utils::packageVersion("sblr")), "\n")
cat("Matrix object sizes (MiB): tiny=",
  format(as.numeric(object.size(tiny$XX)) / 1024^2, digits = 4),
  " moderate=", format(as.numeric(object.size(moderate$XX)) / 1024^2,
    digits = 4), "\n")
cat("completed-fit RSS is not peak RSS\n")
cat("dense matrices may dominate memory\n")
cat("tiny timings are regression signals\n")
cat("no cross-implementation speed ranking\n")
cat("all inputs are in-memory; no MCMC-time file I/O or page-cache dependency\n")
