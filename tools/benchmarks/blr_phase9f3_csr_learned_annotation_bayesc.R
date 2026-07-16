# Post-migration learned-annotation CSR BayesC baseline; run manually from root.
pkgload::load_all(".", compile = FALSE)
source(file.path("tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"))

rss_mib <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2
}

make_work <- function(m = 2000L) {
  ids <- paste0("m", seq_len(m))
  A <- cbind(
    intercept = 1,
    coding = as.numeric(seq_len(m) %% 5L == 0L),
    qtl = as.numeric(seq_len(m) %% 13L == 0L)
  )
  rownames(A) <- ids
  list(
    stats = list(
      wy = list(T1 = stats::setNames(10 * sin(seq_len(m) / 23), ids)),
      ww = list(T1 = stats::setNames(rep(250, m), ids)), yy = c(T1 = 250),
      n = 250L, m = m, marker_names = ids, trait_names = "T1"
    ),
    prefix = phase9a_prefix(m), A = A
  )
}

run_fit <- function(w, mode, chains, cores, keep) {
  flags <- switch(mode,
    fixed = c(FALSE, FALSE), probability = c(TRUE, FALSE),
    multiplier = c(FALSE, TRUE), both = c(TRUE, TRUE)
  )
  sblr::stblr_csr_learn_annot(
    stats = w$stats, ld_prefix = w$prefix, A = w$A,
    add_intercept = FALSE, standardize_annotations = FALSE,
    learn_pi_annot = flags[[1]], learn_vb_annot = flags[[2]],
    eta_pi_init = rep(0, ncol(w$A)), eta_vb_init = rep(0, ncol(w$A)),
    rw_sd_eta_pi = 0.02, rw_sd_eta_vb = 0.02, annot_update_every = 5L,
    pi_min = 1e-8, pi_max = 0.5,
    vb_multiplier_min = 1e-3, vb_multiplier_max = 1e3,
    pi_init = 0.3, pi_prior_mean = 0.3, pi_prior_strength = 2,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 100L, nburn = 25L, nthin = 1L,
    nchains = chains, ncores = cores, keep_chains = keep, seed = 919L
  )
}

bench <- function(markers, mode, chains, cores, keep = FALSE, repetitions = 5L) {
  w <- make_work(markers)
  invisible(run_fit(w, mode, chains, cores, keep))
  elapsed <- rss <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    elapsed[[i]] <- system.time(invisible(run_fit(w, mode, chains, cores, keep)))[["elapsed"]]
    rss[[i]] <- rss_mib()
  }
  data.frame(
    markers, traits = 1L, annotations = ncol(w$A), coefficient_count = ncol(w$A),
    iterations = 100L, burnin = 25L, retained = 100L,
    mode, chains, cores, keep_chains = keep, update_every = 5L,
    pi_bounds = "[1e-8,0.5]", multiplier_bounds = "[1e-3,1e3]",
    elapsed = paste(elapsed, collapse = ","), mean = mean(elapsed),
    median = stats::median(elapsed), minimum = min(elapsed), maximum = max(elapsed),
    range = diff(range(elapsed)), completed_fit_rss_mib = max(rss, na.rm = TRUE)
  )
}

tiny <- bench(40L, "both", 1L, 1L, repetitions = 2L)
primary <- do.call(rbind, list(
  bench(2000L, "fixed", 1L, 1L),
  bench(2000L, "probability", 1L, 1L),
  bench(2000L, "multiplier", 1L, 1L),
  bench(2000L, "both", 1L, 1L),
  bench(2000L, "both", 2L, 1L),
  bench(2000L, "both", 2L, 2L, TRUE)
))
print(rbind(tiny, primary), row.names = FALSE)
cat("warm-up: one untimed fit per configuration; five primary repetitions\n")
cat("memory: whole-process RSS sampled after completed fits, not interval peak\n")
cat("environment:", R.version.string, "; sblr",
    as.character(utils::packageVersion("sblr")), ";", R.version$platform, "\n")
