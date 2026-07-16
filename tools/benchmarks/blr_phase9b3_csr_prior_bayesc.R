# Post-migration fixed-prior CSR BayesC baseline; run from the package root.
pkgload::load_all(".", compile = FALSE)
source(file.path("tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"))

rss_mib <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2
}

make_work <- function(m) {
  ids <- paste0("m", seq_len(m))
  list(
    stats = list(
      wy = list(T1 = stats::setNames(10 * sin(seq_len(m) / 23), ids)),
      ww = list(T1 = stats::setNames(rep(250, m), ids)),
      yy = c(T1 = 250), n = 250L, m = m,
      marker_names = ids, trait_names = "T1"
    ),
    prefix = phase9a_prefix(m),
    pi_marker = list(rep(0.3, m)),
    vb_multiplier = list(ifelse(seq_len(m) %% 5 == 0, 1.5, 0.9))
  )
}

run_prior <- function(work, chains, cores, keep, policy, update_b, update_pi,
                      nit, nburn) {
  args <- list(
    stats = work$stats, ld_prefix = work$prefix,
    pi_init = 0.3, pi_prior_mean = 0.3, pi_prior_strength = 2,
    updateB = update_b, updateE = FALSE, updatePi = update_pi,
    nit = nit, nburn = nburn, nthin = 1L, nchains = chains,
    ncores = cores, keep_chains = keep, seed = 919L,
    fixed_pi_marker = work$pi_marker,
    fixed_vb_multiplier = work$vb_multiplier,
    use_pi_marker = policy %in% c("pi_marker", "both"),
    use_vb_multiplier = policy %in% c("vb_multiplier", "both")
  )
  do.call(sblr::stblr_csr_prior_annot, args)
}

bench <- function(markers, nit, nburn, chains, cores, keep, policy,
                  update_b = FALSE, update_pi = FALSE, repetitions = 5L) {
  work <- make_work(markers)
  invisible(run_prior(work, chains, cores, keep, policy, update_b, update_pi, nit, nburn))
  elapsed <- rss <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    elapsed[i] <- system.time(invisible(
      run_prior(work, chains, cores, keep, policy, update_b, update_pi, nit, nburn)
    ))[["elapsed"]]
    rss[i] <- rss_mib()
  }
  data.frame(
    markers, traits = 1L, iterations = nit, burnin = nburn,
    retained = nit, chains, cores, keep_chains = keep,
    updateB = update_b, updatePi = update_pi, policy,
    elapsed = paste(elapsed, collapse = ","), mean = mean(elapsed),
    median = stats::median(elapsed), minimum = min(elapsed), maximum = max(elapsed),
    iqr = stats::IQR(elapsed), completed_fit_rss_mib = max(rss, na.rm = TRUE)
  )
}

rows <- list(
  bench(4L, 8L, 2L, 1L, 1L, FALSE, "both"),
  bench(2000L, 100L, 25L, 1L, 1L, FALSE, "pi_marker"),
  bench(2000L, 100L, 25L, 1L, 1L, FALSE, "vb_multiplier"),
  bench(2000L, 100L, 25L, 1L, 1L, FALSE, "both"),
  bench(2000L, 100L, 25L, 2L, 1L, FALSE, "both"),
  bench(2000L, 100L, 25L, 2L, 2L, TRUE, "both"),
  bench(2000L, 100L, 25L, 1L, 1L, FALSE, "both", update_b = TRUE),
  bench(2000L, 100L, 25L, 1L, 1L, FALSE, "both", update_pi = TRUE)
)

result <- do.call(rbind, rows)
print(result, row.names = FALSE)
cat("warm-up/order: one untimed fit per row; rows executed in printed order; five repetitions\n")
cat("memory: whole-process RSS after each completed fit, not interval peak\n")
cat("measurement: base system.time elapsed and ps::ps_memory_info RSS\n")
cat("environment:", R.version.string, "; sblr",
    as.character(utils::packageVersion("sblr")), ";", R.version$platform, "\n")
