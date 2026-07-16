# Post-migration group CSR BayesC baseline; run manually from the package root.
pkgload::load_all(".", compile = FALSE)
source(file.path("tests", "testthat", "fixtures", "blr-phase9a-annotation-reference.R"))

rss_mib <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2
}
make_work <- function(m) {
  ids <- paste0("m", seq_len(m))
  list(stats = list(wy = list(T1 = stats::setNames(10 * sin(seq_len(m) / 23), ids)),
                    ww = list(T1 = stats::setNames(rep(250, m), ids)),
                    yy = c(T1 = 250), n = 250L, m = m,
                    marker_names = ids, trait_names = "T1"),
       prefix = phase9a_prefix(m),
       group = stats::setNames(ifelse(seq_len(m) %% 5 == 1, "coding", "background"), ids))
}
run_group <- function(w, chains, cores, keep, update_pi, update_vb, normalize, nit, nburn) {
  sblr::stblr_csr_group_annot(
    stats = w$stats, ld_prefix = w$prefix, group = w$group,
    group_names = c("coding", "background"), group_pi_init = c(.4, .2),
    group_vb_multiplier_init = c(1.4, .8), updateGroupVb = update_vb,
    normalize_group_vb = normalize, pi_init = .3, pi_prior_mean = .3,
    pi_prior_strength = 2, updateB = FALSE, updateE = FALSE, updatePi = update_pi,
    nit = nit, nburn = nburn, nthin = 1L, nchains = chains, ncores = cores,
    keep_chains = keep, seed = 919L, updateLDswap = FALSE
  )
}
bench <- function(markers, chains, cores, keep, update_pi, update_vb, normalize,
                  nit = 100L, nburn = 25L, repetitions = 5L) {
  w <- make_work(markers)
  invisible(run_group(w, chains, cores, keep, update_pi, update_vb, normalize, nit, nburn))
  elapsed <- rss <- numeric(repetitions)
  for (i in seq_len(repetitions)) {
    elapsed[i] <- system.time(invisible(
      run_group(w, chains, cores, keep, update_pi, update_vb, normalize, nit, nburn)
    ))[["elapsed"]]
    rss[i] <- rss_mib()
  }
  data.frame(markers, traits = 1L, groups = 2L, mapping = "20% coding / 80% background",
             iterations = nit, burnin = nburn, retained = nit, chains, cores,
             keep_chains = keep, updatePi = update_pi, updateGroupVb = update_vb,
             normalize_group_vb = normalize, elapsed = paste(elapsed, collapse = ","),
             mean = mean(elapsed), median = median(elapsed), minimum = min(elapsed),
             maximum = max(elapsed), iqr = stats::IQR(elapsed),
             completed_fit_rss_mib = max(rss, na.rm = TRUE))
}
rows <- list(
  bench(4L, 1L, 1L, FALSE, FALSE, TRUE, TRUE, 8L, 2L),
  bench(2000L, 1L, 1L, FALSE, FALSE, FALSE, TRUE),
  bench(2000L, 1L, 1L, FALSE, FALSE, TRUE, TRUE),
  bench(2000L, 1L, 1L, FALSE, FALSE, TRUE, FALSE),
  bench(2000L, 1L, 1L, FALSE, TRUE, TRUE, TRUE),
  bench(2000L, 2L, 1L, FALSE, FALSE, TRUE, TRUE),
  bench(2000L, 2L, 2L, TRUE, FALSE, TRUE, TRUE)
)
print(do.call(rbind, rows), row.names = FALSE)
cat("warm-up/order: one untimed fit per row; rows executed in printed order; five repetitions\n")
cat("memory: whole-process RSS after each completed fit, not interval peak\n")
cat("measurement: base system.time elapsed and ps::ps_memory_info RSS\n")
cat("environment:", R.version.string, "; sblr", as.character(utils::packageVersion("sblr")),
    ";", R.version$platform, "\n")
