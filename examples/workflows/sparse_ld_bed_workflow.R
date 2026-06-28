# Summary-statistics/sparse-LD and individual-level BED BLR workflow
#
# Set SBLR_EXAMPLE_DATA_DIR to a directory containing the qgdata example files.
# Set SBLR_RUN_HEAVY_EXAMPLES=true to run the benchmark-scale computations.

# Packages -----------------------------------------------------------------
library(sblr)

# Data setup ---------------------------------------------------------------
data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR")
if (!nzchar(data_dir)) {
  data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"
}

nthreads <- 4


Glist <- readRDS(file.path(data_dir, "Glist_sparseLD_1k.RDS"))



# Simulate traits ----------------------------------------------------------
sim <- mtsim(
  Glist = Glist,
  chr = 1,
  rsids = Glist$rsidsLD[[1]],
  nt = 3,
  n_shared = 30,
  n_specific = 10,
  h2 = c(0.4, 0.5, 0.3),
  rg = matrix(
    c(
      1.0, 0.7, 0.3,
      0.7, 1.0, 0.5,
      0.3, 0.5, 1.0
    ),
    nrow = 3,
    byrow = TRUE
  ),
  re = 0,
  seed = 1
)
y <- as.matrix(scale(sim$y))


# Compute summary statistics

system.time(
  stats <- make_stats(
    Glist,
    y,
    nthreads = nthreads
  )
)

# Compute sparse LD
ld_prefix <- file.path(data_dir, "ld_test")
system.time(Glist <- make_sparseLD(
  Glist = Glist,
  out_prefix = ld_prefix,
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 1000,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = 4
))


# Fit summary-statistics / sparse-LD BLR model
fit <- stblr_csr(
  stats = stats,
  Glist = Glist,
  ## Conservative sparse architecture.
  pi_init = 0.001,
  pi_prior_a = 1,
  pi_prior_b = 1,
  h2 = 0.3,
  adjE = 0.9,
  nit = 10000,
  nburn = 1000,
  ncores = 3,
  seed = 10,
  scheduled = FALSE
)

fitMH <- stblr_csr(
  stats = stats,
  Glist = Glist,
  pi_init = 0.001,
  pi_prior_a = 1,
  pi_prior_b = 1,
  h2 = 0.3,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  ncores = 3,
  nchains = 3,
  seed = 10,
  scheduled = FALSE,
  updateLDswap = TRUE,
  ld_swap_prob = 0.10,
  ld_swap_r2 = 0.05,
  ld_swap_moves = 5
)

cs_global <- make_stblr_credible_sets(
  fit = fitMH,
  Glist = Glist,
  trait = "D1",
  coverage = 0.95,
  min_r2 = 0.5,
  pip_cutoff = 0.001,
  locus_pip_cutoff = 0.01,
  max_locus_distance = 1e6
)

fm_extract <- extract_stblr_finemap_loci(
  fit = fitMH,
  Glist = Glist,
  locus_sets = cs_global$locus_sets,
  trait = "D1",
  credible_sets = TRUE,
  coverage = 0.95,
  min_r2 = 0.5,
  pip_cutoff = 0.001,
  cs_mode = "multi",
  min_signal_pip = 0.05
)

fm_extract$summary
fm_extract$markers
fm_extract$credible_sets$summary
colMeans(fit$dm)
colMeans(fitMH$dm)


summarise_stblr_posterior(fit)
summarise_stblr_posterior(fitMH)


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

## ------------------------------------------------------------
## Simple ST-BLR backend and fine-mapping checks
## ------------------------------------------------------------

## 1. Basic backend consistency checks
chk_fit <- check_stblr_backend_consistency(
  fit,
  require_chain_summaries = FALSE,
  require_ld_swap = FALSE
)

chk_fitMH <- check_stblr_backend_consistency(
  fitMH,
  require_chain_summaries = TRUE,
  require_ld_swap = TRUE
)

print(chk_fit)
print(chk_fitMH)

stopifnot(isTRUE(chk_fit$ok))
stopifnot(isTRUE(chk_fitMH$ok))


## ------------------------------------------------------------
## 2. Basic dimensions and metadata
## ------------------------------------------------------------

cat("\n--- Basic dimensions ---\n")
cat("dim(fit$dm):   ", paste(dim(fit$dm), collapse = " x "), "\n")
cat("dim(fitMH$dm): ", paste(dim(fitMH$dm), collapse = " x "), "\n")

cat("\n--- Chain metadata ---\n")
cat("fit nchains:   ", fit$input$nchains %||% 1L, "\n")
cat("fitMH nchains: ", fitMH$input$nchains %||% 1L, "\n")

cat("\n--- Scheduled flags ---\n")
cat("fit scheduled:   ", fit$input$scheduled %||% NA, "\n")
cat("fitMH scheduled: ", fitMH$input$scheduled %||% NA, "\n")


## ------------------------------------------------------------
## 3. Check chain summary fields in fitMH
## ------------------------------------------------------------

chain_fields <- c(
  "dm_sd", "dm_min", "dm_max",
  "bm_sd", "bm_min", "bm_max"
)

cat("\n--- Chain summary fields in fitMH ---\n")
print(sapply(chain_fields, function(x) x %in% names(fitMH)))

stopifnot(all(chain_fields %in% names(fitMH)))

stopifnot(identical(dim(fitMH$dm), dim(fitMH$dm_sd)))
stopifnot(identical(dim(fitMH$dm), dim(fitMH$dm_min)))
stopifnot(identical(dim(fitMH$dm), dim(fitMH$dm_max)))

stopifnot(identical(dim(fitMH$bm), dim(fitMH$bm_sd)))
stopifnot(identical(dim(fitMH$bm), dim(fitMH$bm_min)))
stopifnot(identical(dim(fitMH$bm), dim(fitMH$bm_max)))

stopifnot(all(fitMH$dm_sd >= 0, na.rm = TRUE))
stopifnot(all(fitMH$bm_sd >= 0, na.rm = TRUE))

stopifnot(all(fitMH$dm_min <= fitMH$dm + 1e-12, na.rm = TRUE))
stopifnot(all(fitMH$dm <= fitMH$dm_max + 1e-12, na.rm = TRUE))

stopifnot(all(fitMH$bm_min <= fitMH$bm + 1e-12, na.rm = TRUE))
stopifnot(all(fitMH$bm <= fitMH$bm_max + 1e-12, na.rm = TRUE))


## ------------------------------------------------------------
## 4. Check LD-swap diagnostics
## ------------------------------------------------------------

cat("\n--- LD-swap diagnostics ---\n")
print(fitMH$ld_swap)

if (!is.null(fitMH$ld_swap)) {
  stopifnot(all(fitMH$ld_swap$attempted >= fitMH$ld_swap$accepted))
  stopifnot(all(fitMH$ld_swap$accepted >= 0))
  stopifnot(all(fitMH$ld_swap$acceptance_rate >= 0))
  stopifnot(all(fitMH$ld_swap$acceptance_rate <= 1))
}


## ------------------------------------------------------------
## 5. Compare posterior inclusion levels
## ------------------------------------------------------------

cat("\n--- Mean PIP per trait ---\n")

pip_compare <- data.frame(
  trait = colnames(fit$dm),
  scheduled_mean_pip = colMeans(fit$dm),
  exact_mh_mean_pip = colMeans(fitMH$dm)
)

print(pip_compare)


## ------------------------------------------------------------
## 6. Check credible-set object
## ------------------------------------------------------------

cat("\n--- Global credible-set object ---\n")

stopifnot(is.list(cs_global))
stopifnot(all(c("summary", "sets", "locus_sets", "loci", "parameters") %in% names(cs_global)))

cat("Number of detected loci: ", nrow(cs_global$summary), "\n")
cat("Number of locus sets:    ", length(cs_global$locus_sets), "\n")

print(head(cs_global$summary, 10))


## ------------------------------------------------------------
## 7. Check extracted fine-mapping object
## ------------------------------------------------------------

cat("\n--- Fine-mapping extraction object ---\n")

stopifnot(is.list(fm_extract))
stopifnot(all(c("summary", "markers", "credible_sets", "loci", "parameters") %in% names(fm_extract)))

cat("Number of extracted loci:   ", nrow(fm_extract$summary), "\n")
cat("Number of marker rows:      ", nrow(fm_extract$markers), "\n")

print(head(fm_extract$summary, 20))


## ------------------------------------------------------------
## 8. Confirm chain summaries propagated into fine-mapping output
## ------------------------------------------------------------

cat("\n--- Fine-mapping chain-summary propagation ---\n")

required_marker_cols <- c(
  "pip_mean", "pip_sd", "pip_min", "pip_max",
  "bm_mean", "bm_sd", "bm_min", "bm_max"
)

print(sapply(required_marker_cols, function(x) x %in% names(fm_extract$markers)))

stopifnot(all(required_marker_cols %in% names(fm_extract$markers)))

stopifnot(any(!is.na(fm_extract$markers$pip_sd)))
stopifnot(any(!is.na(fm_extract$markers$pip_min)))
stopifnot(any(!is.na(fm_extract$markers$pip_max)))

stopifnot(any(!is.na(fm_extract$markers$bm_sd)))
stopifnot(any(!is.na(fm_extract$markers$bm_min)))
stopifnot(any(!is.na(fm_extract$markers$bm_max)))


## ------------------------------------------------------------
## 9. Summarize locus classes
## ------------------------------------------------------------

fm_extract$summary$class <- with(fm_extract$summary, ifelse(
  lead_pip >= 0.95 & secondary_pip < 0.05,
  "single strong signal",
  ifelse(
    lead_pip >= 0.95 & secondary_pip >= 0.05,
    "strong lead + secondary mass",
    ifelse(
      total_pip >= 0.95,
      "distributed/partial credible signal",
      "weak or exploratory signal"
    )
  )
))

cat("\n--- Locus classes ---\n")
print(table(fm_extract$summary$class))

cat("\n--- Strong or distributed loci ---\n")
print(
  subset(
    fm_extract$summary,
    class != "weak or exploratory signal"
  )
)


## ------------------------------------------------------------
## 10. Credible-set summary
## ------------------------------------------------------------

cat("\n--- Credible-set summary ---\n")

if (!is.null(fm_extract$credible_sets$summary)) {
  print(head(fm_extract$credible_sets$summary, 50))
  
  cat("\nComplete credible sets:\n")
  print(table(fm_extract$credible_sets$summary$complete))
  
  cat("\nNumber of signals per locus:\n")
  print(table(fm_extract$credible_sets$summary$locus))
}


## ------------------------------------------------------------
## 11. Final message
## ------------------------------------------------------------

cat("\nAll basic ST-BLR backend and fine-mapping checks passed.\n")

matplot(fit$pis)
matplot(fit$ves)
matplot(fit$vld)
matplot(fit$vle)
matplot(fit$vgs)



post <- sblr:::summarise_stblr_posterior(fit)

plot_stblr_posterior(
  post,
  parameters = c("vg", "ve", "vb", "varch")
)

plot_stblr_posterior(
  post,
  parameters = c("pi","m_included")
)

cs <- make_stblr_credible_sets(
  fit = fit,
  Glist = Glist,
  trait = "D1",
  coverage = 0.95,
  min_r2 = 0.5,
  pip_cutoff = 0.001,
  locus_pip_cutoff = 0.01,
  max_locus_distance = 1e6,
  method = "pip"
)


cs$loci
cs$summary
cs$sets[[1]]

cs_global <- make_stblr_credible_sets(
  fit = fit,
  Glist = Glist,
  trait = "D1",
  locus_pip_cutoff = 0.01,
  max_locus_distance = 1e6
)

str(cs_global$locus_sets, max.level = 1)
length(cs_global$locus_sets)
head(names(cs_global$locus_sets))
lengths(cs_global$locus_sets)[1:10]

fm <- finemap_stblr_csr(
  fit = fit,
  Glist = Glist,
  stats = stats,
  sets = cs_global$locus_sets,
  trait = "D1",
  nruns = 8,
  use_residual = TRUE,
  verbose = FALSE
)

fm <- finemap_stblr_csr(
  fit = fit,
  Glist = Glist,
  stats = stats,
  sets = cs_global$locus_sets,
  trait = "D1",
  nruns = 8,
  use_residual = TRUE,   # use FALSE first unless residualization is fully implemented
  nit = 20000,
  nburn = 2000,
  seeds = 101:108,
  credible_sets = TRUE,
  coverage = 0.95,
  min_r2 = 0.5,
  pip_cutoff = 0.001
)

fm_sum$secondary_pip <- fm_sum$total_pip - fm_sum$lead_pip

fm_sum$class <- with(fm_sum, ifelse(
  lead_pip >= 0.95 & secondary_pip < 0.1,
  "single strong signal",
  ifelse(
    lead_pip >= 0.95 & secondary_pip >= 0.1,
    "strong lead + secondary mass",
    ifelse(
      lead_pip < 0.95 & total_pip >= 0.95,
      "distributed credible signal",
      "weak/uncertain"
    )
  )
))

table(fm_sum$class)

fm_sum[order(fm_sum$class, -fm_sum$total_pip), c(
  "locus", "chr", "lead_marker",
  "lead_pip", "total_pip", "secondary_pip", "class"
)]



# Individual-level BED BLR models =========================================

# High-level BED wrapper.
#
# Here cls is omitted. The wrapper uses:
#   match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
#
# For subset phenotypes, rownames(y) or names(y) are matched to Glist$ids.

fit_bed <- stblr_bed_marker(
  Glist = Glist,
  y = y,
  ## Same conservative sparse architecture.
  ## Conservative sparse architecture.
  pi_init = 0.001,
  pi_prior_a = 1,
  pi_prior_b = 1,
  h2 = 0.3,
  nit = 1000,
  nburn = 100,
  ncores = nthreads,
  seed = 10
)


matplot(fit_bed$pis)
matplot(fit_bed$ves)
matplot(fit_bed$vld)
matplot(fit_bed$vle)
matplot(fit_bed$vgs)


post <- summarise_stblr_posterior(fit_bed)

plot_stblr_posterior(
  post,
  parameters = c("vg", "ve", "vb", "varch")
)


# Optional whole-genome multi-chain BED example ----------------------------
# This illustrates the new automatic backend selection.

if (run_heavy) {
  fit_bed_wg <- stblr_bed_marker(
    Glist = Glist,
    y = y[, 1],
    pi_init = 0.001,
    pi_prior_mean = 0.001,
    pi_prior_strength = 5e5,
    h2 = 0.3,
    nchains = 2,
    ncores = 2,
    nit = 1000,
    nburn = 100,
    read_block_size = 64,
    progress_every = 100,
    seed = 10
  )
}

# Simple diagnostic plots --------------------------------------------------

plot(
  sim$B[, 1],
  fit_csr$dm[, 1],
  xlab = "True effect",
  ylab = "CSR PIP"
)

plot(
  sim$B[, 1],
  fit_bed$dm[, 1],
  xlab = "True effect",
  ylab = "BED PIP"
)

plot(
  fit_csr$bm[, 1],
  fit_bed$bm[, 1],
  xlab = "CSR effect",
  ylab = "BED effect"
)










# Summary-statistics / sparse-LD and individual-level BED BLR workflow
#
# Set SBLR_EXAMPLE_DATA_DIR to a directory containing the qgdata example files.
# Set SBLR_RUN_HEAVY_EXAMPLES=true to run the benchmark-scale computations.

# Packages -----------------------------------------------------------------

library(qgg)
library(sblr)

# Data setup ---------------------------------------------------------------

data_dir <- "C:/Users/au223366/Documents/GitHub/examples/human"
# Prepare Glist and sparse LD ----------------------------------------------

Glist <- readRDS(file = file.path(data_dir, "Glist_sparseLD_1k.RDS"))
chr <- 1
cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
stopifnot(!anyNA(cls))

ld_prefix <- file.path(data_dir, "ld_test")

# Simulate traits ----------------------------------------------------------

sim <- mtsim(
 Glist = Glist,
 chr = chr,
 rsids = Glist$rsidsLD[[chr]],
 nt = 3,
 n_shared = 30,
 n_specific = 10,
 h2 = c(0.4, 0.5, 0.3),
 rg = matrix(
  c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0),
  nrow = 3,
  byrow = TRUE
 ),
 re = 0,
 seed = 1
)
y <- as.matrix(scale(sim$y))

# qgg comparison using glma() and gbayes() --------------------------------

stat_qgg <- glma(
 y = y[, 1],
 rsids = Glist$rsidsLD[[chr]],
 Glist = Glist
)

# Runtime-heavy benchmark.
fit_qgg <- gbayes(
  stat = stat_qgg,
  Glist = Glist,
  method = "bayesC",
  nit = 1000
)

# Summary-statistics / sparse-LD BLR models ================================

# Compute summary statistics from BED with bed_xtx_xty() ------------------

# Runtime-heavy benchmark.
system.time(stats <- bed_xtx_xty(
  bed_file = Glist$bedfiles[chr],
  n = Glist$n,
  cls = cls,
  af = Glist$af[[chr]][cls],
  y = y,
  scale = TRUE,
  nthreads = 4
))


# Stream sparse LD with sparseLD_stream_CSR() ------------------------------

# Runtime-heavy benchmark. This writes disk-backed CSR files at ld_prefix.
system.time(sparseLD_stream_CSR(
  bed_files = Glist$bedfiles[chr],
  n = Glist$n,
  cls = list(cls),
  out_prefix = ld_prefix,
  rows = NULL,
  af = list(Glist$af[[chr]][cls]),
  pos_bp = NULL,
  max_distance_bp = 0,           # disables bp-distance filtering
  max_distance_variants = 1000,  # local LD window; 0 disables this filter
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = 4
))

ld <- sparseLD_read_CSR(ld_prefix, one_based = FALSE)

# Fit summary-statistics / sparse-LD BLR with stblr_csr() -----------------

# Runtime-heavy benchmark.
fit_csr <- stblr_csr(
  stats = stats,
  ld_prefix = ld_prefix,
  n = Glist$n,
  h2 = 0.5,
  adjE = 0.9,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,
  scheduled = FALSE
)

# Individual-level BED BLR models =========================================

# Scheduled individual-level BED sampler, single chain.
fit_bed_sched_single <- stblr_bed_marker(
  Glist = Glist,
  y = y,
  chr = chr,
  cls = cls,
  block_size = 1000,
  h2 = 0.5,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10
)

# Simple diagnostic plots --------------------------------------------------

plot(sim$B[, 1], fit_csr$dm[, 1], xlab = "True effect", ylab = "CSR PIP")
plot(sim$B[, 1], fit_bed_sched_single$dm[, 1], xlab = "True effect", ylab = "BED PIP")
plot(fit_csr$bm[, 1], fit_bed_sched_single$bm[, 1], xlab = "CSR effect", ylab = "BED effect")


















