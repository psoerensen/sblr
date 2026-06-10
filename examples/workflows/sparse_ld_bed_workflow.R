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
stats <- bed_xtx_xty(
  bed_file = Glist$bedfiles[chr],
  n = Glist$n,
  cls = cls,
  af = Glist$af[[chr]][cls],
  y = y,
  scale = TRUE,
  nthreads = 4
)

# Stream sparse LD with sparseLD_stream_CSR() ------------------------------

# Runtime-heavy benchmark. This writes disk-backed CSR files at ld_prefix.
sparseLD_stream_CSR(
  bed_files = Glist$bedfiles[chr],
  n = Glist$n,
  cls = list(cls),
  out_prefix = ld_prefix,
  rows = NULL,
  af = list(Glist$af[[chr]][cls]),
  pos_bp = NULL,
  max_distance_bp = 0,
  max_distance_variants = 0,
  r2_threshold = 0.001,
  block_size = 1024,
  nthreads = 4
)

ld <- sparseLD_read_CSR(ld_prefix, one_based = FALSE)

# Fit summary-statistics / sparse-LD BLR with stblr_csr() -----------------

# Runtime-heavy benchmark.
fit_csr <- stblr_csr(
  stats = stats,
  ld_prefix = ld_prefix,
  n = Glist$n,
  pi_marker = 0.001,
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
  pi_marker = 0.001,
  h2 = 0.5,
  nit = 1000,
  nburn = 100,
  nthin = 1,
  ncores = 3,
  seed = 10,
  scheduled = TRUE
)

# Simple diagnostic plots --------------------------------------------------

plot(sim$B[, 1], fit_csr$dm[, 1], xlab = "True effect", ylab = "CSR PIP")
plot(sim$B[, 1], fit_bed_sched_single$dm[, 1], xlab = "True effect", ylab = "BED PIP")
plot(fit_csr$bm[, 1], fit_bed_sched_single$bm[, 1], xlab = "CSR effect", ylab = "BED effect")
