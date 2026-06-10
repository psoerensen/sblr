# Core sparse-LD and PLINK BED workflow
#
# Set SBLR_EXAMPLE_DATA_DIR to a directory containing the qgdata example files.
# Set SBLR_RUN_HEAVY_EXAMPLES=true to run the benchmark-scale computations.

# Packages -----------------------------------------------------------------

library(qgg)
library(sblr)

# Data setup ---------------------------------------------------------------

data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR", unset = "path/to/example/data")
glist_file <- file.path(data_dir, "Glist_sparseLD_1k.RDS")
run_heavy <- identical(Sys.getenv("SBLR_RUN_HEAVY_EXAMPLES"), "true")

# Interactive data download and Glist preparation. The qgdata URLs are kept
# unchanged so this block can be run directly when preparing the example data.
if (interactive() && !file.exists(glist_file)) {
 dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
 files <- c("bed", "bim", "fam", "pheno", "covar")

 for (f in files) {
  url <- sprintf(
   "https://github.com/psoerensen/qgdata/raw/main/simulated_human_data/human.%s",
   f
  )
  download.file(
   url,
   destfile = file.path(data_dir, paste0("human.", f)),
   mode = "wb"
  )
 }

 Glist <- gprep(
  study = "Example",
  bedfiles = file.path(data_dir, "human.bed"),
  bimfiles = file.path(data_dir, "human.bim"),
  famfiles = file.path(data_dir, "human.fam")
 )
 rsids <- gfilter(
  Glist = Glist,
  excludeMAF = 0.05,
  excludeMISS = 0.05,
  excludeCGAT = TRUE,
  excludeINDEL = TRUE,
  excludeDUPS = TRUE,
  excludeHWE = 1e-12,
  excludeMHC = FALSE
 )
 Glist <- gprep(
  Glist,
  task = "sparseld",
  msize = 1000,
  rsids = rsids,
  ldfiles = file.path(data_dir, "human.ld"),
  overwrite = TRUE
 )
 saveRDS(Glist, file = glist_file, compress = FALSE)
}

if (!file.exists(glist_file)) {
 stop(
  "Example Glist not found. Set SBLR_EXAMPLE_DATA_DIR and run the ",
  "interactive data setup block first."
 )
}

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
if (run_heavy) {
 fit_qgg <- gbayes(
  stat = stat_qgg,
  Glist = Glist,
  method = "bayesC",
  nit = 1000
 )
}

# Compute BED sufficient statistics with bed_xtx_xty() --------------------

# Runtime-heavy benchmark.
if (run_heavy) {
 stats <- bed_xtx_xty(
  bed_file = Glist$bedfiles[chr],
  n = Glist$n,
  cls = cls,
  af = Glist$af[[chr]][cls],
  y = y,
  scale = TRUE,
  nthreads = 4
 )
}

# Stream sparse LD with sparseLD_stream_CSR() ------------------------------

# Runtime-heavy benchmark. This writes disk-backed CSR files at ld_prefix.
if (run_heavy) {
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
}

# Fit ST-BLR from summary statistics with stblr_csr() ----------------------

# Runtime-heavy benchmark.
if (run_heavy) {
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
}

# Fit ST-BLR directly from BED markers with stblr_bed_marker() -------------

# Runtime-heavy benchmark.
if (run_heavy) {
 fit_bed <- stblr_bed_marker(
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
}

# Simple diagnostic plots --------------------------------------------------

if (run_heavy) {
 plot(sim$B[, 1], fit_csr$dm[, 1], xlab = "True effect", ylab = "CSR PIP")
 plot(sim$B[, 1], fit_bed$dm[, 1], xlab = "True effect", ylab = "BED PIP")
 plot(fit_csr$bm[, 1], fit_bed$bm[, 1], xlab = "CSR effect", ylab = "BED effect")
}
