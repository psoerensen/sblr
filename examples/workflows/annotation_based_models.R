# Annotation-based ST-BLR workflow
#
# Fixed-prior, learned annotation, group annotation, and continuous
# overlapping-annotation SBayesRC models are runnable through package wrappers.
#
# The MCMC settings shown here are demonstration settings. Real analyses need
# longer chains and appropriate convergence and posterior predictive checks.

# Packages -----------------------------------------------------------------

library(qgg)
library(sblr)

# Data directory setup -----------------------------------------------------

data_dir <- Sys.getenv("SBLR_EXAMPLE_DATA_DIR", unset = "path/to/example/data")
glist_file <- file.path(data_dir, "Glist_sparseLD_1k.RDS")
ld_prefix <- file.path(data_dir, "annotation_ld")

# Load or prepare Glist ----------------------------------------------------

# Load a qgg Glist whose BED, BIM, and FAM paths are valid on the current
# machine. Prepare and save this object with qgg::gprep() if needed.
if (!file.exists(glist_file)) {
 stop(
  "Example Glist not found. Set SBLR_EXAMPLE_DATA_DIR to a directory ",
  "containing Glist_sparseLD_1k.RDS."
 )
}

Glist <- readRDS(glist_file)
chr <- 1
cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
stopifnot(!anyNA(cls))

# Simulate or load multi-trait phenotype data ------------------------------

# Replace this simulation with a phenotype matrix from the study if desired.
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

# Compute BED sufficient statistics with bed_xtx_xty() --------------------

# Runtime-heavy example.
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

# Runtime-heavy example. This writes disk-backed CSR files at ld_prefix.
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

# Construct a simple annotation matrix A ----------------------------------

# A is an m x K marker annotation matrix. Rows of A must correspond exactly
# to the marker ordering in stats and the sparse-LD files. These deterministic
# toy annotations demonstrate the interface; real analyses should replace A
# with biological annotations in the same marker order.
m <- length(cls)
marker_id <- Glist$rsidsLD[[chr]]
A <- cbind(
 every_tenth = as.numeric(seq_len(m) %% 10 == 0),
 marker_order = as.numeric(scale(seq_len(m)))
)
rownames(A) <- marker_id
stopifnot(nrow(A) == m, identical(rownames(A), marker_id))

# Fixed marker-specific prior model ---------------------------------------

# stblr_csr_prior_annot() uses fixed marker-specific inclusion and variance
# priors derived from A and supplied coefficients. The small coefficients below
# are conservative demonstration values, not analysis recommendations.
fit_prior_annot <- stblr_csr_prior_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 beta_pi = c(0.25, 0.10),
 beta_vb = c(0.10, 0.10),
 use_pi_marker = TRUE,
 use_vb_multiplier = TRUE,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Learned annotation model ------------------------------------------------

# stblr_csr_learn_annot() learns annotation effects on marker inclusion
# probabilities and, optionally, marker-effect variances. These proposal scales
# are conservative demonstration settings. Real analyses require longer chains,
# convergence checks, and posterior diagnostics.
fit_learn_annot <- stblr_csr_learn_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 learn_pi_annot = TRUE,
 learn_vb_annot = TRUE,
 rw_sd_eta_pi = 0.02,
 rw_sd_eta_vb = 0.02,
 annot_update_every = 10,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Group annotation model --------------------------------------------------

# group is a length-m marker group vector in the same order as stats and the
# sparse LD. These are conservative demonstration settings. Real analyses
# require longer chains, convergence checks, and posterior diagnostics.
group <- ifelse(A[, "every_tenth"] == 1, "annotated", "background")
names(group) <- rownames(A)

fit_group_annot <- stblr_csr_group_annot(
 stats = stats,
 ld_prefix = ld_prefix,
 group = group,
 n = Glist$n,
 group_names = c("annotated", "background"),
 group_pi_init = c(0.002, 0.001),
 group_vb_multiplier_init = c(1.1, 1.0),
 updatePi = TRUE,
 updateGroupVb = TRUE,
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# SBayesRC-style annotation model -----------------------------------------

# stblr_csr_sbayesrc_generic() uses continuous overlapping annotations to
# model SBayesRC-style mixture-component probabilities. These are conservative
# demonstration settings. Real analyses require longer chains, convergence
# checks, and posterior diagnostics.
fit_sbayesrc <- stblr_csr_sbayesrc_generic(
 stats = stats,
 ld_prefix = ld_prefix,
 A = A,
 n = Glist$n,
 gamma = c(0, 0.01, 0.1, 1),
 nit = 1000,
 nburn = 100,
 nthin = 1,
 ncores = 3,
 seed = 10
)

# Inspect posterior summaries ---------------------------------------------

# Fixed-prior posterior summaries:
fit_prior_annot$bm
fit_prior_annot$dm
fit_prior_annot$input$pi_marker
fit_prior_annot$input$vb_multiplier

# Learned-annotation posterior summaries:
fit_learn_annot$bm
fit_learn_annot$dm
fit_learn_annot$eta_pi
fit_learn_annot$eta_vb

# Group-annotation posterior summaries:
fit_group_annot$bm
fit_group_annot$dm
fit_group_annot$group_pi
fit_group_annot$group_vb_multiplier

# SBayesRC posterior summaries:
fit_sbayesrc$bm
fit_sbayesrc$dm
fit_sbayesrc$ncomp
fit_sbayesrc$alpha
