# Annotation-based ST-BLR workflow template
#
# The high-level annotation wrappers used below are currently prototypes in
# examples/mtsim_prior_sparse_ld.R. They must be promoted into R/ and exported
# before the model-fitting calls in this workflow can run through library(sblr).
# The calls are therefore commented templates rather than runnable examples.
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

# Load a qgg Glist whose BED, BIM, and FAM paths are portable or valid on the
# current machine. Prepare and save this object with qgg::gprep() if needed.
#
# Glist <- readRDS(glist_file)
# chr <- 1
# cls <- match(Glist$rsidsLD[[chr]], Glist$rsids[[chr]])
# stopifnot(!anyNA(cls))

# Simulate or load multi-trait phenotype data ------------------------------

# Replace this simulation with a phenotype matrix from the study if desired.
#
# sim <- mtsim(
#  Glist = Glist,
#  chr = chr,
#  rsids = Glist$rsidsLD[[chr]],
#  nt = 3,
#  n_shared = 30,
#  n_specific = 10,
#  h2 = c(0.4, 0.5, 0.3),
#  rg = matrix(
#   c(1.0, 0.7, 0.3, 0.7, 1.0, 0.5, 0.3, 0.5, 1.0),
#   nrow = 3,
#   byrow = TRUE
#  ),
#  re = 0,
#  seed = 1
# )
# y <- as.matrix(scale(sim$y))

# Compute BED sufficient statistics with bed_xtx_xty() --------------------

# Runtime-heavy example.
#
# stats <- bed_xtx_xty(
#  bed_file = Glist$bedfiles[chr],
#  n = Glist$n,
#  cls = cls,
#  af = Glist$af[[chr]][cls],
#  y = y,
#  scale = TRUE,
#  nthreads = 4
# )

# Stream sparse LD with sparseLD_stream_CSR() ------------------------------

# Runtime-heavy example. This writes disk-backed CSR files at ld_prefix.
#
# sparseLD_stream_CSR(
#  bed_files = Glist$bedfiles[chr],
#  n = Glist$n,
#  cls = list(cls),
#  out_prefix = ld_prefix,
#  rows = NULL,
#  af = list(Glist$af[[chr]][cls]),
#  pos_bp = NULL,
#  max_distance_bp = 0,
#  max_distance_variants = 0,
#  r2_threshold = 0.001,
#  block_size = 1024,
#  nthreads = 4
# )

# Construct a simple annotation matrix A ----------------------------------

# A is an m x K marker annotation matrix. Rows of A must correspond exactly
# to the marker ordering in stats and the sparse-LD files.
#
# m <- length(cls)
# marker_id <- Glist$rsidsLD[[chr]]
# A <- cbind(
#  first_half = as.numeric(seq_len(m) <= ceiling(m / 2)),
#  marker_order = as.numeric(scale(seq_len(m)))
# )
# rownames(A) <- marker_id
#
# group is a length-m marker group vector in the same marker order.
# group <- ifelse(A[, "first_half"] == 1, "first_half", "second_half")

# Fixed marker-specific prior model ---------------------------------------

# model = "prior" uses fixed marker-specific inclusion and variance priors
# derived from A and the supplied annotation coefficients.
#
# fit_prior <- stblr_csr_annotation(
#  stats = stats,
#  ld_prefix = ld_prefix,
#  model = "prior",
#  A = A,
#  n = Glist$n,
#  beta_pi = c(0.5, -0.2),
#  beta_vb = c(0.0, 0.3),
#  use_pi_marker = TRUE,
#  use_vb_multiplier = TRUE,
#  nit = 1000,
#  nburn = 100,
#  nthin = 1,
#  ncores = 3,
#  seed = 10
# )

# Learned annotation model ------------------------------------------------

# model = "annot" learns annotation effects on marker inclusion probabilities
# and, optionally, marker-effect variances.
#
# fit_annot <- stblr_csr_annotation(
#  stats = stats,
#  ld_prefix = ld_prefix,
#  model = "annot",
#  A = A,
#  n = Glist$n,
#  learn_pi_annot = TRUE,
#  learn_vb_annot = TRUE,
#  nit = 1000,
#  nburn = 100,
#  nthin = 1,
#  ncores = 3,
#  seed = 10
# )

# Group annotation model --------------------------------------------------

# model = "group" uses group-level priors. group is a length-m marker group
# vector whose entries correspond to the markers in stats and the sparse LD.
#
# fit_group <- stblr_csr_annotation(
#  stats = stats,
#  ld_prefix = ld_prefix,
#  model = "group",
#  group = group,
#  n = Glist$n,
#  updatePi = TRUE,
#  updateGroupVb = TRUE,
#  nit = 1000,
#  nburn = 100,
#  nthin = 1,
#  ncores = 3,
#  seed = 10
# )

# SBayesRC-style annotation model -----------------------------------------

# model = "sbayesrc" uses SBayesRC-style mixture components with annotations.
#
# fit_sbayesrc <- stblr_csr_annotation(
#  stats = stats,
#  ld_prefix = ld_prefix,
#  model = "sbayesrc",
#  A = A,
#  n = Glist$n,
#  gamma = c(0, 0.01, 0.1, 1),
#  nit = 1000,
#  nburn = 100,
#  nthin = 1,
#  ncores = 3,
#  seed = 10
# )

# Inspect posterior summaries ---------------------------------------------

# Common posterior summaries:
#
# fit_prior$bm
# fit_prior$dm
# fit_annot$eta_pi
# fit_annot$eta_vb
# fit_group$group_pi
# fit_group$group_vb_multiplier
# fit_sbayesrc$ncomp
# fit_sbayesrc$alpha
#
# Equivalent convenience wrappers are intended to be:
# stblr_csr_prior_annot(), stblr_csr_learn_annot(),
# stblr_csr_group_annot(), and stblr_csr_sbayesrc_generic().
