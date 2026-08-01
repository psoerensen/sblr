#!/usr/bin/env Rscript

pkgload::load_all(".", quiet = TRUE)
source(file.path("tests", "testthat", "helper-blr-fixtures.R"), local = TRUE)

set.seed(8317)
n <- 80L; m <- 24L
dosage_by_sample <- matrix(rbinom(n * m, 2, rep(seq(0.18, 0.42, length.out = m), each = n)), n, m)
for (marker in seq_len(m)) {
  if (length(unique(dosage_by_sample[, marker])) < 2L)
    dosage_by_sample[seq_len(3L), marker] <- 0:2
}
af <- colMeans(dosage_by_sample) / 2
Z <- sweep(dosage_by_sample, 2L, 2 * af, "-")
Z <- sweep(Z, 2L, sqrt(2 * af * (1 - af)), "/")
truth <- c(rep(0.18, 4L), rep(0, 8L), rep(-0.12, 4L), rep(0, 8L))
y <- drop(Z %*% truth + seq(-0.7, 0.7, length.out = n))
y <- y - mean(y)
score <- drop(crossprod(Z, y))
marker <- paste0("rs", seq_len(m))
bed <- tempfile(fileext = ".bed")
blr_block_fixture_write_bed(bed, t(dosage_by_sample))
Glist <- list(
  n = n, ids = paste0("id", seq_len(n)), bedfiles = bed,
  rsids = list(marker), rsidsLD = list(marker), chr = list(rep(1L, m)),
  pos = list(seq_len(m) * 100), af = list(af), maf = list(pmin(af, 1 - af))
)
stats <- list(
  wy = list(T1 = stats::setNames(score, marker)),
  ww = list(T1 = stats::setNames(colSums(Z * Z), marker)),
  yy = stats::setNames(sum(y * y), "T1"), n = n, m = m,
  bed_files = bed, cls = list(seq_len(m)), rows = seq_len(n), af = list(af),
  marker_names = marker, trait_names = "T1"
)
block_start <- seq.int(1L, m, by = 6L)
annotation <- cbind(intercept = 1, informative = as.numeric(seq_len(m) <= 12L))
common <- list(
  stats = stats, Glist = Glist, block_start = block_start,
  nit = 240L, nburn = 120L, nthin = 2L, nchains = 4L, ncores = 4L,
  keep_chains = TRUE, seed = 9921L, updateB = TRUE, updateE = TRUE
)

fit_pair <- function(method) {
  model <- switch(method,
    sbayesc = list(method = method, pi_init = 0.25, pi_prior_mean = 0.25,
                   pi_prior_strength = 20, updatePi = TRUE),
    sbayesr = list(method = method, updatePi = TRUE),
    sbayesrc = list(method = method, annotation = annotation, updateAlpha = TRUE)
  )
  dense <- do.call(stblr_block_eigen, c(common, model, list(
    representation = "dense_reconstructed", eigen_policy = "absolute_threshold",
    eigen_tau = 1e-12
  )))
  retained <- do.call(stblr_block_eigen, c(common, model, list(
    representation = "low_rank", eigen_policy = "cumulative_positive_mass",
    eigen_prop = 1 - 1e-12
  )))
  list(dense = dense, retained = retained)
}

metric_difference <- function(a, b, field) {
  if (!field %in% names(a) || !field %in% names(b) ||
      is.null(a[[field]]) || is.null(b[[field]])) return(NA_real_)
  av <- as.numeric(unlist(a[[field]], use.names = FALSE))
  bv <- as.numeric(unlist(b[[field]], use.names = FALSE))
  if (length(av) != length(bv)) return(Inf)
  max(abs(av - bv), na.rm = TRUE)
}

limits <- c(bm = 0.12, dm = 0.20, pi_mean = 0.20, vgs = 0.35,
            ves = 0.35, component_probabilities = 0.25)
rows <- list()
for (method in c("sbayesc", "sbayesr", "sbayesrc")) {
  pair <- fit_pair(method)
  differences <- vapply(names(limits), function(field)
    metric_difference(pair$dense, pair$retained, field), numeric(1))
  finite <- is.finite(differences)
  if (any(differences[finite] > limits[finite])) {
    stop(method, " full-rank posterior comparison exceeded its predeclared Monte Carlo envelope: ",
         paste(names(differences)[finite], signif(differences[finite], 4), collapse = ", "))
  }
  rows[[method]] <- data.frame(method = method, metric = names(differences),
                               max_abs_difference = differences,
                               envelope = unname(limits), row.names = NULL)
}
result <- do.call(rbind, rows)
print(result, row.names = FALSE)
cat("four-chain full-rank scalar sampler comparisons passed\n")
