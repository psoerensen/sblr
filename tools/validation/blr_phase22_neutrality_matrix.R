args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}
lib <- get_arg("library")
output <- get_arg("output")
if (is.null(output)) stop("--output is required", call. = FALSE)
if (!is.null(lib)) .libPaths(c(normalizePath(lib, mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(sblr))

write_bed <- function(path, dosage) {
  map <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(j) {
    z <- unname(map[as.character(dosage[j, ])])
    z <- c(z, rep(0L, (-length(z)) %% 4L))
    vapply(seq(1L, length(z), by = 4L), function(i)
      sum(z[i:(i + 3L)] * c(1L, 4L, 16L, 64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

fixture <- function() {
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 1, 2),
    c(0, 0, 1, 1, 2, 2, 1, 0))
  bed <- tempfile(fileext = ".bed")
  write_bed(bed, dosage)
  ids <- paste0("id", seq_len(ncol(dosage)))
  markers <- paste0("m", seq_len(nrow(dosage)))
  glist <- list(
    n = ncol(dosage), ids = ids, bedfiles = bed,
    rsids = list(markers), rsidsLD = list(markers),
    chr = list(rep(1L, length(markers))), pos = list(seq_along(markers) * 100),
    af = list(rowMeans(dosage) / 2),
    maf = list(pmin(rowMeans(dosage) / 2, 1 - rowMeans(dosage) / 2)))
  y <- cbind(T1 = c(-1.2, -.4, .8, .2, 1.1, -.7, .5, -.3))
  rownames(y) <- ids
  stats <- make_summary_stats(glist, y)
  stats$cls <- unname(stats$cls)
  prefix <- tempfile("blr_neutrality_ld_")
  correlation <- stats::cor(t(dosage))
  edges <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edges <- edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
  row_ptr <- c(0, cumsum(tabulate(edges[, 1L], nbins = nrow(dosage))))
  getFromNamespace(".stblr_write_uint64_file", "sblr")(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  getFromNamespace(".stblr_write_uint32_file", "sblr")(
    paste0(prefix, ".col_idx.u32.0based.bin"), edges[, 2L] - 1L)
  writeBin(as.numeric(correlation[edges]), paste0(prefix, ".values.f32.bin"),
           size = 4L, endian = .Platform$endian)
  writeLines(c(sprintf("n_variants=%d", nrow(dosage)),
               sprintf("nnz=%d", nrow(edges))), paste0(prefix, ".meta.txt"))
  metadata <- stats$marker_metadata
  metadata$effect_allele <- "A"; metadata$other_allele <- "C"
  stats$marker_metadata <- metadata
  annotation <- cbind(intercept = 1, coding = c(0, 1, 0))
  rownames(annotation) <- markers
  list(bed = bed, dosage = dosage, Glist = glist, y = y, stats = stats,
       prefix = prefix, annotations = annotation,
       ld_metadata = list(prefix = prefix, marker_ids = markers,
         marker_metadata = metadata, scale = "standardized_genotype",
         source = "make_summary_stats"))
}

cleanup_fixture <- function(x) {
  unlink(x$bed)
  unlink(paste0(x$prefix, c(".row_ptr.u64.bin", ".col_idx.u32.0based.bin",
                           ".values.f32.bin", ".meta.txt")))
}

normalize <- function(x) {
  if (is.environment(x) || is.function(x) || inherits(x, "externalptr")) return(NULL)
  if (is.list(x)) {
    drop <- grepl("(^|_)(time|timing|elapsed|runtime|seconds)(_|$)|^process_id$",
                  names(x), ignore.case = TRUE)
    x <- x[!drop]
    return(lapply(x, normalize))
  }
  if (is.data.frame(x)) {
    drop <- grepl("(^|_)(time|timing|elapsed|runtime|seconds)(_|$)|^process_id$",
                  names(x), ignore.case = TRUE)
    return(x[, !drop, drop = FALSE])
  }
  if (is.character(x)) {
    is_temporary <- grepl("[/\\\\](Temp|tmp)[/\\\\]", x, ignore.case = TRUE) |
      grepl("[.]bed$|[.]row_ptr[.]u64[.]bin$|[.]col_idx[.]u32[.]0based[.]bin$|[.]values[.]f32[.]bin$|[.]meta[.]txt$", x)
    x[is_temporary] <- "<temporary_path>"
  }
  attributes(x) <- attributes(x)[setdiff(names(attributes(x)), c("srcfile", "srcref", "wholeSrcref"))]
  x
}

x <- fixture()
on.exit(cleanup_fixture(x), add = TRUE)
common <- list(updateB = FALSE, updateE = FALSE, updatePi = FALSE,
               nit = 6L, nburn = 2L, nthin = 1L, seed = 22001L,
               convergence = "extended",
               convergence_control = list(warn = FALSE, keep_traces = TRUE,
                 selected_markers = c(2L, 1L),
                 selected_marker_quantities = c("b", "d")))
mix <- list(mixture_var = c(0, .1, 1), pi = c(.7, .2, .1), alpha = c(1, 1, 1))

fits <- list()
fits$st_bed_bayesc <- do.call(stblr_bed, c(list(y = x$y, Glist = x$Glist, method = "bayesc"), common))
fits$st_bed_bayesr <- do.call(stblr_bed, c(list(y = x$y, Glist = x$Glist, method = "bayesr"), common, mix))
fits$st_bed_bayesrc <- do.call(stblr_bed, c(list(y = x$y, Glist = x$Glist, method = "bayesrc",
  annotation = x$annotations, updateAlpha = FALSE, mixture_var = c(0, .1, 1),
  pi = c(.7, .2, .1)), common))
fits$st_csr_sbayesc <- do.call(stblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix, method = "sbayesc"), common))
fits$st_csr_sbayesr <- do.call(stblr_csr, c(list(stats = x$stats, ld_prefix = x$prefix, method = "sbayesr"), common, mix))
fits$st_block_sbayesc <- do.call(stblr_block_eigen, c(list(stats = x$stats, Glist = x$Glist,
  block_start = 1L, eigen_filter = "hard_truncate", eigen_tau = 0, method = "sbayesc"), common))
fits$st_block_sbayesr <- do.call(stblr_block_eigen, c(list(stats = x$stats, Glist = x$Glist,
  block_start = 1L, eigen_filter = "hard_truncate", eigen_tau = 0, method = "sbayesr"), common, mix))
ann_common <- common; ann_common$updatePi <- NULL
fits$st_block_sbayesrc <- do.call(stblr_block_eigen, c(list(stats = x$stats, Glist = x$Glist,
  block_start = 1L, eigen_filter = "hard_truncate", eigen_tau = 0, method = "sbayesrc",
  annotation = x$annotations, updateAlpha = FALSE), ann_common))
fits$st_annot_fixed <- do.call(stblr_csr_annot, c(list(stats = x$stats, Glist = x$Glist,
  ld_prefix = x$prefix, annotations = list(A = x$annotations,
    fixed_pi_marker = list(rep(.3, 3)), fixed_vb_multiplier = list(rep(1, 3))),
  annotation_model = "fixed_marker", use_pi_marker = TRUE, use_vb_multiplier = TRUE), common))
fits$st_annot_group <- do.call(stblr_csr_annot, c(list(stats = x$stats, Glist = x$Glist,
  ld_prefix = x$prefix, annotations = setNames(c("a", "b", "a"), x$stats$marker_names),
  annotation_model = "group", group_names = c("a", "b"), updateGroupVb = FALSE), common))
fits$st_annot_learned <- do.call(stblr_csr_annot, c(list(stats = x$stats, Glist = x$Glist,
  ld_prefix = x$prefix, annotations = x$annotations, annotation_model = "learned_logistic",
  learn_pi_annot = TRUE, learn_vb_annot = TRUE, annot_update_every = 1L), common))
fits$st_annot_sbayesrc <- do.call(stblr_csr_annot, c(list(stats = x$stats, Glist = x$Glist,
  ld_prefix = x$prefix, annotations = x$annotations, annotation_model = "annotation_probit_stick",
  updateAlpha = FALSE), ann_common))

models <- matrix(c(0L, 1L), 2L, 1L)
mt_common <- list(models = models, pimodels = c(.7, .3), vb = matrix(.1), ve = matrix(.5),
  updateB = FALSE, updateE = FALSE, updatePi = FALSE, nit = 6L, nburn = 2L,
  nthin = 1L, seed = 22002L, convergence = "extended",
  convergence_control = list(warn = FALSE, keep_traces = TRUE,
    selected_markers = c(2L, 1L), selected_marker_quantities = c("b", "d")))
mt_mix <- list(mixture_var = c(0, .1, 1))
for (op in c("bed", "csr", "block")) for (kernel in c("bayesc", "bayesr", "bayesrc")) {
  method <- if (op == "bed") kernel else paste0("s", kernel)
  args <- mt_common; args$method <- method
  if (kernel != "bayesc") args <- c(args, mt_mix)
  if (kernel == "bayesrc") {
    args$annotations <- x$annotations; args$add_intercept <- FALSE
    args$standardize_annotations <- FALSE; args$updateAlpha <- FALSE
  }
  if (op == "bed") fit <- do.call(mtblr_bed, c(list(y = x$y, Glist = x$Glist,
    residual_covariance = if (kernel == "bayesc") "full" else "diagonal"), args))
  if (op == "csr") fit <- do.call(mtblr_csr, c(list(stats = x$stats,
    ld_prefix = x$prefix, ld_metadata = x$ld_metadata), args))
  if (op == "block") fit <- do.call(mtblr_block_eigen, c(list(stats = x$stats,
    Glist = x$Glist, block_start = 1L, eigen_filter = "hard_truncate", eigen_tau = 0), args))
  fits[[paste("mt", op, kernel, sep = "_")]] <- fit
}

saveRDS(normalize(fits), output, version = 3L)
cat(sprintf("neutrality_routes=%d output=%s\n", length(fits), normalizePath(output)))
