blr_unified_write_bed <- function(path, dosage) {
  dosage_to_code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    codes <- unname(dosage_to_code[as.character(dosage[marker, ])])
    codes <- c(codes, rep(0L, (-length(codes)) %% 4L))
    vapply(seq(1L, length(codes), by = 4L), function(i) {
      sum(codes[i:(i + 3L)] * c(1L, 4L, 16L, 64L))
    }, integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

blr_unified_fixture <- function() {
  bed <- tempfile(fileext = ".bed")
  dosage <- rbind(
    c(0, 1, 2, 0, 1, 2, 1, 0),
    c(2, 0, 1, 2, 0, 1, 1, 2),
    c(0, 0, 1, 1, 2, 2, 1, 0))
  blr_unified_write_bed(bed, dosage)
  ids <- paste0("id", seq_len(ncol(dosage)))
  markers <- paste0("m", seq_len(nrow(dosage)))
  glist <- list(
    n = ncol(dosage), ids = ids, bedfiles = bed,
    rsids = list(markers), rsidsLD = list(markers),
    chr = list(rep(1L, length(markers))),
    pos = list(seq_along(markers) * 100),
    af = list(rowMeans(dosage) / 2),
    maf = list(pmin(rowMeans(dosage) / 2, 1 - rowMeans(dosage) / 2)))
  y <- cbind(T1 = c(-1.2, -.4, .8, .2, 1.1, -.7, .5, -.3))
  rownames(y) <- ids
  stats <- make_summary_stats(glist, y)
  stats$cls <- unname(stats$cls)
  list(bed = bed, dosage = dosage, Glist = glist, y = y, stats = stats)
}

blr_unified_cleanup <- function(x) {
  unlink(x$bed)
  invisible(NULL)
}

blr_unified_scheduled_csr_fixture <- function(nt = 2L) {
  prefix <- tempfile("blr_unified_scheduled_")
  row_ptr <- c(0, 2, 3, 4, 4, 5, 5)
  col_idx <- c(1L, 2L, 2L, 3L, 5L)
  values <- c(.65, -.25, .45, .55, -.35)
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), col_idx)
  writeBin(values, paste0(prefix, ".values.f32.bin"), size = 4,
           endian = "little")
  writeLines(c(
    "format=sparse_ld_csr", "storage=streamed_upper_triangle",
    "n_bed=NA", "n_used=NA", "n_samples_used=NA", "n_variants=6",
    "nnz=5", "triangle=upper", "diagonal=implicit_1",
    paste0("row_ptr_file=", prefix, ".row_ptr.u64.bin"),
    paste0("col_idx_file=", prefix, ".col_idx.u32.0based.bin"),
    paste0("values_file=", prefix, ".values.f32.bin"),
    "row_ptr_type=uint64", "col_idx_type=uint32", "values_type=float32",
    "index_base=0", "value=r"), paste0(prefix, ".meta.txt"))
  markers <- paste0("m", seq_len(6L))
  trait_names <- paste0("T", seq_len(nt))
  stats <- list(
    wy = setNames(lapply(seq_len(nt), function(tt)
      setNames(c(4, -2, .25, 1.5, -.1, .8) + .05 * (tt - 1L), markers)),
      trait_names),
    ww = setNames(rep(list(setNames(rep(80, 6L), markers)), nt), trait_names),
    yy = setNames(rep(80, nt), trait_names), n = rep(80L, nt), m = 6L,
    marker_names = markers, trait_names = trait_names)
  list(prefix = prefix, stats = stats)
}

blr_unified_cleanup_csr <- function(x) {
  unlink(paste0(x$prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt")))
  invisible(NULL)
}

blr_unified_exact_ld_prefix <- function(dosage) {
  prefix <- tempfile("blr_unified_exact_ld_")
  m <- nrow(dosage)
  correlation <- stats::cor(t(dosage))
  edges <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edges <- edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
  row_ptr <- c(0, cumsum(tabulate(edges[, 1L], nbins = m)))
  sblr:::.stblr_write_uint64_file(
    paste0(prefix, ".row_ptr.u64.bin"), row_ptr)
  sblr:::.stblr_write_uint32_file(
    paste0(prefix, ".col_idx.u32.0based.bin"), edges[, 2L] - 1L)
  writeBin(as.numeric(correlation[edges]), paste0(prefix, ".values.f32.bin"),
           size = 4L, endian = .Platform$endian)
  writeLines(c(sprintf("n_variants=%d", m),
               sprintf("nnz=%d", nrow(edges))), paste0(prefix, ".meta.txt"))
  prefix
}

blr_unified_cleanup_prefix <- function(prefix) {
  unlink(paste0(prefix, c(
    ".row_ptr.u64.bin", ".col_idx.u32.0based.bin", ".values.f32.bin",
    ".meta.txt")))
  invisible(NULL)
}
