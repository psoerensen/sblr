# ---- consolidated from tests/testthat/helper-blr-unified.R ----
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

# ---- consolidated from tests/testthat/helper-mtblr-csr-fixtures.R ----
blr_csr_fixture_src <- function(file) paste(readLines(
  blr_repo_path(file), warn = FALSE), collapse = "\n")

blr_csr_fixture_float <- function(x) {
  con <- rawConnection(raw(), "w+b"); on.exit(close(con))
  writeBin(as.numeric(x), con, size = 4L, endian = .Platform$endian)
  seek(con, 0L)
  readBin(con, "numeric", n = length(x), size = 4L,
          endian = .Platform$endian)
}

blr_csr_fixture_write_operator <- function(prefix, correlation, diagonal) {
  m <- nrow(correlation)
  edges <- which(upper.tri(correlation) & correlation != 0, arr.ind = TRUE)
  edges <- edges[order(edges[, 1L], edges[, 2L]), , drop = FALSE]
  row_ptr <- c(0, cumsum(tabulate(edges[, 1L], nbins = m)))
  writeBin(as.integer(c(rbind(row_ptr, rep(0, length(row_ptr))))),
           paste0(prefix, ".row_ptr.u64.bin"), size = 4L,
           endian = .Platform$endian)
  writeBin(as.integer(edges[, 2L] - 1L),
           paste0(prefix, ".col_idx.u32.0based.bin"), size = 4L,
           endian = .Platform$endian)
  writeBin(as.numeric(correlation[edges]), paste0(prefix, ".values.f32.bin"),
           size = 4L, endian = .Platform$endian)
  writeLines(c(sprintf("n_variants=%d", m), sprintf("nnz=%d", nrow(edges))),
             paste0(prefix, ".meta.txt"))
  indices <- vector("list", m); values <- vector("list", m)
  for (k in seq_len(nrow(edges))) {
    i <- edges[k, 1L]; j <- edges[k, 2L]
    x <- blr_csr_fixture_float(correlation[i, j] * sqrt(diagonal[i] * diagonal[j]))
    indices[[i]] <- c(indices[[i]], j - 1L); values[[i]] <- c(values[[i]], x)
    indices[[j]] <- c(indices[[j]], i - 1L); values[[j]] <- c(values[[j]], x)
  }
  list(indices = indices, values = values)
}

blr_csr_fixture_case <- function(nt = 2L, trait_specific = FALSE,
                          independent = FALSE, nonzero_b = FALSE,
                          multiple_sets = FALSE, updates = TRUE) {
  m <- 4L
  base <- matrix(c(1, .18, 0, .05, .18, 1, .12, 0,
                   0, .12, 1, .22, .05, 0, .22, 1), m, m, byrow = TRUE)
  matrices <- rep(list(base), nt)
  if (trait_specific) {
    for (t in seq_len(nt)) {
      matrices[[t]][upper.tri(matrices[[t]])] <-
        matrices[[t]][upper.tri(matrices[[t]])] * (1 - .12 * (t - 1L))
      matrices[[t]][lower.tri(matrices[[t]])] <-
        t(matrices[[t]])[lower.tri(matrices[[t]])]
    }
  }
  if (independent && nt >= 2L) {
    matrices[[2L]][1, 2] <- matrices[[2L]][2, 1] <- 0
    matrices[[2L]][1, 3] <- matrices[[2L]][3, 1] <- .09
  }
  ww <- lapply(seq_len(nt), function(t)
    rep(1 + if (trait_specific) .1 * (t - 1) else 0, m))
  prefix <- file.path(tempdir(), paste0("blr_csr_fixture_", Sys.getpid(), "_",
                                        sample.int(1e8, 1)))
  prefixes <- paste0(prefix, "_", seq_len(if (trait_specific || independent) nt else 1L))
  operators <- vector("list", nt)
  if (length(prefixes) == 1L) {
    operators <- rep(list(blr_csr_fixture_write_operator(prefixes, matrices[[1]], ww[[1]])), nt)
  } else {
    for (t in seq_len(nt)) operators[[t]] <-
      blr_csr_fixture_write_operator(prefixes[t], matrices[[t]], ww[[t]])
  }
  union_indices <- lapply(seq_len(m), function(i)
    unique(unlist(lapply(operators, function(op) op$indices[[i]]))))
  XXindices <- lapply(seq_len(m), function(i) c(i - 1L, union_indices[[i]]))
  XXvalues <- lapply(seq_len(nt), function(t) lapply(seq_len(m), function(i) {
    positions <- match(union_indices[[i]], operators[[t]]$indices[[i]])
    c(ww[[t]][i], ifelse(is.na(positions), 0,
                         operators[[t]]$values[[i]][positions]))
  }))
  wy <- lapply(seq_len(nt), function(t) c(1.2, -.6, .5, .9) + .08 * (t - 1))
  models <- as.matrix(expand.grid(rep(list(0:1), nt)))
  models <- lapply(seq_len(nrow(models)), function(i) as.integer(models[i, ]))
  pi <- c(.7, rep(.3 / (length(models) - 1), length(models) - 1))
  b <- lapply(seq_len(nt), function(t)
    if (nonzero_b) c(.04 * t, 0, -.02, 0) else rep(0, m))
  sets <- if (multiple_sets) list(c(0L, 2L), c(1L, 3L)) else list(0:(m - 1L))
  mat_list <- function(x) split(x, rep(seq_len(ncol(x)), each = nrow(x)))
  common <- list(wy=wy, ww=ww, yy=rep(50, nt), b=b, sets=sets,
    B=diag(.15, nt), E=diag(.8, nt), ssb_prior=mat_list(diag(.05, nt)),
    sse_prior=mat_list(diag(.3, nt)), models=models, pi=pi, nub=4, nue=4,
    updateB=updates, updateE=updates, updatePi=updates,
    n=as.integer(40 + seq_len(nt)), nit=8L, nburn=3L, nthin=2L,
    seed=17001L, method=4L)
  dense <- c(common, list(XXvalues=XXvalues, XXindices=XXindices))
  dense <- dense[c("wy","ww","yy","b","XXvalues","XXindices","sets","B","E",
    "ssb_prior","sse_prior","models","pi","nub","nue","updateB","updateE",
    "updatePi","n","nit","nburn","nthin","seed","method")]
  csr <- c(common, list(ld_prefixes=prefixes))
  csr <- csr[c("wy","ww","yy","b","ld_prefixes","sets","B","E","ssb_prior",
    "sse_prior","models","pi","nub","nue","updateB","updateE","updatePi",
    "n","nit","nburn","nthin","seed","method")]
  list(dense=dense, csr=csr, prefixes=prefixes)
}

blr_csr_public_public_case <- function(trait_specific = FALSE, independent = FALSE) {
  x <- blr_csr_fixture_case(trait_specific = trait_specific, independent = independent,
                     updates = FALSE)
  nt <- length(x$csr$wy); ids <- paste0("m", seq_along(x$csr$wy[[1L]]))
  alleles <- data.frame(marker_id = ids, effect_allele = rep("A", length(ids)),
                        other_allele = rep("C", length(ids)))
  stats <- list(wy = setNames(lapply(x$csr$wy, function(z) setNames(z, ids)),
                               paste0("T", seq_len(nt))),
    ww = setNames(lapply(x$csr$ww, function(z) setNames(z, ids)),
                  paste0("T", seq_len(nt))),
    yy = setNames(x$csr$yy, paste0("T", seq_len(nt))), n = x$csr$n,
    marker_names = ids, trait_names = paste0("T", seq_len(nt)),
    marker_metadata = alleles, scale = "standardized_genotype", source = "external")
  descriptors <- lapply(rep(x$prefixes, length.out = nt), function(prefix)
    list(prefix = prefix, marker_ids = ids, marker_metadata = alleles,
         scale = "standardized_genotype", source = "external"))
  list(x = x, stats = stats,
       metadata = if (length(x$prefixes) == 1L) descriptors[[1L]] else descriptors)
}

blr_csr_public_call <- function(case, marker_policy = "strict", ...) {
  x <- case$x$csr
  args <- list(stats = case$stats, ld_prefix = x$ld_prefixes,
    ld_metadata = case$metadata, marker_policy = marker_policy,
    b = do.call(cbind, x$b), vb = x$B, ve = x$E,
    ssb_prior = do.call(cbind, x$ssb_prior),
    sse_prior = do.call(cbind, x$sse_prior), models = do.call(rbind, x$models),
    pimodels = x$pi, sets = lapply(x$sets, `+`, 1L), updateB = x$updateB,
    updateE = x$updateE, updatePi = x$updatePi, nub = x$nub, nue = x$nue,
    nit = x$nit, nburn = x$nburn, nthin = x$nthin, seed = x$seed)
  overrides <- list(...); args[names(overrides)] <- overrides
  do.call(mtblr_csr, args)
}

# ---- consolidated from tests/testthat/helper-mtblr-block-eigen-fixtures.R ----
blr_block_fixture_write_bed <- function(path, dosage) {
  code <- function(x) ifelse(is.na(x), 1L, c(`0`=3L, `1`=2L, `2`=0L)[as.character(x)])
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    values <- as.integer(code(dosage[marker, ]))
    values <- c(values, rep(0L, (-length(values)) %% 4L))
    vapply(seq(1L, length(values), by=4L), function(i)
      sum(values[i:(i+3L)]*c(1L,4L,16L,64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c,0x1b,0x01,packed)), path)
}

blr_block_fixture_descriptor <- function(dosage, blocks=c(0L,2L), filter="hard_truncate",
                                tau=.01, eta=0, af=rep(.35,nrow(dosage))) {
  bed <- tempfile(fileext=".bed"); blr_block_fixture_write_bed(bed, dosage)
  list(bed_files=bed, n_bed=as.integer(ncol(dosage)),
       cls=list(seq_len(nrow(dosage))), rows=NULL, af=af,
       block_start=as.integer(blocks), eigen_filter=filter,
       eigen_tau=tau, eigen_eta=eta)
}

blr_block_fixture_inspect <- function(descriptor, wy) {
  sblr:::stblr_block_eigen_contract_internal(
    descriptor$bed_files, descriptor$n_bed, descriptor$cls, descriptor$rows,
    descriptor$af, descriptor$block_start, do.call(rbind, wy),
    rep(0, length(wy[[1L]])), descriptor$eigen_filter,
    descriptor$eigen_tau, descriptor$eigen_eta, "")
}

blr_block_fixture_dense_operator <- function(inspections) {
  nt <- length(inspections); m <- inspections[[1L]]$marker_count
  matrices <- lapply(inspections, function(x) {
    W <- matrix(0, m, m)
    for (g in seq_along(x$block_start)) {
      start <- x$block_start[g]+1L; size <- x$block_size[g]
      p <- x$packed_upper_triangle[[g]]; k <- 1L
      for (i in seq_len(size)) for (j in i:size) {
        W[start+i-1L,start+j-1L] <- W[start+j-1L,start+i-1L] <- p[k]; k <- k+1L
      }
    }
    W
  })
  neighbors <- lapply(seq_len(m), function(i)
    sort(unique(unlist(lapply(matrices, function(W) which(W[i,]!=0L))))))
  indices <- lapply(seq_len(m), function(i) c(i-1L, setdiff(neighbors[[i]],i)-1L))
  values <- lapply(seq_len(nt), function(t) lapply(seq_len(m), function(i)
    matrices[[t]][i, indices[[i]]+1L]))
  list(ww=lapply(matrices, diag), XXindices=indices, XXvalues=values,
       matrices=matrices)
}

blr_block_fixture_case <- function(nt=2L, filters=rep("hard_truncate",nt), shared=TRUE,
                          blocks=rep(list(c(0L,2L)),nt), updates=FALSE,
                          nonzero=FALSE, multiple_sets=FALSE) {
  dosage <- rbind(c(0,1,2,0,1,2,1,0), c(2,1,0,2,1,0,NA,2),
                  c(0,0,1,1,2,2,1,0), c(2,2,1,1,0,0,1,2))
  descriptors <- lapply(seq_len(if(shared) 1L else nt), function(t)
    blr_block_fixture_descriptor(dosage, blocks[[if(shared) 1L else t]],
      filters[[if(shared) 1L else t]], eta=.5, af=c(.25,.35,.4,.3)))
  wy <- lapply(seq_len(nt), function(t) c(1.2,-.6,.5,.9)+.08*(t-1L))
  inspections <- if(shared) {
    one <- blr_block_fixture_inspect(descriptors[[1L]], wy)
    lapply(seq_len(nt), function(t) { z <- one; z$transformed_wy <- one$transformed_wy[t,,drop=FALSE]; z })
  } else lapply(seq_len(nt), function(t) blr_block_fixture_inspect(descriptors[[t]], list(wy[[t]])))
  op <- blr_block_fixture_dense_operator(inspections)
  transformed <- lapply(inspections, function(x) as.numeric(x$transformed_wy[1,]))
  models <- as.matrix(expand.grid(rep(list(0:1),nt)))
  models <- lapply(seq_len(nrow(models)), function(i) as.integer(models[i,]))
  pi <- c(.7, rep(.3/(length(models)-1),length(models)-1))
  b <- lapply(seq_len(nt), function(t) if(nonzero) c(.04*t,0,-.02,0) else rep(0,4))
  sets <- if(multiple_sets) list(c(0L,2L),c(1L,3L)) else list(0:3)
  mat_list <- function(x) split(x, rep(seq_len(ncol(x)),each=nrow(x)))
  common <- list(yy=rep(50,nt), b=b, sets=sets, B=diag(.15,nt), E=diag(.8,nt),
    ssb_prior=mat_list(diag(.05,nt)), sse_prior=mat_list(diag(.3,nt)),
    models=models, pi=pi, nub=4, nue=4, updateB=updates, updateE=updates,
    updatePi=updates, n=as.integer(40+seq_len(nt)), nit=8L, nburn=3L,
    nthin=2L, seed=17012L, method=4L)
  dense <- c(list(wy=transformed,ww=op$ww),common[c("yy","b")],
    list(XXvalues=op$XXvalues,XXindices=op$XXindices),common[setdiff(names(common),c("yy","b"))])
  block <- c(list(wy=wy),common[c("yy","b")],list(operator_descriptors=descriptors),
             common[setdiff(names(common),c("yy","b"))])
  list(dense=dense, block=block, inspections=inspections, matrices=op$matrices,
       transformed=transformed)
}

blr_block_public_public_case <- function(nt = 2L, filters = rep("hard_truncate", nt),
                                 blocks = rep(list(c(1L, 3L)), nt),
                                 rows = NULL) {
  dosage <- rbind(c(0,1,2,0,1,2,1,0), c(2,1,0,2,1,0,NA,2),
                  c(0,0,1,1,2,2,1,0), c(2,2,1,1,0,0,1,2))
  bed <- tempfile(fileext = ".bed")
  blr_block_fixture_write_bed(bed, dosage)
  if (is.null(rows)) rows <- seq_len(ncol(dosage))
  ids <- paste0("m", seq_len(nrow(dosage)))
  af <- c(.25, .35, .4, .3)
  metadata <- data.frame(
    marker_id = ids, chromosome_or_file = 1L,
    bed_column = seq_along(ids), allele_frequency = af)
  wy <- setNames(lapply(seq_len(nt), function(t)
    setNames(c(1.2,-.6,.5,.9) + .08 * (t - 1L), ids)),
    paste0("T", seq_len(nt)))
  ww <- setNames(lapply(seq_len(nt), function(t)
    setNames(rep(length(rows), length(ids)), ids)), names(wy))
  stats <- list(
    wy = wy, ww = ww, yy = setNames(rep(50, nt), names(wy)),
    n = rep(length(rows), nt), n_bed = ncol(dosage), m = length(ids),
    bed_files = bed, cls = list(seq_along(ids)), af = list(af),
    rows = rows, marker_names = ids, trait_names = names(wy),
    marker_metadata = metadata, scale = "standardized_genotype",
    source = "make_summary_stats")
  glist <- list(
    bedfiles = bed, n = ncol(dosage), ids = paste0("i", seq_len(ncol(dosage))),
    rsids = list(ids), rsidsLD = list(ids), af = list(af))
  list(stats = stats, Glist = glist, blocks = blocks, filters = filters)
}

blr_block_public_call <- function(case, operator_sharing = "auto", ...) {
  args <- list(
    stats = case$stats, Glist = case$Glist,
    block_start = case$blocks, eigen_filter = case$filters,
    eigen_tau = .01, eigen_eta = .5,
    operator_sharing = operator_sharing,
    updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 8L, nburn = 3L, nthin = 2L, seed = 17013L)
  overrides <- list(...)
  args[names(overrides)] <- overrides
  do.call(mtblr_block_eigen, args)
}

# ---- consolidated from tests/testthat/helper-mtblr-bayesrc.R ----
.mt_bayesrc_fixture <- function() {
  fixture <- blr_unified_fixture()
  prefix <- blr_unified_exact_ld_prefix(fixture$dosage)
  stats <- fixture$stats
  metadata <- stats$marker_metadata
  metadata$effect_allele <- "A"
  metadata$other_allele <- "C"
  stats$marker_metadata <- metadata
  annotation <- cbind(intercept = 1, coding = c(0, 1, 0))
  rownames(annotation) <- stats$marker_names
  list(
    fixture = fixture,
    prefix = prefix,
    stats = stats,
    ld_metadata = list(
      prefix = prefix,
      marker_ids = stats$marker_names,
      marker_metadata = metadata,
      scale = "standardized_genotype",
      source = "make_summary_stats"
    ),
    annotations = annotation
  )
}

.mt_bayesrc_cleanup <- function(x) {
  blr_unified_cleanup(x$fixture)
  blr_unified_cleanup_prefix(x$prefix)
}

.mt_bayesrc_common <- function(method = "sbayesrc", alpha_init = NULL,
                               maf_effect_s = NULL) {
  list(
    method = method,
    annotations = NULL,
    add_intercept = FALSE,
    standardize_annotations = FALSE,
    mixture_var = c(0, .1, 1),
    models = matrix(c(0L, 1L), 2L, 1L),
    alpha_init = alpha_init,
    maf_effect_s = maf_effect_s,
    vb = matrix(.1),
    ve = matrix(.5),
    updateB = FALSE,
    updateE = FALSE,
    updatePi = FALSE,
    updateAlpha = FALSE,
    nit = 8L,
    nburn = 2L,
    nthin = 1L,
    seed = 42L,
    convergence = "none"
  )
}

