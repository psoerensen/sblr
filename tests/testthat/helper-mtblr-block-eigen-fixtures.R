phase17l_write_bed <- function(path, dosage) {
  code <- function(x) ifelse(is.na(x), 1L, c(`0`=3L, `1`=2L, `2`=0L)[as.character(x)])
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(marker) {
    values <- as.integer(code(dosage[marker, ]))
    values <- c(values, rep(0L, (-length(values)) %% 4L))
    vapply(seq(1L, length(values), by=4L), function(i)
      sum(values[i:(i+3L)]*c(1L,4L,16L,64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c,0x1b,0x01,packed)), path)
}

phase17l_descriptor <- function(dosage, blocks=c(0L,2L), filter="hard_truncate",
                                tau=.01, eta=0, af=rep(.35,nrow(dosage))) {
  bed <- tempfile(fileext=".bed"); phase17l_write_bed(bed, dosage)
  list(bed_files=bed, n_bed=as.integer(ncol(dosage)),
       cls=list(seq_len(nrow(dosage))), rows=NULL, af=af,
       block_start=as.integer(blocks), eigen_filter=filter,
       eigen_tau=tau, eigen_eta=eta)
}

phase17l_inspect <- function(descriptor, wy) {
  sblr:::stblr_block_eigen_contract_internal(
    descriptor$bed_files, descriptor$n_bed, descriptor$cls, descriptor$rows,
    descriptor$af, descriptor$block_start, do.call(rbind, wy),
    rep(0, length(wy[[1L]])), descriptor$eigen_filter,
    descriptor$eigen_tau, descriptor$eigen_eta, "")
}

phase17l_dense_operator <- function(inspections) {
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

phase17l_case <- function(nt=2L, filters=rep("hard_truncate",nt), shared=TRUE,
                          blocks=rep(list(c(0L,2L)),nt), updates=FALSE,
                          nonzero=FALSE, multiple_sets=FALSE) {
  dosage <- rbind(c(0,1,2,0,1,2,1,0), c(2,1,0,2,1,0,NA,2),
                  c(0,0,1,1,2,2,1,0), c(2,2,1,1,0,0,1,2))
  descriptors <- lapply(seq_len(if(shared) 1L else nt), function(t)
    phase17l_descriptor(dosage, blocks[[if(shared) 1L else t]],
      filters[[if(shared) 1L else t]], eta=.5, af=c(.25,.35,.4,.3)))
  wy <- lapply(seq_len(nt), function(t) c(1.2,-.6,.5,.9)+.08*(t-1L))
  inspections <- if(shared) {
    one <- phase17l_inspect(descriptors[[1L]], wy)
    lapply(seq_len(nt), function(t) { z <- one; z$transformed_wy <- one$transformed_wy[t,,drop=FALSE]; z })
  } else lapply(seq_len(nt), function(t) phase17l_inspect(descriptors[[t]], list(wy[[t]])))
  op <- phase17l_dense_operator(inspections)
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

phase17l_compare <- function(case, tolerance=1e-12) {
  dense <- do.call(sblr:::mtblr, case$dense)
  block <- do.call(sblr:::mtblr_block_eigen_internal, case$block)
  testthat::expect_equal(block, dense, tolerance=tolerance)
  invisible(block)
}

phase17l_compare_csr <- function(case, tolerance=1e-10) {
  nt <- length(case$transformed); prefixes <- character(nt)
  diagonals <- lapply(case$matrices, diag)
  for (t in seq_len(nt)) {
    prefixes[t] <- tempfile(paste0("phase17l_csr_",t,"_"))
    W <- case$matrices[[t]]
    correlation <- W / outer(sqrt(diag(W)),sqrt(diag(W)))
    diag(correlation) <- 1
    phase17i_write_operator(prefixes[t],correlation,diagonals[[t]])
  }
  args <- case$block
  args$wy <- case$transformed; args$ww <- diagonals
  args$ld_prefixes <- prefixes; args$operator_descriptors <- NULL
  args <- args[c("wy","ww","yy","b","ld_prefixes","sets","B","E",
    "ssb_prior","sse_prior","models","pi","nub","nue","updateB","updateE",
    "updatePi","n","nit","nburn","nthin","seed","method")]
  block <- do.call(sblr:::mtblr_block_eigen_internal,case$block)
  csr <- do.call(sblr:::mtblr_csr_internal,args)
  testthat::expect_equal(block,csr,tolerance=tolerance)
}

phase17m_public_case <- function(nt = 2L, filters = rep("hard_truncate", nt),
                                 blocks = rep(list(c(1L, 3L)), nt),
                                 rows = NULL) {
  dosage <- rbind(c(0,1,2,0,1,2,1,0), c(2,1,0,2,1,0,NA,2),
                  c(0,0,1,1,2,2,1,0), c(2,2,1,1,0,0,1,2))
  bed <- tempfile(fileext = ".bed")
  phase17l_write_bed(bed, dosage)
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

phase17m_call <- function(case, operator_sharing = "auto", ...) {
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
