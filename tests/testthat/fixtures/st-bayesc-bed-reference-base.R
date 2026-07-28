st_bayesc_bed_base_write_bed <- function(path, dosage) {
  code <- c(`0` = 3L, `1` = 2L, `2` = 0L)
  packed <- unlist(lapply(seq_len(nrow(dosage)), function(j) {
    z <- unname(code[as.character(dosage[j, ])])
    z <- c(z, rep(0L, (-length(z)) %% 4L))
    vapply(seq(1L, length(z), 4L), function(i)
      sum(z[i:(i + 3L)] * c(1L, 4L, 16L, 64L)), integer(1))
  }))
  writeBin(as.raw(c(0x6c, 0x1b, 0x01, packed)), path)
}

st_bayesc_bed_base_fixture <- function() {
  dosage <- rbind(c(0, 1, 2, 0, 1, 2), c(2, 0, 1, 2, 0, 1))
  bed <- tempfile(fileext = ".bed")
  st_bayesc_bed_base_write_bed(bed, dosage)
  list(Glist = list(n = 6L, ids = paste0("id", 1:6), bedfiles = bed,
    rsids = list(c("rs1", "rs2")), rsidsLD = list(c("rs1", "rs2")),
    chr = list(c(1L, 1L)), pos = list(c(100, 200)),
    af = list(rowMeans(dosage) / 2)),
    y = matrix(c(-1, 0, 1, -.5, .5, 1.5), ncol = 1,
      dimnames = list(NULL, "D1")))
}

st_bayesc_bed_base_capture <- function(model, ncores = 1L, nchains = 1L, seed = 71L) {
  x <- st_bayesc_bed_base_fixture()
  ns <- asNamespace("sblr")
  assign(".st_bayesc_bed_base_raw", NULL, envir = .GlobalEnv)
  suppressMessages(trace(".as_stblr_fit",
    tracer = quote(assign(".st_bayesc_bed_base_raw", raw, envir = .GlobalEnv)),
    where = ns, print = FALSE))
  on.exit(suppressMessages(try(untrace(".as_stblr_fit", where = ns), silent = TRUE)))
  args <- list(y = x$y, Glist = x$Glist, method = model, nit = 6L,
    nburn = 2L, nthin = 1L, seed = seed, ncores = ncores,
    nchains = nchains, updateB = FALSE, updateE = FALSE,
    rebuild_every = 2L, read_block_size = 2L, convergence = "none")
  if (model == "bayesc") args <- c(args, list(pi_init = .5,
    pi_prior_mean = .5, pi_prior_strength = 4))
  if (model == "bayesr") args <- c(args, list(updatePi = FALSE,
    mixture_var = c(0, .01, .1, 1), pi = c(.95, .03, .015, .005)))
  if (model == "bayesrc") args <- c(args, list(annotation = matrix(c(0, 1),
    ncol = 1, dimnames = list(c("rs1", "rs2"), "annot1")),
    updateAlpha = FALSE, mixture_var = c(0, .01, .1, 1),
    pi = c(.95, .03, .015, .005)))
  fit <- do.call(sblr::stblr_bed, args)
  list(raw = get(".st_bayesc_bed_base_raw", envir = .GlobalEnv), fit = fit)
}

st_bayesc_bed_base_normalize <- function(x) {
  if (is.list(x) && all(c("raw", "fit") %in% names(x))) {
    x$raw <- st_bayesc_bed_base_normalize(x$raw)
    x$fit <- st_bayesc_bed_base_normalize(x$fit)
  }
  if (is.list(x) && !is.null(x$input)) x$input$ncores <- 0L
  if (is.list(x) && "memory_estimate" %in% names(x)) x$memory_estimate <- NULL
  if (is.list(x)) {
    x <- lapply(x, st_bayesc_bed_base_normalize)
    for (nm in intersect(names(x), c("seconds", "seconds_mean", "seconds_max"))) x[[nm]][] <- 0
  }
  x
}
