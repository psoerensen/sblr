mt_bayesc_reference_config <- function(id = 1L) {
  stopifnot(id %in% 1:3)
  marker_names <- paste0("M", 1:4)
  trait_names <- c("TraitA", "TraitB")
  ld1 <- matrix(c(
    1.00, 0.18, 0.00, 0.05,
    0.18, 1.00, 0.12, 0.00,
    0.00, 0.12, 1.00, 0.22,
    0.05, 0.00, 0.22, 1.00), 4, 4, byrow = TRUE,
    dimnames = list(marker_names, marker_names))
  ld2 <- matrix(c(
    1.00, 0.10, 0.04, 0.00,
    0.10, 1.00, 0.16, 0.03,
    0.04, 0.16, 1.00, 0.14,
    0.00, 0.03, 0.14, 1.00), 4, 4, byrow = TRUE,
    dimnames = list(marker_names, marker_names))
  xy <- list(
    TraitA = setNames(c(1.30, -0.65, 0.55, 0.95), marker_names),
    TraitB = setNames(c(0.85, 0.45, -0.80, 1.10), marker_names))
  common <- list(yy = c(TraitA = 56, TraitB = 61), Xy = xy,
    XX = list(TraitA = ld1, TraitB = ld2), n = c(40L, 42L),
    nit = 18L, nburn = 7L, nthin = 2L, verbose = FALSE,
    method = "bayesC", algorithm = "default", pi = 0.001,
    nub = 4, nue = 4)
  if (id == 1L) return(modifyList(common, list(r_seed = 1701L,
    updateB = TRUE, updateE = TRUE, updatePi = TRUE)))
  if (id == 2L) return(modifyList(common, list(r_seed = 1702L,
    vb = matrix(c(0.18, 0.035, 0.035, 0.16), 2, 2),
    ve = matrix(c(0.75, 0.12, 0.12, 0.82), 2, 2),
    pi = 0.20, updateB = FALSE, updateE = FALSE, updatePi = FALSE,
    nit = 14L, nburn = 5L, nthin = 1L)))
  modifyList(common, list(r_seed = 1703L, sets = list(c(1L, 3L), c(2L, 4L)),
    models = list(c(0L, 0L), c(1L, 0L), c(0L, 1L), c(1L, 1L)),
    pimodels = c(0.55, 0.15, 0.15, 0.15),
    updateB = TRUE, updateE = TRUE, updatePi = TRUE,
    nit = 20L, nburn = 8L, nthin = 2L))
}

mt_bayesc_reference_public_args <- function(config) {
  config$r_seed <- NULL
  config
}

mt_bayesc_reference_native_arguments <- function(config) {
  seed_state <- config$r_seed
  args <- mt_bayesc_reference_public_args(config)
  yy <- args$yy; Xy <- args$Xy; XX <- args$XX; n <- args$n
  wy <- lapply(Xy, as.vector)
  ww <- lapply(XX, diag)
  m <- mean(lengths(wy)); nt <- length(wy)
  b <- args[["b"]]
  if (is.null(b)) b <- lapply(seq_len(nt), function(x) rep(0, m))
  if (is.matrix(b)) b <- split(b, rep(seq_len(ncol(b)), each = nrow(b)))
  models <- args[["models"]]
  if (is.null(models)) {
    models <- rep(list(0:1), nt)
    models <- t(do.call(expand.grid, models))
    models <- split(models, rep(seq_len(ncol(models)), each = nrow(models)))
  }
  pi <- if (is.null(args$pi)) 0.001 else args$pi
  pimodels <- args[["pimodels"]]
  if (is.null(pimodels)) pimodels <- c(1 - pi,
    rep(pi / (length(models) - 1), length(models) - 1))
  vy <- diag(yy / (n - 1), nt)
  h2 <- if (is.null(args[["h2"]])) 0.5 else args[["h2"]]
  vg <- args[["vg"]]; if (is.null(vg)) vg <- diag(diag(vy) * h2)
  ve <- args[["ve"]]; if (is.null(ve)) ve <- diag(diag(vy) * (1 - h2))
  vb <- args[["vb"]]; if (is.null(vb)) vb <- diag((diag(vy) * h2) / (m * pi))
  nub <- args$nub; nue <- args$nue
  ssb <- args[["ssb_prior"]]
  if (is.null(ssb)) ssb <- diag(((nub - 2) / nub) * (diag(vg) / (m * pi)))
  sse <- args[["sse_prior"]]
  if (is.null(sse)) sse <- diag(((nue - 2) / nue) * diag(ve))
  sets <- args[["sets"]]; if (is.null(sets)) sets <- list(seq_len(m))
  sets <- lapply(sets, function(x) x - 1L)
  XXvalues <- lapply(XX, function(x) as.list(as.data.frame(x)))
  XXindices <- lapply(seq_len(m), function(x) seq_len(m) - 1L)
  set.seed(seed_state)
  native_seed <- sample.int(.Machine$integer.max, 1)
  list(wy = wy, ww = ww, yy = yy, b = b, XXvalues = XXvalues,
    XXindices = XXindices, sets = sets, B = vb, E = ve,
    ssb_prior = split(ssb, rep(seq_len(ncol(ssb)), each = nrow(ssb))),
    sse_prior = split(sse, rep(seq_len(ncol(sse)), each = nrow(sse))),
    models = models, pi = pimodels, nub = nub, nue = nue,
    updateB = args$updateB, updateE = args$updateE,
    updatePi = args$updatePi, n = n, nit = args$nit, nburn = args$nburn,
    nthin = args$nthin, seed = native_seed, method = 4L)
}

mt_bayesc_reference_native_raw <- function(config) {
  do.call(sblr:::mtblr, mt_bayesc_reference_native_arguments(config))
}

mt_bayesc_reference_capture <- function(id = 1L, formatted = TRUE) {
  config <- mt_bayesc_reference_config(id)
  if (!formatted) return(mt_bayesc_reference_native_raw(config))
  set.seed(config$r_seed)
  do.call(sblr::sblr, mt_bayesc_reference_public_args(config))
}

mt_bayesc_reference_metadata <- function(id = 1L) {
  x <- mt_bayesc_reference_config(id)
  list(starting_commit = "dc429e1",
    reference_mode = "structure_exact_numeric_tolerance",
    numeric_tolerance = 1e-12, structure_exact = TRUE,
    schema = "legacy positional native / named public fit",
    rng = "R-generated seed; one fit-local std::mt19937",
    samples = x$n, markers = lengths(x$Xy)[1], traits = length(x$Xy),
    marker_names = names(x$Xy[[1]]), trait_names = names(x$Xy),
    r_seed = x$r_seed, iterations = x$nit, burnin = x$nburn,
    thinning = x$nthin, updateB = x$updateB, updateE = x$updateE,
    updatePi = x$updatePi, models = x$models, sets = x$sets)
}
