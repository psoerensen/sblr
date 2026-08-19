# Study 06 final estimability and annotation-contrast analysis.
# This script consumes frozen retained chains; it never runs a sampler.

study06_pi <- function(eta) {
  stopifnot(is.matrix(eta), ncol(eta) == 3L, all(is.finite(eta)))
  continuation <- stats::pnorm(eta)
  remaining <- rep(1, nrow(eta))
  out <- matrix(0, nrow(eta), 4L)
  for (stick in seq_len(3L)) {
    out[, stick] <- remaining * (1 - continuation[, stick])
    remaining <- remaining * continuation[, stick]
  }
  out[, 4L] <- remaining
  colnames(out) <- paste0("pi", 0:3)
  out
}

study06_diag <- function(chains) {
  n <- min(vapply(chains, length, integer(1)))
  mat <- do.call(cbind, lapply(chains, function(x) as.numeric(x)[seq_len(n)]))
  draws <- posterior::as_draws_array(array(mat, dim = c(n, ncol(mat), 1L)))
  pooled <- as.numeric(mat)
  sd_value <- stats::sd(pooled)
  mcse <- posterior::mcse_mean(draws)
  c(rhat = posterior::rhat(draws), bulk_ess = posterior::ess_bulk(draws),
    tail_ess = posterior::ess_tail(draws), mcse = mcse,
    relative_mcse = if (is.finite(sd_value) && sd_value > 0) mcse / sd_value else NA_real_)
}

study06_pair_metrics <- function(x, y, top = c(25L, 50L, 100L)) {
  out <- c(pearson = stats::cor(x, y), spearman = stats::cor(x, y,
    method = "spearman"), rmse = sqrt(mean((x - y)^2)),
    mae = mean(abs(x - y)), max_abs = max(abs(x - y)))
  for (k in top) {
    ix <- head(order(x, decreasing = TRUE), min(k, length(x)))
    iy <- head(order(y, decreasing = TRUE), min(k, length(y)))
    out[paste0("top", k, "_overlap")] <- length(intersect(ix, iy)) / min(k, length(x))
  }
  out
}

study06_auc <- function(score, truth) {
  truth <- as.logical(truth)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - sum(seq_len(sum(truth)))) / (sum(truth) * sum(!truth))
}

study06_auprc <- function(score, truth) {
  truth <- as.logical(truth)
  ord <- order(score, decreasing = TRUE)
  hit <- truth[ord]
  recall <- cumsum(hit) / sum(hit)
  precision <- cumsum(hit) / seq_along(hit)
  sum((recall - c(0, head(recall, -1L))) * precision)
}

study06_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

study06_chain_alpha <- function(fit) {
  convergence <- sblr::extract_diagnostics(fit, "convergence")
  quantities <- convergence$traces$quantities
  descriptor <- quantities[quantities$parameter_name == "alpha" &
    quantities$annotation_index > 0L & quantities$stick_index > 0L, , drop = FALSE]
  descriptor <- descriptor[order(descriptor$quantity_index), , drop = FALSE]
  values <- convergence$traces$values[, , descriptor$quantity_index,
    drop = FALSE]
  traces <- lapply(seq_len(dim(values)[2L]), function(chain) {
    x <- matrix(values[, chain, , drop = FALSE], nrow = dim(values)[1L],
      ncol = length(descriptor$quantity_index))
    colnames(x) <- paste(descriptor$annotation_name, descriptor$stick_name, sep = "::")
    x
  })
  list(traces = traces, descriptor = descriptor)
}

study06_chain_posterior <- function(fit) {
  effects <- sblr::extract_posterior(fit, "realised_effects",
    state = "retained")
  assignments <- sblr::extract_posterior(fit, "component_assignments",
    state = "retained")
  nchains <- dim(effects)[2L]
  lapply(seq_len(nchains), function(chain) {
    summary <- sblr::summarise_posterior(fit, quantity = "realised_effects",
      chains = chain)
    beta <- summary$mean
    names(beta) <- summary$marker
    state <- assignments[, chain, , drop = FALSE]
    dim(state) <- c(dim(assignments)[1L], dim(assignments)[3L])
    pip <- colMeans(state > 0)
    names(pip) <- dimnames(assignments)[[3L]]
    list(beta = beta, pip = pip)
  })
}

study06_raw_alpha <- function(route, alpha, truth) {
  desc <- alpha$descriptor
  rows <- vector("list", nrow(desc))
  chain_rows <- list()
  for (j in seq_len(nrow(desc))) {
    values <- lapply(alpha$traces, function(x) x[, j])
    pooled <- unlist(values, use.names = FALSE)
    target <- truth[desc$annotation_index[j], desc$stick_index[j]]
    dg <- study06_diag(values)
    rows[[j]] <- data.frame(route, stick = desc$stick_name[j],
      annotation = desc$annotation_name[j], posterior_mean = mean(pooled),
      posterior_sd = stats::sd(pooled), median = stats::median(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      truth = target, bias = mean(pooled) - target,
      absolute_error = abs(mean(pooled) - target), t(dg), check.names = FALSE)
    chain_rows[[j]] <- do.call(rbind, lapply(seq_along(values), function(ch) {
      z <- values[[ch]]
      data.frame(route, chain = ch, stick = desc$stick_name[j],
        annotation = desc$annotation_name[j], mean = mean(z), sd = stats::sd(z),
        q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975))
    }))
  }
  list(summary = do.call(rbind, rows), chains = do.call(rbind, chain_rows))
}

study06_design_dependence <- function(route, A, alpha) {
  design <- data.frame(route, statistic = c("rank", "condition_number",
      paste0("singular_value_", seq_len(ncol(A)))),
    value = c(qr(A)$rank, kappa(A), svd(A, nu = 0, nv = 0)$d))
  correlation <- as.data.frame(as.table(stats::cor(A)))
  names(correlation) <- c("annotation_1", "annotation_2", "correlation")
  correlation$route <- route
  covariance <- correlation_rows <- eigen_rows <- list()
  q <- alpha$descriptor
  for (stick in seq_len(3L)) {
    idx <- which(q$stick_index == stick)
    mats <- lapply(alpha$traces, function(x) x[, idx, drop = FALSE])
    names_here <- q$annotation_name[idx]
    for (ch in seq_along(mats)) {
      cv <- stats::cov(mats[[ch]]); cr <- stats::cor(mats[[ch]])
      covariance[[length(covariance) + 1L]] <- transform(as.data.frame(as.table(cv)),
        route = route, chain = ch, stick = stick, matrix = "covariance")
      correlation_rows[[length(correlation_rows) + 1L]] <- transform(as.data.frame(as.table(cr)),
        route = route, chain = ch, stick = stick, matrix = "correlation")
    }
    pooled <- do.call(rbind, mats); cv <- stats::cov(pooled); cr <- stats::cor(pooled)
    covariance[[length(covariance) + 1L]] <- transform(as.data.frame(as.table(cv)),
      route = route, chain = 0L, stick = stick, matrix = "covariance")
    correlation_rows[[length(correlation_rows) + 1L]] <- transform(as.data.frame(as.table(cr)),
      route = route, chain = 0L, stick = stick, matrix = "correlation")
    eg <- eigen(cv, symmetric = TRUE)
    for (direction in seq_along(eg$values)) {
      for (a in seq_along(names_here)) eigen_rows[[length(eigen_rows) + 1L]] <-
        data.frame(route, stick, direction, eigenvalue = eg$values[direction],
          annotation = names_here[a], loading = eg$vectors[a, direction])
    }
  }
  list(design = design, design_correlation = correlation,
    posterior_matrices = rbind(do.call(rbind, covariance), do.call(rbind, correlation_rows)),
    eigen = do.call(rbind, eigen_rows))
}

study06_representatives <- function(A) {
  target <- expand.grid(enriched_binary = 0:1,
    continuous_signal = c(-1, 0, 1), null_annotation = 0)
  unique(vapply(seq_len(nrow(target)), function(i) {
    distance <- rowSums((A[, c("enriched_binary", "continuous_signal", "null_annotation")] -
      as.numeric(target[i, ]))^2)
    which.min(distance)
  }, integer(1)))
}

study06_derive_route <- function(route, fit, alpha, A, truth_alpha, truth_pi,
                                  chunk_size = 500L) {
  m <- nrow(A); q <- alpha$descriptor; representative <- study06_representatives(A)
  mods <- list(enriched_binary = c(0, 1), continuous_signal = c(-1, 1),
    null_annotation = c(-1, 1))
  marker <- vector("list", length(alpha$traces)); scalar <- vector("list", length(alpha$traces))
  contrast <- vector("list", length(alpha$traces)); observed <- vector("list", length(alpha$traces))
  normalization <- list()
  for (ch in seq_along(alpha$traces)) {
    message(sprintf("%s: deriving chain %d/%d", route, ch, length(alpha$traces)))
    tr <- alpha$traces[[ch]]; nd <- nrow(tr)
    eta_sum <- matrix(0, m, 3L); eta_sq <- matrix(0, m, 3L)
    pi_sum <- matrix(0, m, 4L); pi_sq <- matrix(0, m, 4L)
    eta_selected <- array(NA_real_, c(nd, length(representative), 3L))
    contrast_draw <- array(NA_real_, c(nd, length(mods), 5L),
      dimnames = list(NULL, names(mods), c("active", paste0("pi", 0:3))))
    observed_draw <- matrix(NA_real_, nd, 5L,
      dimnames = list(NULL, c("active", paste0("pi", 0:3))))
    max_row_error <- 0; min_probability <- 1; max_probability <- 0
    for (lo in seq.int(1L, nd, by = chunk_size)) {
      ix <- lo:min(nd, lo + chunk_size - 1L); nc <- length(ix)
      eta <- array(NA_real_, c(m, nc, 3L)); coef <- vector("list", 3L)
      for (stick in seq_len(3L)) {
        idx <- which(q$stick_index == stick)
        coef[[stick]] <- tr[ix, idx, drop = FALSE]
        eta[, , stick] <- A %*% t(coef[[stick]])
        eta_sum[, stick] <- eta_sum[, stick] + rowSums(eta[, , stick])
        eta_sq[, stick] <- eta_sq[, stick] + rowSums(eta[, , stick]^2)
        eta_selected[ix, , stick] <- t(eta[representative, , stick, drop = FALSE][, , 1L])
      }
      cont <- lapply(seq_len(3L), function(stick) stats::pnorm(eta[, , stick]))
      pi <- array(0, c(m, nc, 4L)); rem <- matrix(1, m, nc)
      for (stick in seq_len(3L)) { pi[, , stick] <- rem * (1 - cont[[stick]]); rem <- rem * cont[[stick]] }
      pi[, , 4L] <- rem
      for (k in seq_len(4L)) { pi_sum[, k] <- pi_sum[, k] + rowSums(pi[, , k]); pi_sq[, k] <- pi_sq[, k] + rowSums(pi[, , k]^2) }
      row_error <- abs(apply(pi, c(1, 2), sum) - 1)
      max_row_error <- max(max_row_error, row_error); min_probability <- min(min_probability, pi)
      max_probability <- max(max_probability, pi)
      enriched <- A[, "enriched_binary"] == 1
      for (k in seq_len(4L)) observed_draw[ix, k + 1L] <- colMeans(pi[enriched, , k]) - colMeans(pi[!enriched, , k])
      observed_draw[ix, 1L] <- -observed_draw[ix, 2L]
      for (a in seq_along(mods)) {
        annotation <- names(mods)[a]; pair <- mods[[a]]; pi_cf <- vector("list", 2L)
        annotation_index <- match(annotation, colnames(A))
        for (side in 1:2) {
          cfs <- lapply(seq_len(3L), function(stick) stats::pnorm(eta[, , stick] +
            outer(pair[side] - A[, annotation], coef[[stick]][, annotation_index])))
          z <- array(0, c(m, nc, 4L)); rr <- matrix(1, m, nc)
          for (stick in seq_len(3L)) { z[, , stick] <- rr * (1 - cfs[[stick]]); rr <- rr * cfs[[stick]] }
          z[, , 4L] <- rr; pi_cf[[side]] <- z
        }
        for (k in seq_len(4L)) contrast_draw[ix, a, k + 1L] <-
          colMeans(pi_cf[[2L]][, , k]) - colMeans(pi_cf[[1L]][, , k])
        contrast_draw[ix, a, 1L] <- -contrast_draw[ix, a, 2L]
      }
    }
    colnames(eta_sum) <- paste0("stick_", 1:3); colnames(pi_sum) <- paste0("pi", 0:3)
    marker[[ch]] <- list(eta_mean = eta_sum / nd,
      eta_sd = sqrt(pmax(0, (eta_sq - eta_sum^2 / nd) / (nd - 1))),
      pi_mean = pi_sum / nd, pi_sd = sqrt(pmax(0, (pi_sq - pi_sum^2 / nd) / (nd - 1))))
    scalar[[ch]] <- eta_selected; contrast[[ch]] <- contrast_draw; observed[[ch]] <- observed_draw
    normalization[[ch]] <- data.frame(route, chain = ch, min_probability,
      max_probability, max_row_sum_error = max_row_error, all_finite = TRUE)
    message(sprintf("%s: completed chain %d/%d", route, ch, length(alpha$traces)))
  }
  eta_stability <- pi_stability <- list()
  for (x in seq_len(length(marker) - 1L)) for (y in (x + 1L):length(marker)) {
    for (stick in seq_len(3L)) eta_stability[[length(eta_stability) + 1L]] <-
      data.frame(route, chain_1 = x, chain_2 = y, quantity = paste0("eta_stick_", stick),
        t(study06_pair_metrics(marker[[x]]$eta_mean[, stick], marker[[y]]$eta_mean[, stick])))
    for (k in seq_len(4L)) pi_stability[[length(pi_stability) + 1L]] <-
      data.frame(route, chain_1 = x, chain_2 = y, quantity = paste0("pi", k - 1L),
        t(study06_pair_metrics(marker[[x]]$pi_mean[, k], marker[[y]]$pi_mean[, k])))
    pi_stability[[length(pi_stability) + 1L]] <- data.frame(route, chain_1 = x,
      chain_2 = y, quantity = "p_active", t(study06_pair_metrics(
        1 - marker[[x]]$pi_mean[, 1L], 1 - marker[[y]]$pi_mean[, 1L])))
  }
  eta_scalar <- list()
  for (r in seq_along(representative)) for (stick in seq_len(3L)) {
    vals <- lapply(scalar, function(z) z[, r, stick]); dg <- study06_diag(vals)
    target <- sum(A[representative[r], ] * truth_alpha[, stick])
    eta_scalar[[length(eta_scalar) + 1L]] <- data.frame(route,
      marker_id = rownames(A)[representative[r]], stick, truth = target,
      posterior_mean = mean(unlist(vals)), bias = mean(unlist(vals)) - target, t(dg))
  }
  contrast_table <- list(); observed_table <- list()
  truth_contrast <- function(annotation, target) {
    pair <- mods[[annotation]]; aa <- A; aa[, annotation] <- pair[1]; p0 <- study06_pi(aa %*% truth_alpha)
    aa[, annotation] <- pair[2]; p1 <- study06_pi(aa %*% truth_alpha)
    values <- c(active = mean(1 - p1[, 1]) - mean(1 - p0[, 1]), colMeans(p1 - p0))
    values[target]
  }
  for (a in seq_along(mods)) for (target in dimnames(contrast[[1]])[[3]]) {
    vals <- lapply(contrast, function(z) z[, a, target]); pooled <- unlist(vals)
    dg <- study06_diag(vals); tv <- truth_contrast(names(mods)[a], target)
    contrast_table[[length(contrast_table) + 1L]] <- data.frame(route,
      annotation = names(mods)[a], target, comparison = if (a == 1L) "1_vs_0" else "+1SD_vs_-1SD",
      posterior_mean = mean(pooled), posterior_sd = stats::sd(pooled), median = stats::median(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      p_gt_0 = mean(pooled > 0), p_lt_0 = mean(pooled < 0), t(dg), truth = tv,
      bias = mean(pooled) - tv, coverage = tv >= stats::quantile(pooled, .025) && tv <= stats::quantile(pooled, .975))
    ovals <- lapply(observed, function(z) z[, target]); op <- unlist(ovals); od <- study06_diag(ovals)
    observed_table[[length(observed_table) + 1L]] <- data.frame(route,
      annotation = "enriched_binary", target, posterior_mean = mean(op),
      q025 = stats::quantile(op, .025), q975 = stats::quantile(op, .975), t(od))
  }
  pooled_eta <- Reduce(`+`, lapply(marker, `[[`, "eta_mean")) / length(marker)
  pooled_pi <- Reduce(`+`, lapply(marker, `[[`, "pi_mean")) / length(marker)
  eta_truth <- A %*% truth_alpha
  recovery <- rbind(data.frame(route, level = "eta", quantity = paste0("stick_", 1:3),
    rmse = sqrt(colMeans((pooled_eta - eta_truth)^2)), mae = colMeans(abs(pooled_eta - eta_truth))),
    data.frame(route, level = "pi", quantity = paste0("pi", 0:3),
      rmse = sqrt(colMeans((pooled_pi - truth_pi)^2)), mae = colMeans(abs(pooled_pi - truth_pi))))
  list(marker = marker, eta_stability = do.call(rbind, eta_stability),
    pi_stability = do.call(rbind, pi_stability), eta_scalar = do.call(rbind, eta_scalar),
    contrasts = do.call(rbind, contrast_table), observed = unique(do.call(rbind, observed_table)),
    normalization = do.call(rbind, normalization), recovery = recovery,
    representatives = data.frame(route, marker_index = representative, marker_id = rownames(A)[representative]))
}

study06_snp_stability <- function(route, fit, truth) {
  rows <- list(); chain_metrics <- list()
  chains <- study06_chain_posterior(fit)
  for (ch in seq_along(chains)) {
    pip <- chains[[ch]]$pip; beta <- chains[[ch]]$beta
    gv <- drop(truth$validation_x %*% beta)
    chain_metrics[[ch]] <- data.frame(route, chain = ch,
      auprc = study06_auprc(pip, truth$marker_truth$true_nonnull),
      auroc = study06_auc(pip, truth$marker_truth$true_nonnull),
      validation_genetic_correlation = stats::cor(gv, truth$validation_genetic),
      phenotype_prediction_correlation = stats::cor(gv, truth$validation_phenotype))
  }
  for (x in seq_len(length(chains) - 1L)) for (y in (x + 1L):length(chains)) {
    for (quantity in c("pip", "abs_beta")) {
      a <- if (quantity == "pip") chains[[x]]$pip else abs(chains[[x]]$beta)
      b <- if (quantity == "pip") chains[[y]]$pip else abs(chains[[y]]$beta)
      rows[[length(rows) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y, quantity,
        t(study06_pair_metrics(a, b, top = c(50L, 100L))))
    }
    bx <- chains[[x]]$beta; by <- chains[[y]]$beta
    rows[[length(rows) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y,
      quantity = "beta_signed", pearson = stats::cor(bx, by),
      spearman = stats::cor(bx, by, method = "spearman"), rmse = sqrt(mean((bx - by)^2)),
      mae = mean(abs(bx - by)), max_abs = max(abs(bx - by)), top50_overlap = NA, top100_overlap = NA)
  }
  list(stability = do.call(rbind, rows), performance = do.call(rbind, chain_metrics))
}

study06_official <- function(A, truth_alpha, base) {
  files <- file.path(base, paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))
  if (!all(file.exists(files))) return(list(status = "alpha traces unavailable"))
  traces <- lapply(files, function(f) as.matrix(utils::read.table(f, header = TRUE,
    check.names = FALSE)))
  n <- min(vapply(traces, nrow, integer(1))); traces <- lapply(traces, function(x) tail(x, n))
  raw <- list()
  for (stick in seq_len(3L)) for (a in seq_len(4L)) {
    z <- traces[[stick]][, a]; first <- head(z, floor(n / 2)); last <- tail(z, floor(n / 2))
    raw[[length(raw) + 1L]] <- data.frame(stick, annotation = colnames(traces[[stick]])[a],
      mean = mean(z), sd = stats::sd(z), q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975),
      truth = truth_alpha[a, stick], bias = mean(z) - truth_alpha[a, stick],
      ess = coda::effectiveSize(z), relative_mcse = stats::sd(z) / sqrt(coda::effectiveSize(z)) / stats::sd(z),
      half_mean_shift_sd = (mean(last) - mean(first)) / stats::sd(z))
  }
  contrast <- list(); chunk <- 500L
  mods <- list(enriched_binary = c(0, 1), continuous_signal = c(-1, 1), null_annotation = c(-1, 1))
  draws <- array(NA_real_, c(n, length(mods), 5L), dimnames = list(NULL, names(mods), c("active", paste0("pi", 0:3))))
  for (lo in seq.int(1L, n, by = chunk)) {
    ix <- lo:min(n, lo + chunk - 1L)
    for (a in seq_along(mods)) {
      ps <- vector("list", 2L)
      for (side in 1:2) {
        aa <- A; aa[, names(mods)[a]] <- mods[[a]][side]
        eta <- array(NA_real_, c(nrow(A), length(ix), 3L))
        for (stick in 1:3) eta[, , stick] <- aa %*% t(traces[[stick]][ix, , drop = FALSE])
        cc <- lapply(1:3, function(s) stats::pnorm(eta[, , s])); p <- array(0, c(nrow(A), length(ix), 4L)); rem <- matrix(1, nrow(A), length(ix))
        for (s in 1:3) { p[, , s] <- rem * (1 - cc[[s]]); rem <- rem * cc[[s]] }; p[, , 4] <- rem; ps[[side]] <- p
      }
      for (k in 1:4) draws[ix, a, k + 1L] <- colMeans(ps[[2]][, , k]) - colMeans(ps[[1]][, , k])
      draws[ix, a, 1L] <- -draws[ix, a, 2L]
    }
  }
  for (a in seq_along(mods)) for (target in dimnames(draws)[[3]]) {
    z <- draws[, a, target]; first <- head(z, floor(n / 2)); last <- tail(z, floor(n / 2))
    contrast[[length(contrast) + 1L]] <- data.frame(annotation = names(mods)[a], target,
      mean = mean(z), sd = stats::sd(z), q025 = stats::quantile(z, .025), q975 = stats::quantile(z, .975),
      p_gt_0 = mean(z > 0), ess = coda::effectiveSize(z),
      half_mean_shift_sd = (mean(last) - mean(first)) / stats::sd(z))
  }
  list(status = "single deterministic native trajectory; multitrajectory inference unavailable",
    raw = do.call(rbind, raw), contrasts = do.call(rbind, contrast))
}

study06_estimability_validate <- function(root = ".") {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  required <- c(
    file.path(root, "studies/06_annotation_models/06a_estimability/spec.R"),
    file.path(root, "studies/06_annotation_models/shared/annotation-design.R"),
    file.path(root, "results/reference/06_annotation_models/final_decision.json"),
    file.path(root, "results/reference/06_annotation_models/final_hierarchy_of_evidence.csv"))
  if (basename(root) != "sblrbench" || any(!file.exists(required)))
    stop("Study 06A maintained specification or accepted evidence is incomplete.",
      call. = FALSE)
  invisible(list(root = root, required = required,
    status = "CLOSED — EST-R2", writes = FALSE))
}

study06_reproduction_arguments <- function(root, input_dir, output_dir) {
  if (is.null(input_dir) || !nzchar(input_dir))
    stop(paste(
      "A full offline reproduction requires an explicit --input path",
      "containing the historical Study 06A truth, retained sblr fits, and",
      "official D1 annotation-effect traces."), call. = FALSE)
  if (is.null(output_dir) || !nzchar(output_dir))
    stop("A full offline reproduction requires an explicit --output path.",
      call. = FALSE)

  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  reference_dir <- normalizePath(file.path(root, "results/reference"),
    winslash = "/", mustWork = TRUE)
  if (startsWith(tolower(output_dir), paste0(tolower(reference_dir), "/")) ||
      identical(tolower(output_dir), tolower(reference_dir)))
    stop("Study 06A reproduction output cannot be results/reference.",
      call. = FALSE)

  required_input <- c(
    truth = file.path(input_dir, "gctb_parity/export/truth.rds"),
    bed_fit = file.path(input_dir,
      "v2_identifiable_qualification/qualification/checkpoints",
      "informative_annotations--r1--st_bed_bayesrc.rds"),
    block_fit = file.path(input_dir,
      "v2_identifiable_qualification/qualification/checkpoints",
      "informative_annotations--r1--st_block_eigen_sbayesrc.rds"),
    official_p1 = file.path(input_dir,
      "gctb_single_trajectory/runs/D1/D1.mcmcsamples",
      "D1.mcmcsamples.AnnoEffects_p1"),
    official_p2 = file.path(input_dir,
      "gctb_single_trajectory/runs/D1/D1.mcmcsamples",
      "D1.mcmcsamples.AnnoEffects_p2"),
    official_p3 = file.path(input_dir,
      "gctb_single_trajectory/runs/D1/D1.mcmcsamples",
      "D1.mcmcsamples.AnnoEffects_p3"))
  if (!dir.exists(input_dir) || any(!file.exists(required_input)))
    stop(paste(
      "The explicit --input path does not contain the complete historical",
      "Study 06A reproduction inputs:", input_dir), call. = FALSE)

  list(input_dir = input_dir, output_dir = output_dir,
    required_input = required_input)
}

run_study06_estimability <- function(root = ".", input_dir = NULL,
                                     output_dir = NULL, dry_run = TRUE) {
  validation <- study06_estimability_validate(root)
  root <- validation$root
  if (isTRUE(dry_run)) return(validation)
  arguments <- study06_reproduction_arguments(root, input_dir, output_dir)
  input_dir <- arguments$input_dir
  output_dir <- arguments$output_dir
  suppressPackageStartupMessages({ library(posterior); library(ggplot2) })
  source(file.path(root, "studies/06_annotation_models/06a_estimability/spec.R"), local = TRUE)
  source(file.path(root, "studies/06_annotation_models/shared/annotation-design.R"), local = TRUE)
  out <- output_dir
  dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
  truth_path <- file.path(input_dir, "gctb_parity/export/truth.rds")
  truth <- readRDS(truth_path); A <- truth$annotations
  stopifnot(identical(rownames(A), truth$marker_ids), identical(colnames(A), spec$annotation_design$columns))
  truth_object <- construct_annotation_truth(A, spec); truth_alpha <- truth_object$informative_annotations
  truth_pi <- study06_pi(A %*% truth_alpha)
  stopifnot(max(abs(truth_pi - as.matrix(truth$marker_truth[paste0("true_prior_component_", 0:3)]))) < 1e-10)
  fit_paths <- c(bed = file.path(input_dir, "v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_bed_bayesrc.rds"),
    block = file.path(input_dir, "v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_block_eigen_sbayesrc.rds"))
  fits <- lapply(fit_paths, function(p) readRDS(p)$result$native_fit)
  all_raw <- all_raw_chain <- all_design <- all_design_cor <- all_postmat <- all_eigen <- list()
  all_eta <- all_pi <- all_eta_scalar <- all_contrast <- all_observed <- all_norm <- all_recovery <- all_snp <- all_perf <- list(); derived <- list()
  for (route in names(fits)) {
    message(sprintf("Starting retained-draw analysis for %s", route))
    alpha <- study06_chain_alpha(fits[[route]])
    raw <- study06_raw_alpha(route, alpha, truth_alpha); dep <- study06_design_dependence(route, A, alpha)
    drv <- study06_derive_route(route, fits[[route]], alpha, A, truth_alpha, truth_pi); derived[[route]] <- drv
    snp <- study06_snp_stability(route, fits[[route]], truth)
    all_raw[[route]] <- raw$summary; all_raw_chain[[route]] <- raw$chains
    all_design[[route]] <- dep$design; all_design_cor[[route]] <- dep$design_correlation
    all_postmat[[route]] <- dep$posterior_matrices; all_eigen[[route]] <- dep$eigen
    all_eta[[route]] <- drv$eta_stability; all_pi[[route]] <- drv$pi_stability; all_eta_scalar[[route]] <- drv$eta_scalar
    all_contrast[[route]] <- drv$contrasts; all_observed[[route]] <- drv$observed; all_norm[[route]] <- drv$normalization
    all_recovery[[route]] <- drv$recovery; all_snp[[route]] <- snp$stability; all_perf[[route]] <- snp$performance
    message(sprintf("Completed retained-draw analysis for %s", route))
  }
  tables <- list(raw_alpha_convergence = do.call(rbind, all_raw), raw_alpha_chain = do.call(rbind, all_raw_chain),
    annotation_design = do.call(rbind, all_design), annotation_design_correlations = do.call(rbind, all_design_cor),
    alpha_posterior_matrices = do.call(rbind, all_postmat), alpha_eigendirections = do.call(rbind, all_eigen),
    A_alpha_stability = do.call(rbind, all_eta), selected_A_alpha_diagnostics = do.call(rbind, all_eta_scalar),
    prior_probability_stability = do.call(rbind, all_pi), probability_normalization = do.call(rbind, all_norm),
    derived_truth_recovery = do.call(rbind, all_recovery), counterfactual_annotation_contrasts = do.call(rbind, all_contrast),
    observed_group_contrasts = do.call(rbind, all_observed), SNP_chain_stability = do.call(rbind, all_snp),
    SNP_chain_performance = do.call(rbind, all_perf))
  official_base <- file.path(input_dir, "gctb_single_trajectory/runs/D1/D1.mcmcsamples")
  official <- study06_official(A, truth_alpha, official_base)
  message("Completed official single-trajectory derived summaries")
  if (!is.null(official$raw)) tables$official_alpha_summary <- official$raw
  if (!is.null(official$contrasts)) tables$official_contrast_summary <- official$contrasts
  for (nm in names(tables)) study06_write_csv(tables[[nm]], file.path(out, "tables", paste0(nm, ".csv")))
  saveRDS(lapply(derived, `[[`, "marker"), file.path(out, "marker_derived_chain_summaries.rds"))
  manifest <- data.frame(role = c("truth", names(fit_paths), "official_D1_alpha_p1", "official_D1_alpha_p2", "official_D1_alpha_p3"),
    path = c(truth_path, fit_paths, file.path(official_base, paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))), stringsAsFactors = FALSE)
  manifest$bytes <- file.info(manifest$path)$size
  manifest$sha256 <- vapply(manifest$path, digest::digest, character(1), algo = "sha256", file = TRUE)
  manifest$specification_hash <- "241c15afab8fefc571e38e625130de6e4ab58b958c30cff369e26998ee30fa56"
  manifest$truth_hash <- "169d52bef390022a9106d7e61b200493869b40cbb76f1c8d911eebbc80fea1eb"
  study06_write_csv(manifest, file.path(out, "analysis_input_manifest.csv"))
  p1 <- ggplot(tables$raw_alpha_chain, aes(chain, mean, colour = factor(chain))) + geom_point() +
    facet_grid(annotation ~ route + stick, scales = "free_y") + theme_bw() + guides(colour = "none")
  ggsave(file.path(out, "figures/alpha_chain_disagreement.png"), p1, width = 12, height = 8, dpi = 160)
  p2 <- ggplot(tables$A_alpha_stability, aes(pearson, rmse, colour = route)) + geom_point() + facet_wrap(~quantity, scales = "free") + theme_bw()
  ggsave(file.path(out, "figures/A_alpha_chain_agreement.png"), p2, width = 9, height = 5, dpi = 160)
  p3 <- ggplot(tables$prior_probability_stability, aes(pearson, rmse, colour = route)) + geom_point() + facet_wrap(~quantity, scales = "free") + theme_bw()
  ggsave(file.path(out, "figures/marker_prior_probability_agreement.png"), p3, width = 9, height = 5, dpi = 160)
  p4 <- ggplot(tables$counterfactual_annotation_contrasts[tables$counterfactual_annotation_contrasts$target == "active", ],
    aes(annotation, posterior_mean, ymin = q025, ymax = q975, colour = route)) + geom_pointrange(position = position_dodge(.4)) + geom_hline(yintercept = 0, linetype = 2) + theme_bw()
  ggsave(file.path(out, "figures/annotation_contrast_posteriors.png"), p4, width = 8, height = 5, dpi = 160)
  invisible(list(tables = tables, derived = derived, official = official, manifest = manifest))
}


# Fast final aggregation for Study 06. Requires estimability-and-contrasts.R.

finalize_study06_estimability <- function(root = ".", input_dir = NULL,
                                          output_dir = NULL) {
  root <- study06_estimability_validate(root)$root
  arguments <- study06_reproduction_arguments(root, input_dir, output_dir)
  input_dir <- arguments$input_dir
  output_dir <- arguments$output_dir
  out <- output_dir
  td <- file.path(out, "tables")
  truth <- readRDS(file.path(input_dir, "gctb_parity/export/truth.rds"))
  A <- truth$annotations
  fit_paths <- c(bed = file.path(input_dir, "v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_bed_bayesrc.rds"),
    block = file.path(input_dir, "v2_identifiable_qualification/qualification/checkpoints/informative_annotations--r1--st_block_eigen_sbayesrc.rds"))
  fits <- lapply(fit_paths, function(p) readRDS(p)$result$native_fit)
  source(file.path(root, "studies/06_annotation_models/06a_estimability/spec.R"), local = TRUE)
  source(file.path(root, "studies/06_annotation_models/shared/annotation-design.R"), local = TRUE)
  truth_alpha <- construct_annotation_truth(A, spec)$informative_annotations
  truth_pi <- study06_pi(A %*% truth_alpha)
  marker <- readRDS(file.path(out, "marker_derived_chain_summaries.rds"))

  raw <- raw_chain <- expected <- representative <- ranking <- gv <- eta_chain <- probability_crosscheck <- list()
  reps <- study06_representatives(A)
  for (route in names(fits)) {
    alpha <- study06_chain_alpha(fits[[route]])
    chains <- study06_chain_posterior(fits[[route]])
    rr <- study06_raw_alpha(route, alpha, truth_alpha)
    raw[[route]] <- rr$summary; raw_chain[[route]] <- rr$chains
    q <- alpha$descriptor; active_draws <- vector("list", length(alpha$traces))
    rep_draws <- lapply(seq_along(alpha$traces), function(i) array(NA_real_, c(nrow(alpha$traces[[i]]), length(reps), 5L),
      dimnames = list(NULL, NULL, c("active", paste0("pi", 0:3)))))
    for (ch in seq_along(alpha$traces)) {
      tr <- alpha$traces[[ch]]
      idx1 <- which(q$stick_index == 1L)
      active_draws[[ch]] <- rowSums(stats::pnorm(A %*% t(tr[, idx1, drop = FALSE])))
      eta <- lapply(seq_len(3L), function(s) tr[, which(q$stick_index == s), drop = FALSE] %*% t(A[reps, , drop = FALSE]))
      c1 <- stats::pnorm(eta[[1]]); c2 <- stats::pnorm(eta[[2]]); c3 <- stats::pnorm(eta[[3]])
      rep_draws[[ch]][, , 2L] <- 1 - c1
      rep_draws[[ch]][, , 3L] <- c1 * (1 - c2)
      rep_draws[[ch]][, , 4L] <- c1 * c2 * (1 - c3)
      rep_draws[[ch]][, , 5L] <- c1 * c2 * c3
      rep_draws[[ch]][, , 1L] <- c1
    }
    dg <- study06_diag(active_draws); pooled <- unlist(active_draws)
    expected[[route]] <- data.frame(route, posterior_mean = mean(pooled), posterior_sd = stats::sd(pooled),
      q025 = stats::quantile(pooled, .025), q975 = stats::quantile(pooled, .975),
      truth = sum(1 - truth_pi[, 1]), bias = mean(pooled) - sum(1 - truth_pi[, 1]), t(dg))
    for (r in seq_along(reps)) for (target in dimnames(rep_draws[[1]])[[3]]) {
      vals <- lapply(rep_draws, function(z) z[, r, target]); z <- unlist(vals); d <- study06_diag(vals)
      tv <- if (target == "active") 1 - truth_pi[reps[r], 1] else truth_pi[reps[r], match(target, paste0("pi", 0:3))]
      representative[[length(representative) + 1L]] <- data.frame(route, marker_id = rownames(A)[reps[r]], target,
        posterior_mean = mean(z), truth = tv, bias = mean(z) - tv, t(d))
    }
    mm <- marker[[route]]
    for (ch in seq_along(alpha$traces)) {
      eta_sd <- matrix(mm[[ch]]$eta_sd, nrow = nrow(A), ncol = 3L)
      for (stick in seq_len(3L)) eta_chain[[length(eta_chain) + 1L]] <- data.frame(route, chain = ch, stick,
        posterior_mean_across_markers = mean(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_marker_sd = stats::sd(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_q025 = stats::quantile(mm[[ch]]$eta_mean[, stick], .025),
        posterior_mean_median = stats::median(mm[[ch]]$eta_mean[, stick]),
        posterior_mean_q975 = stats::quantile(mm[[ch]]$eta_mean[, stick], .975),
        posterior_sd_across_markers = mean(eta_sd[, stick]),
        posterior_sd_q025 = stats::quantile(eta_sd[, stick], .025),
        posterior_sd_median = stats::median(eta_sd[, stick]),
        posterior_sd_q975 = stats::quantile(eta_sd[, stick], .975))
      checked_draws <- unique(c(1L, ceiling(nrow(alpha$traces[[ch]]) / 2), nrow(alpha$traces[[ch]])))
      for (draw in checked_draws) {
        alpha_draw <- matrix(alpha$traces[[ch]][draw, ], nrow = ncol(A), ncol = 3L)
        offline <- study06_pi(A %*% alpha_draw)
        package <- sblr::sbayesrc_marker_pi(A, alpha_draw, gamma = c(0, .01, .1, 1))
        probability_crosscheck[[length(probability_crosscheck) + 1L]] <- data.frame(route, chain = ch,
          draw, maximum_absolute_difference = max(abs(offline - package)),
          comparison = "offline exact transform versus sblr::sbayesrc_marker_pi")
      }
    }
    for (x in seq_len(length(chains) - 1L)) for (y in (x + 1L):length(chains)) {
      quantities <- list(prior_active = 1 - mm[[x]]$pi_mean[, 1],
        largest_component_prior = mm[[x]]$pi_mean[, 4], pip = chains[[x]]$pip,
        absolute_posterior_beta = abs(chains[[x]]$beta))
      quantities_y <- list(prior_active = 1 - mm[[y]]$pi_mean[, 1],
        largest_component_prior = mm[[y]]$pi_mean[, 4], pip = chains[[y]]$pip,
        absolute_posterior_beta = abs(chains[[y]]$beta))
      for (nm in names(quantities)) ranking[[length(ranking) + 1L]] <- data.frame(route,
        chain_1 = x, chain_2 = y, quantity = nm,
        t(study06_pair_metrics(quantities[[nm]], quantities_y[[nm]], top = c(50L, 100L))))
      gx <- drop(truth$validation_x %*% chains[[x]]$beta)
      gy <- drop(truth$validation_x %*% chains[[y]]$beta)
      gv[[length(gv) + 1L]] <- data.frame(route, chain_1 = x, chain_2 = y,
        correlation = stats::cor(gx, gy), rmse = sqrt(mean((gx - gy)^2)))
    }
  }
  raw <- do.call(rbind, raw); raw_chain <- do.call(rbind, raw_chain)
  expected <- do.call(rbind, expected); representative <- do.call(rbind, representative)
  ranking <- do.call(rbind, ranking); gv <- do.call(rbind, gv)
  eta_chain <- do.call(rbind, eta_chain); probability_crosscheck <- do.call(rbind, probability_crosscheck)
  study06_write_csv(raw, file.path(td, "raw_alpha_convergence.csv"))
  study06_write_csv(raw_chain, file.path(td, "raw_alpha_chain.csv"))
  study06_write_csv(expected, file.path(td, "expected_active_count.csv"))
  study06_write_csv(representative, file.path(td, "representative_prior_probability_diagnostics.csv"))
  study06_write_csv(ranking, file.path(td, "ranking_stability.csv"))
  study06_write_csv(gv, file.path(td, "genetic_value_stability.csv"))
  study06_write_csv(eta_chain, file.path(td, "A_alpha_chain_summary.csv"))
  study06_write_csv(probability_crosscheck, file.path(td, "probability_package_crosscheck.csv"))
  study06_write_csv(read.csv(file.path(td, "selected_A_alpha_diagnostics.csv"), check.names = FALSE),
    file.path(td, "selected_alpha_linear_contrasts.csv"))

  eta <- read.csv(file.path(td, "A_alpha_stability.csv"), check.names = FALSE)
  pi_stab <- read.csv(file.path(td, "prior_probability_stability.csv"), check.names = FALSE)
  recovery <- read.csv(file.path(td, "derived_truth_recovery.csv"), check.names = FALSE)
  contrast <- read.csv(file.path(td, "counterfactual_annotation_contrasts.csv"), check.names = FALSE)
  snp <- read.csv(file.path(td, "SNP_chain_stability.csv"), check.names = FALSE)
  perf <- read.csv(file.path(td, "SNP_chain_performance.csv"), check.names = FALSE)
  selected <- read.csv(file.path(td, "selected_A_alpha_diagnostics.csv"), check.names = FALSE)
  hierarchy_row <- function(route, level, d = NULL, corr = NA_real_, truth_error = NA_real_) {
    data.frame(route, level,
      worst_rhat = if (is.null(d) || !"rhat" %in% names(d)) NA else max(d$rhat, na.rm = TRUE),
      median_rhat = if (is.null(d) || !"rhat" %in% names(d)) NA else stats::median(d$rhat, na.rm = TRUE),
      minimum_bulk_ess = if (is.null(d) || !"bulk_ess" %in% names(d)) NA else min(d$bulk_ess, na.rm = TRUE),
      maximum_relative_mcse = if (is.null(d) || !"relative_mcse" %in% names(d)) NA else max(d$relative_mcse, na.rm = TRUE),
      minimum_between_chain_correlation = corr, truth_error = truth_error)
  }
  hierarchy <- list()
  for (route in names(fits)) {
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "raw alpha", raw[raw$route == route, ], truth_error = max(raw$absolute_error[raw$route == route]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "selected alpha linear contrasts", selected[selected$route == route, ], truth_error = max(abs(selected$bias[selected$route == route])))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "A alpha", selected[selected$route == route, ], min(eta$pearson[eta$route == route]), max(recovery$rmse[recovery$route == route & recovery$level == "eta"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "prior active probability", expected[expected$route == route, ], min(pi_stab$pearson[pi_stab$route == route & pi_stab$quantity == "p_active"]), abs(expected$bias[expected$route == route]) / nrow(A))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "component probabilities", representative[representative$route == route, ], min(pi_stab$pearson[pi_stab$route == route & grepl("^pi", pi_stab$quantity)]), max(recovery$rmse[recovery$route == route & recovery$level == "pi"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "counterfactual annotation contrasts", contrast[contrast$route == route, ], truth_error = max(abs(contrast$bias[contrast$route == route])))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "SNP PIP", corr = min(snp$pearson[snp$route == route & snp$quantity == "pip"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "posterior beta", corr = min(snp$pearson[snp$route == route & snp$quantity == "beta_signed"]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "genetic value", corr = min(gv$correlation[gv$route == route]))
    hierarchy[[length(hierarchy)+1L]] <- hierarchy_row(route, "prediction", truth_error = 1 - min(perf$phenotype_prediction_correlation[perf$route == route]))
  }
  hierarchy <- do.call(rbind, hierarchy)
  study06_write_csv(hierarchy, file.path(td, "hierarchy_of_stability.csv"))

  official_raw <- read.csv(file.path(td, "official_alpha_summary.csv"), check.names = FALSE)
  official_contrast <- read.csv(file.path(td, "official_contrast_summary.csv"), check.names = FALSE)
  contrast_truth <- function(annotation, target) {
    pair <- if (annotation == "enriched_binary") c(0, 1) else c(-1, 1)
    aa <- A; aa[, annotation] <- pair[1]; p0 <- study06_pi(aa %*% truth_alpha)
    aa[, annotation] <- pair[2]; p1 <- study06_pi(aa %*% truth_alpha)
    values <- c(active = mean(1 - p1[, 1]) - mean(1 - p0[, 1]), colMeans(p1 - p0))
    unname(values[target])
  }
  official_contrast$truth <- mapply(contrast_truth, official_contrast$annotation, official_contrast$target)
  official_contrast$bias <- official_contrast$mean - official_contrast$truth
  official_contrast$coverage <- official_contrast$truth >= official_contrast$q025 & official_contrast$truth <= official_contrast$q975
  official_contrast$p_lt_0 <- 1 - official_contrast$p_gt_0
  official_contrast$relative_mcse <- 1 / sqrt(official_contrast$ess)
  study06_write_csv(official_contrast, file.path(td, "official_contrast_summary.csv"))
  ofiles <- file.path(input_dir, "gctb_single_trajectory/runs/D1/D1.mcmcsamples",
    paste0("D1.mcmcsamples.AnnoEffects_p", 1:3))
  otr <- lapply(ofiles, function(f) as.matrix(read.table(f, header = TRUE, check.names = FALSE)))
  n <- min(vapply(otr, nrow, integer(1))); half <- floor(n/2); official_level <- list()
  for (part in list(first = seq_len(half), second = (n-half+1L):n)) {
    am <- sapply(otr, function(z) colMeans(z[part, , drop = FALSE])); rownames(am) <- colnames(otr[[1]])
    ep <- A %*% am; pp <- study06_pi(ep)
    official_level[[length(official_level)+1L]] <- list(eta = ep, pi = pp)
  }
  official_compare <- data.frame(level = c("raw alpha", "A alpha", "component probabilities", "annotation contrasts"),
    official_evidence = c(sprintf("single-chain ESS %.1f-%.1f; max half shift %.2f SD", min(official_raw$ess), max(official_raw$ess), max(abs(official_raw$half_mean_shift_sd))),
      sprintf("half-vector minimum Pearson %.3f", min(vapply(1:3, function(k) cor(official_level[[1]]$eta[,k], official_level[[2]]$eta[,k]), numeric(1)))),
      sprintf("half-vector minimum Pearson %.3f", min(vapply(1:4, function(k) suppressWarnings(cor(official_level[[1]]$pi[,k], official_level[[2]]$pi[,k])), numeric(1)), na.rm = TRUE)),
      sprintf("active contrast max half shift %.2f SD", max(abs(official_contrast$half_mean_shift_sd[official_contrast$target == "active"])))),
    sblr_evidence = c(sprintf("max R-hat %.3f", max(raw$rhat)),
      sprintf("selected eta max R-hat %.3f; min vector Pearson %.3f", max(selected$rhat), min(eta$pearson)),
      sprintf("representative max R-hat %.3f; min vector Pearson %.3f", max(representative$rhat), min(pi_stab$pearson)),
      sprintf("max R-hat %.3f", max(contrast$rhat))),
    limitation = c("official trajectories cannot be independently seeded", rep("single official trajectory only", 3)))
  study06_write_csv(official_compare, file.path(td, "official_sbayesrc_comparison.csv"))

  p1 <- ggplot2::ggplot(raw_chain, ggplot2::aes(chain, mean, colour = factor(chain))) + ggplot2::geom_point() +
    ggplot2::facet_grid(annotation ~ route + stick, scales = "free_y") + ggplot2::theme_bw() + ggplot2::guides(colour = "none")
  ggplot2::ggsave(file.path(out, "figures/alpha_chain_disagreement.png"), p1, width = 12, height = 8, dpi = 160)
  hp <- rbind(data.frame(route = hierarchy$route, level = hierarchy$level, metric = "R-hat", value = hierarchy$worst_rhat),
    data.frame(route = hierarchy$route, level = hierarchy$level, metric = "between-chain correlation", value = hierarchy$minimum_between_chain_correlation))
  p5 <- ggplot2::ggplot(hp[is.finite(hp$value), ], ggplot2::aes(level, value, fill = route)) + ggplot2::geom_col(position = "dodge") +
    ggplot2::facet_wrap(~metric, scales = "free_y") + ggplot2::coord_flip() + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(out, "figures/hierarchy_of_stability.png"), p5, width = 10, height = 8, dpi = 160)
  invisible(list(hierarchy = hierarchy, official = official_compare))
}
