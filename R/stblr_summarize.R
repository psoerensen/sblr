# Convergence diagnostics are owned by the shared rank-normalized engine.
# This internal convenience returns the validated modern result already owned
# by the fit and never calculates a competing statistic.
check_stblr_convergence <- function(fit) {
  if (!inherits(fit, "stblr_fit") || is.null(fit$convergence)) {
    stop("fit must be an stblr_fit with modern convergence diagnostics.",
         call. = FALSE)
  }
  fit$convergence
}
#' Plot ST-BLR Posterior Summaries
#'
#' Plots posterior means or medians with HPD intervals, equal-tail intervals,
#' or both interval types from the data frame returned by
#' [summarise_posterior()]. The `parameters` argument uses internal
#' parameter names such as `vg`, `ve`, `h2`, `pi`, `vb`, and `varch`, while plot
#' labels use the `label` column or built-in plotmath labels.
#'
#' The function uses base R graphics.
#'
#' @param post Posterior summary data frame from
#'   [summarise_posterior()].
#' @param parameters Optional character vector of internal parameter names to
#'   plot.
#' @param traits Optional character vector of traits to plot.
#' @param interval Interval type to plot: `"hpd"`, `"quantile"`, or `"both"`.
#' @param point Point estimate to plot: `"mean"` or `"median"`.
#' @param facet_by Faceting variable: `"parameter"` or `"trait"`.
#' @param log_scale Logical; use a log-scaled x-axis when possible.
#' @param xlab Optional x-axis label.
#' @param main Optional panel title override.
#' @param cex_axis Axis-label character expansion.
#' @param pch Point plotting character.
#' @param lwd_hpd,lwd_quantile Line widths for HPD and equal-tail intervals.
#' @param col_point,col_hpd,col_quantile Colors for points, HPD intervals, and
#'   equal-tail intervals.
#' @param parse_labels Logical; parse parameter labels as plotmath
#'   expressions.
#'
#' @return Invisibly returns the filtered data frame used for plotting.
#'
#' @examples
#' fake_post <- data.frame(
#'   parameter = c("vg", "ve", "h2", "pi", "varch"),
#'   label = c("V[g]", "V[e]", "h^2", "pi", "V[arch]"),
#'   trait = "trait1",
#'   n = 20L,
#'   mean = c(0.6, 0.4, 0.6, 0.02, 0.1),
#'   median = c(0.6, 0.4, 0.6, 0.02, 0.1),
#'   sd = 0.01,
#'   mcse = 0.002,
#'   q_lower = c(0.55, 0.35, 0.55, 0.01, 0.08),
#'   q_upper = c(0.65, 0.45, 0.65, 0.03, 0.12),
#'   hpd_lower = c(0.56, 0.36, 0.56, 0.01, 0.08),
#'   hpd_upper = c(0.64, 0.44, 0.64, 0.03, 0.12),
#'   ess = NA_real_,
#'   autocorr_lag1 = NA_real_,
#'   min = c(0.54, 0.34, 0.54, 0.01, 0.07),
#'   max = c(0.66, 0.46, 0.66, 0.03, 0.13)
#' )
#' if (interactive()) {
#'   plot_posterior(fake_post, parameters = c("vg", "ve", "h2"))
#' }
#'
#' @export
plot_posterior <- function(
    post,
    parameters = NULL,
    traits = NULL,
    interval = c("hpd", "quantile", "both"),
    point = c("mean", "median"),
    facet_by = c("parameter", "trait"),
    log_scale = FALSE,
    xlab = NULL,
    main = NULL,
    cex_axis = 0.8,
    pch = 19,
    lwd_hpd = 4,
    lwd_quantile = 1,
    col_point = "black",
    col_hpd = "gray35",
    col_quantile = "gray70",
    parse_labels = TRUE
) {
  interval <- match.arg(interval)
  point <- match.arg(point)
  facet_by <- match.arg(facet_by)
  
  required <- c(
    "parameter", "trait", "mean", "median",
    "q_lower", "q_upper", "hpd_lower", "hpd_upper"
  )
  
  missing <- setdiff(required, names(post))
  if (length(missing) > 0L) {
    stop(
      "post is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
  
  parameter_labels <- c(
    vb = "V[b]",
    vg = "V[g]",
    ve = "V[e]",
    pi = "pi",
    vle = "V[LE]",
    vld = "V[LD]",
    h2 = "h^2",
    ve_ratio = "V[e] / (V[g] + V[e])",
    le_ratio = "V[LE] / V[g]",
    ld_ratio = "V[LD] / V[g]",
    varch = "V[arch]",
    varch_ratio = "V[arch] / V[g]",
    h2_arch = "h[arch]^2",
    m_included = "pi * m"
  )
  label_source <- parameter_labels
  
  label_for <- function(x) {
    x <- as.character(x)
    out <- unname(label_source[x])
    miss <- is.na(out)
    out[miss] <- x[miss]
    out
  }
  
  make_labels <- function(x) {
    labs <- label_for(x)
    
    if (!parse_labels) {
      return(labs)
    }
    
    parsed <- tryCatch(
      parse(text = labs),
      error = function(e) labs
    )
    
    parsed
  }
  
  panel_title <- function(x) {
    lab <- label_for(x)
    
    if (!parse_labels) {
      return(lab)
    }
    
    tryCatch(
      parse(text = lab)[[1]],
      error = function(e) lab
    )
  }
  
  dat <- post
  
  if (!is.null(parameters)) {
    dat <- dat[dat$parameter %in% parameters, , drop = FALSE]
  }
  
  if (!is.null(traits)) {
    dat <- dat[dat$trait %in% traits, , drop = FALSE]
  }
  
  if (nrow(dat) < 1L) {
    stop("No rows left after filtering parameters/traits.")
  }
  
  dat$parameter <- as.character(dat$parameter)
  dat$trait <- as.character(dat$trait)
  
  if ("label" %in% names(dat)) {
    supplied <- dat[!is.na(dat$label), c("parameter", "label"), drop = FALSE]
    supplied$label <- as.character(supplied$label)
    supplied <- supplied[nzchar(supplied$label), , drop = FALSE]
    if (nrow(supplied) > 0L) {
      supplied_labels <- supplied$label
      names(supplied_labels) <- supplied$parameter
      label_source[names(supplied_labels)] <- supplied_labels
    }
  }
  
  group_var <- facet_by
  groups <- split(dat, dat[[group_var]], drop = TRUE)
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  ng <- length(groups)
  nr <- ceiling(sqrt(ng))
  nc <- ceiling(ng / nr)
  
  graphics::par(
    mfrow = c(nr, nc),
    mar = c(4, 8, 3, 1),
    oma = c(0, 0, 2, 0)
  )
  
  point_col <- point
  
  for (g in names(groups)) {
    d <- groups[[g]]
    
    if (facet_by == "parameter") {
      d$plot_label <- d$trait
      y_labels <- d$plot_label
      this_main <- panel_title(g)
    } else {
      d$plot_label <- d$parameter
      y_labels <- make_labels(d$plot_label)
      this_main <- g
    }
    
    d <- d[order(d$plot_label), , drop = FALSE]
    
    if (facet_by == "trait") {
      y_labels <- make_labels(d$plot_label)
    } else {
      y_labels <- d$plot_label
    }
    
    y <- seq_len(nrow(d))
    x <- d[[point_col]]
    
    qlo <- d$q_lower
    qhi <- d$q_upper
    hlo <- d$hpd_lower
    hhi <- d$hpd_upper
    
    if (interval == "hpd") {
      lo <- hlo
      hi <- hhi
    } else if (interval == "quantile") {
      lo <- qlo
      hi <- qhi
    } else {
      lo <- pmin(qlo, hlo, na.rm = TRUE)
      hi <- pmax(qhi, hhi, na.rm = TRUE)
    }
    
    xlim <- range(c(lo, hi, x), finite = TRUE)
    
    if (log_scale) {
      positive <- c(lo, hi, x)
      positive <- positive[is.finite(positive) & positive > 0]
      
      if (length(positive) < 1L) {
        warning("Cannot use log scale for group ", g, "; no positive values.")
        log_now <- FALSE
      } else {
        xlim <- range(positive)
        log_now <- TRUE
      }
    } else {
      log_now <- FALSE
    }
    
    graphics::plot(
      x,
      y,
      type = "n",
      yaxt = "n",
      xlim = xlim,
      ylim = c(0.5, length(y) + 0.5),
      xlab = if (is.null(xlab)) point_col else xlab,
      ylab = "",
      main = if (is.null(main)) this_main else main,
      log = if (log_now) "x" else ""
    )
    
    graphics::axis(
      side = 2,
      at = y,
      labels = y_labels,
      las = 2,
      cex.axis = cex_axis
    )
    
    if (interval %in% c("quantile", "both")) {
      graphics::segments(
        qlo,
        y,
        qhi,
        y,
        lwd = lwd_quantile,
        col = col_quantile
      )
    }
    
    if (interval %in% c("hpd", "both")) {
      graphics::segments(
        hlo,
        y,
        hhi,
        y,
        lwd = lwd_hpd,
        col = col_hpd
      )
    }
    
    graphics::points(
      x,
      y,
      pch = pch,
      col = col_point
    )
    
    graphics::grid(nx = NA, ny = NULL)
    graphics::box()
  }
  
  interval_text <- switch(
    interval,
    hpd = "HPD intervals",
    quantile = "equal-tail intervals",
    both = "equal-tail and HPD intervals"
  )
  
  graphics::mtext(
    paste(point_col, "with", interval_text),
    outer = TRUE,
    line = 0.5
  )
  
  invisible(dat)
}


#' Summarise ST-BLR Posterior Traces
#'
#' Computes posterior summaries for global ST-BLR trace parameters. The
#' function summarizes traces such as `vbs`, `vgs`, `ves`, `pi_trace`, `vle`, and
#' `vld`.
#'
#' For each trace it computes posterior mean, median, posterior standard
#' deviation, Monte Carlo standard error, equal-tail credible intervals, HPD
#' intervals, effective sample size, lag-1 autocorrelation, minimum, and
#' maximum. HPD intervals use `coda::HPDinterval()` when `coda` is available and
#' fall back to equal-tail intervals otherwise.
#'
#' When requested, derived quantities are computed from post-burn-in samples:
#' `h2 = vg / (vg + ve)`, `ve_ratio = ve / (vg + ve)`,
#' `le_ratio = vle / vg`, `ld_ratio = vld / vg`,
#' `varch = vb * pi * m`, `varch_ratio = varch / vg`,
#' `h2_arch = varch / (varch + ve)`, and `m_included = pi * m`.
#' `m_included` is the expected number of included markers and is derived from
#' `pi` and `fit$input$m`.
#'
#' @param fit Fitted ST-BLR object.
#' @param nburn Number of initial iterations to discard. If `NULL`, uses
#'   `fit$input$nburn` when available and otherwise zero.
#' @param traces Character vector of trace components to summarize.
#' @param prob Credible interval probability.
#' @param derived Logical; include derived posterior quantities when the needed
#'   traces are available.
#' @param include_m_included Logical; include `m_included = pi * m` when
#'   `derived = TRUE`, `pi_trace` is available, and `fit$input$m` is present.
#' @param include_diagnostics Logical; include effective sample size and lag-1
#'   autocorrelation diagnostics when possible.
#'
#' @return A data frame with columns `parameter`, `label`, `trait`, `n`,
#'   `mean`, `median`, `sd`, `mcse`, `q_lower`, `q_upper`, `hpd_lower`,
#'   `hpd_upper`, `ess`, `autocorr_lag1`, `min`, and `max`.
#'
#' @examples
#' set.seed(1)
#' fake_fit <- list(
#'   input = list(nburn = 5, m = 100),
#'   vbs = matrix(rnorm(40, mean = 0.005, sd = 0.001), ncol = 1,
#'                dimnames = list(NULL, "trait1")),
#'   vgs = matrix(rnorm(40, mean = 0.6, sd = 0.02), ncol = 1,
#'                dimnames = list(NULL, "trait1")),
#'   ves = matrix(rnorm(40, mean = 0.4, sd = 0.02), ncol = 1,
#'                dimnames = list(NULL, "trait1")),
#'   pi_trace = matrix(runif(40, min = 0.01, max = 0.03), ncol = 1,
#'                dimnames = list(NULL, "trait1")),
#'   vle = matrix(rnorm(40, mean = 0.1, sd = 0.01), ncol = 1,
#'                dimnames = list(NULL, "trait1")),
#'   vld = matrix(rnorm(40, mean = 0.2, sd = 0.01), ncol = 1,
#'                dimnames = list(NULL, "trait1"))
#' )
#' post <- summarise_posterior(fake_fit)
#' post[post$parameter %in% c("vg", "ve", "h2", "pi", "varch"), ]
#'
#' @export
summarise_posterior <- function(
    fit,
    nburn = NULL,
    traces = c("vbs", "vgs", "ves", "pi_trace", "vle", "vld"),
    prob = 0.95,
    derived = TRUE,
    include_m_included = TRUE,
    include_diagnostics = TRUE
) {
  if (is.null(nburn)) {
    nburn <- if (!is.null(fit$input$nburn)) fit$input$nburn else 0L
  }
  nburn <- as.integer(nburn)
  
  if (!is.numeric(prob) || length(prob) != 1L ||
      !is.finite(prob) || prob <= 0 || prob >= 1) {
    stop("prob must be a finite scalar in (0, 1).")
  }
  
  has_coda <- requireNamespace("coda", quietly = TRUE)
  
  trace_parameter <- c(
    vbs = "vb",
    vgs = "vg",
    ves = "ve",
    pi_trace = "pi",
    vle = "vle",
    vld = "vld"
  )
  
  parameter_label <- c(
    vb = "V[b]",
    vg = "V[g]",
    ve = "V[e]",
    pi = "pi",
    vle = "V[LE]",
    vld = "V[LD]",
    h2 = "h^2",
    ve_ratio = "V[e] / (V[g] + V[e])",
    le_ratio = "V[LE] / V[g]",
    ld_ratio = "V[LD] / V[g]",
    varch = "V[arch]",
    varch_ratio = "V[arch] / V[g]",
    h2_arch = "h[arch]^2",
    m_included = "pi * m"
  )
  
  get_label <- function(parameter) {
    label <- unname(parameter_label[parameter])
    ifelse(is.na(label), parameter, label)
  }
  
  as_trace_matrix <- function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  }
  
  post_trace <- function(x) {
    x <- as_trace_matrix(x)
    
    if (nrow(x) < 1L) {
      stop("Trace has zero rows.")
    }
    
    if (nburn >= nrow(x)) {
      stop("nburn must be smaller than the number of trace rows.")
    }
    
    x[seq.int(nburn + 1L, nrow(x)), , drop = FALSE]
  }
  
  hpd_interval <- function(x) {
    alpha <- (1 - prob) / 2
    q <- stats::quantile(
      x,
      probs = c(alpha, 1 - alpha),
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )
    
    fallback <- c(lower = q[1], upper = q[2])
    
    if (!has_coda || length(x) < 2L || stats::sd(x, na.rm = TRUE) == 0) {
      return(fallback)
    }
    
    out <- tryCatch(
      coda::HPDinterval(coda::mcmc(x), prob = prob),
      error = function(e) NULL
    )
    
    if (is.null(out)) {
      return(fallback)
    }
    
    c(lower = unname(out[1, "lower"]), upper = unname(out[1, "upper"]))
  }
  
  ess_value <- function(x) {
    if (!include_diagnostics || !has_coda ||
        length(x) < 2L || stats::sd(x, na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    
    tryCatch(
      as.numeric(coda::effectiveSize(coda::mcmc(x))),
      error = function(e) NA_real_
    )
  }
  
  autocorr_lag1_value <- function(x) {
    if (!include_diagnostics ||
        length(x) < 2L || stats::sd(x, na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    
    tryCatch(
      as.numeric(stats::acf(x, lag.max = 1, plot = FALSE)$acf[2]),
      error = function(e) NA_real_
    )
  }
  
  summarise_matrix <- function(mat, parameter) {
    mat <- as_trace_matrix(mat)
    
    trait_names <- colnames(mat)
    if (is.null(trait_names)) {
      trait_names <- paste0("T", seq_len(ncol(mat)))
    }
    
    alpha <- (1 - prob) / 2
    label <- get_label(parameter)
    
    rows <- lapply(seq_len(ncol(mat)), function(j) {
      x <- mat[, j]
      x <- x[is.finite(x)]
      
      if (length(x) < 1L) {
        return(data.frame(
          parameter = parameter,
          label = label,
          trait = trait_names[j],
          n = 0L,
          mean = NA_real_,
          median = NA_real_,
          sd = NA_real_,
          mcse = NA_real_,
          q_lower = NA_real_,
          q_upper = NA_real_,
          hpd_lower = NA_real_,
          hpd_upper = NA_real_,
          ess = NA_real_,
          autocorr_lag1 = NA_real_,
          min = NA_real_,
          max = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      
      sx <- stats::sd(x)
      ess <- ess_value(x)
      autocorr_lag1 <- autocorr_lag1_value(x)
      
      if (is.na(ess) || ess <= 0) {
        mcse <- sx / sqrt(length(x))
      } else {
        mcse <- sx / sqrt(ess)
      }
      
      q <- stats::quantile(
        x,
        probs = c(alpha, 1 - alpha),
        na.rm = TRUE,
        names = FALSE,
        type = 8
      )
      
      hpd <- hpd_interval(x)
      
      data.frame(
        parameter = parameter,
        label = label,
        trait = trait_names[j],
        n = length(x),
        mean = mean(x),
        median = stats::median(x),
        sd = sx,
        mcse = mcse,
        q_lower = q[1],
        q_upper = q[2],
        hpd_lower = hpd[1],
        hpd_upper = hpd[2],
        ess = ess,
        autocorr_lag1 = autocorr_lag1,
        min = min(x),
        max = max(x),
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, rows)
  }
  
  safe_divide <- function(num, den) {
    out <- num / den
    out[!is.finite(out) | abs(den) <= .Machine$double.eps] <- NA_real_
    out
  }
  
  out <- list()
  post <- list()
  
  for (nm in traces) {
    if (is.null(fit[[nm]])) next
    
    mat <- post_trace(fit[[nm]])
    post[[nm]] <- mat
    
    parameter <- unname(trace_parameter[nm])
    if (is.na(parameter)) parameter <- nm
    
    out[[parameter]] <- summarise_matrix(
      mat = mat,
      parameter = parameter
    )
  }
  
  if (isTRUE(derived)) {
    if (!is.null(post$vgs) && !is.null(post$ves)) {
      denom <- post$vgs + post$ves
      
      out$h2 <- summarise_matrix(
        safe_divide(post$vgs, denom),
        parameter = "h2"
      )
      
      out$ve_ratio <- summarise_matrix(
        safe_divide(post$ves, denom),
        parameter = "ve_ratio"
      )
    }
    
    if (!is.null(post$vle) && !is.null(post$vgs)) {
      out$le_ratio <- summarise_matrix(
        safe_divide(post$vle, post$vgs),
        parameter = "le_ratio"
      )
    }
    
    if (!is.null(post$vld) && !is.null(post$vgs)) {
      out$ld_ratio <- summarise_matrix(
        safe_divide(post$vld, post$vgs),
        parameter = "ld_ratio"
      )
    }
    
    if (!is.null(post$vbs) &&
        !is.null(post$pi_trace) &&
        !is.null(fit$input$m)) {
      m <- as.numeric(fit$input$m)
      varch <- post$vbs * post$pi_trace * m
      
      out$varch <- summarise_matrix(
        varch,
        parameter = "varch"
      )
      
      if (!is.null(post$vgs)) {
        out$varch_ratio <- summarise_matrix(
          safe_divide(varch, post$vgs),
          parameter = "varch_ratio"
        )
      }
      
      if (!is.null(post$ves)) {
        out$h2_arch <- summarise_matrix(
          safe_divide(varch, varch + post$ves),
          parameter = "h2_arch"
        )
      }
      
      if (isTRUE(include_m_included)) {
        out$m_included <- summarise_matrix(
          post$pi_trace * m,
          parameter = "m_included"
        )
      }
    }
  }
  
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  
  ans
}


plot_stblr_traces <- function(
    fit,
    traces = c("vgs", "ves", "vbs", "pi_trace", "vle", "vld"),
    traits = NULL,
    nburn = NULL,
    max_cols = 2
) {
  if (is.null(nburn)) {
    nburn <- if (!is.null(fit$input$nburn)) fit$input$nburn else 0L
  }
  
  available <- traces[!vapply(traces, function(x) is.null(fit[[x]]), logical(1))]
  
  if (length(available) < 1L) {
    stop("No requested traces were found in fit.")
  }
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  
  nr <- ceiling(length(available) / max_cols)
  nc <- min(max_cols, length(available))
  
  graphics::par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
  
  for (nm in available) {
    x <- as.matrix(fit[[nm]])
    
    if (!is.null(traits)) {
      keep <- colnames(x) %in% traits
      x <- x[, keep, drop = FALSE]
    }
    
    graphics::matplot(
      x,
      type = "l",
      lty = 1,
      xlab = "Iteration",
      ylab = nm,
      main = nm
    )
    
    if (!is.null(nburn) && nburn > 0) {
      graphics::abline(v = nburn, lty = 2)
    }
    
    if (ncol(x) > 1L) {
      graphics::legend(
        "topright",
        legend = colnames(x),
        lty = 1,
        col = seq_len(ncol(x)),
        bty = "n",
        cex = 0.8
      )
    }
  }
  
  invisible(available)
}
