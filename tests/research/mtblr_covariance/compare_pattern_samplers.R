script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
research_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L]), mustWork = TRUE))
} else normalizePath("tests/research/mtblr_covariance", mustWork = TRUE)

for (file in c("mtblr_reference_model.R", "mtblr_exact_reference.R",
               "mtblr_pattern_reference.R", "mtblr_pattern_samplers.R",
               "mtblr_regional.R", "compare_samplers.R")) {
  source(file.path(research_dir, file), local = FALSE)
}

run_pattern_comparison <- function(n_iter = 7000L, burn = 1200L,
                                   seeds = c(611L, 612L)) {
  f <- fixture(2L)
  patterns <- mt_pattern_space()
  prior <- c(6, 1, 1, 3)
  nu0 <- 7
  Psi0 <- mt_iw_scale_from_mean(matrix(c(.30, .08, .08, .25), 2L), nu0)
  exact <- mt_pattern_configuration_reference(
    f$X, f$Y, f$Ve, f$Vb, patterns, dirichlet_prior = prior
  )
  rows <- list()
  fits <- list()
  q <- 0L
  for (seed in seeds) for (augmentation in c("full", "completed_active")) {
    for (transition in c("joint", "coordinate")) {
      fit <- mtblr_pattern_sampler(
        f$X, f$Y, f$Ve, patterns, prior, nu0, Psi0,
        n_iter, burn, seed, augmentation, transition
      )
      summary <- mt_pattern_summary(fit)
      pleio <- fit$state[, 1L] == which(rowSums(patterns) == 2L)
      diagnostic <- mt_chain_diagnostics(as.numeric(pleio))
      q <- q + 1L
      rows[[q]] <- data.frame(
        augmentation = augmentation, transition = transition, seed = seed,
        pip1_trait1 = summary$trait_pip[1L, 1L],
        pip1_trait2 = summary$trait_pip[1L, 2L],
        pleiotropic1 = summary$pleiotropic_probability[1L],
        alpha11 = summary$alpha_mean[1L, 1L],
        alpha12 = summary$alpha_mean[1L, 2L],
        Vb11 = summary$Vb_mean[1L, 1L], Vb12 = summary$Vb_mean[1L, 2L],
        pi_null = summary$Pi_mean[1L], pi_both = summary$Pi_mean[4L],
        pleio_lag1 = diagnostic$lag1, pleio_ess = diagnostic$ess
      )
      fits[[paste(augmentation, transition, seed, sep = "_")]] <- fit
    }
  }
  rows <- do.call(rbind, rows)

  # A fixed highly imbalanced Pi isolates state-transition mixing.
  imbalanced <- c(.49, .01, .01, .49)
  mixing <- lapply(c("joint", "coordinate"), function(transition) {
    fit <- mtblr_pattern_sampler(
      f$X[, 1L, drop = FALSE], f$Y, f$Ve, patterns, rep(1, 4),
      nu0, Psi0, 9000L, 1500L, 901L,
      "completed_active", transition, f$Vb, imbalanced,
      update_vb = FALSE, update_pi = FALSE
    )
    pleio <- as.numeric(fit$state[, 1L] == 4L)
    d <- mt_chain_diagnostics(pleio)
    data.frame(transition = transition, pleio = mean(pleio),
               lag1 = d$lag1, ess = d$ess,
               null_both_moves = fit$transitions[1L, 4L] + fit$transitions[4L, 1L])
  })

  list(exact = exact, chains = rows,
       means = aggregate(. ~ augmentation + transition,
                         rows[, setdiff(names(rows), "seed")], mean),
       mixing = do.call(rbind, mixing), fits = fits)
}

if (sys.nframe() == 0L) {
  result <- run_pattern_comparison()
  cat("Exact posterior mean Pi\n")
  print(signif(result$exact$pi_mean, 5))
  cat("\nSampler means\n")
  print(result$means, digits = 5, row.names = FALSE)
  cat("\nJoint versus coordinate mixing under imbalanced Pi\n")
  print(result$mixing, digits = 5, row.names = FALSE)
}
