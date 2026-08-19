# Deterministic, Study-06-independent calibration of the default probit-stick
# intercept scale. This script reads no benchmark data and writes no files.

candidates <- c(0.5, 1.0, 1.5, 2.5)
architectures <- list(
  balanced = rep(0.25, 4),
  sparse_equal = c(0.999, rep(0.001 / 3, 3)),
  moderate = c(0.95, 0.02, 0.02, 0.01),
  skewed_active = c(0.98, 0.001, 0.004, 0.015)
)
practical_boundary <- 1e-8

component_to_stick <- function(pi) {
  remaining <- rev(cumsum(rev(pi / sum(pi))))
  remaining[-1] / remaining[-length(remaining)]
}

rows <- lapply(names(architectures), function(name) {
  p <- component_to_stick(architectures[[name]])
  mu <- qnorm(p)
  do.call(rbind, lapply(candidates, function(s) {
    lower <- pnorm(qnorm(0.025, mu, s))
    upper <- pnorm(qnorm(0.975, mu, s))
    degenerate <- pnorm(qnorm(practical_boundary), mu, s) +
      1 - pnorm(qnorm(1 - practical_boundary), mu, s)
    data.frame(
      architecture = name,
      sd = s,
      maximum_degenerate_mass = max(degenerate),
      minimum_95_width = min(upper - lower),
      maximum_95_width = max(upper - lower)
    )
  }))
})
summary <- do.call(rbind, rows)

# Preregistered rule: choose the largest (least informative) candidate whose
# practically-degenerate mass is at most 1% for every representative stick.
by_sd <- aggregate(maximum_degenerate_mass ~ sd, summary, max)
eligible <- by_sd$sd[by_sd$maximum_degenerate_mass <= 0.01]
selected <- if (length(eligible)) max(eligible) else NA_real_

print(summary, row.names = FALSE)
cat("\nSelection rule: largest candidate with maximum degenerate mass <= 0.01\n")
cat("Selected SD:", selected, "\n")

