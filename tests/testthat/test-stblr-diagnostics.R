fake_stblr_fit <- function() {
  set.seed(1)
  trait <- "trait1"
  list(
    input = list(nburn = 5L, m = 100L),
    vbs = matrix(rnorm(40, 0.005, 0.001), ncol = 1,
                 dimnames = list(NULL, trait)),
    vgs = matrix(rnorm(40, 0.6, 0.02), ncol = 1,
                 dimnames = list(NULL, trait)),
    ves = matrix(rnorm(40, 0.4, 0.02), ncol = 1,
                 dimnames = list(NULL, trait)),
    pis = matrix(runif(40, 0.01, 0.03), ncol = 1,
                 dimnames = list(NULL, trait)),
    vle = matrix(rnorm(40, 0.1, 0.01), ncol = 1,
                 dimnames = list(NULL, trait)),
    vld = matrix(rnorm(40, 0.2, 0.01), ncol = 1,
                 dimnames = list(NULL, trait))
  )
}

test_that("summarise_posterior returns expected columns", {
  post <- sblr::summarise_posterior(
    fake_stblr_fit(),
    include_diagnostics = FALSE
  )

  expect_named(
    post,
    c(
      "parameter", "label", "trait", "n", "mean", "median", "sd", "mcse",
      "q_lower", "q_upper", "hpd_lower", "hpd_upper", "ess",
      "autocorr_lag1", "min", "max"
    )
  )
})

test_that("summarise_posterior includes expected marker count", {
  post <- sblr::summarise_posterior(
    fake_stblr_fit(),
    include_m_included = TRUE,
    include_diagnostics = FALSE
  )

  expect_true("m_included" %in% post$parameter)
})

test_that("summarise_posterior uses internal parameter names", {
  post <- sblr::summarise_posterior(
    fake_stblr_fit(),
    include_diagnostics = FALSE
  )

  expect_true("vg" %in% post$parameter)
  expect_false("V_g" %in% post$parameter)
})

test_that("check_stblr_convergence returns diagnostics list", {
  conv <- sblr::check_stblr_convergence(
    fake_stblr_fit(),
    traces = c("vgs", "ves", "pis"),
    require = c("vgs", "ves"),
    use_coda = FALSE
  )

  expect_type(conv, "list")
  expect_named(conv, c("passed", "diagnostics", "nburn", "crit", "required"))
  expect_true("diagnostics" %in% names(conv))
  expect_true("passed" %in% names(conv))
  expect_s3_class(conv$diagnostics, "data.frame")
})

test_that("plot_posterior invisibly returns filtered data", {
  post <- sblr::summarise_posterior(
    fake_stblr_fit(),
    include_diagnostics = FALSE
  )
  post <- post[post$parameter %in% c("vg", "ve", "h2"), , drop = FALSE]

  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted <- sblr::plot_posterior(
    post,
    parameters = c("vg", "ve"),
    facet_by = "trait"
  )

  expect_s3_class(plotted, "data.frame")
  expect_setequal(plotted$parameter, c("vg", "ve"))
})
