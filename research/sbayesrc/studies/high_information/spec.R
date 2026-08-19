# Study 08: one-replicate high-information SBayesRC annotation positive control.
# Scientific inputs are frozen before any method fit. Local outputs are
# review-pending evidence and are not promoted by this study.

spec <- list(
  study = "06c_high_information",
  former_study = "08_high_information_sbayesrc",
  task = "high_information_annotation_positive_control",
  title = "Study 06C: High-information SBayesRC annotation positive control",
  study_version = "v1.1_one_replicate_operator_qualified",
  status = "design_frozen_execution_pending",
  design_amendments = c(
    calibration = paste(
      "Before any allocation or effect draw, replace ill-conditioned joint",
      "BFGS calibration with the exact triangular sequential uniroot solution",
      "implied by continuation-stick probabilities."),
    operator_qualification = paste(
      "Before any successful fit or inspection of method output, align the",
      "executable residual flag with the pre-specified fixVe block contract:",
      "updateE=FALSE for all methods. The rejected updateE=TRUE qualification",
      "failed at SBayesR iteration zero with a negative CSR residual scale.")),
  scope = list(
    simulated_data_replicates = 1L,
    factorial_benchmark = FALSE,
    population_calibration_claim = FALSE,
    causal_N_or_M_claim = FALSE),
  source = list(
    study06 = "06_annotation_models",
    study07 = "07_joint_em_sbayesrc",
    qgdata_sha = "6cca5819e711d326cfb2614d7e9d9f34942612cd",
    glist = file.path("results", "local", "06_annotation_models",
      "checkpoints", "data", "human_glist.rds"),
    qc_marker_glist = file.path("results", "local",
      "06_annotation_models", "checkpoints", "ld",
      "training_ld_train1400-test600-seed3101_glist.rds")),
  data = list(
    chromosome = 1L, sample_count = 5000L,
    training_count = 4500L, validation_count = 500L,
    trait = "trait1",
    split = list(train_fraction = 0.90, seed = 3801L),
    marker_target = 35000L,
    marker_selection = paste(
      "first 35 complete source-order 1000-marker blocks from the",
      "37,991 QC-qualified chromosome-1 marker order"),
    block_size = 1000L, expected_block_count = 35L,
    sparse_ld = list(r2_threshold = 0.001,
      rule = "training correlation within registered blocks; upper triangle; implicit unit diagonal; no cross-block edges")),
  annotation = list(
    columns = c("Intercept", "enriched_binary", "continuous_signal",
      "null_annotation"),
    seed = 8201L, enriched_fraction = 0.15,
    standardization = "sample mean zero and sample SD one",
    maximum_absolute_pairwise_correlation = 0.10,
    informative_nonintercept_alpha = matrix(c(
      1.60, .30, .20,
       .30, .15, .10,
       .00, .00, .00), nrow = 3L, byrow = TRUE,
      dimnames = list(c("enriched_binary", "continuous_signal",
        "null_annotation"), paste0("step_", 1:3)))),
  simulation = list(
    h2 = 0.50,
    gamma = c(0, 0.01, 0.1, 1),
    target_expected_component_counts = c(null = 34820, small = 90,
      medium = 54, large = 36),
    study06_expected_active_component_counts = c(small = 90,
      medium = 54, large = 36),
    study06_realized_component_counts = c(null = 1329, small = 84,
      medium = 50, large = 37),
    calibration = list(method =
      "deterministic sequential continuation-stick uniroot",
      tolerance = 1e-8, root_interval = c(-20, 20)),
    realized_active_sanity = c(150L, 215L)),
  methods = c("SBayesR", "SBayesRC", "SBayesRC-S"),
  operator = list(
    representation = "sparse_ld_csr",
    residual_policy = "global_projected_legacy",
    block_ve_mode = "fixVe",
    updateB = TRUE, updateE = FALSE,
    common_summary_statistics = TRUE,
    common_marker_order = TRUE,
    common_residual_contract = TRUE),
  priors = list(
    nub = 4, nue = 4,
    sigmaSqAlpha_a = 2, sigmaSqAlpha_b = 2,
    sigmaSqAlpha_init = c(step_1 = 1, step_2 = 1, step_3 = 1),
    intercept_prior = list(distribution = "normal", sd = 1),
    selection = list(
      delta_init = c(enriched_binary = 1L, continuous_signal = 1L,
        null_annotation = 1L),
      pi_A_init = 0.30, tau2_init = c(step_1 = 1, step_2 = 1,
        step_3 = 1),
      pi_A_prior = c(a = 1, b = 1),
      tau2_prior = c(a = 3, b = 1.6))),
  mcmc = list(
    nit = 6000L, nburn = 3000L, nthin = 1L,
    nchains = 4L, ncores = 4L,
    alpha_update_every = 1L,
    convergence = list(rhat = 1.01, ess_bulk = 400,
      ess_tail = 400, relative_mcse = 0.05, chain_count = 4L)),
  seeds = list(
    marker_selection = 8301L,
    component = 805412L, effect = 805424L, residual = 805438L,
    baseline = 806040L,
    baseline_chains = c(806141L, 806242L, 806343L, 806444L),
    continuous = 807040L,
    continuous_chains = c(807141L, 807242L, 807343L, 807444L),
    selection = 808040L,
    selection_chains = c(808141L, 808242L, 808343L, 808444L)),
  packages = list(
    sblr = list(version = "0.2.0",
      sha = "01196e3a9ca65e7fbf60650fd285576d5d2b25d8",
      installed_tree_sha256 =
        "86ab17e019534ef851be7d2b4f097f55e12327f65cfdaa62da0257ac868a9d93")),
  outputs = list(
    local_dir = file.path("results", "local",
      "08_high_information_sbayesrc"),
    capsule_promoted = FALSE),
  decision_screen = list(
    maximum_nonintercept_rhat = 1.01,
    minimum_nonintercept_bulk_ess = 400,
    maximum_nonintercept_alpha_rmse = 0.30,
    minimum_component_prior_correlation = 0.98,
    minimum_active_prior_correlation = 0.95,
    minimum_chain_active_prior_pearson = 0.98,
    signal_annotation_pip = 0.80,
    null_annotation_pip = 0.20,
    maximum_genomic_beta_degradation = 0.05),
  decisions = c("STUDY08-R1", "STUDY08-R2", "STUDY08-R3", "STUDY08-R4")
)
