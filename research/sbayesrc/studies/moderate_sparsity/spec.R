# Study 10: one-replicate moderate-sparsity SBayesRC annotation positive
# control. This specification is a controlled layer over frozen Study 08.

.study10_root <- if (file.exists("DESCRIPTION")) "." else
  file.path("..", "..", "..")
.study10_study08_spec <- source(file.path(.study10_root, "studies",
  "06_annotation_models", "06c_high_information", "spec.R"),
  local = TRUE)$value

spec <- .study10_study08_spec
spec$study <- "06d_moderate_sparsity"
spec$former_study <- "10_moderate_sparsity_sbayesrc"
spec$task <- "moderate_sparsity_annotation_positive_control"
spec$title <- "Study 06D: Moderate-sparsity SBayesRC annotation positive control"
spec$study_version <- "v1_one_replicate_matched_study08"
spec$status <- "design_frozen_execution_pending"
spec$design_amendments <- NULL
spec$scope <- list(
  simulated_data_replicates = 1L,
  factorial_benchmark = FALSE,
  calibration_claim = FALSE,
  sample_size_claim = FALSE,
  marker_count_claim = FALSE,
  active_count_isolated_causal_claim = FALSE,
  fixed_h2_effect_size_confounding = paste(
    "At fixed h2, increasing active count necessarily reduces average",
    "per-active SNP effect magnitude."))
spec$source$study08 <- "08_high_information_sbayesrc"
spec$source$study08_truth <- file.path("results", "local",
  "08_high_information_sbayesrc", "checkpoints", "truth.rds")
spec$shared_input_pins <- list(
  study08_specification_hash =
    "0eabd79cf8f780bb00868177254d68109aa8a56dec98127d0d71abca87bc8d32",
  study08_truth_hash =
    "f8ff0fa329a0776af09c84b504a201c295624b775dbe01b9fd72edd81f013528",
  study08_truth_file_sha256 =
    "17b7e04e547683eed4cc9dec04abb3c65c56091d2bffe5917f2b4e64aaf469ce",
  marker_hash =
    "0dfaa465f0ed80517552fa0c73410f0e7d40983c2ea7f8082828652fa4dbb155",
  block_hash =
    "96a95ee80cb6150617acad65caa5a6553af18732d4b3b0ac71c476b63d70c5ed",
  training_sample_hash =
    "77aa28dee26d71574148f4805cf902e62d33d96741f41b951dea9086e359e78c",
  validation_sample_hash =
    "31e482448eff22dd5894b48fc253006c9eb4664d3ce5f71939592389763ddc23",
  annotation_hash =
    "cdf7f3edd04bfcb3b6c3e5943e976bc3ce378d76e13f8827786125da2a59f7a1",
  training_af_hash =
    "d8f65f3eaed2d25ddc7037fa4e0b9f6c5601ada24c6a6680ac7e307b203d021e",
  working_glist_hash =
    "e6a76b829e493f74fae537f58cb039eb06ef6d5c9e554341f327924f44039985",
  csr_hashes = c(
    "1df5e86fd1236b4d1a17ae798a00e63d85ea430c915db203b43a6ca4c77d0c7f",
    "09cac8c046bdcb1f0b80945f5048236c8aa8a9dd627a642db653acc7b5da48a5",
    "e1781f1de06c2709ad74b040335cdf21d8a10da49071f3f6161323645e85798f",
    "4a875cad2341c393612ce026fb4b368b582817debbb545e7921f1d10a728fb9f",
    "cca657897532aaade9353f960da78f8ec942080b6c16f89fe771446b1bf932b2"))
spec$study08_evidence_pins <- c(
  alpha_recovery.csv =
    "0be785c3c7418b8efce5132f04c6653a4b35b5d68bb7412010d71b6ab748833e",
  prior_recovery.csv =
    "af8105d9d6e012b7de3a4a95d113ec9b1514bdac344c00cc63562859473e04cc",
  annotation_selection.csv =
    "96a40200e2a977eb65ba9c067082ea4251aa0d6908d916d225953fb75b1953fd",
  annotation_selection_chains.csv =
    "f4bd9c6ab647a8a6b78ba1a10460bb64fe6db8d176e75cff6e83409fc7ca85a9",
  genomic_summary.csv =
    "48441ec0cdb262738abac7fd0d2a6cc67238d3c0d86a92493b86f45bf4d7138f")
spec$simulation$target_expected_component_counts <- c(
  null = 34000, small = 500, medium = 300, large = 200)
spec$simulation$study08_target_expected_component_counts <- c(
  null = 34820, small = 90, medium = 54, large = 36)
spec$simulation$expected_active_fraction <- 1000 / 35000
spec$simulation$realized_active_sanity <- c(850L, 1150L)
spec$simulation$minimum_realized_active_counts <- c(
  small = 350L, medium = 200L, large = 120L)
spec$seeds <- list(
  marker_selection = 8301L,
  component = 1005412L, effect = 1005424L, residual = 1005438L,
  baseline = 1006040L,
  baseline_chains = c(1006141L, 1006242L, 1006343L, 1006444L),
  continuous = 1007040L,
  continuous_chains = c(1007141L, 1007242L, 1007343L, 1007444L),
  selection = 1008040L,
  selection_chains = c(1008141L, 1008242L, 1008343L, 1008444L))
spec$outputs <- list(local_dir = file.path("results", "local",
  "10_moderate_sparsity_sbayesrc"), capsule_promoted = FALSE)
spec$decision_screen <- list(
  maximum_nonintercept_rhat = 1.10,
  minimum_nonintercept_bulk_ess = 100,
  minimum_nonintercept_tail_ess = 100,
  maximum_nonintercept_alpha_rmse = 0.30,
  minimum_nonintercept_alpha_correlation = 0.85,
  minimum_component_prior_correlation = 0.98,
  minimum_active_prior_correlation = 0.95,
  minimum_chain_active_prior_pearson = 0.90,
  minimum_chain_active_prior_spearman = 0.85,
  signal_annotation_pip = 0.80,
  null_annotation_pip = 0.20,
  maximum_null_annotation_chain_pip = 0.20,
  maximum_null_annotation_chain_range = 0.10,
  minimum_signal_annotation_chain_pip = 0.80,
  maximum_genomic_beta_degradation = 0.05)
spec$decisions <- c("STUDY10-R1", "STUDY10-R2", "STUDY10-R3",
  "STUDY10-R4")

spec
