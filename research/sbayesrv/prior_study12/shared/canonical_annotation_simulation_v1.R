# Canonical annotation simulation manifest.
#
# This object points to the immutable Study 10 semantic checkpoint. It does
# not copy or regenerate the large genotype, phenotype, LD, or truth objects.
canonical_annotation_simulation_v1 <- list(
  id = "canonical_annotation_simulation_v1",
  schema = "sblrbench-canonical-simulation-v1",
  owner_study = "10_moderate_sparsity_sbayesrc",
  checkpoint = file.path("results", "local",
    "10_moderate_sparsity_sbayesrc", "checkpoints", "truth.rds"),
  checkpoint_sha256 =
    "97786396a3e5e297eaed1dacfda8fd58f2c4284a999c324397fac98b3f690d07",
  specification_hash =
    "cb01b88e4e7ebf87b5ace27b7b80d64e6414498f1242dff340e27a89874dadc8",
  truth_hash =
    "44352063d8b0527a19dcff3fefd97b194098d3fc00d7272cbb122720f67d2d4c",
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
  N_train = 4500L,
  N_validation = 500L,
  M = 35000L,
  block_count = 35L,
  block_size = 1000L,
  target_h2 = 0.5,
  target_component_counts = c(null = 34000L, small = 500L,
    medium = 300L, large = 200L),
  realized_component_counts = c(null = 33989L, small = 502L,
    medium = 305L, large = 204L),
  realized_active_count = 1011L,
  realized_h2 = 0.5,
  annotations = c("Intercept", "enriched_binary", "continuous_signal",
    "null_annotation"),
  seeds = c(component = 1005412L, effect = 1005424L,
    residual = 1005438L),
  immutable = TRUE,
  regeneration_allowed = FALSE
)
