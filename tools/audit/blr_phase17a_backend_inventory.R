root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)

inventory <- data.frame(
  route = c(
    "CSR BayesC", "scheduled CSR BayesC", "CSR BayesR", "CSR SBayesRC",
    "fixed-prior CSR BayesC", "group CSR BayesC", "learned-annotation CSR BayesC",
    "block-eigen CSR BayesC", "block-eigen CSR BayesR", "block-eigen CSR SBayesRC",
    "sblr default mtblr", "sblr cpg", "sblr cpg_arma", "sblr cpg_omp",
    "sblr eigen", "mtblr_hybrid", "mtblr_cpg_omp_csr"),
  family = c(rep("scalar CSR", 7), rep("block-eigen", 3), rep("multivariate", 7)),
  r_route = c("stblr_csr", "stblr_csr(scheduled=TRUE)", "stblr_csr(method='bayesr')",
    "stblr_csr_annot(annotation_model='sbayesrc')", "stblr_csr_prior_annot",
    "stblr_csr_group_annot", "stblr_csr_learn_annot",
    ".stblr_csr_bayesc_block_eigen", ".stblr_csr_bayesr_block_eigen",
    ".stblr_csr_sbayesrc_block_eigen", "sblr(algorithm='default')",
    "sblr(algorithm='cpg')", "sblr(algorithm='cpg_arma')",
    "sblr(algorithm='cpg_omp')", "sblr(algorithm='eigen')", NA, NA),
  native = c("stblr_cpg_omp_csr", "stblr_cpg_omp_csr_scheduled",
    "stblr_cpg_omp_csr_bayesr", "stblr_cpg_omp_csr_sbayesrc",
    "stblr_cpg_omp_csr_prior", "stblr_cpg_omp_csr_group_annot",
    "stblr_cpg_omp_csr_annot", "stblr_cpg_omp_csr_block_eigen",
    "stblr_cpg_omp_csr_bayesr_block_eigen", "stblr_cpg_omp_csr_sbayesrc_block_eigen",
    "mtblr", "mtblr_cpg", "mtblr_cpg_arma", "mtblr_cpg_omp", "mtblr_eigen",
    "mtblr_hybrid", "mtblr_cpg_omp_csr"),
  reachability = c(rep("supported public", 7), rep("internal research", 3),
    rep("supported public legacy", 5), "native-only", "native-only research"),
  maturity = c(rep("canonical architecture", 7), rep("audited but noncanonical", 3),
    rep("legacy architecture", 7)),
  rng = c(rep("logical-chain safe", 7), rep("inherits canonical CSR core", 3),
    "single-engine; R-generated seed", "single-engine; R-generated seed",
    "single-engine; R-generated seed", "worker-sensitive risk",
    "single-engine; R-generated seed", "single-engine; native-only",
    "single-engine; native-only"),
  references = c("strong permanent references", "strong permanent references",
    "strong permanent references", "strong permanent references",
    "strong permanent references", "strong permanent references",
    "strong permanent references", rep("smoke test only", 3),
    rep("smoke test only", 5), "no deterministic reference", "no deterministic reference"),
  disposition = c(rep("canonical - no action", 7), rep("retain experimental", 3),
    "retain supported - migrate", "retain experimental", "retain experimental",
    "correct defect before migration", "retain experimental",
    "retire/remove candidate", "retain experimental"),
  priority = c(rep("P2", 7), rep("P2", 3), "P1", "P3", "P3", "P0", "P3", "P4", "P2"),
  stringsAsFactors = FALSE)

exports <- readLines(file.path(root, "R", "RcppExports.R"), warn = FALSE)
inventory$generated_wrapper <- vapply(inventory$native, function(x)
  any(grepl(paste0("^", x, " <- function"), exports)), logical(1))

print(inventory, row.names = FALSE)
if (length(commandArgs(trailingOnly = TRUE))) {
  utils::write.csv(inventory, commandArgs(trailingOnly = TRUE)[1], row.names = FALSE)
}
