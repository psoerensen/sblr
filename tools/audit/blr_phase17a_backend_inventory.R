root <- normalizePath(if (file.exists("DESCRIPTION")) "." else file.path("..", ".."),
  winslash = "/", mustWork = TRUE)

inventory <- data.frame(
  route = c(
    "CSR BayesC", "scheduled CSR BayesC", "CSR BayesR", "CSR SBayesRC",
    "fixed-prior CSR BayesC", "group CSR BayesC", "learned-annotation CSR BayesC",
    "block-eigen CSR BayesC", "block-eigen CSR BayesR", "block-eigen CSR SBayesRC",
    "sblr default mtblr", "mtblr eigen research", "mtblr CSR research"),
  family = c(rep("scalar CSR", 7), rep("block-eigen", 3), rep("multivariate", 3)),
  r_route = c("stblr_csr", "stblr_csr(scheduled=TRUE)", "stblr_csr(method='bayesr')",
    "stblr_csr_annot(annotation_model='sbayesrc')", "stblr_csr_prior_annot",
    "stblr_csr_group_annot", "stblr_csr_learn_annot",
    ".stblr_csr_bayesc_block_eigen", ".stblr_csr_bayesr_block_eigen",
    ".stblr_csr_sbayesrc_block_eigen", "sblr(algorithm='default')", NA, NA),
  native = c("stblr_cpg_omp_csr", "stblr_cpg_omp_csr_scheduled",
    "stblr_cpg_omp_csr_bayesr", "stblr_cpg_omp_csr_sbayesrc",
    "stblr_cpg_omp_csr_prior", "stblr_cpg_omp_csr_group_annot",
    "stblr_cpg_omp_csr_annot", "stblr_cpg_omp_csr_block_eigen",
    "stblr_cpg_omp_csr_bayesr_block_eigen", "stblr_cpg_omp_csr_sbayesrc_block_eigen",
    "mtblr", "mtblr_eigen", "mtblr_cpg_omp_csr"),
  reachability = c(rep("supported public", 7), rep("internal research", 3),
    "supported public legacy", "native-only research", "native-only research"),
  maturity = c(rep("canonical architecture", 7), rep("audited but noncanonical", 3),
    "typed noncanonical architecture", "legacy research", "legacy research"),
  references = c(rep("strong permanent references", 7), rep("smoke test only", 3),
    "strong permanent references", "no deterministic reference", "no deterministic reference"),
  disposition = c(rep("canonical - no action", 7), rep("retain experimental", 3),
    "authoritative supported public", "retain internal research", "retain internal research"),
  stringsAsFactors = FALSE)

exports <- readLines(file.path(root, "R", "RcppExports.R"), warn = FALSE)
inventory$generated_wrapper <- vapply(inventory$native, function(x)
  any(grepl(paste0("^", x, " <- function"), exports)), logical(1))

historical_disposition <- data.frame(
  route = c("mtblr_cpg", "mtblr_cpg_arma", "mtblr_cpg_omp", "mtblr_hybrid"),
  phase17g_action = rep("retired and removed", 4),
  stringsAsFactors = FALSE)

print(inventory, row.names = FALSE)
if (length(commandArgs(trailingOnly = TRUE)))
  utils::write.csv(inventory, commandArgs(trailingOnly = TRUE)[1], row.names = FALSE)
