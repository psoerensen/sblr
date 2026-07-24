args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], mustWork = TRUE) else
  normalizePath(".", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse = "\n")
count <- function(pattern, text) lengths(regmatches(text, gregexpr(pattern, text, perl = TRUE)))
flag <- function(name, value) cat(name, "=", toupper(as.character(isTRUE(value))), "\n", sep = "")

packed <- read("src/packed_bed.h")
bayesc <- read("src/st_cpg_omp_individual_scheduled_chains.cpp")
bayesr <- read("src/st_bed_bayesr_common.h")
family <- read("src/blr_bed_family_types.h")
bayesrc_types <- read("src/blr_bed_bayesrc_types.h")
contract <- read("docs/dev/blr_mt_bed_internal_contract.md")
namespace <- read("NAMESPACE")

cat("OWNERS=PackedBedMatrix,FastPackedBedMatrix,FastPackedBedMatrixBR\n")
cat("READERS=read_bedfiles_to_packed_matrix,read_bedfiles_to_fast_packed_matrix_blocked,br_read_bed_blocked\n")
cat("VIEWS=BedPackedGenotypeView,BedBayesRCPackedGenotypeView\n")
cat("PACKED_OWNER_FORMULA=m*64*ceiling(ceiling(n/4)/64)\n")
cat("PHENOTYPE_FORMULA=8*n*nt\n")
cat("SAMPLE_RESIDUAL_FORMULA=8*n*nt\n")
cat("DECODED_MARKER_FORMULA=8*n\n")
cat("EFFECTIVE_EFFECT_FORMULA=8*m*nt\n")
cat("LATENT_EFFECT_FORMULA=8*m*nt\n")
cat("STATE_FORMULA=4*m*nt\n")
cat("MARKER_MAP_FORMULA=40*m\n")

examples <- rbind(
  small = c(n = 1000, m = 10000, nt = 2),
  moderate = c(n = 10000, m = 100000, nt = 4)
)
for (label in rownames(examples)) {
  n <- examples[label, "n"]; m <- examples[label, "m"]; nt <- examples[label, "nt"]
  stride <- 64 * ceiling(ceiling(n / 4) / 64)
  cat(sprintf(
    "MEMORY_%s packed=%g Y=%g R=%g effective=%g latent=%g state=%g workspace=%g\n",
    toupper(label), m * stride, 8 * n * nt, 8 * n * nt,
    8 * m * nt, 8 * m * nt, 4 * m * nt, 8 * n
  ))
}

flag("THREE_PACKED_OWNERS_INVENTORIED",
     all(vapply(c("PackedBedMatrix", "FastPackedBedMatrix", "FastPackedBedMatrixBR"),
                function(x) grepl(x, paste(packed, bayesc, bayesr)), logical(1))))
flag("THREE_READERS_INVENTORIED",
     all(vapply(c("read_bedfiles_to_packed_matrix",
                  "read_bedfiles_to_fast_packed_matrix_blocked",
                  "br_read_bed_blocked"),
                function(x) grepl(x, contract, fixed = TRUE), logical(1))))
flag("COMMON_VIEW_RETAINED",
     count("struct BedPackedGenotypeView", family) == 1L &&
       grepl("reuses `BedPackedGenotypeView` unchanged", contract, fixed = TRUE))
flag("BAYESRC_VIEW_RETAINED",
     count("struct BedBayesRCPackedGenotypeView", bayesrc_types) == 1L)
flag("BED_CODE_CONTRACT",
     all(vapply(c("`00` / 0 | 2", "`01` / 1 | missing",
                  "`10` / 2 | 1", "`11` / 3 | 0"),
                function(x) grepl(x, contract, fixed = TRUE), logical(1))))
flag("PADDING_CONTRACT",
     grepl("partial-byte nor stride padding is treated as a sample", contract, fixed = TRUE))
flag("STANDARDIZED_MISSING_ZERO",
     grepl("map.val[1] = 0.0", bayesc, fixed = TRUE) &&
       grepl("missing maps to zero", contract, fixed = TRUE))
flag("RAW_MISSING_2P",
     grepl("map.val[1] = 2.0 * p", bayesc, fixed = TRUE) &&
       grepl("missing maps to `2p_j`", contract, fixed = TRUE))
flag("FULL_E_CONDITIONAL",
     all(vapply(c("C_k   = P + w_j D_k Omega D_k", "rhs_k = D_k Omega s_j",
                  "log q_k = log(pi_k)"),
                function(x) grepl(x, contract, fixed = TRUE), logical(1))))
flag("DIAGONAL_E_REDUCTION",
     grepl("C_k = P + diag(w_j d_kt/e_t)", contract, fixed = TRUE))
flag("SAMPLE_RESIDUAL_DISTINCT",
     grepl("r_jt = x_j'R_t", contract, fixed = TRUE) &&
       grepl("R = Y - X B_eff", contract, fixed = TRUE))
flag("RAW_V1_REUSED", grepl("reuses `mtblr_raw` version 1", contract, fixed = TRUE))
flag("NO_PUBLIC_MTBLR_BED", !grepl("export(mtblr_bed)", namespace, fixed = TRUE))
flag("NO_PHASE17O_SAMPLER",
     !any(file.exists(file.path(root, c("R/mtblr-bed.R", "src/blr_mt_bed_core_impl.h")))))
flag("NO_PRODUCTION_OWNER_CONSOLIDATION",
     grepl("Phase 17N does\\s+not consolidate them", contract, perl = TRUE))

required <- c(
  "THREE_PACKED_OWNERS_INVENTORIED", "THREE_READERS_INVENTORIED",
  "COMMON_VIEW_RETAINED", "BAYESRC_VIEW_RETAINED", "BED_CODE_CONTRACT",
  "PADDING_CONTRACT", "STANDARDIZED_MISSING_ZERO", "RAW_MISSING_2P",
  "FULL_E_CONDITIONAL", "DIAGONAL_E_REDUCTION", "SAMPLE_RESIDUAL_DISTINCT",
  "RAW_V1_REUSED", "NO_PUBLIC_MTBLR_BED", "NO_PHASE17O_SAMPLER",
  "NO_PRODUCTION_OWNER_CONSOLIDATION"
)
values <- vapply(required, function(name) {
  switch(name,
    THREE_PACKED_OWNERS_INVENTORIED =
      all(vapply(c("PackedBedMatrix", "FastPackedBedMatrix", "FastPackedBedMatrixBR"),
                 function(x) grepl(x, paste(packed, bayesc, bayesr)), logical(1))),
    THREE_READERS_INVENTORIED =
      all(vapply(c("read_bedfiles_to_packed_matrix",
                   "read_bedfiles_to_fast_packed_matrix_blocked", "br_read_bed_blocked"),
                 function(x) grepl(x, contract, fixed = TRUE), logical(1))),
    COMMON_VIEW_RETAINED = grepl("reuses `BedPackedGenotypeView` unchanged", contract, fixed = TRUE),
    BAYESRC_VIEW_RETAINED = grepl("BedBayesRCPackedGenotypeView", contract, fixed = TRUE),
    BED_CODE_CONTRACT = grepl("`00` / 0 | 2", contract, fixed = TRUE),
    PADDING_CONTRACT = grepl("partial-byte nor stride padding", contract, fixed = TRUE),
    STANDARDIZED_MISSING_ZERO = grepl("missing maps to zero", contract, fixed = TRUE),
    RAW_MISSING_2P = grepl("missing maps to `2p_j`", contract, fixed = TRUE),
    FULL_E_CONDITIONAL = grepl("C_k   = P + w_j D_k Omega D_k", contract, fixed = TRUE),
    DIAGONAL_E_REDUCTION = grepl("C_k = P + diag(w_j d_kt/e_t)", contract, fixed = TRUE),
    SAMPLE_RESIDUAL_DISTINCT = grepl("r_jt = x_j'R_t", contract, fixed = TRUE),
    RAW_V1_REUSED = grepl("reuses `mtblr_raw` version 1", contract, fixed = TRUE),
    NO_PUBLIC_MTBLR_BED = !grepl("export(mtblr_bed)", namespace, fixed = TRUE),
    NO_PHASE17O_SAMPLER = !file.exists(file.path(root, "src/blr_mt_bed_core_impl.h")),
    NO_PRODUCTION_OWNER_CONSOLIDATION = grepl("does\\s+not consolidate", contract, perl = TRUE)
  )
}, logical(1))
if (!all(values)) stop("Phase 17N contract audit failed: ", paste(required[!values], collapse = ", "))
