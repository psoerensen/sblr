root <- normalizePath(".", winslash="/")
read <- function(...) paste(readLines(file.path(root,...),warn=FALSE),collapse="\n")
contract <- read("src","blr_block_eigen.h"); builder <- read("src","st_block_eigen.cpp")
bridge <- read("src","st_block_eigen_rcpp.h")
scalar <- paste(read("src","st_cpg_omp_csr.cpp"),read("src","st_cpg_omp_csr_bayesr.cpp"),read("src","st_sbayesrc_omp_csr.cpp"))
mt <- read("src","mtblr.cpp")
checks <- c(
 float_packing=grepl("std::vector<float> upper_triangle",contract,fixed=TRUE),
 hard_wy_projection=grepl("wy_mat.cols(c0, c1) = (W * Lk) * Rk",builder,fixed=TRUE),
 ridge_does_not_project=grepl("bool do_project = false",builder,fixed=TRUE),
 threshold_floor=grepl("std::max(tau, mu_floor)",builder,fixed=TRUE),
 fallback=grepl("mu.index_max()",builder,fixed=TRUE),
 packed_diagonal=grepl("stored.sym_at(i, i)",builder,fixed=TRUE),
 no_cross_block_update=grepl("block.start + j",contract,fixed=TRUE),
 mappings_validated=grepl("mappings do not identify",contract,fixed=TRUE),
 vector_order=length(gregexpr("for (int j = 0; j < block.size; ++j)",contract,fixed=TRUE)[[1]])>=4,
 scalar_uses_storage=length(gregexpr("BlockEigenOperator",scalar,fixed=TRUE)[[1]])>=6,
 public_csr_unchanged=!grepl("(?s)stblr_csr\\s*<-\\s*function\\([^)]*ld_backend",read("R","sparse_ld_bed_helper.R"),perl=TRUE),
 mt_research_not_canonical=grepl("mtblr_eigen",mt,fixed=TRUE)&&!grepl("BlockEigenView",mt,fixed=TRUE),
 one_storage_definition=length(gregexpr("struct BlockEigenStorage",paste(contract,read("src","st_ld_operator.h")),fixed=TRUE)[[1]])==1,
 inspection_reuses_builder=grepl("build_block_eigen(",builder,fixed=TRUE)&&!grepl("eig_sym|tildeA",bridge)
)
print(data.frame(mutation=names(checks),detected=unname(checks)),row.names=FALSE)
if(!all(checks)) stop("Undetected: ",paste(names(checks)[!checks],collapse=", "))
