read <- function(path) paste(readLines(path,warn=FALSE),collapse="\n")
types <- read("src/blr_mt_block_eigen_types.h")
access <- read("src/blr_mt_ld_access.h")
core <- read("src/blr_mt_default_core_impl.h")
native <- read("src/mtblr.cpp")
namespace <- read("NAMESPACE")
checks <- c(
  own_trait_view=grepl("trait_operator\\[static_cast<std::size_t>\\(trait\\)\\]",access),
  transformed_wy_retained=grepl("transformed_wy=wy",native,fixed=TRUE),
  filter_specific_transform=grepl("build_block_eigen",native,fixed=TRUE),
  transformed_output=grepl("final_result, transformed_wy",native,fixed=TRUE),
  diagonal_subtraction=grepl("mt_diagonal(data, trait, marker)*difference",access,fixed=TRUE),
  one_gibbs_loop=length(gregexpr("for \\( int it =",core)[[1]])==1L,
  one_marker_call=length(gregexpr("sampleBetaCPG_Mt_latent\\(i",core)[[1]])==1L,
  build_before_core=regexpr("build_block_eigen",native)<regexpr("run_mt_block_eigen_core",native),
  descriptor_count=grepl("length one or trait count",native,fixed=TRUE),
  independent_blocks=!grepl("trait-specific blocks must match",native,fixed=TRUE),
  no_chain_copy=!grepl("chain.*BlockEigenStorage",native),
  no_research_call=!grepl("run_mt_block_eigen_core[\\s\\S]*mtblr_eigen\\(",native),
  no_public_export=!grepl("export\\(mtblr_block_eigen",namespace),
  all_wrappers=all(vapply(c("run_mt_default_core","run_mt_csr_core","run_mt_block_eigen_core"),
    grepl,logical(1),x=core,fixed=TRUE)),
  no_ww_argument=!grepl("mtblr_block_eigen_internal\\(\\s*std::vector<std::vector<double>> wy,\\s*std::vector<std::vector<double>> ww",native)
)
print(data.frame(mutation=names(checks),detected=unname(checks)),row.names=FALSE)
if(!all(checks)) stop("Phase 17L critical mutation sensitivity failed")
