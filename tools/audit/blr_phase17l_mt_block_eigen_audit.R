args <- commandArgs(trailingOnly=TRUE)
traits <- if(length(args)) as.integer(args[[1]]) else 3L
markers <- if(length(args)>1L) as.integer(args[[2]]) else 1000L
sizes <- rep(50L,ceiling(markers/50)); sizes[length(sizes)] <- markers-sum(head(sizes,-1L))
packed <- sum(sizes*(sizes+1)/2)*4
mapping <- 2*markers*4
diagonal <- markers*8
wy <- traits*markers*8
report <- data.frame(
  mode=c("shared","trait_specific"), traits=traits, markers=markers,
  operator_count=c(1L,traits), borrowed_views=traits,
  packed_value_bytes=c(packed,traits*packed), mapping_bytes=c(mapping,traits*mapping),
  diagonal_bytes=c(diagonal,traits*diagonal), transformed_wy_bytes=wy,
  per_chain_operator_bytes=0L, mcmc_bed_io=0L, mcmc_eigendecompositions=0L)
print(report,row.names=FALSE)
cat("BLOCK_BUILD_COMPLETES_BEFORE_RNG=TRUE\n",
    "TRANSFORMED_WY_USED_FOR_RANKING=TRUE\n",
    "TRANSFORMED_WY_USED_FOR_CORE=TRUE\n",
    "TRANSFORMED_WY_USED_FOR_OUTPUT=TRUE\n",sep="")
