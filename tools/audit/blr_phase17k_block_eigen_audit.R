args <- commandArgs(trailingOnly = TRUE)
sizes <- if (length(args)) as.integer(strsplit(args[[1]], ",", fixed = TRUE)[[1]]) else c(1L, 4L, 8L, 16L)
stopifnot(all(is.finite(sizes)), all(sizes > 0L))
m <- sum(sizes); packed <- sum(sizes * (sizes + 1) / 2)
print(data.frame(marker_count=m, blocks=length(sizes), block_sizes=paste(sizes,collapse=","),
 packed_value_bytes=packed*4, mapping_bytes=2*m*4, diagonal_bytes=m*8,
 dense_full_bytes=m*m*8, max_Z_float_bytes_per_sample=max(sizes)*4,
 max_quadratic_double_bytes_per_matrix=max(sizes)^2*8,
 full_sweep_value_visits=sum(sizes^2)), row.names=FALSE)
cat("runtime_representation=float-packed filtered dense blocks\n")
cat("hard_truncation_reduces_runtime_rank=FALSE\n")
cat("per_chain_operator_bytes=0\nMCMC_time_BED_IO=0\nMCMC_time_eigendecomposition=0\n")
