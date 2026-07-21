cases <- data.frame(mode=c("shared", "trait_specific"), markers=10000,
  traits=4, input_nnz_per_trait=500000)
cases$symmetric_nnz_per_trait <- 2 * cases$input_nnz_per_trait
operators <- c(1, cases$traits[2])
cases$storage_owners <- operators
cases$diagonal_owners <- operators
cases$borrowed_views <- cases$traits
cases$row_pointer_bytes <- operators * (cases$markers + 1) * 8
cases$index_bytes <- operators * cases$symmetric_nnz_per_trait * 4
cases$offdiag_value_bytes <- operators * cases$symmetric_nnz_per_trait * 4
cases$diagonal_bytes <- operators * cases$markers * 8
cases$csr_bytes <- with(cases, row_pointer_bytes + index_bytes +
  offdiag_value_bytes + diagonal_bytes)
cases$dense_oracle_bytes <- cases$traits * cases$markers^2 * 8
cases$per_chain_csr_bytes <- 0
print(cases, row.names=FALSE)
cat("All disk I/O completes before the Gibbs core; completed-fit RSS is not peak RSS.\n")
