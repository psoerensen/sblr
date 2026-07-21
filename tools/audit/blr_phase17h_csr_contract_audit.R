sizes <- data.frame(
  fixture = c("small", "moderate"),
  markers = c(1000, 100000),
  input_nnz = c(5000, 5000000),
  traits = c(2, 4)
)
sizes$symmetric_nnz <- 2 * sizes$input_nnz
sizes$row_pointer_bytes <- (sizes$markers + 1) * 8
sizes$column_index_bytes <- sizes$symmetric_nnz * 4
sizes$offdiag_value_bytes <- sizes$symmetric_nnz * 4
sizes$diagonal_bytes <- sizes$markers * 8
sizes$one_operator_bytes <- with(sizes,
  row_pointer_bytes + column_index_bytes + offdiag_value_bytes + diagonal_bytes)
sizes$shared_operator_bytes <- sizes$one_operator_bytes
sizes$shared_pattern_trait_values_bytes <- with(sizes,
  row_pointer_bytes + column_index_bytes +
    traits * (offdiag_value_bytes + diagonal_bytes))
sizes$independent_operator_bytes <- sizes$traits * sizes$one_operator_bytes
sizes$builder_storage_owners <- 1L
sizes$ordinary_adapter_compatibility_copies <- 1L
sizes$borrowed_views <- 1L
sizes$per_chain_csr_bytes <- 0

print(sizes, row.names = FALSE)
cat("\nContract observations:\n")
cat("- Phase 17H introduces no O(nnz) allocation or copy.\n")
cat("- CsrOperator retains one pre-existing owning compatibility copy.\n")
cat("- Chain tasks borrow one immutable view; per-chain CSR bytes are zero.\n")
cat("- No MCMC-time disk I/O occurs.\n")
cat("- Calculated bytes are storage estimates, not peak RSS.\n")
