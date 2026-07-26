phase17q_uint32 <- function(x) {
  modulus <- 2^32
  ((as.numeric(x) %% modulus) + modulus) %% modulus
}

phase17q_tasks <- function(nchains, nt = 1L) {
  stopifnot(length(nchains) == 1L, nchains >= 1L, nt >= 1L)
  data.frame(chain = seq.int(0L, as.integer(nchains) - 1L),
             result_slot = seq.int(0L, as.integer(nchains) - 1L))
}

phase17q_seeds <- function(seed, nchains, chain_seeds = NULL) {
  if (!is.null(chain_seeds)) {
    stopifnot(length(chain_seeds) == nchains, all(is.finite(chain_seeds)))
    return(vapply(chain_seeds, phase17q_uint32, numeric(1)))
  }
  vapply(seq.int(0, nchains - 1), function(chain) {
    phase17q_uint32(as.numeric(seed) + 9176 * chain)
  }, numeric(1))
}

phase17q_chain <- function(index, retained, bm, dm, final_state,
                           trace_offset = 0, failed = FALSE, error = "") {
  bm <- as.matrix(bm); dm <- as.matrix(dm)
  nt <- ncol(bm)
  list(
    chain = as.integer(index), failed = failed, error = error,
    retained = as.integer(retained),
    bm_acc = bm * retained, dm_acc = dm * retained,
    covb_acc = diag(nt) * retained * index,
    covg_acc = diag(nt) * retained * (index + 1),
    cove_acc = diag(nt) * retained * (index + 2),
    pi_acc = c(.25, .75) * retained,
    vbs = matrix(seq_len(3 * nt) + trace_offset, nt),
    vgs = matrix(seq_len(3 * nt) + trace_offset + 10, nt),
    ves = matrix(seq_len(3 * nt) + trace_offset + 20, nt),
    b = bm + index, state = as.matrix(final_state), r = bm - index,
    B = diag(nt) * index, G = diag(nt) * (index + 1),
    E = diag(nt) * (index + 2), pi_final = c(.1, .9),
    seed = phase17q_uint32(index), seconds = index / 10
  )
}

phase17q_stability <- function(values) {
  arr <- simplify2array(values)
  if (length(values) == 1L) {
    return(list(sd = array(0, dim(values[[1]])),
                min = values[[1]], max = values[[1]]))
  }
  margin <- seq_len(length(dim(arr)) - 1L)
  list(sd = apply(arr, margin, stats::sd),
       min = apply(arr, margin, min),
       max = apply(arr, margin, max))
}

phase17q_aggregate <- function(chains, keep_chains = FALSE) {
  failures <- chains[vapply(chains, `[[`, logical(1), "failed")]
  if (length(failures)) {
    failures <- failures[order(vapply(failures, `[[`, integer(1), "chain"))]
    stop(paste(vapply(failures, function(x) {
      sprintf("chain %d: %s", x$chain, x$error)
    }, character(1)), collapse = "; "), call. = FALSE)
  }
  chains <- chains[order(vapply(chains, `[[`, integer(1), "chain"))]
  retained <- sum(vapply(chains, `[[`, integer(1), "retained"))
  sum_field <- function(field) Reduce(`+`, lapply(chains, `[[`, field))
  mean_trace <- function(field) sum_field(field) / length(chains)
  bm_chain <- lapply(chains, function(x) x$bm_acc / x$retained)
  dm_chain <- lapply(chains, function(x) x$dm_acc / x$retained)
  compact <- function(x) x[c("chain", "seed", "bm_acc", "dm_acc", "b",
                              "state", "vbs", "vgs", "ves", "B", "G", "E",
                              "pi_final", "seconds")]
  list(
    bm = sum_field("bm_acc") / retained,
    dm = sum_field("dm_acc") / retained,
    covb = sum_field("covb_acc") / retained,
    covg = sum_field("covg_acc") / retained,
    cove = sum_field("cove_acc") / retained,
    pi_mean = sum_field("pi_acc") / retained,
    vbs = mean_trace("vbs"), vgs = mean_trace("vgs"),
    ves = mean_trace("ves"),
    b = chains[[1]]$b, state = chains[[1]]$state,
    r = chains[[1]]$r, B = chains[[1]]$B, G = chains[[1]]$G,
    E = chains[[1]]$E, pi_final = chains[[1]]$pi_final,
    bm_stability = phase17q_stability(bm_chain),
    dm_stability = phase17q_stability(dm_chain),
    chains = if (keep_chains) setNames(lapply(chains, compact),
                                      paste0("chain", seq_along(chains))) else NULL,
    retained = retained
  )
}

phase17q_memory <- function(n, m, nt, nmodels, trace_length,
                            nchains, ncores, keep_chains = FALSE) {
  packed <- m * ceiling(n / 4)
  shared <- packed + 8 * (n * nt + 5 * m + m * nt) + 4 * m
  private <- 8 * (n * nt + 2 * m * nt + nt^2 * 4 + n) + 4 * m * nt
  result <- 8 * (4 * m * nt + 3 * nt * trace_length + 4 * nt^2 + nmodels)
  retained <- 8 * (3 * m * nt + 3 * nt * trace_length + 3 * nt^2 + nmodels)
  workers <- min(as.integer(nchains), as.integer(ncores))
  retained_total <- if (keep_chains) nchains * retained else 0
  total <- shared + workers * private + nchains * result + retained_total
  list(shared_bytes = shared, private_state_bytes_per_chain = private,
       result_bytes_per_chain = result,
       retained_chain_bytes_per_chain = retained,
       worker_count = workers, nchains = nchains,
       estimated_concurrent_bytes = shared + workers * private,
       estimated_retained_output_bytes = retained_total,
       estimated_total_bytes = total, estimated_total_gib = total / 2^30,
       estimate_kind = "analytical upper-bound estimate; not measured RSS or peak RSS")
}

phase17q_future_controls <- function() {
  list(nchains = 1L, ncores = 1L, chain_seeds = NULL, keep_chains = FALSE,
       task_topology = "one complete joint MT chain",
       openmp_unavailable = "warn once and run serial",
       raw_schema = "mtblr_raw version 1",
       final_state_policy = "primary_chain",
       primary_chain = 1L,
       posterior_summary_policy = "pooled_retained_samples",
       trace_policy = "iterationwise_chain_mean")
}

phase17q_future_internal_signature <- function() {
  "mtblr_bed_chains_internal(current arguments, int nchains, int ncores, std::vector<int> chain_seeds, bool keep_chains)"
}
