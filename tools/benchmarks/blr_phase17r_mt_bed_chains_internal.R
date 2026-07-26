pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source("tests/testthat/helper-mtblr-bed-contract.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-internal.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-chains-internal.R", local = TRUE)
source("tests/testthat/helper-mtblr-bed-multichain-contract.R", local = TRUE)
phase17r_benchmark_case <- function(size, mode) {
  if (size == "small") return(phase17o_case(
    nt=2L,residual_covariance=mode,updates=TRUE,multiple_sets=TRUE))
  n <- 28L; m <- 24L
  dosage <- outer(seq_len(n), seq_len(m), function(i, j) (i + 2*j) %% 3)
  dosage[(row(dosage) + 3*col(dosage)) %% 17 == 0] <- NA_real_
  path <- tempfile("phase17r-moderate-", fileext = ".bed")
  phase17n_write_bed(path, dosage)
  af <- seq(.18, .42, length.out = m)
  X <- phase17n_transform_genotypes(dosage, af, TRUE)
  case <- phase17o_case(nt=3L,residual_covariance=mode,updates=TRUE,
                        multiple_sets=TRUE)
  unlink(case$fixture$paths)
  case$fixture <- list(paths=path,n_bed=n,rows=NULL,cls=list(seq_len(m)),af=af,X=X)
  case$Y <- vapply(seq_len(3L), function(t) {
    y <- sin(seq_len(n)*(.11+t/19))+cos(seq_len(n)*(.07+t/23)); y-mean(y)
  }, numeric(n)); colnames(case$Y) <- paste0("T",1:3)
  case$beta <- replicate(3L, numeric(m), simplify=FALSE)
  case$b <- replicate(3L, numeric(m), simplify=FALSE)
  case$state <- replicate(3L, integer(m), simplify=FALSE)
  case$sets <- list(as.integer(seq(0L,m-1L,2L)),as.integer(seq(1L,m-1L,2L)))
  case$nit <- 8L; case$nburn <- 2L
  case
}
rows <- list(); index <- 0L
for (size in c("small", "moderate")) for (mode in c("diagonal", "full")) for (nchains in c(1L,2L,4L))
 for (ncores in c(1L,2L,4L)) for (keep in c(FALSE,TRUE)) {
  case <- phase17r_benchmark_case(size,mode)
  elapsed <- system.time(raw <- phase17r_call(case,nchains,ncores,keep_chains=keep))["elapsed"]
  bed <- raw$diagnostics$mt_bed; index <- index+1L
  mem <- phase17q_memory(nrow(case$Y),ncol(case$fixture$X),ncol(case$Y),
    length(case$models),case$nit+case$nburn,nchains,ncores,keep)
  rows[[index]] <- data.frame(size=size,n=nrow(case$Y),m=ncol(case$fixture$X),
   nt=ncol(case$Y),mode=mode,nchains=nchains,requested_cores=ncores,
   used_workers=bed$used_workers,keep_chains=keep,
   dispatch_seconds=bed$dispatch_seconds,chain_seconds_mean=bed$seconds_mean,
   chain_seconds_max=bed$seconds_max,
   preparation_and_aggregation_seconds=max(0,unname(elapsed)-bed$dispatch_seconds),
   total_seconds=unname(elapsed),shared_bytes=mem$shared_bytes,
   private_bytes_per_worker=mem$private_state_bytes_per_chain,
   result_bytes_per_chain=mem$result_bytes_per_chain,
   retained_output_bytes=mem$estimated_retained_output_bytes,
   raw_bytes=as.numeric(object.size(raw)))
  phase17o_cleanup(case)
 }
print(do.call(rbind,rows),row.names=FALSE)
cat("Phase 17R regression signals only; no linear-speedup claim.\n")
