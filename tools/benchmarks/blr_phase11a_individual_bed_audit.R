# Phase 11A packed-BED audit baseline. Completed-fit RSS is not peak memory.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11a-bed-reference.R"))

rss_mb <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2
}
bench <- function(model, ncores, nchains, reps = 5L, markers = 2L, samples = 6L) {
  warm <- phase11a_capture(model, ncores, nchains, 900L)
  before <- rss_mb()
  times <- vapply(seq_len(reps), function(i) system.time(
    phase11a_capture(model, ncores, nchains, 900L + i))[["elapsed"]], numeric(1))
  data.frame(model = model, markers = markers, samples = samples, chains = nchains,
    cores = ncores, times = paste(times, collapse = ","), mean = mean(times),
    median = median(times), minimum = min(times), maximum = max(times),
    range = diff(range(times)), rss_before_mb = before, rss_after_mb = rss_mb(),
    memory_method = "ps whole-process RSS after completed fits; not peak RSS")
}

results <- do.call(rbind, list(
  bench("bayesc", 1L, 1L), bench("bayesc", 1L, 2L), bench("bayesc", 2L, 2L),
  bench("bayesr", 1L, 1L), bench("bayesr", 2L, 2L),
  bench("bayesrc", 1L, 1L), bench("bayesrc", 2L, 2L)))

# Moderate on-disk workload. This exercises decoding and page-cache effects;
# it is deliberately separate from the tiny exact-reference fixtures.
original_fixture <- phase11a_fixture
set.seed(1101)
moderate_dosage <- matrix(sample(0:2, 2000L * 200L, replace = TRUE), 2000L, 200L)
moderate_bed <- tempfile(fileext = ".bed")
phase11a_write_bed(moderate_bed, moderate_dosage)
moderate <- list(Glist = list(n = 200L, ids = paste0("id", 1:200),
  bedfiles = moderate_bed, rsids = list(paste0("rs", 1:2000)),
  rsidsLD = list(paste0("rs", 1:2000)), chr = list(rep(1L, 2000)),
  pos = list(seq_len(2000)), af = list(rowMeans(moderate_dosage) / 2)),
  y = matrix(scale(seq_len(200)), ncol = 1, dimnames = list(NULL, "D1")))
phase11a_fixture <- function() moderate
moderate_results <- rbind(
  bench("bayesc", 1L, 1L, markers = 2000L, samples = 200L),
  bench("bayesr", 1L, 1L, markers = 2000L, samples = 200L))
phase11a_fixture <- original_fixture
results <- rbind(results, moderate_results)
cat("Phase 11A individual/packed-BED backend audit baseline\n")
cat("All production paths decode SNP-major packed BED into fit-local native storage.\n")
cat("Tiny fixture file size:", file.info(phase11a_fixture()$Glist$bedfiles)$size, "bytes\n")
cat("Moderate BED file size:", file.info(moderate_bed)$size, "bytes\n")
cat("Cold/warm distinction: first call is warm-up; later calls may use OS page cache.\n")
cat("R:", R.version.string, "\n")
cat("sblr:", as.character(utils::packageVersion("sblr")), "\n")
print(results, row.names = FALSE)
