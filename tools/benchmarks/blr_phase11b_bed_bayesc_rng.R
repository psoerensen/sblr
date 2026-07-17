# Phase 11B corrected packed-BED BayesC baseline. Completed-fit RSS is not peak.
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
source(file.path("tests", "testthat", "fixtures", "blr-phase11b-bed-bayesc-reference.R"))
rss_mb <- function() if (requireNamespace("ps", quietly = TRUE))
  as.numeric(ps::ps_memory_info()[["rss"]]) / 1024^2 else NA_real_
bench <- function(route, cores, chains, markers = 2L, samples = 6L, reps = 5L) {
  invisible(phase11b_capture(route, cores, chains, 901L)); before <- rss_mb()
  times <- vapply(seq_len(reps), function(i) system.time(
    phase11b_capture(route, cores, chains, 901L + i))[["elapsed"]], numeric(1))
  data.frame(route, markers, samples, chains, cores,
    scheduler = "full=10,skip=50:200,candidate=.001/life=20",
    times = paste(times, collapse = ","), mean = mean(times), median = median(times),
    minimum = min(times), maximum = max(times), range = diff(range(times)),
    rss_before_mb = before, rss_after_mb = rss_mb(),
    memory_method = "ps whole-process RSS after completed fits; not peak RSS")
}
results <- rbind(bench("single", 1L, 1L), bench("multichain", 1L, 2L),
  bench("multichain", 2L, 2L))
original_fixture <- phase11a_fixture; set.seed(1112)
dosage <- matrix(sample(0:2, 2000L * 200L, replace = TRUE), 2000L, 200L)
bed <- tempfile(fileext = ".bed"); phase11a_write_bed(bed, dosage)
moderate <- list(Glist = list(n = 200L, ids = paste0("id", 1:200), bedfiles = bed,
  rsids = list(paste0("rs", 1:2000)), rsidsLD = list(paste0("rs", 1:2000)),
  chr = list(rep(1L, 2000)), pos = list(seq_len(2000)),
  af = list(rowMeans(dosage) / 2)), y = matrix(scale(seq_len(200)), ncol = 1,
  dimnames = list(NULL, "D1")))
phase11a_fixture <- function() moderate
results <- rbind(results, bench("single", 1L, 1L, 2000L, 200L),
  bench("multichain", 1L, 2L, 2000L, 200L),
  bench("multichain", 2L, 2L, 2000L, 200L)); phase11a_fixture <- original_fixture
cat("Phase 11B corrected packed-BED BayesC RNG baseline\n")
cat("BED bytes:", file.info(bed)$size, "; repeated fits may use OS page cache.\n")
cat("R:", R.version.string, "; sblr:", as.character(packageVersion("sblr")), "\n")
print(results, row.names = FALSE)
