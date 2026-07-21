root <- normalizePath(".", winslash = "/", mustWork = TRUE)
read <- function(path) paste(readLines(file.path(root, path), warn = FALSE), collapse="\n")
tests <- read("tests/testthat/test-blr-framework-phase17j2.R")
helper <- read("tests/testthat/helper-blr-test-contracts.R")
workflow <- read(".github/workflows/blr-framework.yml")
checks <- c(
  fixture_independent = grepl("blr_fixture_path", helper, fixed=TRUE),
  root_structural = all(vapply(c("DESCRIPTION", '"R"', '"src"'), grepl,
                               logical(1), x=helper, fixed=TRUE)),
  no_whole_file_skip = !grepl("skip_on_cran|NOT_CRAN", tests),
  no_production_source = grepl("source_sblr_test_file", tests, fixed=TRUE),
  no_cross_test_parse = grepl("eval", tests, fixed=TRUE),
  fixture_md5_installed = grepl("blr_expect_fixture_md5", helper, fixed=TRUE),
  runtime_hash_policy = grepl("runtime-object MD5", read("docs/dev/blr_installed_check_contract.md"), fixed=TRUE),
  warnings_enforced = grepl("tools/check/check_package.R", workflow, fixed=TRUE),
  source_architecture_active = grepl("blr_repo_path", helper, fixed=TRUE),
  scientific_not_source_only = grepl("No scientific or public-contract owner", read("docs/dev/blr_source_only_test_inventory.md"), fixed=TRUE)
)
print(data.frame(mutation=names(checks), detected=unname(checks)), row.names=FALSE)
if (!all(checks)) stop("Phase 17J2 mutation-sensitivity audit failed: ",
                       paste(names(checks)[!checks], collapse=", "))
