test_that("test and fixture roots are separate portable concepts", {
  expect_true(dir.exists(blr_test_dir))
  expect_true(dir.exists(blr_fixture_path()))
  expect_false(identical(blr_test_dir, blr_fixture_path()))
  expect_true(file.exists(blr_fixture_path("blr-phase2-reference.R")))
})

test_that("source-root discovery validates an sblr checkout", {
  blr_skip_if_no_source_tree()
  expect_true(blr_is_source_root(blr_source_root))
  fake <- tempfile("not-an-sblr-root-"); dir.create(fake)
  on.exit(unlink(fake, recursive = TRUE), add = TRUE)
  expect_false(blr_is_source_root(fake))
  expect_true(is.na(blr_find_source_root(explicit = fake, start = fake,
                                        test_dir = fake)))
})

test_that("fixture access survives simulated source-tree absence", {
  fake <- tempfile("installed-test-dir-"); dir.create(fake)
  on.exit(unlink(fake, recursive = TRUE), add = TRUE)
  expect_true(file.exists(blr_fixture_path("blr-phase2-reference.R")))
  expect_true(is.na(blr_find_source_root(explicit = fake, start = fake,
                                        test_dir = fake)))
})

test_that("test sources obey installed-package boundaries", {
  files <- list.files(blr_test_dir, "[.]R$", full.names = TRUE)
  text <- paste(vapply(files, function(x) paste(readLines(x, warn = FALSE),
                                               collapse = "\n"), character(1)),
                collapse = "\n")
  expect_false(grepl('source_sblr_test_file\\(', text))
  expect_false(grepl('eval\\(parse\\(.*test-blr-framework', text))
  expect_false(grepl('test_path\\("\\.\\.",\\s*"\\.\\."', text))
  expect_false(grepl('source\\([^\n]*"R/', text))
})

test_that("portable fixtures and installed owners remain declared", {
  inventory <- blr_source_text("docs/dev/blr_source_only_test_inventory.md")
  contract <- blr_source_text("docs/dev/blr_installed_check_contract.md")
  expect_match(contract, "fixtures are installed-check requirements", fixed = TRUE)
  expect_match(contract, "runtime-object MD5", fixed = TRUE)
  expect_match(inventory, "No scientific or public-contract owner", fixed = TRUE)
})

test_that("CI checks a built tarball and does not ignore warnings", {
  workflow <- blr_source_text(".github/workflows/blr-framework.yml")
  driver <- blr_source_text("tools/check/check_package.R")
  expect_match(workflow, "tools/check/check_package.R", fixed = TRUE)
  expect_match(driver, "R CMD build", fixed = TRUE)
  expect_match(driver, "WARNING", fixed = TRUE)
  expect_false(grepl("ignore.*WARNING", driver, ignore.case = TRUE))
})
