# BLR installed-check contract

## Purpose

The package is validated both from a checkout and from the namespace installed from a source tarball.

## Test categories

Scientific/public tests and portable fixture tests run in both contexts. Repository architecture tests run in a checkout and skip only their individual block when repository assets are absent.

## Path model

`blr_test_dir` identifies `tests/testthat`. `blr_fixture_path()` resolves only beneath `tests/testthat/fixtures`. `blr_find_source_root()` returns a root only after validating `DESCRIPTION`, `R/`, `src/`, `tests/testthat/`, and `Package: sblr`. `blr_repo_path()` issues the narrow source-architecture skip if no root exists.

## Fixture policy

Portable fixtures are installed-check requirements. Frozen file hashes remain exact and never depend on a source root.

## Source architecture policy

Source text, workflow, developer-document, audit-tool, benchmark, and active-source hash assertions are source-only. The complete inventory is in `blr_source_only_test_inventory.md`.

## Namespace policy

Tests call exported functions normally and obtain internals with `getFromNamespace()` or `sblr:::`. They never source production `R/*.R` files.

## Numerical reference policy

Structure, names, dimensions, discrete states, seed mapping, retained counts, and same-build determinism remain exact. Portable continuous fixture comparisons use their model-owned fieldwise tolerance, normally `1e-12`. The learned-annotation chain-summary standard deviation alone uses `1e-8` after a measured GCC-versus-Rtools difference of `1.756119e-9` around an exact-zero reference. The runtime-object MD5 is toolchain-specific developer evidence, not a portable installed reference. No tolerance exceeds `1e-8`.

## Documentation policy

Roxygen comments are authoritative; generated Rd and NAMESPACE are regenerated and checked.

## Local and CI check

Run `Rscript tools/check/check_package.R .`. It builds a source tarball in a temporary directory, checks that tarball, prints the log, and fails on ERROR or WARNING. CI uses the same command.

## Expected skips

Installed checks skip the source-only blocks listed in the inventory plus explicitly disabled extended reproducibility and peak-RSS blocks. Scientific/public and fixture owners never skip for lack of a checkout.

## Completion criterion

The source suite has no failures or warnings, and the built-tarball check has no ERROR or WARNING.
