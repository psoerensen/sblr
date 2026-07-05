# `make_credible_sets()` Internals

This note documents how `make_credible_sets()` currently works, based on
inspection of `R/credible_sets.R`. No R or C++ behaviour was changed to
produce this note.

`make_credible_sets()` is R-level post-processing over an already-fitted
ST-BLR object. It does not run MCMC and does not call any native code
directly; native code is only reached indirectly through `sparseLD_read_CSR()`
when sparse LD is used as the LD source.

## Function Family

`R/credible_sets.R` defines three related exported functions:

- `make_credible_sets_from_ld(pip, LD, ...)` — single-region credible-set
  construction from a raw named PIP vector and a dense LD matrix. No fit
  object or `Glist` involved.
- `make_multisignal_credible_sets_from_ld(pip, LD, ...)` — same inputs, but
  with explicit signal-level stopping rules (`signal_coverage`,
  `min_signal_pip`) intended for regions where the total regional PIP is much
  larger than the lead marker's PIP.
- `make_credible_sets(fit, Glist, ...)` — the fit-level entry point. It
  resolves PIPs and loci from a fitted object, resolves regional LD (dense or
  sparse), and calls `make_credible_sets_from_ld()` once per locus.

`make_credible_sets()` calls the **single-signal** `make_credible_sets_from_ld()`
per locus, not `make_multisignal_credible_sets_from_ld()`. See "Sets per
locus" below for what this implies in practice.

## Formatted Fit Fields It Depends On

- `fit$dm` (required) — marker PIPs. Read by the internal helper
  `.stblr_extract_pip()`. If `fit$dm` is a matrix/data frame, `trait` selects
  a column (by name or index) and marker names come from `rownames(fit$dm)`.
  If `fit$dm` is a plain vector, marker names come from `names(fit$dm)`, and
  `trait` is only used to label the output.
- No other stable fit fields (`bm`, `vbs`, `chains`, `ld_swap`, etc.) are read
  by `make_credible_sets()` itself. Those are consumed downstream by
  `extract_stblr_finemap_loci()` (see
  `docs/dev/stblr_csr_finemap_overview.md`), not here.
- `Glist` is optional and is not a formatted-fit field, but it supplies marker
  positions (`Glist$rsids`, `Glist$pos`, optionally `Glist$chr`) for automatic
  locus definition, and `Glist$sparseLD$prefix` (plus optional
  `Glist$sparseLD$r2_threshold`) when dense `LD` is not supplied directly.

## How It Works

1. **Extract PIPs.** `.stblr_extract_pip(fit, trait)` pulls a named PIP vector
   for one trait out of `fit$dm`.
2. **Resolve loci** (mutually exclusive paths):
   - **No `sets` supplied:** `.stblr_marker_map_from_Glist(Glist, fit)` builds
     a marker/chromosome/position table from `Glist$rsids` and `Glist$pos`
     (using `Glist$sparseLD$chr`/`Glist$sparseLD$cls` component indexing when
     present, otherwise treating every `Glist` component as one chromosome).
     `.stblr_define_loci()` then greedily picks the highest-remaining-PIP
     marker above `locus_pip_cutoff` as a lead, takes all markers on the same
     chromosome within `max_locus_distance` base pairs, removes them from
     consideration, and repeats until no marker clears `locus_pip_cutoff` or
     `max_loci` is reached. This requires `Glist` — it errors without one.
   - **`sets` supplied:** `.stblr_clean_marker_sets()` deduplicates and drops
     markers not present in `fit$dm` (dropping empty sets with a warning, and
     erroring if every set ends up empty). `.stblr_loci_from_sets()` then
     builds one locus row per supplied set, with the highest-PIP member as
     `lead_marker`; `chr`/`start`/`end` are only filled in when a marker map is
     available (via `Glist`) and only when all markers in the set share one
     chromosome, otherwise they are `NA`.
3. **Resolve regional LD**, per locus:
   - If `LD` is supplied, it is used directly (a single dense matrix requires
     row/column names when there is more than one locus; a named list of
     matrices is indexed by locus name).
   - Otherwise (`Glist$sparseLD$prefix` required), `sparseLD_read_CSR()` reads
     the sparse LD once, and `.extract_sparseLD_region_dense()` densifies just
     the locus's rows/columns into a dense regional matrix. A warning is
     raised up front if `Glist$sparseLD$r2_threshold` is higher than the
     requested `min_r2`, since needed LD pairs may then be missing from the
     sparse structure. `.extract_sparseLD_region_dense()` separately warns if
     a locus has more than 10 markers but every densified off-diagonal entry
     is exactly zero, as a signal that indexing is likely wrong rather than
     that the region is genuinely LD-independent.
4. **Build credible set(s) per locus** by calling `make_credible_sets_from_ld()`
   with the locus's regional PIPs and regional LD (see "Ranking and Grouping"
   below). `max_sets` is not exposed as a `make_credible_sets()` argument, so
   it is left at its `make_credible_sets_from_ld()` default of `Inf`.
5. **Assemble output** by row-binding per-locus summaries (tagging each with
   `locus`, `trait`, `chr`, `start`, `end`) and collecting per-locus set lists,
   loci metadata, and the resolved parameters.

## Ranking and Grouping (`make_credible_sets_from_ld()`)

Within one region:

1. PIPs below `pip_cutoff` are zeroed; `r2 = LD^2` with the diagonal forced
   to `1`.
2. While any marker still has nonzero PIP (`remaining`) and the per-locus set
   count is below `max_sets`:
   - **Lead selection** (`method`): `"pip"` picks the remaining marker with
     the highest PIP; `"ld_pip"` picks the remaining marker with the highest
     LD-smoothed score `r2 %*% pip` (restricted to `remaining` positions).
   - **Candidate membership**: all remaining markers with `r2` to the lead
     `>= min_r2`; if none qualify, the candidate set falls back to the lead
     alone.
   - **Ranking within the candidate set**: sorted by decreasing PIP.
   - **Coverage cutoff**: markers are taken in that order until the
     cumulative PIP first reaches `coverage`. If it never reaches `coverage`
     using all candidates:
     - `allow_incomplete = FALSE` (default): the whole candidate set is
       dropped from further lead selection (`remaining[candidate_idx] <- FALSE`)
       and the loop tries a new lead — this locus's summary/sets simply omit
       that attempted set rather than recording an incomplete one.
     - `allow_incomplete = TRUE`: all candidates are taken as the set even
       though `cs_pip < coverage`.
   - **Marker removal for the next iteration** (`remove`): `"ld_neighborhood"`
     (default) removes every candidate (not just the selected members) from
     future lead consideration; `"credible_set"` removes only the markers
     actually selected into the set, leaving other same-neighborhood markers
     eligible to become a later lead.
3. Each accepted set is recorded with lead marker, lead PIP, total `cs_pip`,
   member count, and min/mean `r2` to the lead.

**Sets per locus in practice:** because `max_sets` defaults to `Inf` and the
loop only stops when `remaining` is empty, a single locus from
`make_credible_sets()` can yield more than one credible set (`CS1`, `CS2`, ...)
if it contains multiple LD-separated high-PIP peaks that each independently
reach `coverage`. This is a real (if implicit) form of multi-signal detection
within `make_credible_sets()`, distinct from — and simpler than —
`make_multisignal_credible_sets_from_ld()`, which additionally tracks
`signal_coverage`/`min_signal_pip` stopping rules and a `complete` flag per
signal. `make_credible_sets()` never surfaces that `complete` flag; an
under-coverage attempt is just silently dropped (per point 2 above) rather
than returned and marked incomplete, unless `allow_incomplete = TRUE`.

## What It Returns

`make_credible_sets()` returns a list with:

- `summary` — one row per accepted credible set across all loci:
  `locus`, `trait`, `chr`, `start`, `end`, `cs`, `lead_marker`, `lead_pip`,
  `cs_pip`, `n_markers`, `min_r2_to_lead`, `mean_r2_to_lead`, plus the
  parameters used (`coverage`, `min_r2`, `pip_cutoff`, `method`). Empty
  (`data.frame()`) if no locus produced any accepted set.
- `sets` — named list (by locus) of named lists (by `CS1`, `CS2`, ...) of
  per-marker data frames (`marker`, `pip`, `r2_to_lead`, `rank`).
- `locus_sets` — named list (by locus) of the full marker vectors used as
  input to each locus's credible-set construction (i.e. all markers in the
  locus, not just those selected into a credible set). Documented as the
  intended input to `finemap_stblr_csr(..., sets = ...)` and used as
  `extract_stblr_finemap_loci(..., locus_sets = ...)` in the recommended
  workflow (`docs/dev/stblr_csr_finemap_overview.md`).
- `loci` — the locus table (`locus`, `chr`, `start`, `end`, `lead_marker`,
  `lead_pip`, `n_markers`), with the internal `markers`/`indices` list-columns
  stripped out.
- `parameters` — the resolved call parameters (`coverage`, `min_r2`,
  `pip_cutoff`, `locus_pip_cutoff`, `max_locus_distance`, `max_loci`, `method`,
  `allow_incomplete`, `remove`).
- `trait` — the resolved trait label from `.stblr_extract_pip()`.

## Limitations

- **Marginal PIPs, not posterior configurations.** Both the single- and
  multi-signal builders operate on marginal per-marker PIPs and pairwise LD.
  This is not equivalent to SuSiE-style per-effect credible sets, which
  condition on a fitted number of causal effects; the roxygen docs for
  `make_multisignal_credible_sets_from_ld()` state this explicitly, and the
  same caveat applies to the (simpler) per-locus loop used inside
  `make_credible_sets()`.
- **Sparse LD thresholding.** When sparse `Glist$sparseLD` is used, any pair
  excluded at build time (`r2_threshold`, window) is unavailable here. A
  warning fires when `min_r2` is looser than `r2_threshold`, but there is no
  automatic correction — the resulting credible set can silently omit
  markers that are actually in LD with the lead.
- **Locus definition is positional and greedy, not model-based.** Automatic
  loci are chromosome + base-pair windows around successive highest-PIP
  markers; they do not account for LD structure when defining locus
  boundaries (only when building credible sets within a locus), so a locus
  boundary can split or merge signals in ways that depend only on
  `max_locus_distance`.
- **Single trait per call.** `trait` selects one column of `fit$dm`; multi-trait
  runs require looping over traits at the call site.
- **`max_sets` is not user-adjustable from `make_credible_sets()`.** It is
  always `Inf` inside the per-locus `make_credible_sets_from_ld()` call, so
  locus-level runaway multi-set behavior (see above) cannot be capped except
  indirectly through `locus_pip_cutoff`/`max_locus_distance` at the locus-
  definition stage.
- **Incomplete sets are dropped, not flagged, by default.** With the default
  `allow_incomplete = FALSE`, an attempted set that cannot reach `coverage`
  is silently excluded from `summary`/`sets` for that locus rather than
  appearing with a `complete = FALSE` marker (contrast with
  `make_multisignal_credible_sets_from_ld()`, which does carry a `complete`
  column).

## Documentation Gaps / Unclear Behaviour Discovered

- `.stblr_loci_from_sets()` silently sets `chr`/`start`/`end` to `NA` when a
  user-supplied set spans more than one chromosome, with no warning. Callers
  passing hand-built marker sets that cross chromosomes will get loci rows
  with missing position metadata and no diagnostic.
- Whether `make_credible_sets()` can produce more than one `CS*` per locus is
  not mentioned in its own roxygen documentation or in
  `docs/dev/stblr_csr_finemap_overview.md`; both describe it as producing
  "credible sets" per locus without stating the possible one-to-many
  relationship documented above.
- No existing developer doc previously covered `R/credible_sets.R`'s internal
  helpers (`.stblr_extract_pip`, `.stblr_marker_map_from_Glist`,
  `.stblr_define_loci`, `.stblr_loci_from_sets`, `.extract_sparseLD_region_dense`)
  at all; this note is the first.

## Related Documents

- `docs/dev/stblr_csr_finemap_overview.md` — the broader CSR fine-mapping
  architecture, including how `make_credible_sets()`'s `locus_sets` output
  feeds into `extract_stblr_finemap_loci()` and `finemap_stblr_csr()`.
