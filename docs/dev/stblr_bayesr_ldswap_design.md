# BayesR CSR LD-Swap Design

## Executive Summary

Exact CSR BayesR should not copy the BayesC LD-swap move verbatim. BayesC
swaps an active binary inclusion/effect state:

```text
d_j in {0, 1}
b_j
```

BayesR marker state is richer:

```text
component_j in {0, 1, ..., K - 1}
b_j
prior variance vb * mixture_var[component_j]
```

The first statistically valid BayesR LD-swap should move the full marker state
`(component, b)` between LD neighbors, not the effect alone. For the plain exact
CSR BayesR backend with global mixture probabilities and shared component
variance multipliers, component and effect prior terms cancel under a full-state
swap. The Metropolis-Hastings ratio then uses the same summary-statistic
likelihood change and the same asymmetric proposal correction as BayesC. This
statement is specific to the current plain CSR BayesR model; marker-specific
priors or annotation-informed mixture probabilities would add non-cancelling
prior terms.

Recommended first implementation scope is **Scope B: active/null relocation**
using a full BayesR state move:

```text
j active: (component_j > 0, b_j)
k null:   (0, 0)

after proposal:
j null:   (0, 0)
k active: (old component_j, old b_j)
```

This matches the current BayesC proposal topology, keeps diagnostics and
proposal correction reusable, and avoids introducing active/active state swaps
before the active/null move is tested.

## Current BayesC LD-Swap Mechanism

Implementation file:

```text
src/st_cpg_omp_csr.cpp
```

The native arguments controlling LD-swap are:

- `updateLDswap`
- `ld_swap_prob`
- `ld_swap_r2`
- `ld_swap_max_friends`
- `ld_swap_moves`

The R wrapper validates the same controls in `.validate_ld_swap_args()` in
`R/sparse_ld_bed_helper.R`. Scheduled CSR BayesC rejects `updateLDswap = TRUE`.

LD friends are identified in `build_ld_swap_friends_st_csr()`. The function
walks the exact CSR LD rows, computes:

```text
r2(i, j) = xij^2 / (ww_i * ww_j)
```

keeps pairs with `r2 >= ld_swap_r2`, stores both directions, sorts each friend
list by descending `r2`, removes adjacent duplicate marker ids, and truncates
each row to `ld_swap_max_friends`.

Candidate active markers are selected by `collect_ld_swap_candidates()`. A
marker is a candidate when:

```text
d_i > 0
and at least one stored LD friend has d_friend == 0
```

`attempt_ld_swap_st_csr()` then samples:

1. one candidate active marker `j` uniformly among current candidates;
2. one excluded LD friend `k` uniformly among the excluded friends of `j`.

The current BayesC move is active/null relocation, not a general active/active
exchange. It requires:

```text
d_j > 0
d_k == 0
b_j != 0
```

and proposes:

```text
j: (d_j, b_j) -> (0, 0)
k: (0, 0)     -> (1, b_j)
```

Thus BayesC moves inclusion and effect together. It does not resample the
effect during the swap and does not change residual variance, effect variance,
or mixture probability parameters inside the move.

The residual state is updated by `set_marker_effect_st_csr()`. For a marker
effect change `diff = b_new - b_old`, it applies:

```text
r_i -= ww_i * diff
r_neighbor -= xij * diff
```

for the marker and all stored LD neighbors, preserving:

```text
r = X'y - X'Xb
```

BayesC computes the old and new residual SSE with:

```text
SSE = y'y - sum_i b_i * (r_i + wy_i)
```

which is equivalent to:

```text
y'y - 2 b'X'y + b'X'Xb
```

because `r = X'y - X'Xb`.

The BayesC log acceptance ratio is:

```text
log_alpha =
  -0.5 * (SSE_new - SSE_old) / vei
  + log_q_reverse
  - log_q_forward
```

where:

```text
log_q_forward = -log(n_forward_candidates) - log(n_forward_excluded_friends)
log_q_reverse = -log(n_reverse_candidates) - log(n_reverse_excluded_friends)
```

No explicit BayesC prior terms enter the current ratio. Under the implemented
active/null relocation with marker-exchangeable `pi` and a shared active normal
effect prior, the inclusion prior and effect prior cancel between old and new
states. No `updateB`, `updateE`, or `updatePi` step is part of the proposal.

Diagnostics are accumulated per trait-chain task:

```text
ld_swap_attempted_task
ld_swap_accepted_task
```

and aggregated per trait:

```text
ld_swap_attempted_vec
ld_swap_accepted_vec
```

The raw C++ output stores trait-level diagnostics in slot 22 as:

```text
attempted
accepted
acceptance_rate
sentinel = 1
```

When `keep_chains = TRUE`, slot 31 stores chain-major blocks of the same four
values. `.format_stblr_fit()` exposes:

```r
fit$ld_swap
fit$ld_swap_chains
fit$chains[[trait]][[chain]]$ld_swap
```

with user-visible columns:

```text
attempted
accepted
acceptance_rate
```

## Current BayesR CSR State

Implementation file:

```text
src/st_cpg_omp_csr_bayesr.cpp
```

Component states are stored in an `arma::Row<int> comp_t` parallel to the
effect vector `b_t`.

Component `0` is the null component. The backend validates:

```text
mixture_var[0] == 0
mixture_var[k] > 0 for k > 0
```

Null-component effects are forced to zero by construction in
`sampleBetaR_ST_csr()`: when the sampled component has `mixture_var[k] <= 0`,
`b_new` remains `0.0`. `ensure_null_effects_bayesr_ST_csr()` verifies before
residual-variance updates that no marker with `comp == 0` has a nonzero effect.

For each marker update, BayesR computes a current conditional score:

```text
score = r_i + ww_i * b_i
```

and a log posterior mass for every component. The null component uses only
`log(pi_0)`. A non-null component `k > 0` uses prior variance:

```text
vbk = vb_t * mixture_var[k]
```

and log Bayes factor:

```text
0.5 * log(vei / (vei + ww_i * vbk))
+ 0.5 * score^2 * vbk / (vei * (vei + ww_i * vbk))
```

After sampling the component from these log probabilities, the effect is sampled
from its conditional normal when `k > 0`:

```text
lhs = ww_i + vei / vbk
mean = score / lhs
sd = sqrt(vei / lhs)
```

The residual update uses the same sign convention as BayesC:

```text
diff = b_new - b_old
r_i -= ww_i * diff
r_neighbor -= xij * diff
```

`comp_prob` is accumulated after burn-in/thinning by incrementing
`comp_prob_t[i, comp_i]`. The standard BayesR `dm` is accumulated as:

```text
dm_i += 1{comp_i > 0}
```

The formatter then enforces:

```r
fit$dm = 1 - fit$comp_prob[[trait]][, "component_0"]
```

`updateB` uses active components only:

```text
sum_{i: comp_i > 0} b_i^2 / mixture_var[comp_i]
```

`updateE` rebuilds `r`, verifies null effects, computes strict SSE diagnostics,
and then calls shared `sampleE_ST_csr()`. `updatePi` applies a global
Dirichlet update:

```text
pi ~ Dirichlet(alpha + component_counts)
```

When `keep_chains = TRUE`, BayesR currently returns compact per-chain:

- `dm`
- `bm`
- `comp_prob`
- `dm_component_mean`
- `final_pi`
- `mean_pi`
- `vbs`, `vgs`, `ves`, `vle`, `vld`
- `updateE_diagnostics`

Trait-level chain stability summaries include:

- `dm_sd`, `dm_min`, `dm_max`
- `bm_sd`, `bm_min`, `bm_max`

## BayesR Swap State

The candidate options are:

### Option A: Effect Only

```text
j: b_j
k: b_k
components unchanged
```

This is not a valid first design. The BayesR effect prior depends on the
component. Moving an effect without its component can leave a nonzero effect in
component 0 or can evaluate an effect under the wrong component variance.

### Option B: Component Plus Effect

```text
j: (component_j, b_j)
k: (component_k, b_k)
```

This is the natural BayesR analogue of the current BayesC move. It preserves the
effect with the component that defines its prior variance. For an active/null
move, it relocates the active BayesR effect and its component assignment to an
LD-correlated null marker.

### Option C: New Component/Effect at Target

```text
remove j
sample or propose component/effect for k
```

This could mix better, but it requires a more complex proposal density and an
acceptance ratio that accounts for the sampled component and effect density.
It is not the safest first implementation.

### Option D: Relocation Plus Local Refresh

```text
move component/effect from j to k
then optionally Gibbs refresh one or both markers
```

This may be useful later, but the refresh step must be treated either as part of
the proposal kernel or as a separate Gibbs step with its own invariance
argument. It should not be bundled into the first LD-swap change.

Recommended first state is **Option B**, restricted initially to active/null
relocation.

## MH Ratio for a Full-State Swap

Let the current BayesR state for one trait-chain be:

```text
b
z = component
r = X'y - X'Xb
```

A full-state proposal swaps `(z_j, b_j)` and `(z_k, b_k)`. For the recommended
first active/null implementation:

```text
z_j > 0, b_j != 0
z_k = 0, b_k = 0
```

and the proposal is:

```text
z'_j = 0,   b'_j = 0
z'_k = z_j, b'_k = b_j
```

The generic MH ratio is:

```text
log_accept =
  log p(y | b', ve, other fixed parameters)
- log p(y | b,  ve, other fixed parameters)
+ log p(b', z' | vb, pi, mixture_var)
- log p(b,  z  | vb, pi, mixture_var)
+ log q(b, z | b', z')
- log q(b', z' | b, z)
```

### Likelihood Term

Use the same summary-statistic identity as BayesC:

```text
SSE(b) = y'y - sum_i b_i * (r_i + wy_i)
       = y'y - 2 b'X'y + b'X'Xb
```

with `r` updated to match each proposed `b`. Conditional on the current
adjusted residual variance `vei`, the likelihood contribution is:

```text
-0.5 * (SSE_new - SSE_old) / vei
```

This matches the current BayesC code.

### Effect Prior Term

For `k > 0`:

```text
b_i | z_i = k, vb ~ N(0, vb * mixture_var[k])
```

The null component has:

```text
z_i = 0 => b_i = 0
```

Under a full-state swap with global `mixture_var`, the same non-null pair
`(z_j, b_j)` is moved from marker `j` to marker `k`. The normal prior density
is unchanged, so it cancels. The null point mass also cancels.

If a future BayesR variant uses marker-specific variance multipliers, annotation
classes that alter effect variance, or marker-specific scaling beyond the
current shared `ww` likelihood scaling, this cancellation must be revisited.

### Component Prior Term

Current plain CSR BayesR uses a global trait-chain vector `pi`. Under a
full-state swap, the multiset of component assignments is unchanged:

```text
{z_j, z_k} = {z'_j, z'_k}
```

so:

```text
log pi[z'_j] + log pi[z'_k]
- log pi[z_j] - log pi[z_k] = 0
```

If marker-specific component probabilities are introduced, as in
annotation-informed SBayesRC-style models, this term does not cancel. Then the
ratio must include:

```text
log pi_k[z_j] + log pi_j[z_k]
- log pi_j[z_j] - log pi_k[z_k]
```

using the appropriate marker-specific probability model.

### Proposal Ratio

The current BayesC proposal is not automatically symmetric because candidate
sets and excluded-friend counts can differ after relocation. BayesC therefore
computes:

```text
log_q_forward = -log(n_candidates_old) - log(n_excluded_friends_of_j_old)
log_q_reverse = -log(n_candidates_new) - log(n_excluded_friends_of_k_new)
```

and adds:

```text
log_q_reverse - log_q_forward
```

BayesR should reuse this correction for the active/null first scope, replacing
`d_i > 0` with `component_i > 0`.

### Plain CSR BayesR First-Scope Ratio

For active/null full-state relocation in current plain CSR BayesR, the ratio is:

```text
log_alpha =
  -0.5 * (SSE_new - SSE_old) / vei
  + log_q_reverse
  - log_q_forward
```

This is algebraically the same as BayesC only because the BayesR component and
effect are moved together and current BayesR priors are marker-exchangeable.

## Reusable BayesC Code

Reusable with light generalization:

- `LDLDFriends`
- `build_ld_swap_friends_st_csr()`
- LD friend thresholding by `ld_swap_r2`
- per-row truncation by `ld_swap_max_friends`
- trigger scheduling with `ld_swap_prob` and `ld_swap_moves`
- proposal diagnostics and aggregation
- residual update machinery and SSE identity
- asymmetric proposal correction
- chain-level diagnostic layout and R formatting pattern

Not directly reusable:

- logic assuming binary `d`
- checks such as `d_j_old <= 0 || d_k_old != 0`
- setting `d(k) = 1` without carrying a component id
- any future BayesR variant with marker-specific priors, because prior terms
  would need explicit evaluation

Recommended helper split for implementation:

- Generalize candidate collection around an "active" predicate
  `component_i > 0`.
- Add a BayesR-specific setter that updates `b` and `component` together.
- Keep the BayesC move unchanged; do not alter `src/st_cpg_omp_csr.cpp` except
  by extracting shared helpers only if that extraction is mechanically small
  and separately tested.

## API Design

Future public arguments for `stblr_csr_bayesr()` should match BayesC:

```r
updateLDswap = FALSE
ld_swap_prob = 0.05
ld_swap_r2 = 0.8
ld_swap_max_friends = 50L
ld_swap_moves = 1L
```

Rules:

- Default remains `updateLDswap = FALSE`.
- `scheduled = TRUE, updateLDswap = TRUE` remains unsupported.
- LD-swap tuning arguments should be validated with the same rules as BayesC.
- `updateLDswap = TRUE` requires exact CSR LD availability through `ld_prefix`
  or `Glist$sparseLD$prefix`.
- BayesC behavior and BayesC defaults must not change.

`updateE` can be allowed in principle because the swap conditions on the current
`vei`, updates `r` exactly, and the BayesR backend already rebuilds `r` before
`sampleE_ST_csr()`. For the first implementation pass, the safer scope is to
support and test:

```r
updateLDswap = TRUE, updateE = FALSE
```

first. If `updateE = TRUE` is accepted in the same pass, it should require an
explicit tiny-fixture test with strict `updateE_diagnostics` and should rebuild
`r` before `updateE` exactly as the current backend does. If that validation is
not included, reject:

```r
updateLDswap = TRUE, updateE = TRUE
```

with a clear message until a follow-up test covers the combination.

## Diagnostics Design

Mirror BayesC for the first implementation:

```r
fit$ld_swap
```

with columns:

```text
attempted
accepted
acceptance_rate
```

When `keep_chains = TRUE`, also expose:

```r
fit$ld_swap_chains
fit$chains[[trait]][[chain]]$ld_swap
```

No BayesR-specific diagnostic columns are required for the first pass.
Additional counters may be useful later:

```text
attempted_active_null
accepted_active_null
attempted_active_active
accepted_active_active
```

but they should not be added until there is a non-active/null proposal or a
debugging need.

Because CSR BayesR returns a named `Rcpp::List`, the cleanest native layout is a
named `ld_swap` matrix or data frame equivalent plus optional chain entries in
the existing `chains` list. It does not need to mimic BayesC positional slot 22.

## Testing Plan

Future implementation tests should be written before enabling the API.

1. API validation:
   - invalid `updateLDswap`
   - invalid `ld_swap_prob`
   - invalid `ld_swap_r2`
   - invalid `ld_swap_max_friends`
   - invalid `ld_swap_moves`
   - unsupported `scheduled = TRUE, updateLDswap = TRUE`
   - unsupported `updateLDswap = TRUE, updateE = TRUE` if first implementation
     does not test that combination

2. Tiny CSR fixture:
   - `updateLDswap = FALSE` still works
   - `updateLDswap = TRUE` works with `updateE = FALSE`
   - `fit$ld_swap` is present
   - `attempted >= accepted >= 0`
   - `acceptance_rate` is in `[0, 1]`

3. Multi-chain:
   - `nchains = 2`
   - aggregated attempted/accepted counts equal the sum of chain counts

4. `keep_chains`:
   - per-chain LD-swap diagnostics are present
   - existing `dm`, `bm`, `comp_prob`, `final_pi`, and `mean_pi` chain
     aggregation remains valid

5. BayesR conventions:
   - `dm = 1 - component_0`
   - component probabilities sum to 1
   - `check_stblr_backend_consistency()` passes
   - fine-mapping extractor compatibility remains intact

6. Controlled LD toy example:
   - construct two highly correlated markers
   - verify attempted swaps occur with `ld_swap_prob = 1` and
     `ld_swap_moves > 0`
   - if possible, initialize one active and one null marker using
     `use_comp_init = TRUE`, `comp_init`, `use_r_init = FALSE`, and verify that
     accepted moves can occur under a simple state

7. Regression:
   - existing BayesC LD-swap tests remain unchanged and pass
   - existing BayesR CSR backend tests remain unchanged and pass with
     `updateLDswap = FALSE`

Avoid expensive MCMC. Use the existing tiny CSR fixtures as the template.

## Recommended Implementation Scope

Recommended first scope is **Scope B: active/null relocation only**, implemented
as a full BayesR state move:

```text
(component, b) moves from an active marker to a null LD friend
```

Acceptance ratio:

```text
log_alpha =
  -0.5 * (SSE_new - SSE_old) / vei
  + log_q_reverse
  - log_q_forward
```

with the explicit implementation invariant that prior terms cancel only because
the full component/effect pair is moved and current plain CSR BayesR has global
`pi` and global `mixture_var`.

Scope A, full active/active state swaps, can follow after Scope B. It is not
much harder algebraically under global priors, but it changes proposal topology
and requires more careful testing of reverse candidates. Scope C is unnecessary
unless the implementation pass needs exploratory diagnostics. Scope D is not
recommended because the active/null derivation is sufficiently clear.

## Suggested Implementation Prompt

Implement BayesR CSR LD-swap active/null relocation only.

- Do not change BayesC behavior.
- Do not implement scheduled CSR BayesR.
- Add BayesR LD-swap arguments to `stblr_csr_bayesr()` with BayesC defaults.
- Validate arguments with BayesC-compatible rules.
- In `src/st_cpg_omp_csr_bayesr.cpp`, add a BayesR-specific active/null
  full-state relocation:
  - candidate marker has `component > 0`;
  - friend marker has `component == 0`;
  - move `(component, b)` together;
  - update `r` through the exact CSR residual update;
  - compute old/new SSE with the existing BayesR/BayesC identity;
  - include BayesC-style forward/reverse proposal correction;
  - condition on current `vei`.
- Return trait-level `ld_swap` diagnostics and per-chain diagnostics when
  `keep_chains = TRUE`.
- Add focused tests for API validation, tiny CSR execution, multi-chain
  aggregation, chain diagnostics, BayesR component conventions, and BayesC
  regression.
- Initially support `updateLDswap = TRUE, updateE = FALSE`; only allow
  `updateE = TRUE` in the same pass if a dedicated strict diagnostic test is
  included.
