# Permanent BLR test ownership

| Contract | Permanent owners | Tier |
|---|---|---|
| Public API, model semantics, output and raw schema | `test-blr-unified-public-contract.R`, `test-blr-model-semantics.R`, `test-stblr-raw-schema.R` | fast |
| BayesC/BayesR/BayesRC scientific behavior | model/backend tests named by model or operator | fast/full |
| Operator reductions | `test-blr-operator-reductions.R` | fast/extended |
| Selection-S and annotations | `test-stblr-selection-s.R`, annotation backend tests, MT BayesRC tests | fast/full |
| Multichain and RNG | `test-blr-unified-reproducibility.R`, `test-blr-extended-reproducibility.R` | fast/extended |
| Core and extended convergence | `test-blr-unified-convergence.R`, extended diagnostic owners | fast/full |
| Selected markers | `test-blr-selected-marker-diagnostics.R` | fast/full |
| Generated/native architecture | permanent scripts under `tools/audit` | fast/extended |

Installed tests own scientific and public behavior. Source-only checks are
limited to architecture that cannot be observed from an installed package.
Phase-numbered ownership and source-hash migration assertions are retired.
