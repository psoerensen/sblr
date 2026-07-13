# Frozen before editing src/st_cpg_omp_csr.cpp at commit
# 7df2d1ddec1c5f45be03e399dc7afb419afd56d8. Hashes cover complete raw and
# formatted objects after normalization of runtime seconds and the temporary
# CSR fixture path. Fixture metadata and construction live beside this file.
phase2_reference_expected_hashes <- list(
  one_trait_one_chain_one_core = c(
    raw = "1bdbead1dbdc2f3b5f9c3576d2f78669",
    fit = "b2abcc834a085970048880aaefc9ee98"
  ),
  one_trait_two_chains_one_core = c(
    raw = "6c5667d2ab9a1aa9fdd1e4d3bd289122",
    fit = "068d142a01eb2ce9d01e48eefcda4e09"
  ),
  one_trait_two_chains_two_cores = c(
    raw = "6c5667d2ab9a1aa9fdd1e4d3bd289122",
    fit = "11d89e8e46d4908f088a3099e64ad39b"
  ),
  multiple_traits = c(
    raw = "e63936437a33fd3f2b479401a2520fac",
    fit = "b91f69491bf599b232252132166e8a45"
  ),
  explicit_chain_seeds = c(
    raw = "fc9847ed816c0cc2e41b1ded38dccc0a",
    fit = "4794db2fb163d1cc99c03d6719f2cd91"
  ),
  keep_chains = c(
    raw = "871871a5687ee5412b4d68b7b62fe86d",
    fit = "aa7c7b183748f14020c65f2ef4e4701d"
  ),
  fixed_selection_s = c(
    raw = "1c50036f9bca4b12456d97d333bbedc0",
    fit = "10bb707927928dd8b8cba3594a89ff03"
  )
)
