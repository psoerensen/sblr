#ifndef SBLR_ST_CHAIN_UTILS_H
#define SBLR_ST_CHAIN_UTILS_H

#include <algorithm>
#include <cstdint>

inline unsigned int stblr_trait_seed(int seed, int trait) {
 return static_cast<unsigned int>(
  seed + 1000003 * (trait + 1)
 );
}

inline unsigned int stblr_chain_seed(int seed, int trait, int chain) {
 return static_cast<unsigned int>(
  seed + 1000003 * (trait + 1) + 9176 * (chain + 1)
 );
}

inline unsigned int stblr_seed_with_chain_base(int chain_seed, int trait) {
 return static_cast<unsigned int>(
  chain_seed + 1000003 * (trait + 1)
 );
}

inline int stblr_task_trait(int task, int nchains) {
 return task / nchains;
}

inline int stblr_task_chain(int task, int nchains) {
 return task % nchains;
}

inline int stblr_num_chain_tasks(int ntraits, int nchains) {
 return ntraits * nchains;
}

inline int stblr_num_threads_for_tasks(int ncores, int ntasks) {
 return std::max(1, std::min(ncores, ntasks));
}

#endif
