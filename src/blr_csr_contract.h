#ifndef SBLR_CORE_BLR_CSR_CONTRACT_H
#define SBLR_CORE_BLR_CSR_CONTRACT_H

#include <cstddef>
#include <stdexcept>
#include <string>

namespace sblr {
namespace core {

// Describes an existing CSR resource; it does not own or read the CSR arrays.
// The storage owner must outlive every chain using the resource. All chains
// share the same immutable storage instance. Mutable effects, residuals,
// parameters, accumulators, random engines, and workspaces remain chain-owned.
struct CsrResourceSpec {
  std::string resource_id;
  std::size_t marker_count = 0;
  bool shared_read_only = true;
  bool per_chain_data = false;
  bool lifetime_exceeds_chains = true;
};

inline void validate_csr_resource(const CsrResourceSpec& resource) {
  if (resource.resource_id.empty()) {
    throw std::invalid_argument("data$csr$resource_id must not be empty");
  }
  if (resource.marker_count == 0) {
    throw std::invalid_argument("data$csr$marker_count must be positive");
  }
  if (!resource.shared_read_only) {
    throw std::invalid_argument("data$csr$shared_read_only must be true");
  }
  if (resource.per_chain_data) {
    throw std::invalid_argument("data$csr$per_chain_data must be false");
  }
  if (!resource.lifetime_exceeds_chains) {
    throw std::invalid_argument(
      "data$csr$lifetime_exceeds_chains must be true"
    );
  }
}

}  // namespace core
}  // namespace sblr

#endif
