#pragma once

#include <cstdint>
#include <vector>

// Simplified Li–Stephens transition scoring utilities.
// This module deliberately avoids tight coupling to internal HMM code
// and uses a standard approximation until tighter integration is added.

namespace shapeit5 {
namespace modules {

struct EdgeInfo {
  int neighbor_idx;
  double distance_cM;
  int m;
};

double transScore(bool stay, const EdgeInfo& e);

double scoreCandidateFlip(bool flip_k1,
                          int p,
                          int k1,
                          int k2,
                          int n,
                          const EdgeInfo& L,
                          const EdgeInfo& R,
                          const std::vector<uint8_t>& hapBit);

}  // namespace modules
}  // namespace shapeit5
