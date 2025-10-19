#include "transition_scorer.h"

#include <algorithm>
#include <cmath>

namespace shapeit5 {
namespace modules {

double transScore(bool stay, const EdgeInfo& e) {
  if (e.neighbor_idx < 0) return 0.0;
  // Simple Li–Stephens approximation with Ne=10,000
  // Convert cM to recombination fraction proxy then to probability
  const double Ne = 10000.0;
  const double dist_cm = std::max(e.distance_cM, 1e-8);
  const int m = std::max(e.m, 2); // avoid division by zero

  // Approximate rho from genetic distance. Keep within (0,1)
  double rho = 4.0 * Ne * (dist_cm / 100.0);
  rho = rho / (1.0 + rho);
  rho = std::min(std::max(rho, 1e-12), 1.0 - 1e-12);

  if (stay) {
    return std::log(1.0 - rho);
  } else {
    return std::log(rho / static_cast<double>(m - 1));
  }
}

double scoreCandidateFlip(bool flip_k1,
                          int p,
                          int k1,
                          int k2,
                          int n,
                          const EdgeInfo& L,
                          const EdgeInfo& R,
                          const std::vector<uint8_t>& hapBit) {
  const bool h_p = (p >= 0) ? static_cast<bool>(hapBit[static_cast<std::size_t>(p)]) : false;
  const bool h_k1_new = flip_k1 ? !static_cast<bool>(hapBit[static_cast<std::size_t>(k1)])
                                : static_cast<bool>(hapBit[static_cast<std::size_t>(k1)]);
  const bool h_k2_new = flip_k1 ? static_cast<bool>(hapBit[static_cast<std::size_t>(k2)])
                                : !static_cast<bool>(hapBit[static_cast<std::size_t>(k2)]);
  const bool h_n = (n >= 0) ? static_cast<bool>(hapBit[static_cast<std::size_t>(n)]) : false;

  const bool stayL = (p >= 0) ? (h_p == h_k1_new) : true;
  const bool stayR = (n >= 0) ? (h_k2_new == h_n) : true;

  return transScore(stayL, L) + transScore(stayR, R);
}

}  // namespace modules
}  // namespace shapeit5
