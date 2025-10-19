#include "multiallelic_position_map.h"

#include <algorithm>
#include <unordered_map>

namespace shapeit5 {
namespace modules {

namespace {
struct PositionKey {
  int32_t chrom_id;
  int32_t pos;
  bool operator==(const PositionKey& o) const noexcept {
    return chrom_id == o.chrom_id && pos == o.pos;
  }
};

struct PositionKeyHash {
  std::size_t operator()(const PositionKey& k) const noexcept {
    return (static_cast<std::size_t>(k.chrom_id) << 32) ^ static_cast<std::size_t>(k.pos);
  }
};
}  // namespace

MultiallelicPositionMap::MultiallelicPositionMap() = default;

void MultiallelicPositionMap::clear() { groups_.clear(); }

void MultiallelicPositionMap::build(const std::vector<int32_t>& chrom_ids,
                                    const std::vector<int32_t>& positions) {
  groups_.clear();
  const std::size_t n = positions.size();
  std::unordered_map<PositionKey, std::size_t, PositionKeyHash> idx;
  idx.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    PositionKey key{chrom_ids[i], positions[i]};
    auto it = idx.find(key);
    if (it == idx.end()) {
      PositionGroup g;
      g.chrom_id = key.chrom_id;
      g.pos = key.pos;
      g.variant_indices.push_back(static_cast<int>(i));
      groups_.push_back(std::move(g));
      idx.emplace(key, groups_.size() - 1);
    } else {
      groups_[it->second].variant_indices.push_back(static_cast<int>(i));
    }
  }
  std::sort(groups_.begin(), groups_.end(), [](const PositionGroup& a, const PositionGroup& b) {
    if (a.chrom_id != b.chrom_id) return a.chrom_id < b.chrom_id;
    return a.pos < b.pos;
  });
}

const std::vector<PositionGroup>& MultiallelicPositionMap::groups() const { return groups_; }

std::size_t MultiallelicPositionMap::size() const { return groups_.size(); }

}  // namespace modules
}  // namespace shapeit5

