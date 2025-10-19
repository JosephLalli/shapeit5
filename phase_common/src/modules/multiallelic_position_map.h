#pragma once

#include <cstdint>
#include <vector>

namespace shapeit5 {
namespace modules {

struct PositionGroup {
  int32_t chrom_id;
  int32_t pos;
  std::vector<int> variant_indices;
};

class MultiallelicPositionMap {
 public:
  MultiallelicPositionMap();
  void clear();
  void build(const std::vector<int32_t>& chrom_ids,
             const std::vector<int32_t>& positions);
  const std::vector<PositionGroup>& groups() const;
  std::size_t size() const;

 private:
  std::vector<PositionGroup> groups_;
};

}  // namespace modules
}  // namespace shapeit5
