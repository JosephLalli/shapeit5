/*******************************************************************************
 * Supersite utility helpers.
 *
 * These helpers expose compact per-supersite allele codes for a chosen set of
 * conditioning haplotypes. They rely on the metadata curated in
 * variant_map::buildSupersites and keep temporary buffers aligned for SIMD
 * consumers.
 ******************************************************************************/

#ifndef SHAPEIT5_SUPERSITE_ACCESSOR_H
#define SHAPEIT5_SUPERSITE_ACCESSOR_H

#include <containers/variant_map.h>

#include <boost/align/aligned_allocator.hpp>

#include <cassert>
#include <cstdint>
#include <vector>

class supersite_accessor {
public:
	using aligned_u8 = std::vector<
		uint8_t,
		boost::alignment::aligned_allocator<uint8_t, 32>
	>;

	supersite_accessor() = default;

	void reset(const variant_map * map, const std::vector<unsigned int> * cond_haps) {
		map_ = map;
		cond_haps_ = cond_haps;
		buffer_.clear();
	}

	bool ready() const {
		return map_ != nullptr && cond_haps_ != nullptr;
	}

	bool is_super_site(uint32_t site_idx) const {
		return ready() && site_idx < map_->supersites.size() &&
			map_->supersites[site_idx].is_super_site;
	}

	const uint8_t * load(uint32_t site_idx) {
		assert(ready());
		const supersite_desc & desc = map_->supersites[site_idx];
		assert(desc.is_super_site);
		buffer_.resize(cond_haps_->size());
		const uint8_t * codes = map_->supersite_codes.data() + desc.panel_offset;
		for (size_t i = 0; i < cond_haps_->size(); ++i) {
			const unsigned int hap_index = (*cond_haps_)[i];
			buffer_[i] = codes[hap_index];
		}
		return buffer_.data();
	}

	const aligned_u8 & buffer() const { return buffer_; }

private:
	const variant_map * map_ = nullptr;
	const std::vector<unsigned int> * cond_haps_ = nullptr;
	aligned_u8 buffer_;
};

#endif // SHAPEIT5_SUPERSITE_ACCESSOR_H

