#include "../../phase_common/src/models/super_site_accessor.h"
#include <cassert>
#include <iostream>
#include <cstring>

// Simple test to verify super_site_accessor.h compiles and functions work
int main() {
	std::cout << "Testing super_site_accessor.h..." << std::endl;
	
	// Test 1: Create a SuperSite structure
	SuperSite site;
	site.global_site_id = 0;
	site.chr = 1;
	site.bp = 12345;
	site.n_alts = 3;
	site.panel_offset = 0;
	
	std::cout << "  SuperSite structure created successfully" << std::endl;
	
	// Test 2: Test unpacking functions with sample data
	// Create a packed buffer with some test data
	// Example: 4 haplotypes with codes: 0, 1, 2, 3 (packed as 0x10, 0x32)
	uint8_t packed_buffer[2];
	packed_buffer[0] = 0x10; // hap 0 = 0x0, hap 1 = 0x1
	packed_buffer[1] = 0x32; // hap 2 = 0x2, hap 3 = 0x3
	
	// Test single unpacking
	uint8_t code0 = unpackSuperSiteCode(packed_buffer, 0, 0);
	uint8_t code1 = unpackSuperSiteCode(packed_buffer, 0, 1);
	uint8_t code2 = unpackSuperSiteCode(packed_buffer, 0, 2);
	uint8_t code3 = unpackSuperSiteCode(packed_buffer, 0, 3);
	
	assert(code0 == 0x0);
	assert(code1 == 0x1);
	assert(code2 == 0x2);
	assert(code3 == 0x3);
	
	std::cout << "  Single code unpacking works correctly" << std::endl;
	
	// Test 3: Test batch unpacking
	uint8_t codes_out[4];
	unpackSuperSiteCodesBatch(packed_buffer, 0, 0, 4, codes_out);
	
	assert(codes_out[0] == 0x0);
	assert(codes_out[1] == 0x1);
	assert(codes_out[2] == 0x2);
	assert(codes_out[3] == 0x3);
	
	std::cout << "  Batch code unpacking works correctly" << std::endl;
	
	// Test 4: Test vectorized unpacking (currently just calls batch)
	uint8_t codes_vec[4];
	unpackSuperSiteCodesVectorized_PEXT(packed_buffer, 0, 4, codes_vec);
	
	assert(codes_vec[0] == 0x0);
	assert(codes_vec[1] == 0x1);
	assert(codes_vec[2] == 0x2);
	assert(codes_vec[3] == 0x3);
	
	std::cout << "  Vectorized code unpacking works correctly" << std::endl;
	
	// Test 5: Test aligned_vector32
	aligned_vector32<uint8_t> vec;
	vec.push_back(1);
	vec.push_back(2);
	vec.push_back(3);
	
	assert(vec.size() == 3);
	assert(vec[0] == 1);
	assert(vec[1] == 2);
	assert(vec[2] == 3);
	
	std::cout << "  aligned_vector32 works correctly" << std::endl;
	
	std::cout << "All tests passed!" << std::endl;
	return 0;
}
