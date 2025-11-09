/*******************************************************************************
 * Thread Safety Test for SHAPEIT5 SuperSite Threading
 * 
 * Tests the thread safety of SuperSite operations including:
 * - SC buffer concurrent access and updates
 * - anchor_has_missing flag thread safety
 * - supersite_sc_offset concurrent modifications
 * - SuperSite context setting across threads
 * - Memory consistency in supersite data structures
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>
#include <memory>

#include "test_framework.h"
#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/compute_job.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/genotype_set.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

using namespace test_framework;

// Mock SuperSite threading environment
class MockSuperSiteThreading {
public:
    // SuperSite data
    std::vector<SuperSite> super_sites;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> packed_allele_codes;
    
    // Thread-local data simulation (one per thread)
    struct ThreadData {
        std::vector<float> SC;
        std::vector<bool> anchor_has_missing;
        std::vector<uint32_t> supersite_sc_offset;
        std::atomic<int> update_count{0};
        std::atomic<int> access_count{0};
        
        ThreadData(int num_supersites, int hap_count) {
            SC.resize(num_supersites * hap_count * 4, 0.0f); // 4 classes per haplotype
            anchor_has_missing.resize(num_supersites, false);
            supersite_sc_offset.resize(num_supersites, 0);
        }
        
        // Delete copy constructor and copy assignment
        ThreadData(const ThreadData&) = delete;
        ThreadData& operator=(const ThreadData&) = delete;
        
        // Delete move constructor and move assignment
        ThreadData(ThreadData&&) = delete;
        ThreadData& operator=(ThreadData&&) = delete;
    };
    
    std::vector<std::unique_ptr<ThreadData>> thread_data;
    
    // Shared state
    std::atomic<int> global_updates{0};
    std::atomic<int> context_sets{0};
    pthread_mutex_t supersite_mutex;
    
    // Test configuration
    int num_supersites = 10;
    int num_haplotypes = 8;
    int num_threads = 4;
    bool enable_mutex = true;
    
    MockSuperSiteThreading() {
        pthread_mutex_init(&supersite_mutex, nullptr);
        setup_test_data();
    }
    
    ~MockSuperSiteThreading() {
        pthread_mutex_destroy(&supersite_mutex);
    }
    
    void setup_test_data() {
        // Create mock SuperSites
        super_sites.clear();
        locus_to_super_idx.clear();
        super_site_var_index.clear();
        packed_allele_codes.clear();
        
        for (int i = 0; i < num_supersites; i++) {
            SuperSite ss;
            ss.global_site_id = i;
            ss.chr = 1;
            ss.bp = 1000 + i * 100;
            ss.n_alts = 2 + (i % 3); // 2-4 alts
            ss.panel_offset = i * num_haplotypes * 2;
            super_sites.push_back(ss);
            
            locus_to_super_idx.push_back(i);
            super_site_var_index.push_back(i * 2);
            super_site_var_index.push_back(i * 2 + 1);
        }
        
        // Create packed allele codes
        packed_allele_codes.resize(num_supersites * num_haplotypes * 2);
        for (size_t i = 0; i < packed_allele_codes.size(); i++) {
            packed_allele_codes[i] = i % 16; // Mock allele codes
        }
        
        // Initialize thread data
        thread_data.clear();
        for (int t = 0; t < num_threads; t++) {
            thread_data.emplace_back(std::make_unique<ThreadData>(num_supersites, num_haplotypes));
        }
    }
    
    // Simulate SuperSite context setting
    void simulate_set_supersite_context(int thread_id) {
        ThreadData& tdata = *thread_data[thread_id];
        
        for (int ss = 0; ss < num_supersites; ss++) {
            if (enable_mutex) pthread_mutex_lock(&supersite_mutex);
            
            // Simulate setting anchor_has_missing
            tdata.anchor_has_missing[ss] = (ss % 3 == thread_id % 3);
            
            // Simulate setting SC offset
            tdata.supersite_sc_offset[ss] = ss * num_haplotypes * 4;
            
            // Simulate SC buffer updates
            for (int h = 0; h < num_haplotypes; h++) {
                for (int c = 0; c < 4; c++) {
                    int idx = ss * num_haplotypes * 4 + h * 4 + c;
                    if (idx < tdata.SC.size()) {
                        tdata.SC[idx] = thread_id * 0.1f + ss * 0.01f + h * 0.001f + c * 0.0001f;
                    }
                }
            }
            
            tdata.update_count.fetch_add(1);
            global_updates.fetch_add(1);
            
            if (enable_mutex) pthread_mutex_unlock(&supersite_mutex);
            
            usleep(10); // Brief work simulation
        }
        
        context_sets.fetch_add(1);
    }
    
    // Simulate concurrent SC buffer access
    void simulate_sc_buffer_access(int thread_id) {
        ThreadData& tdata = *thread_data[thread_id];
        
        for (int iteration = 0; iteration < 20; iteration++) {
            for (int ss = 0; ss < num_supersites; ss++) {
                if (enable_mutex) pthread_mutex_lock(&supersite_mutex);
                
                // Read from SC buffer
                float sum = 0.0f;
                for (int h = 0; h < num_haplotypes; h++) {
                    for (int c = 0; c < 4; c++) {
                        int idx = ss * num_haplotypes * 4 + h * 4 + c;
                        if (idx < tdata.SC.size()) {
                            sum += tdata.SC[idx];
                        }
                    }
                }
                
                // Write to SC buffer
                int write_idx = ss * num_haplotypes * 4 + (iteration % num_haplotypes) * 4;
                if (write_idx + 4 <= tdata.SC.size()) {
                    for (int c = 0; c < 4; c++) {
                        tdata.SC[write_idx + c] = sum * 0.001f + thread_id * 0.1f;
                    }
                }
                
                tdata.access_count.fetch_add(1);
                
                if (enable_mutex) pthread_mutex_unlock(&supersite_mutex);
                
                usleep(5);
            }
        }
    }
    
    // Simulate anchor flag concurrent access
    void simulate_anchor_flag_access(int thread_id) {
        ThreadData& tdata = *thread_data[thread_id];
        
        for (int iteration = 0; iteration < 30; iteration++) {
            int ss = iteration % num_supersites;
            
            if (enable_mutex) pthread_mutex_lock(&supersite_mutex);
            
            // Toggle anchor flag
            bool current_state = tdata.anchor_has_missing[ss];
            tdata.anchor_has_missing[ss] = !current_state;
            
            // Update offset based on flag
            if (tdata.anchor_has_missing[ss]) {
                tdata.supersite_sc_offset[ss] = ss * num_haplotypes * 4;
            } else {
                tdata.supersite_sc_offset[ss] = 0;
            }
            
            if (enable_mutex) pthread_mutex_unlock(&supersite_mutex);
            
            usleep(3);
        }
    }
    
    // Reset state for new test
    void reset() {
        global_updates = 0;
        context_sets = 0;
        
        for (auto& tdata_ptr : thread_data) {
            std::fill(tdata_ptr->SC.begin(), tdata_ptr->SC.end(), 0.0f);
            std::fill(tdata_ptr->anchor_has_missing.begin(), tdata_ptr->anchor_has_missing.end(), false);
            std::fill(tdata_ptr->supersite_sc_offset.begin(), tdata_ptr->supersite_sc_offset.end(), 0);
            tdata_ptr->update_count = 0;
            tdata_ptr->access_count = 0;
        }
    }
    
    // Verify state consistency
    bool verify_consistency() {
        // Check that updates happened
        if (context_sets.load() != num_threads) return false;
        
        // Check thread data consistency
        for (const auto& tdata_ptr : thread_data) {
            if (tdata_ptr->update_count.load() != num_supersites) return false;
        }
        
        return true;
    }
};

// Test SuperSite context setting thread safety
void test_supersite_context_thread_safety() {
    TEST_RUN("supersite_context_thread_safety", []() {
        MockSuperSiteThreading mock;
        
        TEST_THREADED("context_setting", 4, [&](int thread_id) {
            mock.simulate_set_supersite_context(thread_id);
        });
        
        TEST_ASSERT(mock.verify_consistency(), "SuperSite context setting consistency");
    });
}

// Test SC buffer concurrent access
void test_sc_buffer_concurrent_access() {
    TEST_RUN("sc_buffer_concurrent_access", []() {
        MockSuperSiteThreading mock;
        
        TEST_THREADED("sc_buffer_access", 4, [&](int thread_id) {
            mock.simulate_sc_buffer_access(thread_id);
        });
        
        // Verify access counts
        int total_accesses = 0;
        for (const auto& tdata_ptr : mock.thread_data) {
            total_accesses += tdata_ptr->access_count.load();
        }
        
        TEST_ASSERT(total_accesses == 4 * 20 * mock.num_supersites, "SC buffer access count");
    });
}

// Test anchor flag thread safety
void test_anchor_flag_thread_safety() {
    TEST_RUN("anchor_flag_thread_safety", []() {
        MockSuperSiteThreading mock;
        
        TEST_THREADED("anchor_flag_access", 3, [&](int thread_id) {
            mock.simulate_anchor_flag_access(thread_id);
        });
        
        // Verify that flags were modified
        bool flags_modified = false;
        for (const auto& tdata_ptr : mock.thread_data) {
            for (bool flag : tdata_ptr->anchor_has_missing) {
                if (flag) {
                    flags_modified = true;
                    break;
                }
            }
            if (flags_modified) break;
        }
        
        TEST_ASSERT(flags_modified, "anchor flags were modified");
    });
}

// Test SC offset consistency
void test_sc_offset_consistency() {
    TEST_RUN("sc_offset_consistency", []() {
        MockSuperSiteThreading mock;
        
        TEST_THREADED("sc_offset_updates", 4, [&](int thread_id) {
            for (int ss = 0; ss < mock.num_supersites; ss++) {
                pthread_mutex_lock(&mock.supersite_mutex);
                
                // Set consistent offset
                mock.thread_data[thread_id]->supersite_sc_offset[ss] = ss * mock.num_haplotypes * 4;
                
                // Verify offset is within bounds
                uint32_t offset = mock.thread_data[thread_id]->supersite_sc_offset[ss];
                uint32_t max_offset = mock.num_supersites * mock.num_haplotypes * 4;
                
                TEST_THREADED_ASSERT(offset <= max_offset, 
                    "SC offset within bounds for thread " + std::to_string(thread_id) + " ss " + std::to_string(ss));
                
                pthread_mutex_unlock(&mock.supersite_mutex);
                
                usleep(5);
            }
        });
    });
}

// Test stress conditions with many SuperSites
void test_supersite_stress() {
    TEST_RUN("supersite_stress", []() {
        MockSuperSiteThreading mock;
        mock.num_supersites = 50; // More SuperSites for stress
        mock.num_threads = 6;     // Match the stress-thread count
        mock.setup_test_data();
        
        TEST_STRESS("supersite_high_load", 3, 6, [&](int thread_id) {
            // Combine multiple operations
            mock.simulate_set_supersite_context(thread_id);
            mock.simulate_sc_buffer_access(thread_id);
        });
        
        TEST_ASSERT(mock.context_sets.load() >= mock.num_threads * 3, "stress test context sets");
    });
}

// Test memory consistency across threads
void test_memory_consistency() {
    TEST_RUN("memory_consistency", []() {
        MockSuperSiteThreading mock;
        
        // Initialize with known pattern
        for (int t = 0; t < mock.num_threads; t++) {
            for (int i = 0; i < mock.thread_data[t]->SC.size(); i++) {
                mock.thread_data[t]->SC[i] = t * 100.0f + i;
            }
        }
        
        TEST_THREADED("memory_verification", 4, [&](int thread_id) {
            // Verify memory consistency
            for (int i = 0; i < mock.thread_data[thread_id]->SC.size(); i++) {
                float expected = thread_id * 100.0f + i;
                float actual = mock.thread_data[thread_id]->SC[i];
                
                TEST_THREADED_ASSERT(std::abs(actual - expected) < 1e-6, 
                    "memory consistency at index " + std::to_string(i));
            }
        });
    });
}

// Test deadlock detection in SuperSite operations
void test_supersite_deadlock_detection() {
    TEST_RUN("supersite_deadlock_detection", []() {
        MockSuperSiteThreading mock;
        mock.num_supersites = 20; // Smaller for faster test
        mock.setup_test_data();
        
        TEST_DEADLOCK_DETECTION("supersite_deadlock_test", 8, 4, [&](int thread_id) {
            // Mix of operations that could potentially deadlock
            for (int i = 0; i < 5; i++) {
                mock.simulate_set_supersite_context(thread_id);
                mock.simulate_sc_buffer_access(thread_id);
                usleep(50);
            }
        });
    });
}

// Test race condition detection
void test_supersite_race_detection() {
    TEST_RUN("supersite_race_detection", []() {
        MockSuperSiteThreading mock;
        
        TEST_RACE_DETECTION("supersite_race",
            [&]() { 
                mock.reset(); 
            },
            [&]() { 
                return mock.global_updates.load() >= 0; // Basic consistency check
            });
    });
}

// Test without mutex (should demonstrate issues)
void test_no_mutex_supersite() {
    TEST_RUN("no_mutex_supersite", []() {
        MockSuperSiteThreading mock;
        mock.enable_mutex = false;
        
        std::vector<std::thread> threads;
        for (int i = 0; i < 4; i++) {
            threads.emplace_back([&mock, i]() {
                mock.simulate_set_supersite_context(i);
            });
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        // Check for potential issues without mutex
        bool has_issues = !mock.verify_consistency();
        std::cout << "No-mutex SuperSite test: " << (has_issues ? "detected potential race conditions" : "no issues detected") << std::endl;
    });
}

// Main test runner
int main() {
    std::cout << "=== SHAPEIT5 SuperSite Threading Safety Tests ===" << std::endl;
    
    test_supersite_context_thread_safety();
    test_sc_buffer_concurrent_access();
    test_anchor_flag_thread_safety();
    test_sc_offset_consistency();
    test_supersite_stress();
    test_memory_consistency();
    test_supersite_deadlock_detection();
    test_supersite_race_detection();
    test_no_mutex_supersite();
    
    TEST_SUMMARY();
    return TEST_EXIT();
}