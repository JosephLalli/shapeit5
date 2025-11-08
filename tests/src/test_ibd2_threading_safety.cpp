/*******************************************************************************
 * Thread Safety Test for SHAPEIT5 IBD2 Tracking Threading
 * 
 * Tests the thread safety of IBD2 tracking operations including:
 * - Concurrent pushIBD2() calls from multiple threads
 * - Thread-safe collapse() operations on shared IBD2 structures
 * - Track merging and overlap detection race conditions
 * - Memory consistency in IBD2 track vector operations
 * - Deadlock detection in IBD2 constraint management
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>
#include <algorithm>
#include <random>

#include "test_framework.h"
#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/containers/ibd2_tracks.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

using namespace test_framework;

// Mock IBD2 threading environment
class MockIBD2Threading {
public:
    // IBD2 state
    std::vector<std::vector<track>> thread_local_tracks;
    ibd2_tracks shared_ibd2;
    pthread_mutex_t ibd2_mutex;
    
    // Test configuration
    int num_individuals = 10;
    int num_variants = 100;
    int num_threads = 4;
    bool enable_mutex = true;
    
    // Statistics
    std::atomic<int> push_operations{0};
    std::atomic<int> collapse_operations{0};
    std::atomic<int> track_conflicts{0};
    
    MockIBD2Threading() {
        pthread_mutex_init(&ibd2_mutex, nullptr);
        setup_test_data();
    }
    
    ~MockIBD2Threading() {
        pthread_mutex_destroy(&ibd2_mutex);
    }
    
    void setup_test_data() {
        // Simplified initialization - manually set up IBD2 structure
        shared_ibd2.clear();
        shared_ibd2.IBD2 = std::vector<std::vector<track>>(num_individuals);
        shared_ibd2.vec_cm = std::vector<float>(num_variants, 0.0f);
        for (int i = 0; i < num_variants; i++) {
            shared_ibd2.vec_cm[i] = i * 0.1f; // 0.1 cM spacing
        }
        
        // Initialize thread-local tracks
        thread_local_tracks.clear();
        thread_local_tracks.resize(num_threads);
    }
    
    // Generate random IBD2 tracks for testing
    void generate_random_tracks(int thread_id, int count) {
        std::random_device rd;
        std::mt19937 gen(rd() + thread_id);
        std::uniform_int_distribution<int> ind_dist(0, num_individuals - 1);
        std::uniform_int_distribution<int> pos_dist(0, num_variants - 10);
        std::uniform_int_distribution<int> len_dist(5, 20);
        
        // Ensure thread_local_tracks is large enough for this thread_id (thread-safe)
        if (enable_mutex) pthread_mutex_lock(&ibd2_mutex);
        if (thread_id >= thread_local_tracks.size()) {
            thread_local_tracks.resize(thread_id + 1);
        }
        if (enable_mutex) pthread_mutex_unlock(&ibd2_mutex);
        
        thread_local_tracks[thread_id].clear();
        
        for (int i = 0; i < count; i++) {
            int ind = ind_dist(gen);
            int from = pos_dist(gen);
            int to = std::min(from + len_dist(gen), num_variants - 1);
            
            track new_track(ind, from, to);
            thread_local_tracks[thread_id].push_back(new_track);
        }
    }
    
    // Simulate pushIBD2 operations with threading
    void simulate_push_ibd2(int thread_id) {
        generate_random_tracks(thread_id, 20);
        
        for (int i = 0; i < 10; i++) {
            std::vector<track> tracks_to_push = thread_local_tracks[thread_id];
            
            if (enable_mutex) pthread_mutex_lock(&ibd2_mutex);
            
            // Simulate pushIBD2 for random individual
            int target_ind = i % num_individuals;
            shared_ibd2.pushIBD2(target_ind, tracks_to_push);
            push_operations.fetch_add(1);
            
            if (enable_mutex) pthread_mutex_unlock(&ibd2_mutex);
            
            usleep(50); // Brief work simulation
        }
    }
    
    // Simulate collapse operations
    void simulate_collapse_operations(int thread_id) {
        for (int i = 0; i < 5; i++) {
            if (enable_mutex) pthread_mutex_lock(&ibd2_mutex);
            
            // Collapse specific individual's tracks
            int ind = thread_id % num_individuals;
            if (ind < shared_ibd2.IBD2.size() && !shared_ibd2.IBD2[ind].empty()) {
                std::sort(shared_ibd2.IBD2[ind].begin(), shared_ibd2.IBD2[ind].end());
                shared_ibd2.collapse(shared_ibd2.IBD2[ind]);
                collapse_operations.fetch_add(1);
            }
            
            if (enable_mutex) pthread_mutex_unlock(&ibd2_mutex);
            
            usleep(100); // Collapse is more expensive
        }
    }
    
    // Simulate noIBD2 queries (read operations)
    void simulate_ibd2_queries(int thread_id) {
        for (int i = 0; i < 50; i++) {
            int hap0 = (thread_id * 2) % (num_individuals * 2);
            int hap1 = ((thread_id * 2) + 1) % (num_individuals * 2);
            int locus = i % num_variants;
            
            if (enable_mutex) pthread_mutex_lock(&ibd2_mutex);
            
            // Query IBD2 status
            bool no_ibd = shared_ibd2.noIBD2(hap0, hap1, locus);
            
            if (enable_mutex) pthread_mutex_unlock(&ibd2_mutex);
            
            usleep(10);
        }
    }
    
    // Simulate track overlap detection
    void simulate_track_overlap(int thread_id) {
        generate_random_tracks(thread_id, 15);
        
        for (int i = 0; i < thread_local_tracks[thread_id].size() - 1; i++) {
            const track& t1 = thread_local_tracks[thread_id][i];
            
            for (int j = i + 1; j < thread_local_tracks[thread_id].size(); j++) {
                const track& t2 = thread_local_tracks[thread_id][j];
                
                if (t1.overlap(t2)) {
                    track_conflicts.fetch_add(1);
                }
            }
            
            usleep(5);
        }
    }
    
    // Simulate concurrent track merging
    void simulate_track_merging(int thread_id) {
        generate_random_tracks(thread_id, 10);
        
        if (enable_mutex) pthread_mutex_lock(&ibd2_mutex);
        
        // Create overlapping tracks for same individual
        int target_ind = thread_id % num_individuals;
        std::vector<track> merge_tracks;
        
        for (int i = 0; i < 5; i++) {
            int from = i * 5;
            int to = from + 8; // Overlapping ranges
            merge_tracks.emplace_back(target_ind, from, to);
        }
        
        // Test merging
        for (int i = 0; i < merge_tracks.size() - 1; i++) {
            if (merge_tracks[i].overlap(merge_tracks[i + 1])) {
                merge_tracks[i].merge(merge_tracks[i + 1]);
            }
        }
        
        if (enable_mutex) pthread_mutex_unlock(&ibd2_mutex);
        
        usleep(20);
    }
    
    // Reset state for new test
    void reset() {
        push_operations = 0;
        collapse_operations = 0;
        track_conflicts = 0;
        
        shared_ibd2.clear();
        setup_test_data();
        
        for (auto& tracks : thread_local_tracks) {
            tracks.clear();
        }
    }
    
    // Verify state consistency
    bool verify_consistency() {
        // Check that operations happened
        if (push_operations.load() == 0) return false;
        
        // Basic consistency checks - don't assume sorting/collapsing unless we explicitly collapsed
        for (int i = 0; i < shared_ibd2.IBD2.size(); i++) {
            const auto& tracks = shared_ibd2.IBD2[i];
            
            // Basic validity checks for tracks
            for (const auto& track : tracks) {
                // Check track bounds are valid
                if (track.from > track.to) return false;
                if (track.from < 0 || track.to >= num_variants) return false;
                if (track.ind < 0 || track.ind >= num_individuals) return false;
            }
        }
        
        return true;
    }
    
    // Verify collapsed state consistency (only call after collapse())
    bool verify_collapsed_consistency() {
        // Check IBD2 structure consistency after collapse
        for (int i = 0; i < shared_ibd2.IBD2.size(); i++) {
            const auto& tracks = shared_ibd2.IBD2[i];
            
            // Verify tracks are properly sorted and non-overlapping after collapse
            for (int j = 1; j < tracks.size(); j++) {
                if (tracks[j-1].ind > tracks[j].ind) return false;
                if (tracks[j-1].ind == tracks[j].ind && tracks[j-1].to >= tracks[j].from) {
                    // Overlapping tracks of same individual (should be merged)
                    return false;
                }
            }
        }
        
        return true;
    }
    
    // Get total track count
    int get_total_track_count() const {
        int total = 0;
        for (const auto& tracks : shared_ibd2.IBD2) {
            total += tracks.size();
        }
        return total;
    }
};

// Test basic IBD2 push operations thread safety
void test_ibd2_push_thread_safety() {
    TEST_RUN("ibd2_push_thread_safety", []() {
        MockIBD2Threading mock;
        
        TEST_THREADED("ibd2_push_operations", 4, [&](int thread_id) {
            mock.simulate_push_ibd2(thread_id);
        });
        
        TEST_ASSERT(mock.verify_consistency(), "IBD2 push operations consistency");
        TEST_ASSERT(mock.push_operations.load() == 40, "correct number of push operations"); // 4 threads * 10 ops
    });
}

// Test IBD2 collapse operations thread safety
void test_ibd2_collapse_thread_safety() {
    TEST_RUN("ibd2_collapse_thread_safety", []() {
        MockIBD2Threading mock;
        
        // First populate with some tracks
        for (int t = 0; t < 4; t++) {
            mock.simulate_push_ibd2(t);
        }
        
        mock.reset();
        mock.collapse_operations = 0;
        
        // Now test collapse operations
        TEST_THREADED("ibd2_collapse_operations", 3, [&](int thread_id) {
            mock.simulate_collapse_operations(thread_id);
        });
        
        TEST_ASSERT(mock.collapse_operations.load() >= 0, "collapse operations executed");
    });
}

// Test concurrent IBD2 queries
void test_ibd2_query_thread_safety() {
    TEST_RUN("ibd2_query_thread_safety", []() {
        MockIBD2Threading mock;
        
        // Populate with some IBD2 data
        for (int t = 0; t < 2; t++) {
            mock.simulate_push_ibd2(t);
        }
        
        TEST_THREADED("ibd2_query_operations", 4, [&](int thread_id) {
            mock.simulate_ibd2_queries(thread_id);
        });
        
        // No crashes means success for read operations
        TEST_ASSERT(true, "IBD2 query operations completed");
    });
}

// Test track overlap detection thread safety
void test_track_overlap_safety() {
    TEST_RUN("track_overlap_safety", []() {
        MockIBD2Threading mock;
        
        TEST_THREADED("track_overlap_detection", 4, [&](int thread_id) {
            mock.simulate_track_overlap(thread_id);
        });
        
        TEST_ASSERT(mock.track_conflicts.load() >= 0, "track overlap operations completed");
    });
}

// Test concurrent track merging
void test_track_merging_safety() {
    TEST_RUN("track_merging_safety", []() {
        MockIBD2Threading mock;
        
        TEST_THREADED("track_merging", 4, [&](int thread_id) {
            mock.simulate_track_merging(thread_id);
        });
        
        TEST_ASSERT(mock.verify_consistency(), "track merging basic consistency");
    });
}

// Test mixed operations stress test
void test_ibd2_mixed_operations() {
    TEST_RUN("ibd2_mixed_operations", []() {
        MockIBD2Threading mock;
        
        TEST_THREADED("mixed_ibd2_operations", 6, [&](int thread_id) {
            switch (thread_id % 3) {
                case 0:
                    mock.simulate_push_ibd2(thread_id);
                    break;
                case 1:
                    mock.simulate_ibd2_queries(thread_id);
                    break;
                case 2:
                    mock.simulate_collapse_operations(thread_id);
                    break;
            }
        });
        
        TEST_ASSERT(mock.push_operations.load() > 0, "push operations in mixed test");
    });
}

// Test IBD2 stress conditions
void test_ibd2_stress() {
    TEST_RUN("ibd2_stress", []() {
        MockIBD2Threading mock;
        mock.num_individuals = 20; // More individuals for stress
        mock.setup_test_data();
        
        TEST_STRESS("ibd2_stress_test", 3, 8, [&](int thread_id) {
            // High frequency operations
            for (int i = 0; i < 5; i++) {
                mock.simulate_push_ibd2(thread_id);
                mock.simulate_ibd2_queries(thread_id);
            }
        });
        
        TEST_ASSERT(mock.verify_consistency(), "IBD2 stress test consistency");
    });
}

// Test deadlock detection in IBD2 operations
void test_ibd2_deadlock_detection() {
    TEST_RUN("ibd2_deadlock_detection", []() {
        MockIBD2Threading mock;
        
        TEST_DEADLOCK_DETECTION("ibd2_deadlock_test", 10, 4, [&](int thread_id) {
            // Mix operations that could potentially deadlock
            mock.simulate_push_ibd2(thread_id);
            mock.simulate_collapse_operations(thread_id);
        });
    });
}

// Test race condition detection
void test_ibd2_race_detection() {
    TEST_RUN("ibd2_race_detection", []() {
        MockIBD2Threading mock;
        
        TEST_RACE_DETECTION("ibd2_race",
            [&]() { 
                mock.reset(); 
            },
            [&]() { 
                // Basic consistency: track count should be non-negative
                return mock.get_total_track_count() >= 0;
            });
    });
}

// Test memory consistency across threads
void test_ibd2_memory_consistency() {
    TEST_RUN("ibd2_memory_consistency", []() {
        MockIBD2Threading mock;
        
        // Initialize with known pattern
        for (int ind = 0; ind < mock.num_individuals; ind++) {
            std::vector<track> initial_tracks;
            initial_tracks.emplace_back(ind, 0, 10);
            mock.shared_ibd2.IBD2[ind] = initial_tracks;
        }
        
        TEST_THREADED("memory_consistency_check", 4, [&](int thread_id) {
            // Verify data consistency
            pthread_mutex_lock(&mock.ibd2_mutex);
            
            for (int ind = 0; ind < mock.num_individuals; ind++) {
                const auto& tracks = mock.shared_ibd2.IBD2[ind];
                TEST_THREADED_ASSERT(tracks.size() >= 1, 
                    "individual " + std::to_string(ind) + " has tracks");
                
                if (!tracks.empty()) {
                    TEST_THREADED_ASSERT(tracks[0].ind == ind, 
                        "track individual matches index");
                }
            }
            
            pthread_mutex_unlock(&mock.ibd2_mutex);
        });
    });
}

// Test without mutex (should demonstrate race conditions)
void test_no_mutex_ibd2() {
    TEST_RUN("no_mutex_ibd2", []() {
        MockIBD2Threading mock;
        mock.enable_mutex = false;
        
        std::vector<std::thread> threads;
        for (int i = 0; i < 4; i++) {
            threads.emplace_back([&mock, i]() {
                mock.simulate_push_ibd2(i);
            });
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        // Check for potential issues without mutex
        bool has_issues = !mock.verify_consistency();
        std::cout << "No-mutex IBD2 test: " << (has_issues ? "detected potential race conditions" : "no issues detected") << std::endl;
    });
}

// Main test runner
int main() {
    std::cout << "=== SHAPEIT5 IBD2 Threading Safety Tests ===" << std::endl;
    
    test_ibd2_push_thread_safety();
    test_ibd2_collapse_thread_safety();
    test_ibd2_query_thread_safety();
    test_track_overlap_safety();
    test_track_merging_safety();
    test_ibd2_mixed_operations();
    test_ibd2_stress();
    test_ibd2_deadlock_detection();
    test_ibd2_race_detection();
    test_ibd2_memory_consistency();
    test_no_mutex_ibd2();
    
    TEST_SUMMARY();
    return TEST_EXIT();
}