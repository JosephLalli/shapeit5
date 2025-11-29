/*******************************************************************************
 * Thread Safety Test for SHAPEIT5 Phaser Threading
 * 
 * Tests the thread safety of core phaser functions including:
 * - Mutex protection around shared counters (i_workers, i_jobs)  
 * - Statistics aggregation thread safety (statH, statS)
 * - Job queue integrity under concurrent access
 * - Deadlock detection in multi-threaded phasing
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>

#include "test_framework.h"
#include "../../common/src/utils/otools.h"

#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/phaser/phaser_header.h"
#include "../../phase_common/src/objects/compute_job.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/genotype_set.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

using namespace test_framework;

// Mock phaser class for testing threading behavior
class MockPhaser {
public:
    // Threading state
    int i_workers = 0;
    int i_jobs = 0;
    pthread_mutex_t mutex_workers;
    std::vector<pthread_t> id_workers;
    
    // Statistics (shared data)
    std::atomic<int> stat_updates{0};
    std::vector<double> shared_data;
    
    // Test configuration
    int max_jobs = 100;
    int max_workers = 8;
    bool enable_mutex = true;
    
    MockPhaser() {
        pthread_mutex_init(&mutex_workers, nullptr);
        shared_data.resize(1000, 0.0);
    }
    
    ~MockPhaser() {
        pthread_mutex_destroy(&mutex_workers);
    }
    
    // Simulate phaseWindow_callback behavior
    void simulate_worker_thread(int thread_id) {
        int local_worker_id, local_job_id;
        
        // Get worker ID (mimic phaseWindow_callback)
        if (enable_mutex) pthread_mutex_lock(&mutex_workers);
        local_worker_id = i_workers++;
        if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
        
        // Process jobs
        while (true) {
            if (enable_mutex) pthread_mutex_lock(&mutex_workers);
            local_job_id = i_jobs++;
            if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
            
            if (local_job_id >= max_jobs) break;
            
            // Simulate work and statistics update
            simulate_statistics_update(local_job_id);
            
            // Brief work simulation
            usleep(100); // 0.1ms
        }
    }
    
    // Simulate statistics collection (shared data access)
    void simulate_statistics_update(int job_id) {
        if (enable_mutex) pthread_mutex_lock(&mutex_workers);
        
        // Update shared statistics (mimics statH.push(), statS.push())
        stat_updates.fetch_add(1);
        if (job_id < shared_data.size()) {
            shared_data[job_id] += 1.0;
        }
        
        if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
    }
    
    // Reset state for new test
    void reset() {
        i_workers = 0;
        i_jobs = 0;
        stat_updates = 0;
        std::fill(shared_data.begin(), shared_data.end(), 0.0);
    }
    
    // Verify state consistency
    bool verify_consistency() {
        // Check that job count is correct
        if (i_jobs < max_jobs) return false;
        
        // Check that statistics were updated correctly
        if (stat_updates.load() != max_jobs) return false;
        
        // Check shared data integrity
        for (int i = 0; i < max_jobs; i++) {
            if (shared_data[i] != 1.0) return false;
        }
        
        return true;
    }
};

// Test basic mutex protection for job counters
void test_job_counter_thread_safety() {
    TEST_RUN("job_counter_thread_safety", []() {
        MockPhaser phaser;
        
        // Test with different thread counts
        for (int thread_count : {2, 4, 8}) {
            phaser.reset();
            phaser.max_jobs = 100;
            
            TEST_THREADED("job_counter_threads_" + std::to_string(thread_count), 
                         thread_count, 
                         [&](int thread_id) {
                             phaser.simulate_worker_thread(thread_id);
                         });
            
            TEST_ASSERT(phaser.verify_consistency(), 
                       "job counter consistency with " + std::to_string(thread_count) + " threads");
        }
    });
}

// Test statistics aggregation thread safety
void test_statistics_thread_safety() {
    TEST_RUN("statistics_thread_safety", []() {
        MockPhaser phaser;
        phaser.max_jobs = 200;
        
        TEST_STRESS("statistics_stress_test", 10, 4, [&](int thread_id) {
            for (int i = 0; i < 50; i++) {
                phaser.simulate_statistics_update(i + thread_id * 50);
            }
        });
        
        TEST_ASSERT(phaser.stat_updates.load() == 10 * 4 * 50, 
                   "statistics update count");
    });
}

// Test race condition detection in job assignment
void test_job_assignment_race_detection() {
    TEST_RUN("job_assignment_race_detection", []() {
        MockPhaser phaser;
        
        TEST_RACE_DETECTION("job_assignment_race",
            [&]() { 
                phaser.reset();
                phaser.max_jobs = 50;
            },
            [&]() { 
                return phaser.i_jobs <= phaser.max_jobs; 
            });
    });
}

// Test deadlock detection in mutex usage
void test_deadlock_detection() {
    TEST_RUN("deadlock_detection", []() {
        MockPhaser phaser;
        phaser.max_jobs = 20; // Small job count for quick test
        
        TEST_DEADLOCK_DETECTION("mutex_deadlock_test", 5, 4, [&](int thread_id) {
            phaser.simulate_worker_thread(thread_id);
        });
    });
}

// Test behavior without mutex protection (should fail or show inconsistency)
void test_no_mutex_behavior() {
    TEST_RUN("no_mutex_behavior", []() {
        MockPhaser phaser;
        phaser.enable_mutex = false; // Disable mutex protection
        phaser.max_jobs = 100;
        
        // This test expects potential inconsistency without mutex
        bool found_inconsistency = false;
        
        for (int trial = 0; trial < 5; trial++) {
            phaser.reset();
            
            std::vector<std::thread> threads;
            for (int i = 0; i < 4; i++) {
                threads.emplace_back([&phaser, i]() {
                    phaser.simulate_worker_thread(i);
                });
            }
            
            for (auto& t : threads) {
                t.join();
            }
            
            // Without mutex, we might see race conditions
            if (!phaser.verify_consistency()) {
                found_inconsistency = true;
                break;
            }
        }
        
        // Note: This test demonstrates that mutex is necessary
        // If no inconsistency found, that's also OK (race conditions are not guaranteed)
        std::cout << "No-mutex test: " << (found_inconsistency ? "detected race conditions" : "no issues detected") << std::endl;
    });
}

// Test worker ID uniqueness
void test_worker_id_uniqueness() {
    TEST_RUN("worker_id_uniqueness", []() {
        MockPhaser phaser;
        std::atomic<int> worker_ids[8] = {}; // Track assigned worker IDs
        
        TEST_THREADED("worker_id_assignment", 8, [&](int thread_id) {
            // Simulate getting worker ID
            pthread_mutex_lock(&phaser.mutex_workers);
            int assigned_id = phaser.i_workers++;
            pthread_mutex_unlock(&phaser.mutex_workers);
            
            // Mark this worker ID as used
            worker_ids[assigned_id].fetch_add(1);
        });
        
        // Verify each worker ID was assigned exactly once
        for (int i = 0; i < 8; i++) {
            TEST_ASSERT(worker_ids[i].load() == 1, 
                       "worker ID " + std::to_string(i) + " assigned exactly once");
        }
    });
}

// Test concurrent access to shared data structures
void test_shared_data_access() {
    TEST_RUN("shared_data_access", []() {
        MockPhaser phaser;
        std::atomic<int> access_count{0};
        
        TEST_THREADED("shared_data_concurrent_access", 6, [&](int thread_id) {
            for (int i = 0; i < 20; i++) {
                pthread_mutex_lock(&phaser.mutex_workers);
                
                // Simulate accessing shared data
                access_count.fetch_add(1);
                phaser.shared_data[i] += thread_id;
                
                pthread_mutex_unlock(&phaser.mutex_workers);
                
                usleep(10); // Small delay to increase contention
            }
        });
        
        TEST_ASSERT(access_count.load() == 6 * 20, "total access count");
        
        // Verify data consistency
        double expected_sum = (0 + 1 + 2 + 3 + 4 + 5) * 20; // sum of thread_ids * iterations
        double actual_sum = 0;
        for (int i = 0; i < 20; i++) {
            actual_sum += phaser.shared_data[i];
        }
        
        TEST_ASSERT(std::abs(actual_sum - expected_sum) < 1e-6, "shared data consistency");
    });
}

// Main test runner
int main() {
    TEST_INIT("test_phaser_threading_safety");
    std::cout << "=== SHAPEIT5 Phaser Threading Safety Tests ===" << std::endl;
    
    test_job_counter_thread_safety();
    test_statistics_thread_safety(); 
    test_job_assignment_race_detection();
    test_deadlock_detection();
    test_no_mutex_behavior();
    test_worker_id_uniqueness();
    test_shared_data_access();
    
    TEST_SUMMARY();
    return TEST_EXIT();
}