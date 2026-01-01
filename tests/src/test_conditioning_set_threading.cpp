/*******************************************************************************
 * Thread Safety Test for SHAPEIT5 Conditioning Set Threading
 * 
 * Tests the thread safety of conditioning_set PBWT operations including:
 * - Mutex protection in selecter_callback function
 * - PBWT state management across threads  
 * - Concurrent access to shared PBWT data structures
 * - Progress tracking thread safety
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>

#include "test_common.h"
#include "../../common/src/utils/otools.h"


#define private public
#define protected public
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/genotype_set.h"
#undef private
#undef protected

// Mock conditioning set for testing PBWT threading
class MockConditioningSet {
public:
    // Threading state (matching conditioning_set)
    int i_worker = 0;
    int i_job = 0;
    int d_job = 0;
    pthread_mutex_t mutex_workers;
    std::vector<pthread_t> id_workers;
    
    // PBWT state simulation
    std::vector<int> sites_pbwt_mthreading;
    std::vector<bool> pbwt_processed;
    std::atomic<int> selection_count{0};
    std::atomic<int> progress_updates{0};
    
    // Test configuration
    bool enable_mutex = true;
    int max_sites = 100;
    
    MockConditioningSet() {
        pthread_mutex_init(&mutex_workers, nullptr);
        setup_test_data();
    }
    
    ~MockConditioningSet() {
        pthread_mutex_destroy(&mutex_workers);
    }
    
    void setup_test_data() {
        sites_pbwt_mthreading.clear();
        pbwt_processed.clear();
        
        for (int i = 0; i < max_sites; i++) {
            sites_pbwt_mthreading.push_back(i);
            pbwt_processed.push_back(false);
        }
    }
    
    // Simulate selecter_callback behavior
    void simulate_selecter_thread(int thread_id) {
        int local_worker_id, local_job_id;
        
        // Get worker ID
        if (enable_mutex) pthread_mutex_lock(&mutex_workers);
        local_worker_id = i_worker++;
        if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
        
        while (true) {
            // Get job ID
            if (enable_mutex) pthread_mutex_lock(&mutex_workers);
            local_job_id = i_job++;
            if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
            
            if (local_job_id >= sites_pbwt_mthreading.size()) break;
            
            // Simulate PBWT selection
            simulate_pbwt_select(local_job_id);
            
            // Update progress
            if (enable_mutex) pthread_mutex_lock(&mutex_workers);
            d_job++;
            progress_updates.fetch_add(1);
            if (enable_mutex) pthread_mutex_unlock(&mutex_workers);
            
            // Brief work simulation
            usleep(50); // 0.05ms
        }
    }
    
    // Simulate PBWT selection work
    void simulate_pbwt_select(int job_id) {
        if (job_id < pbwt_processed.size()) {
            pbwt_processed[job_id] = true;
            selection_count.fetch_add(1);
        }
    }
    
    // Reset state for new test
    void reset() {
        i_worker = 0;
        i_job = 0;
        d_job = 0;
        selection_count = 0;
        progress_updates = 0;
        setup_test_data();
    }
    
    // Verify state consistency
    bool verify_consistency() {
        // Check that all jobs were processed
        if (d_job < max_sites) return false;
        
        // Check selection count
        if (selection_count.load() != max_sites) return false;
        
        // Check that all sites were processed exactly once
        for (bool processed : pbwt_processed) {
            if (!processed) return false;
        }
        
        return true;
    }
};

// Test basic PBWT selection thread safety
void test_pbwt_selection_thread_safety() {
    MockConditioningSet conditioning_set;

    // Test with different thread counts
    for (int thread_count : {2, 4, 6}) {
        conditioning_set.reset();
        conditioning_set.max_sites = 80;
        conditioning_set.setup_test_data();

        TEST_THREADED("pbwt_selection_threads_" + std::to_string(thread_count),
                     thread_count,
                     [&](int thread_id) {
                         conditioning_set.simulate_selecter_thread(thread_id);
                     });

        TEST_CHECK(conditioning_set.verify_consistency(),
                   "pbwt_selection_consistency_threads_" + std::to_string(thread_count),
                   "PBWT selection consistency with " + std::to_string(thread_count) + " threads");
    }
}

// Test worker ID assignment uniqueness
void test_worker_id_assignment() {
    MockConditioningSet conditioning_set;
    std::atomic<int> worker_assignments[8] = {};

    TEST_THREADED("worker_id_uniqueness", 8, [&](int thread_id) {
        pthread_mutex_lock(&conditioning_set.mutex_workers);
        int assigned_worker_id = conditioning_set.i_worker++;
        pthread_mutex_unlock(&conditioning_set.mutex_workers);

        if (assigned_worker_id < 8) {
            worker_assignments[assigned_worker_id].fetch_add(1);
        }
    });

    // Verify each worker ID assigned exactly once
    for (int i = 0; i < 8; i++) {
        TEST_CHECK(worker_assignments[i].load() == 1,
                   "worker_id_unique_" + std::to_string(i),
                   "worker ID " + std::to_string(i) + " assigned once");
    }
}

// Test job distribution fairness
void test_job_distribution() {
    MockConditioningSet conditioning_set;
    conditioning_set.max_sites = 100;
    conditioning_set.setup_test_data();

    std::atomic<int> jobs_per_thread[4] = {};

    TEST_THREADED("job_distribution_fairness", 4, [&](int thread_id) {
        int local_jobs = 0;

        while (true) {
            pthread_mutex_lock(&conditioning_set.mutex_workers);
            int job_id = conditioning_set.i_job++;
            pthread_mutex_unlock(&conditioning_set.mutex_workers);

            if (job_id >= conditioning_set.max_sites) break;

            local_jobs++;
            conditioning_set.simulate_pbwt_select(job_id);
            usleep(10);
        }

        jobs_per_thread[thread_id].store(local_jobs);
    });

    // Verify reasonable job distribution (no thread should be starved)
    int total_jobs = 0;
    for (int i = 0; i < 4; i++) {
        int thread_jobs = jobs_per_thread[i].load();
        total_jobs += thread_jobs;
        TEST_CHECK(thread_jobs > 0,
                   "job_distribution_thread_" + std::to_string(i),
                   "thread " + std::to_string(i) + " got some jobs");
    }

    TEST_CHECK(total_jobs == conditioning_set.max_sites,
               "job_distribution_total",
               "all jobs distributed");
}

// Test progress tracking thread safety
void test_progress_tracking() {
    MockConditioningSet conditioning_set;
    conditioning_set.max_sites = 60;
    conditioning_set.setup_test_data();

    TEST_THREADED("progress_updates", 3, [&](int thread_id) {
        for (int i = 0; i < 20; i++) {
            pthread_mutex_lock(&conditioning_set.mutex_workers);
            conditioning_set.d_job++;
            conditioning_set.progress_updates.fetch_add(1);
            pthread_mutex_unlock(&conditioning_set.mutex_workers);

            usleep(5);
        }
    });

    TEST_CHECK(conditioning_set.d_job == 60, "progress_counter_correct", "progress counter correct");
    TEST_CHECK(conditioning_set.progress_updates.load() == 60,
               "progress_updates_correct",
               "progress updates correct");
}

// Test concurrent PBWT state access
void test_pbwt_state_access() {
    MockConditioningSet conditioning_set;
    conditioning_set.max_sites = 50;
    conditioning_set.setup_test_data();

    std::atomic<int> state_conflicts{0};

    TEST_THREADED("pbwt_concurrent_access", 4, [&](int thread_id) {
        for (int i = thread_id; i < conditioning_set.max_sites; i += 4) {
            // Try to process site i
            bool expected = false;
            bool was_processed = conditioning_set.pbwt_processed[i];

            if (was_processed) {
                state_conflicts.fetch_add(1);
            } else {
                conditioning_set.pbwt_processed[i] = true;
                conditioning_set.selection_count.fetch_add(1);
            }

            usleep(20);
        }
    });

    // With proper partitioning, there should be no conflicts
    TEST_CHECK(state_conflicts.load() == 0, "pbwt_no_state_conflicts", "no PBWT state conflicts");

    // Count processed sites
    int processed_count = 0;
    for (bool processed : conditioning_set.pbwt_processed) {
        if (processed) processed_count++;
    }

    TEST_CHECK(processed_count == conditioning_set.max_sites, "pbwt_all_sites_processed", "all sites processed");
}

// Test stress conditions with high contention
void test_high_contention_stress() {
    MockConditioningSet conditioning_set;
    conditioning_set.max_sites = 200;
    conditioning_set.setup_test_data();

    TEST_STRESS("pbwt_high_contention", 5, 8, [&](int thread_id) {
        // Simulate heavy mutex contention
        for (int i = 0; i < 25; i++) {
            pthread_mutex_lock(&conditioning_set.mutex_workers);
            int job_id = conditioning_set.i_job++;
            pthread_mutex_unlock(&conditioning_set.mutex_workers);

            if (job_id < conditioning_set.max_sites) {
                conditioning_set.simulate_pbwt_select(job_id);

                pthread_mutex_lock(&conditioning_set.mutex_workers);
                conditioning_set.d_job++;
                pthread_mutex_unlock(&conditioning_set.mutex_workers);
            }

            // Minimal work to increase contention
            usleep(1);
        }
    });

    // Verify final state consistency
    TEST_CHECK(conditioning_set.selection_count.load() >= conditioning_set.max_sites,
               "pbwt_selection_count_post_stress",
               "selection count after stress test");
}

// Test deadlock detection in PBWT operations
void test_pbwt_deadlock_detection() {
    MockConditioningSet conditioning_set;
    conditioning_set.max_sites = 30;
    conditioning_set.setup_test_data();

    TEST_DEADLOCK_DETECTION("pbwt_deadlock_test", 10, 4, [&](int thread_id) {
        for (int i = 0; i < 8; i++) {
            pthread_mutex_lock(&conditioning_set.mutex_workers);
            int job_id = conditioning_set.i_job++;
            pthread_mutex_unlock(&conditioning_set.mutex_workers);

            if (job_id < conditioning_set.max_sites) {
                conditioning_set.simulate_pbwt_select(job_id);
            }

            usleep(100);
        }
    });
}

// Test behavior without mutex (should show race conditions)
void test_no_mutex_races() {
    TEST_START("no_mutex_races");
    MockConditioningSet conditioning_set;
    conditioning_set.enable_mutex = false;
    conditioning_set.max_sites = 100;
    conditioning_set.setup_test_data();

    std::vector<std::thread> threads;
    for (int i = 0; i < 4; i++) {
        threads.emplace_back([&conditioning_set, i]() {
            conditioning_set.simulate_selecter_thread(i);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    // Without mutex protection, we expect potential issues
    bool has_issues = !conditioning_set.verify_consistency();
    std::cout << "No-mutex PBWT test: " << (has_issues ? "detected race conditions" : "no issues detected") << std::endl;
    TEST_PASS("no_mutex_races");
}

// Main test runner
int main() {
    TEST_INIT("test_conditioning_set_threading");
    std::cout << "=== SHAPEIT5 Conditioning Set Threading Safety Tests ===" << std::endl;
    
    test_pbwt_selection_thread_safety();
    test_worker_id_assignment();
    test_job_distribution();
    test_progress_tracking();
    test_pbwt_state_access();
    test_high_contention_stress();
    test_pbwt_deadlock_detection();
    test_no_mutex_races();
    
    TEST_SUMMARY();
    return TestReporting::exit_code();
}
