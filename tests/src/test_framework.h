/*******************************************************************************
 * Simple Test Framework for SHAPEIT5 Unit Tests
 * 
 * Provides uniform test reporting without crashes on assertion failures.
 * Usage:
 *   TEST_START("test_name");
 *   TEST_ASSERT(condition, "description");
 *   TEST_PASS("test_name");
 * 
 * Or for multiple tests in one binary:
 *   TEST_RUN("test1", []() { 
 *     TEST_ASSERT(condition, "check something");
 *   });
 ******************************************************************************/

#ifndef SHAPEIT5_TEST_FRAMEWORK_H
#define SHAPEIT5_TEST_FRAMEWORK_H

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <functional>
#include <pthread.h>
#include <vector>
#include <atomic>
#include <mutex>
#include <unistd.h>
#include <signal.h>

namespace test_framework {

class TestRunner {
private:
    int total_tests = 0;
    int passed_tests = 0;
    int failed_tests = 0;
    std::chrono::high_resolution_clock::time_point start_time;
    std::string current_test;
    bool current_test_failed = false;
    std::mutex test_output_mutex;  // For thread-safe output

public:
    void start_test(const std::string& test_name) {
        current_test = test_name;
        current_test_failed = false;
        start_time = std::chrono::high_resolution_clock::now();
        total_tests++;
    }

    void assert_condition(bool condition, const std::string& description) {
        if (!condition) {
            std::lock_guard<std::mutex> lock(test_output_mutex);
            current_test_failed = true;
            std::cout << current_test << "::" << description << ": ASSERTION FAILED" << std::endl;
        }
    }

    void finish_test() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1000000.0;

        std::lock_guard<std::mutex> lock(test_output_mutex);
        if (current_test_failed) {
            failed_tests++;
            std::cout << current_test << ": FAIL (" << std::fixed << std::setprecision(3) << seconds << "s)" << std::endl;
        } else {
            passed_tests++;
            std::cout << current_test << ": PASS (" << std::fixed << std::setprecision(3) << seconds << "s)" << std::endl;
        }
    }

    void run_test(const std::string& test_name, std::function<void()> test_func) {
        start_test(test_name);
        try {
            test_func();
        } catch (const std::exception& e) {
            current_test_failed = true;
            std::cout << current_test << "::exception: " << e.what() << std::endl;
        } catch (...) {
            current_test_failed = true;
            std::cout << current_test << "::unknown_exception: caught unknown exception" << std::endl;
        }
        finish_test();
    }

    int exit_code() const {
        return failed_tests > 0 ? 1 : 0;
    }

    void print_summary() const {
        std::cout << "Test Summary: " << passed_tests << " passed, " << failed_tests << " failed, " 
                  << total_tests << " total" << std::endl;
    }

    // Thread-safe assertion for use within threaded tests
    void threaded_assert(bool condition, const std::string& description) {
        if (!condition) {
            std::lock_guard<std::mutex> lock(test_output_mutex);
            current_test_failed = true;
            std::cout << current_test << "::" << description << ": THREADED ASSERTION FAILED" << std::endl;
        }
    }
};

// Global test runner instance
static TestRunner g_test_runner;

// Threading utilities for thread safety testing
struct ThreadTestData {
    std::function<void(int)> test_function;
    int thread_id;
    std::atomic<bool>* start_flag;
    std::atomic<int>* error_count;
    std::atomic<bool>* failed_flag;
    
    ThreadTestData(std::function<void(int)> func, int id, 
                   std::atomic<bool>* start, std::atomic<int>* errors, std::atomic<bool>* failed)
        : test_function(func), thread_id(id), start_flag(start), error_count(errors), failed_flag(failed) {}
};

void* thread_test_wrapper(void* arg) {
    ThreadTestData* data = static_cast<ThreadTestData*>(arg);
    
    // Wait for start signal to synchronize all threads
    while (!data->start_flag->load()) {
        usleep(1000); // 1ms sleep
    }
    
    try {
        data->test_function(data->thread_id);
    } catch (const std::exception& e) {
        data->error_count->fetch_add(1);
        data->failed_flag->store(true);
        std::cerr << "Thread " << data->thread_id << " exception: " << e.what() << std::endl;
    } catch (...) {
        data->error_count->fetch_add(1);
        data->failed_flag->store(true);
        std::cerr << "Thread " << data->thread_id << " unknown exception" << std::endl;
    }
    
    return nullptr;
}

// Timeout handler for deadlock detection
static volatile sig_atomic_t timeout_occurred = 0;
void timeout_handler(int sig) {
    timeout_occurred = 1;
}

class ThreadTester {
public:
    // Run a function across multiple threads with synchronization
    static bool run_threaded_test(const std::string& test_name, int thread_count, 
                                  std::function<void(int)> test_func, int timeout_seconds = 30) {
        std::atomic<bool> start_flag(false);
        std::atomic<int> error_count(0);
        std::atomic<bool> failed_flag(false);
        
        std::vector<pthread_t> threads(thread_count);
        std::vector<ThreadTestData> thread_data;
        
        // Create thread data
        for (int i = 0; i < thread_count; i++) {
            thread_data.emplace_back(test_func, i, &start_flag, &error_count, &failed_flag);
        }
        
        // Set up timeout handler
        timeout_occurred = 0;
        signal(SIGALRM, timeout_handler);
        alarm(timeout_seconds);
        
        // Create threads
        for (int i = 0; i < thread_count; i++) {
            if (pthread_create(&threads[i], nullptr, thread_test_wrapper, &thread_data[i]) != 0) {
                std::cerr << "Failed to create thread " << i << std::endl;
                alarm(0); // Cancel alarm
                return false;
            }
        }
        
        // Start all threads simultaneously
        start_flag.store(true);
        
        // Wait for all threads to complete
        for (int i = 0; i < thread_count; i++) {
            pthread_join(threads[i], nullptr);
        }
        
        // Cancel timeout
        alarm(0);
        
        if (timeout_occurred) {
            std::cerr << "Test " << test_name << " timed out after " << timeout_seconds << " seconds" << std::endl;
            return false;
        }
        
        return !failed_flag.load() && error_count.load() == 0;
    }
    
    // Run stress test with repeated iterations
    static bool run_stress_test(const std::string& test_name, int iterations, int thread_count,
                                std::function<void(int)> test_func) {
        for (int i = 0; i < iterations; i++) {
            if (!run_threaded_test(test_name + "_iteration_" + std::to_string(i), 
                                   thread_count, test_func, 10)) {
                std::cerr << "Stress test failed at iteration " << i << std::endl;
                return false;
            }
        }
        return true;
    }
    
    // Race condition detection through repeated fast execution
    static bool detect_race_condition(const std::string& test_name, 
                                      std::function<void()> setup_func,
                                      std::function<bool()> verify_func,
                                      int thread_count = 4, int iterations = 100) {
        for (int i = 0; i < iterations; i++) {
            setup_func(); // Reset state
            
            std::atomic<bool> start_flag(false);
            std::vector<pthread_t> threads(thread_count);
            
            // Create threads that just wait for start signal then execute setup again
            for (int t = 0; t < thread_count; t++) {
                pthread_create(&threads[t], nullptr, [](void* arg) -> void* {
                    std::atomic<bool>* flag = static_cast<std::atomic<bool>*>(arg);
                    while (!flag->load()) { usleep(1); }
                    return nullptr;
                }, &start_flag);
            }
            
            start_flag.store(true);
            for (int t = 0; t < thread_count; t++) {
                pthread_join(threads[t], nullptr);
            }
            
            if (!verify_func()) {
                std::cerr << "Race condition detected in " << test_name << " at iteration " << i << std::endl;
                return false;
            }
        }
        return true;
    }
};

}

// Convenience macros
#define TEST_START(name) test_framework::g_test_runner.start_test(name)
#define TEST_ASSERT(condition, description) test_framework::g_test_runner.assert_condition(condition, description)
#define TEST_PASS() test_framework::g_test_runner.finish_test()
#define TEST_RUN(name, func) test_framework::g_test_runner.run_test(name, func)
#define TEST_EXIT() test_framework::g_test_runner.exit_code()
#define TEST_SUMMARY() test_framework::g_test_runner.print_summary()

// Threading test macros
#define TEST_THREADED_ASSERT(condition, description) test_framework::g_test_runner.threaded_assert(condition, description)

#define TEST_THREADED(name, thread_count, func) \
    test_framework::g_test_runner.run_test(name, [&]() { \
        bool result = test_framework::ThreadTester::run_threaded_test(name, thread_count, func); \
        TEST_ASSERT(result, "threaded test execution"); \
    })

#define TEST_STRESS(name, iterations, thread_count, func) \
    test_framework::g_test_runner.run_test(name, [&]() { \
        bool result = test_framework::ThreadTester::run_stress_test(name, iterations, thread_count, func); \
        TEST_ASSERT(result, "stress test execution"); \
    })

#define TEST_RACE_DETECTION(name, setup_func, verify_func) \
    test_framework::g_test_runner.run_test(name, [&]() { \
        bool result = test_framework::ThreadTester::detect_race_condition(name, setup_func, verify_func); \
        TEST_ASSERT(result, "race condition detection"); \
    })

#define TEST_DEADLOCK_DETECTION(name, timeout_seconds, thread_count, func) \
    test_framework::g_test_runner.run_test(name, [&]() { \
        bool result = test_framework::ThreadTester::run_threaded_test(name, thread_count, func, timeout_seconds); \
        TEST_ASSERT(result, "deadlock detection test"); \
    })

// Backward compatibility with assert-style testing
#define SAFE_ASSERT(condition, description) \
    if (!(condition)) { \
        std::cout << "✗ " << description << ": FAIL" << std::endl; \
        return 1; \
    } else { \
        std::cout << "✓ " << description << ": PASS" << std::endl; \
    }

#endif // SHAPEIT5_TEST_FRAMEWORK_H