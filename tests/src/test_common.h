#ifndef SHAPEIT5_TEST_COMMON_H
#define SHAPEIT5_TEST_COMMON_H

#include "test_reporting.h"

#include <atomic>
#include <exception>
#include <functional>
#include <mutex>
#include <pthread.h>
#include <signal.h>
#include <string>
#include <unistd.h>
#include <vector>

namespace TestThreading {

struct ThreadTestData {
    std::function<void(int)> test_function;
    int thread_id;
    std::atomic<bool>* start_flag;
    std::atomic<int>* error_count;
    std::atomic<bool>* failed_flag;
    std::mutex* output_mutex;

    ThreadTestData(std::function<void(int)> func, int id,
                   std::atomic<bool>* start, std::atomic<int>* errors,
                   std::atomic<bool>* failed, std::mutex* mutex)
        : test_function(func),
          thread_id(id),
          start_flag(start),
          error_count(errors),
          failed_flag(failed),
          output_mutex(mutex) {}
};

struct ThreadAssertionContext {
    std::atomic<bool>* failed_flag;
    std::mutex* output_mutex;
    std::string test_name;
};

static ThreadAssertionContext* current_ctx = nullptr;

inline void set_context(ThreadAssertionContext* ctx) {
    current_ctx = ctx;
}

inline void threaded_assert(bool condition, const std::string& description) {
    if (condition) return;
    if (current_ctx && current_ctx->failed_flag) {
        current_ctx->failed_flag->store(true);
        if (current_ctx->output_mutex) {
            std::lock_guard<std::mutex> lock(*current_ctx->output_mutex);
            std::cerr << current_ctx->test_name << "::" << description
                      << ": THREADED ASSERTION FAILED" << std::endl;
        } else {
            std::cerr << current_ctx->test_name << "::" << description
                      << ": THREADED ASSERTION FAILED" << std::endl;
        }
    } else {
        std::cerr << "threaded_assert: " << description << ": FAILED" << std::endl;
    }
}

void* thread_test_wrapper(void* arg) {
    ThreadTestData* data = static_cast<ThreadTestData*>(arg);

    while (!data->start_flag->load()) {
        usleep(1000);
    }

    try {
        data->test_function(data->thread_id);
    } catch (const std::exception& e) {
        data->error_count->fetch_add(1);
        data->failed_flag->store(true);
        if (data->output_mutex) {
            std::lock_guard<std::mutex> lock(*data->output_mutex);
            std::cerr << "Thread " << data->thread_id << " exception: " << e.what() << std::endl;
        } else {
            std::cerr << "Thread " << data->thread_id << " exception: " << e.what() << std::endl;
        }
    } catch (...) {
        data->error_count->fetch_add(1);
        data->failed_flag->store(true);
        if (data->output_mutex) {
            std::lock_guard<std::mutex> lock(*data->output_mutex);
            std::cerr << "Thread " << data->thread_id << " unknown exception" << std::endl;
        } else {
            std::cerr << "Thread " << data->thread_id << " unknown exception" << std::endl;
        }
    }

    return nullptr;
}

static volatile sig_atomic_t timeout_occurred = 0;
inline void timeout_handler(int) {
    timeout_occurred = 1;
}

class ThreadTester {
public:
    static bool run_threaded_test_with_state(const std::string& test_name,
                                             int thread_count,
                                             std::function<void(int)> test_func,
                                             std::atomic<bool>& failed_flag,
                                             std::atomic<int>& error_count,
                                             std::mutex& output_mutex,
                                             int timeout_seconds = 30) {
        std::atomic<bool> start_flag(false);
        std::vector<pthread_t> threads(thread_count);
        std::vector<ThreadTestData> thread_data;
        thread_data.reserve(thread_count);

        for (int i = 0; i < thread_count; i++) {
            thread_data.emplace_back(test_func, i, &start_flag, &error_count, &failed_flag, &output_mutex);
        }

        timeout_occurred = 0;
        signal(SIGALRM, timeout_handler);
        alarm(timeout_seconds);

        for (int i = 0; i < thread_count; i++) {
            if (pthread_create(&threads[i], nullptr, thread_test_wrapper, &thread_data[i]) != 0) {
                alarm(0);
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "Failed to create thread " << i << " in " << test_name << std::endl;
                return false;
            }
        }

        start_flag.store(true);

        for (int i = 0; i < thread_count; i++) {
            pthread_join(threads[i], nullptr);
        }

        alarm(0);

        if (timeout_occurred) {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cerr << "Test " << test_name << " timed out after " << timeout_seconds << " seconds" << std::endl;
            return false;
        }

        return !failed_flag.load() && error_count.load() == 0;
    }

    static bool run_stress_test_with_state(const std::string& test_name,
                                           int iterations,
                                           int thread_count,
                                           std::function<void(int)> test_func,
                                           std::atomic<bool>& failed_flag,
                                           std::atomic<int>& error_count,
                                           std::mutex& output_mutex) {
        for (int i = 0; i < iterations; i++) {
            failed_flag.store(false);
            error_count.store(0);
            if (!run_threaded_test_with_state(test_name + "_iteration_" + std::to_string(i),
                                              thread_count, test_func, failed_flag,
                                              error_count, output_mutex, 10)) {
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "Stress test failed at iteration " << i << std::endl;
                return false;
            }
        }
        return true;
    }

    static bool detect_race_condition(const std::string& test_name,
                                      std::function<void()> setup_func,
                                      std::function<bool()> verify_func,
                                      int thread_count = 4,
                                      int iterations = 100) {
        for (int i = 0; i < iterations; i++) {
            setup_func();

            std::atomic<bool> start_flag(false);
            std::vector<pthread_t> threads(thread_count);

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
                std::cerr << "Race condition detected in " << test_name
                          << " at iteration " << i << std::endl;
                return false;
            }
        }
        return true;
    }
};

} // namespace TestThreading

#define TEST_THREADED_ASSERT(condition, description) \
    TestThreading::threaded_assert(condition, description)

#define TEST_THREADED(name, thread_count, func) \
    do { \
        TEST_START(name); \
        std::atomic<bool> _failed(false); \
        std::atomic<int> _errors(0); \
        std::mutex _output_mutex; \
        TestThreading::ThreadAssertionContext _ctx{&_failed, &_output_mutex, name}; \
        TestThreading::set_context(&_ctx); \
        bool _result = TestThreading::ThreadTester::run_threaded_test_with_state( \
            name, thread_count, func, _failed, _errors, _output_mutex); \
        TestThreading::set_context(nullptr); \
        if (_result) TEST_PASS(name); else TEST_FAIL(name, "threaded test execution"); \
    } while (0)

#define TEST_STRESS(name, iterations, thread_count, func) \
    do { \
        TEST_START(name); \
        std::atomic<bool> _failed(false); \
        std::atomic<int> _errors(0); \
        std::mutex _output_mutex; \
        TestThreading::ThreadAssertionContext _ctx{&_failed, &_output_mutex, name}; \
        TestThreading::set_context(&_ctx); \
        bool _result = TestThreading::ThreadTester::run_stress_test_with_state( \
            name, iterations, thread_count, func, _failed, _errors, _output_mutex); \
        TestThreading::set_context(nullptr); \
        if (_result) TEST_PASS(name); else TEST_FAIL(name, "stress test execution"); \
    } while (0)

#define TEST_RACE_DETECTION(name, setup_func, verify_func) \
    do { \
        TEST_START(name); \
        bool _result = TestThreading::ThreadTester::detect_race_condition( \
            name, setup_func, verify_func); \
        if (_result) TEST_PASS(name); else TEST_FAIL(name, "race condition detection"); \
    } while (0)

#define TEST_DEADLOCK_DETECTION(name, timeout_seconds, thread_count, func) \
    do { \
        TEST_START(name); \
        std::atomic<bool> _failed(false); \
        std::atomic<int> _errors(0); \
        std::mutex _output_mutex; \
        TestThreading::ThreadAssertionContext _ctx{&_failed, &_output_mutex, name}; \
        TestThreading::set_context(&_ctx); \
        bool _result = TestThreading::ThreadTester::run_threaded_test_with_state( \
            name, thread_count, func, _failed, _errors, _output_mutex, timeout_seconds); \
        TestThreading::set_context(nullptr); \
        if (_result) TEST_PASS(name); else TEST_FAIL(name, "deadlock detection test"); \
    } while (0)

#endif // SHAPEIT5_TEST_COMMON_H
