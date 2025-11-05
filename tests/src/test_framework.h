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

namespace test_framework {

class TestRunner {
private:
    int total_tests = 0;
    int passed_tests = 0;
    int failed_tests = 0;
    std::chrono::high_resolution_clock::time_point start_time;
    std::string current_test;
    bool current_test_failed = false;

public:
    void start_test(const std::string& test_name) {
        current_test = test_name;
        current_test_failed = false;
        start_time = std::chrono::high_resolution_clock::now();
        total_tests++;
    }

    void assert_condition(bool condition, const std::string& description) {
        if (!condition) {
            current_test_failed = true;
            std::cout << current_test << "::" << description << ": ASSERTION FAILED" << std::endl;
        }
    }

    void finish_test() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1000000.0;

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
};

// Global test runner instance
static TestRunner g_test_runner;

}

// Convenience macros
#define TEST_START(name) test_framework::g_test_runner.start_test(name)
#define TEST_ASSERT(condition, description) test_framework::g_test_runner.assert_condition(condition, description)
#define TEST_PASS() test_framework::g_test_runner.finish_test()
#define TEST_RUN(name, func) test_framework::g_test_runner.run_test(name, func)
#define TEST_EXIT() test_framework::g_test_runner.exit_code()
#define TEST_SUMMARY() test_framework::g_test_runner.print_summary()

// Backward compatibility with assert-style testing
#define SAFE_ASSERT(condition, description) \
    if (!(condition)) { \
        std::cout << "✗ " << description << ": FAIL" << std::endl; \
        return 1; \
    } else { \
        std::cout << "✓ " << description << ": PASS" << std::endl; \
    }

#endif // SHAPEIT5_TEST_FRAMEWORK_H