/*******************************************************************************
 * Standardized Test Reporting Utilities
 *
 * Provides consistent test output format and crash diagnostics across all tests.
 *
 * Usage:
 *   #include "test_reporting.h"
 *
 *   int main() {
 *       TEST_INIT("my_test_name");  // Sets up signal handlers
 *
 *       TEST_START("scenario_1", "Testing basic functionality");
 *       // ... test code ...
 *       TEST_PASS("scenario_1");
 *
 *       TEST_START("scenario_2", "Testing edge cases");
 *       // ... test code ...
 *       if (condition) {
 *           TEST_FAIL("scenario_2", "Expected X but got Y");
 *           return 1;
 *       }
 *       TEST_PASS("scenario_2");
 *
 *       TEST_SUMMARY();
 *       return 0;
 *   }
 ******************************************************************************/

#ifndef TEST_REPORTING_H
#define TEST_REPORTING_H

#include <iostream>
#include <string>
#include <vector>
#include <ctime>
#include <csignal>
#include <cstdlib>
#include <sstream>
#include <iomanip>

namespace TestReporting {

// Global state for tracking test context and results
struct TestState {
    std::string binary_name;
    std::string current_test_context;
    std::vector<std::string> passed_tests;
    std::vector<std::pair<std::string, std::string>> failed_tests;  // name, reason
    time_t start_time;
    bool initialized = false;
} static test_state;

// Signal handler for crash diagnostics
static void crash_handler(int signal) {
    const char* signal_name = "UNKNOWN";
    switch (signal) {
        case SIGSEGV: signal_name = "SIGSEGV (segmentation fault)"; break;
        case SIGABRT: signal_name = "SIGABRT (abort)"; break;
        case SIGFPE:  signal_name = "SIGFPE (floating point exception)"; break;
        case SIGILL:  signal_name = "SIGILL (illegal instruction)"; break;
    }

    std::cerr << "\n[TEST_CRASH] Signal: " << signal_name << std::endl;
    std::cerr << "[TEST_CRASH] Binary: " << test_state.binary_name << std::endl;
    std::cerr << "[TEST_CRASH] Context: " << test_state.current_test_context << std::endl;
    std::cerr << "[TEST_CRASH] Tests passed before crash: " << test_state.passed_tests.size() << std::endl;

    // Re-raise signal to get core dump
    std::signal(signal, SIG_DFL);
    std::raise(signal);
}

// Initialize test reporting system
inline void init(const std::string& binary_name) {
    test_state.binary_name = binary_name;
    test_state.current_test_context = "initialization";
    test_state.start_time = std::time(nullptr);
    test_state.initialized = true;

    // Install signal handlers
    std::signal(SIGSEGV, crash_handler);
    std::signal(SIGABRT, crash_handler);
    std::signal(SIGFPE, crash_handler);
    std::signal(SIGILL, crash_handler);

    std::cout << "[TEST_BINARY] " << binary_name << std::endl;
}

// Update current test context (for crash diagnostics)
inline void set_context(const std::string& context) {
    test_state.current_test_context = context;
}

// Start a test
inline void start_test(const std::string& test_name, const std::string& description = "") {
    test_state.current_test_context = test_name;
    if (description.empty()) {
        std::cout << "[TEST_START] " << test_name << std::endl;
    } else {
        std::cout << "[TEST_START] " << test_name << ": " << description << std::endl;
    }
}

// Mark test as passed
inline void pass_test(const std::string& test_name) {
    test_state.passed_tests.push_back(test_name);
    std::cout << "[TEST_PASS] " << test_name << std::endl;
}

// Mark test as failed
inline void fail_test(const std::string& test_name, const std::string& reason = "") {
    test_state.failed_tests.push_back({test_name, reason});
    if (reason.empty()) {
        std::cerr << "[TEST_FAIL] " << test_name << std::endl;
    } else {
        std::cerr << "[TEST_FAIL] " << test_name << ": " << reason << std::endl;
    }
}

// Scenario reporting (for multi-scenario tests)
inline void start_scenario(const std::string& scenario_name, const std::string& description = "") {
    test_state.current_test_context = "scenario: " + scenario_name;
    if (description.empty()) {
        std::cout << "[SCENARIO_START] " << scenario_name << std::endl;
    } else {
        std::cout << "[SCENARIO_START] " << scenario_name << ": " << description << std::endl;
    }
}

inline void pass_scenario(const std::string& scenario_name, const std::string& details = "") {
    if (details.empty()) {
        std::cout << "[SCENARIO_PASS] " << scenario_name << std::endl;
    } else {
        std::cout << "[SCENARIO_PASS] " << scenario_name << ": " << details << std::endl;
    }
}

inline void fail_scenario(const std::string& scenario_name, const std::string& reason) {
    std::cerr << "[SCENARIO_FAIL] " << scenario_name << ": " << reason << std::endl;
}

// Iteration reporting (for iterative tests)
inline void iteration(int current, int total, const std::string& stage_name = "") {
    std::ostringstream ctx;
    ctx << "iteration " << current << "/" << total;
    if (!stage_name.empty()) {
        ctx << " (" << stage_name << ")";
    }
    test_state.current_test_context = ctx.str();
    std::cout << "[ITERATION] " << current << "/" << total;
    if (!stage_name.empty()) {
        std::cout << ": " << stage_name;
    }
    std::cout << std::endl;
}

// Print test summary
inline void print_summary() {
    if (!test_state.initialized) return;

    time_t end_time = std::time(nullptr);
    double duration = std::difftime(end_time, test_state.start_time);

    size_t total = test_state.passed_tests.size() + test_state.failed_tests.size();

    std::cout << "\n[TEST_SUMMARY]" << std::endl;
    std::cout << "  Binary: " << test_state.binary_name << std::endl;
    std::cout << "  Total tests: " << total << std::endl;
    std::cout << "  Passed: " << test_state.passed_tests.size() << std::endl;
    std::cout << "  Failed: " << test_state.failed_tests.size() << std::endl;
    std::cout << "  Duration: " << std::fixed << std::setprecision(2) << duration << "s" << std::endl;

    if (!test_state.failed_tests.empty()) {
        std::cout << "\n[FAILED_TESTS]" << std::endl;
        for (const auto& fail : test_state.failed_tests) {
            std::cout << "  - " << fail.first;
            if (!fail.second.empty()) {
                std::cout << ": " << fail.second;
            }
            std::cout << std::endl;
        }
    }
}

// Return non-zero when any test failed
inline int exit_code() {
    return test_state.failed_tests.empty() ? 0 : 1;
}

// Helper to check condition and report
inline bool check(bool condition, const std::string& test_name, const std::string& failure_msg = "") {
    if (condition) {
        pass_test(test_name);
        return true;
    } else {
        fail_test(test_name, failure_msg);
        return false;
    }
}

} // namespace TestReporting

// Convenience macros for cleaner test code
#define TEST_INIT(name) TestReporting::init(name)
#define TEST_CONTEXT(ctx) TestReporting::set_context(ctx)
#define TEST_START(name, ...) TestReporting::start_test(name, ##__VA_ARGS__)
#define TEST_PASS(name) TestReporting::pass_test(name)
#define TEST_FAIL(name, ...) TestReporting::fail_test(name, ##__VA_ARGS__)
#define TEST_CHECK(cond, name, ...) TestReporting::check(cond, name, ##__VA_ARGS__)
#define TEST_SUMMARY() TestReporting::print_summary()

#define SCENARIO_START(name, ...) TestReporting::start_scenario(name, ##__VA_ARGS__)
#define SCENARIO_PASS(name, ...) TestReporting::pass_scenario(name, ##__VA_ARGS__)
#define SCENARIO_FAIL(name, reason) TestReporting::fail_scenario(name, reason)

#define ITERATION(current, total, ...) TestReporting::iteration(current, total, ##__VA_ARGS__)

#endif // TEST_REPORTING_H
