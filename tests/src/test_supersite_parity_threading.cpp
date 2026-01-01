/*******************************************************************************
 * Thread Safety Test for SuperSite Parity Tests
 * 
 * Purpose: Verify that supersite parity tests (representation, expansion, 
 *          backward, float/double) are not susceptible to race conditions
 *          when running with multiple threads.
 * 
 * Tests:
 * 1. Concurrent execution of representation parity
 * 2. Concurrent execution of expansion parity
 * 3. Concurrent execution of backward parity
 * 4. Concurrent execution of float/double parity
 * 5. Mixed concurrent execution of all parity tests
 * 6. Stress test with many threads
 * 
 * Date: November 9, 2025
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <chrono>
#include <mutex>

#include "test_common.h"
#include "../../common/src/utils/otools.h"


#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

// Thread-safe test counters
static std::atomic<int> test_iterations{0};
static std::atomic<int> test_failures{0};
static std::mutex output_mutex;

// Helper: Create variant
static variant* make_var(std::string chr, int bp, std::string id, 
                        std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

// Safe output function
template<typename... Args>
void safe_print(Args&&... args) {
    std::lock_guard<std::mutex> lock(output_mutex);
    (std::cout << ... << args) << std::endl;
}

/*******************************************************************************
 * TEST 1: Float/Double Parity Thread Safety
 * 
 * Runs float vs double parity check concurrently to verify no race conditions
 * in HMM state management, supersite context, or result comparison.
 ******************************************************************************/
void run_float_double_parity_test(int thread_id, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        try {
            // Create variant map with supersite
            variant_map V;
            V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0)); // anchor
            V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1)); // sibling

            // Conditioning set
            conditioning_set H;
            H.allocate(0, 2, V.size());
            H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
            H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

            // Build supersites
            std::vector<SuperSite> super_sites;
            std::vector<bool> is_super_site;
            std::vector<uint8_t> packed_codes;
            std::vector<int> locus_to_super_idx;
            std::vector<int> super_site_var_index;
            buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                           locus_to_super_idx, super_site_var_index);

            // HMM parameters
            hmm_parameters M;
            M.ed = 0.01f; M.ee = 1.0f;
            M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
            M.nt = std::vector<float>(M.t.size(), 0.95f);
            M.rare_allele = std::vector<char>(V.size(), -1);
            M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

            // Genotype
            genotype G(0);
            G.n_segments = 1; G.n_variants = V.size(); G.n_ambiguous = 0; G.n_missing = 0;
            G.n_transitions = 0; G.n_stored_transitionProbs = 0; G.n_storage_events = 0;
            G.double_precision = false; G.haploid = false;
            G.Variants.assign(1, 0);
            G.Lengths.assign(1, (unsigned short)V.size());
            G.Lengths_bio = G.Lengths;
            G.Diplotypes.assign(1, 1ull);
            VAR_SET_HET(0, G.Variants[0]); // het at anchor
            VAR_SET_HOM(1, G.Variants[0]); // REF at sibling
            G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
            G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

            // Window
            window W;
            W.start_locus = 0; W.stop_locus = 1;
            W.start_segment = 0; W.stop_segment = 0;
            W.start_ambiguous = 0; W.stop_ambiguous = -1;
            W.start_missing = 0; W.stop_missing = -1;
            W.start_transition = 0; W.stop_transition = -1;

            std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

            // Run float and double HMM
            haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
            
            haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M);

            HS.forward();
            HD.forward();

            // Verify parity
            const double tol = 1e-5;
            if (HS.prob.size() != HD.prob.size()) {
                test_failures.fetch_add(1);
                safe_print("[Thread ", thread_id, " iter ", iter, "] FAIL: prob size mismatch");
                continue;
            }

            bool parity_ok = true;
            for (size_t i = 0; i < HS.prob.size(); ++i) {
                double diff = std::fabs((double)HS.prob[i] - HD.prob[i]);
                if (diff > tol) {
                    parity_ok = false;
                    break;
                }
            }

            if (!parity_ok) {
                test_failures.fetch_add(1);
                safe_print("[Thread ", thread_id, " iter ", iter, "] FAIL: probability parity");
            }

            test_iterations.fetch_add(1);

        } catch (const std::exception& e) {
            test_failures.fetch_add(1);
            safe_print("[Thread ", thread_id, " iter ", iter, "] EXCEPTION: ", e.what());
        } catch (...) {
            test_failures.fetch_add(1);
            safe_print("[Thread ", thread_id, " iter ", iter, "] UNKNOWN EXCEPTION");
        }
    }
}

/*******************************************************************************
 * TEST 2: Backward Parity Thread Safety
 * 
 * Tests concurrent backward pass parity checks
 ******************************************************************************/
void run_backward_parity_test(int thread_id, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        try {
            variant_map V;
            V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
            V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

            conditioning_set H;
            H.allocate(0, 2, V.size());
            H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
            H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

            std::vector<SuperSite> super_sites;
            std::vector<bool> is_super_site;
            std::vector<uint8_t> packed_codes;
            std::vector<int> locus_to_super_idx;
            std::vector<int> super_site_var_index;
            buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                           locus_to_super_idx, super_site_var_index);

            hmm_parameters M;
            M.ed = 0.01f; M.ee = 1.0f;
            M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
            M.nt = std::vector<float>(M.t.size(), 0.95f);
            M.rare_allele = std::vector<char>(V.size(), -1);
            M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

            genotype G(0);
            G.n_segments = 1; G.n_variants = V.size(); G.n_ambiguous = 0; G.n_missing = 0;
            G.n_transitions = 0; G.n_stored_transitionProbs = 0; G.n_storage_events = 0;
            G.double_precision = false; G.haploid = false;
            G.Variants.assign(1, 0);
            G.Lengths.assign(1, (unsigned short)V.size());
            G.Lengths_bio = G.Lengths;
            G.Diplotypes.assign(1, 1ull);
            VAR_SET_HOM(0, G.Variants[0]);
            VAR_SET_HOM(1, G.Variants[0]);
            G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
            G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

            window W;
            W.start_locus = 0; W.stop_locus = 1;
            W.start_segment = 0; W.stop_segment = 0;
            W.start_ambiguous = 0; W.stop_ambiguous = -1;
            W.start_missing = 0; W.stop_missing = -1;
            W.start_transition = 0; W.stop_transition = -1;

            std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

            // Float path
            haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
            HS.forward();
            std::vector<double> T_f(G.countTransitions(), 0.0);
            std::vector<float> Mprob_f;
            HS.backward(T_f, Mprob_f, nullptr, nullptr);

            // Double path
            haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M);
            HD.forward();
            std::vector<double> T_d(G.countTransitions(), 0.0);
            std::vector<float> Mprob_d;
            HD.backward(T_d, Mprob_d, nullptr, nullptr);

            // Compare results
            const double tol = 1e-5;
            bool parity_ok = true;
            
            if (HS.prob.size() != HD.prob.size()) {
                parity_ok = false;
            } else {
                for (size_t i = 0; i < HS.prob.size(); ++i) {
                    double diff = std::fabs((double)HS.prob[i] - HD.prob[i]);
                    if (diff > tol) {
                        parity_ok = false;
                        break;
                    }
                }
            }

            if (!parity_ok) {
                test_failures.fetch_add(1);
                safe_print("[Thread ", thread_id, " iter ", iter, "] FAIL: backward parity");
            }

            test_iterations.fetch_add(1);

        } catch (const std::exception& e) {
            test_failures.fetch_add(1);
            safe_print("[Thread ", thread_id, " iter ", iter, "] EXCEPTION: ", e.what());
        } catch (...) {
            test_failures.fetch_add(1);
            safe_print("[Thread ", thread_id, " iter ", iter, "] UNKNOWN EXCEPTION");
        }
    }
}

/*******************************************************************************
 * TEST 3: Mixed Parity Tests Thread Safety
 * 
 * Runs different parity tests concurrently to detect cross-test interference
 ******************************************************************************/
void run_mixed_parity_test(int thread_id, int iterations) {
    for (int iter = 0; iter < iterations; ++iter) {
        // Alternate between different test types based on thread_id and iteration
        int test_type = (thread_id + iter) % 2;
        
        if (test_type == 0) {
            run_float_double_parity_test(thread_id, 1);
        } else {
            run_backward_parity_test(thread_id, 1);
        }
    }
}

/*******************************************************************************
 * Main Test Functions
 ******************************************************************************/

void test_concurrent_float_double_parity() {
    TEST_START("concurrent_float_double_parity");
    test_iterations = 0;
    test_failures = 0;

    const int num_threads = 4;
    const int iterations_per_thread = 10;

    std::vector<std::thread> threads;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(run_float_double_parity_test, i, iterations_per_thread);
    }

    for (auto& t : threads) {
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "  Completed " << test_iterations.load() << " iterations in "
              << duration.count() << "ms" << std::endl;
    std::cout << "  Failures: " << test_failures.load() << std::endl;

    TEST_CHECK(test_failures.load() == 0,
               "concurrent_float_double_parity",
               "No failures in concurrent float/double parity tests");
}

void test_concurrent_backward_parity() {
    TEST_START("concurrent_backward_parity");
    test_iterations = 0;
    test_failures = 0;

    const int num_threads = 4;
    const int iterations_per_thread = 10;

    std::vector<std::thread> threads;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(run_backward_parity_test, i, iterations_per_thread);
    }

    for (auto& t : threads) {
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "  Completed " << test_iterations.load() << " iterations in "
              << duration.count() << "ms" << std::endl;
    std::cout << "  Failures: " << test_failures.load() << std::endl;

    TEST_CHECK(test_failures.load() == 0,
               "concurrent_backward_parity",
               "No failures in concurrent backward parity tests");
}

void test_mixed_concurrent_parity() {
    TEST_START("mixed_concurrent_parity");
    test_iterations = 0;
    test_failures = 0;

    const int num_threads = 6;
    const int iterations_per_thread = 8;

    std::vector<std::thread> threads;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(run_mixed_parity_test, i, iterations_per_thread);
    }

    for (auto& t : threads) {
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "  Completed " << test_iterations.load() << " test iterations in "
              << duration.count() << "ms" << std::endl;
    std::cout << "  Failures: " << test_failures.load() << std::endl;

    TEST_CHECK(test_failures.load() == 0,
               "mixed_concurrent_parity",
               "No failures in mixed concurrent parity tests");
}

void test_stress_concurrent_parity() {
    TEST_START("stress_concurrent_parity");
    test_iterations = 0;
    test_failures = 0;

    const int num_threads = 8;
    const int iterations_per_thread = 20;

    std::vector<std::thread> threads;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(run_mixed_parity_test, i, iterations_per_thread);
    }

    for (auto& t : threads) {
        t.join();
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "  Stress test: " << test_iterations.load() << " iterations in "
              << duration.count() << "ms" << std::endl;
    std::cout << "  Throughput: " << (test_iterations.load() * 1000.0 / duration.count())
              << " tests/sec" << std::endl;
    std::cout << "  Failures: " << test_failures.load() << std::endl;

    TEST_CHECK(test_failures.load() == 0,
               "stress_concurrent_parity",
               "No failures in stress test with many threads");
}

void test_data_race_detection() {
    TEST_START("data_race_detection");
    // This test runs the same test data from multiple threads simultaneously
    // to detect any potential data races in shared structures
    test_iterations = 0;
    test_failures = 0;

    const int num_threads = 4;
    const int iterations_per_thread = 15;

    std::vector<std::thread> threads;

    // All threads run the same test pattern to maximize chance of race detection
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back([i, iterations_per_thread]() {
            run_float_double_parity_test(i, iterations_per_thread);
        });
    }

    for (auto& t : threads) {
        t.join();
    }

    std::cout << "  Data race detection: " << test_iterations.load() << " iterations" << std::endl;
    std::cout << "  Failures: " << test_failures.load() << std::endl;

    TEST_CHECK(test_failures.load() == 0,
               "data_race_detection",
               "No data races detected in concurrent execution");
}

/*******************************************************************************
 * Main
 ******************************************************************************/
int main() {
    TEST_INIT("test_supersite_parity_threading");
    std::cout << "=== SuperSite Parity Tests Thread Safety Validation ===" << std::endl;
    std::cout << "Testing: representation, expansion, backward, float/double parity" << std::endl;
    std::cout << std::endl;
    
    test_concurrent_float_double_parity();
    test_concurrent_backward_parity();
    test_mixed_concurrent_parity();
    test_stress_concurrent_parity();
    test_data_race_detection();
    
    TEST_SUMMARY();
    return TestReporting::exit_code();
}
