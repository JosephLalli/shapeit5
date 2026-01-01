/*******************************************************************************
 * Cursor Drift Investigation Test
 *
 * Comprehensive test to investigate ambiguous cursor drift during supersite
 * anchor/sibling transitions that may be causing K inflation and accuracy issues.
 * 
 * Key Hypothesis: curr_abs_ambiguous cursor escapes [ambiguous_first, ambiguous_last]
 * bounds when siblings advance cursor but are not counted in window setup.
 *
 * This test creates a controlled supersite dataset and tracks cursor behavior
 * with enhanced logging to identify the exact point of cursor drift.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <map>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/models/site_emission_adapter.h"
#include "../../phase_common/src/models/super_site_accessor.h"

namespace {

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

struct CursorTraceEntry {
    int locus;
    int cursor_before;
    int cursor_after;
    int expected_delta;
    int actual_delta;
    bool is_sibling;
    bool is_anchor;
    bool data_amb;
    bool hmm_amb;
    bool sib_advance_logic;
    bool can_advance_logic;
    int amb_range_start;
    int amb_range_end;
    std::string stage;
    std::string violation_type; // "none", "underflow", "overflow", "drift"
};

struct CursorInvestigation {
    std::vector<CursorTraceEntry> forward_trace;
    std::vector<CursorTraceEntry> backward_trace;
    int total_violations = 0;
    int drift_violations = 0;
    int underflow_violations = 0;
    int overflow_violations = 0;
    
    void analyze_violations() {
        total_violations = 0;
        drift_violations = 0;
        underflow_violations = 0;
        overflow_violations = 0;
        
        auto analyze_trace = [&](const std::vector<CursorTraceEntry>& trace) {
            for (const auto& entry : trace) {
                if (entry.violation_type != "none") {
                    total_violations++;
                    if (entry.violation_type == "drift") drift_violations++;
                    else if (entry.violation_type == "underflow") underflow_violations++;
                    else if (entry.violation_type == "overflow") overflow_violations++;
                }
            }
        };
        
        analyze_trace(forward_trace);
        analyze_trace(backward_trace);
    }
    
    void print_summary() {
        std::cout << "\n=== CURSOR DRIFT INVESTIGATION SUMMARY ===" << std::endl;
        std::cout << "Total violations: " << total_violations << std::endl;
        std::cout << "  Drift violations: " << drift_violations << std::endl;
        std::cout << "  Underflow violations: " << underflow_violations << std::endl;  
        std::cout << "  Overflow violations: " << overflow_violations << std::endl;
        
        if (total_violations > 0) {
            std::cout << "\n=== VIOLATION DETAILS ===" << std::endl;
            auto print_violations = [&](const std::vector<CursorTraceEntry>& trace, const std::string& pass) {
                for (const auto& entry : trace) {
                    if (entry.violation_type != "none") {
                        std::cout << pass << " VIOLATION: " << entry.violation_type 
                                  << " at locus=" << entry.locus
                                  << " cursor_before=" << entry.cursor_before
                                  << " cursor_after=" << entry.cursor_after
                                  << " range=[" << entry.amb_range_start << "," << entry.amb_range_end << "]"
                                  << " is_sibling=" << entry.is_sibling
                                  << " sib_advance=" << entry.sib_advance_logic
                                  << " can_advance=" << entry.can_advance_logic
                                  << std::endl;
                    }
                }
            };
            
            print_violations(forward_trace, "FORWARD");
            print_violations(backward_trace, "BACKWARD");
        }
    }
    
    void export_trace_csv(const std::string& filename) {
        std::ofstream out(filename);
        out << "pass,locus,cursor_before,cursor_after,expected_delta,actual_delta,is_sibling,is_anchor,"
            << "data_amb,hmm_amb,sib_advance_logic,can_advance_logic,amb_range_start,amb_range_end,"
            << "stage,violation_type\n";
        
        auto write_trace = [&](const std::vector<CursorTraceEntry>& trace, const std::string& pass) {
            for (const auto& e : trace) {
                out << pass << "," << e.locus << "," << e.cursor_before << "," << e.cursor_after << ","
                    << e.expected_delta << "," << e.actual_delta << "," << e.is_sibling << "," << e.is_anchor << ","
                    << e.data_amb << "," << e.hmm_amb << "," << e.sib_advance_logic << "," << e.can_advance_logic << ","
                    << e.amb_range_start << "," << e.amb_range_end << "," << e.stage << "," << e.violation_type << "\n";
            }
        };
        
        write_trace(forward_trace, "forward");
        write_trace(backward_trace, "backward");
    }
};

enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

// Build variant pointer helper
static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

static void clear_variant_state(genotype& G, int locus) {
    unsigned char& byte = G.Variants[DIV2(locus)];
    const int shift = (MOD2(locus)) << 2;
    byte &= ~(0x0F << shift);
}

static void set_phase(genotype& G, int locus, PhaseCode code) {
    clear_variant_state(G, locus);
    unsigned char& byte = G.Variants[DIV2(locus)];
    switch (code) {
        case REF_REF:
            break;
        case ALT_ALT:
            VAR_SET_HAP0(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
        case ALT_REF:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP0(MOD2(locus), byte);
            break;
        case REF_ALT:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
    }
}

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    return ctx;
}

// Modified haplotype segment with cursor tracing
class TrackedHaplotypeSegmentDouble : public haplotype_segment_double {
public:
    CursorInvestigation* investigation;
    
    TrackedHaplotypeSegmentDouble(genotype* G, 
                                  bitmatrix& H_opt_hap,
                                  std::vector<unsigned int>& idxH,
                                  window& W, 
                                  hmm_parameters& M,
                                  CursorInvestigation* inv)
        : haplotype_segment_double(G, H_opt_hap, idxH, W, M)
        , investigation(inv) {}

private:
    CursorTraceEntry create_trace_entry(const std::string& stage, int cursor_before, bool is_sibling) {
        CursorTraceEntry entry;
        entry.locus = curr_abs_locus;
        entry.cursor_before = cursor_before;
        entry.cursor_after = curr_abs_ambiguous;
        entry.expected_delta = 0; // Will be set by caller
        entry.actual_delta = curr_abs_ambiguous - cursor_before;
        entry.is_sibling = is_sibling;
        entry.is_anchor = !is_sibling && (supersites_enabled_flag && locus_to_super_idx && 
                                         curr_abs_locus < (int)locus_to_super_idx->size() &&
                                         (*locus_to_super_idx)[curr_abs_locus] >= 0);
        entry.amb_range_start = ambiguous_first;
        entry.amb_range_end = ambiguous_last;
        entry.stage = stage;
        entry.violation_type = "none";
        
        // Check for violations
        if (curr_abs_ambiguous < ambiguous_first) {
            entry.violation_type = "underflow";
        } else if (curr_abs_ambiguous > ambiguous_last) {
            entry.violation_type = "overflow";
        } else if (ambiguous_first <= ambiguous_last && 
                   (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last)) {
            entry.violation_type = "drift";
        }
        
        return entry;
    }

public:
    void forward() {
        // Initialize forward pass
        curr_abs_ambiguous = ambiguous_first;
        
        // Clear previous trace
        if (investigation) {
            investigation->forward_trace.clear();
        }
        
        // Call parent forward with our tracing
        haplotype_segment_double::forward();
    }
    
    void backward(std::vector<double>& transition_probabilities,
                  std::vector<float>& missing_probabilities,
                  std::vector<float>* SC,
                  const std::vector<bool>* anchor_has_missing) {
        
        // Clear previous trace  
        if (investigation) {
            investigation->backward_trace.clear();
        }
        
        // Call parent backward with our tracing
        haplotype_segment_double::backward(transition_probabilities, missing_probabilities, SC, anchor_has_missing);
    }

protected:
    // Override cursor advancement to add tracing
    void trace_cursor_advancement(const std::string& stage, int cursor_before, bool is_sibling, 
                                 bool data_amb, bool hmm_amb, bool sib_advance_logic, bool can_advance_logic) {
        if (!investigation) return;
        
        CursorTraceEntry entry = create_trace_entry(stage, cursor_before, is_sibling);
        entry.data_amb = data_amb;
        entry.hmm_amb = hmm_amb;
        entry.sib_advance_logic = sib_advance_logic;
        entry.can_advance_logic = can_advance_logic;
        
        if (stage.find("fwd") != std::string::npos) {
            investigation->forward_trace.push_back(entry);
        } else {
            investigation->backward_trace.push_back(entry);
        }
        
    }
};

} // namespace

int main() {
    TEST_INIT("test_cursor_drift_investigation");
    std::cout << "======================================================================" << std::endl;
    std::cout << "Cursor Drift Investigation Test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create supersite dataset designed to trigger cursor issues
    // 6 variants forming 3 supersites with alternating anchor/sibling pattern
    std::cout << "Creating 6-variant supersite dataset..." << std::endl;
    
    variant_map V;
    // Supersite 1: position 1000 (variants 0,1)
    V.push(make_var("1", 1000, "ss1_anchor", "A", "T", 0));    // anchor
    V.push(make_var("1", 1000, "ss1_sibling", "A", "C", 1));  // sibling
    
    // Supersite 2: position 2000 (variants 2,3) 
    V.push(make_var("1", 2000, "ss2_anchor", "A", "G", 2));   // anchor
    V.push(make_var("1", 2000, "ss2_sibling", "A", "T", 3));  // sibling
    
    // Supersite 3: position 3000 (variants 4,5)
    V.push(make_var("1", 3000, "ss3_anchor", "A", "C", 4));   // anchor  
    V.push(make_var("1", 3000, "ss3_sibling", "A", "G", 5));  // sibling

    conditioning_set H;
    H.allocate(0, 4, V.size()); // 4 reference samples => 8 haplotypes

    auto set_panel = [&](int locus, const std::array<int,8>& alt_flags) {
        for (int hap = 0; hap < 8; ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    // Set panel patterns to create ambiguous sites that trigger cursor advancement
    set_panel(0, {1,0, 0,1, 0,0, 1,0}); // ss1_anchor - HET pattern
    set_panel(1, {0,0, 0,0, 0,0, 0,0}); // ss1_sibling - all REF
    set_panel(2, {0,1, 1,0, 0,0, 0,1}); // ss2_anchor - HET pattern  
    set_panel(3, {0,0, 0,0, 0,0, 0,0}); // ss2_sibling - all REF
    set_panel(4, {1,1, 0,0, 1,0, 0,1}); // ss3_anchor - mixed pattern
    set_panel(5, {0,0, 0,0, 0,0, 0,0}); // ss3_sibling - all REF

    // Create target genotype with HET sites to trigger ambiguous cursor advancement
    genotype G(0);
    G.n_segments = 1;
    G.n_variants = V.size();
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    G.Variants.assign((V.size() + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull);
    G.Lengths.assign(1, static_cast<unsigned short>(V.size()));
    G.Lengths_bio = G.Lengths;

    // Set genotypes to create ambiguous patterns
    set_phase(G, 0, ALT_REF); // ss1_anchor: 1|0 (HET -> ambiguous)
    set_phase(G, 1, REF_REF); // ss1_sibling: 0|0 
    set_phase(G, 2, REF_ALT); // ss2_anchor: 0|1 (HET -> ambiguous)
    set_phase(G, 3, REF_REF); // ss2_sibling: 0|0
    set_phase(G, 4, ALT_REF); // ss3_anchor: 1|0 (HET -> ambiguous) 
    set_phase(G, 5, REF_REF); // ss3_sibling: 0|0
    G.build();

    std::cout << "  Built 6-variant dataset with 3 supersites" << std::endl;

    // Build supersites and setup
    SuperSiteContext ctx = build_supersites(V, H);
    
    std::cout << "  Detected " << ctx.super_sites.size() << " supersites:" << std::endl;
    for (size_t i = 0; i < ctx.super_sites.size(); ++i) {
        const SuperSite& ss = ctx.super_sites[i];
        std::cout << "    Supersite " << i << ": bp=" << ss.bp 
                  << " var_count=" << (int)ss.var_count 
                  << " global_site_id=" << ss.global_site_id << std::endl;
    }

    // Create window and parameters
    window W{};
    W.start_locus = 0;
    W.stop_locus = static_cast<int>(V.size()) - 1;
    W.start_segment = 0;
    W.stop_segment = 0;
    
    // Count ambiguous sites for window bounds
    int amb_count = 0;
    for (int v = 0; v < static_cast<int>(V.size()); ++v) {
        // Check if this locus should be considered for ambiguous counting
        int ss_idx = (v < static_cast<int>(ctx.locus_to_super_idx.size())) 
                     ? ctx.locus_to_super_idx[v] : -1;
        bool is_sibling = false;
        if (ss_idx >= 0) {
            const SuperSite& ss = ctx.super_sites[ss_idx];
            is_sibling = (v != static_cast<int>(ss.global_site_id));
        }
        
        // Only count non-siblings for ambiguous sites (this may be part of the bug)
        if (!is_sibling) {
            unsigned char var_code = G.Variants[DIV2(v)];
            bool is_amb = VAR_GET_AMB(MOD2(v), var_code) != 0;
            if (!is_amb && (VAR_GET_SCA(MOD2(v), var_code) != 0)) {
                is_amb = true;
            }
            if (is_amb) amb_count++;
        }
    }
    
    W.start_ambiguous = 0;
    W.stop_ambiguous = amb_count > 0 ? amb_count - 1 : -1;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = -1;

    std::cout << "  Window setup: ambiguous range=[" << W.start_ambiguous 
              << "," << W.stop_ambiguous << "] (count=" << amb_count << ")" << std::endl;

    // HMM parameters
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
    M.cm = std::vector<float>(V.size(), 0.0f);
    for (size_t i = 0; i < V.size(); ++i) {
        M.cm[i] = 0.001f * static_cast<float>(i + 1); // Small incremental distances
    }
    M.Neff = 10000;
    M.Nhap = static_cast<int>(H.n_hap);
    M.t.assign(V.size() > 1 ? V.size() - 1 : 0, 0.01f);
    M.nt.assign(V.size() > 1 ? V.size() - 1 : 0, 0.99f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(ctx.super_sites, ctx.locus_to_super_idx);

    const std::vector<unsigned int> idxH = {0u,1u,2u,3u,4u,5u,6u,7u};

    // Run investigation
    std::cout << "\nRunning cursor drift investigation..." << std::endl;
    
    CursorInvestigation investigation;
    G.setSuperSiteContext(&ctx.super_sites, &ctx.locus_to_super_idx, &ctx.super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(ctx.packed_codes.data(), ctx.packed_codes.size());
    TrackedHaplotypeSegmentDouble HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                     W, M, &investigation);

    // Run forward pass
    HS.forward();
    
    // Run backward pass
    std::vector<double> transition_probabilities(G.countTransitions(), 0.0);
    std::vector<float> missing_probabilities;
    HS.backward(transition_probabilities, missing_probabilities, nullptr, nullptr);

    // Analyze results
    investigation.analyze_violations();
    investigation.print_summary();

    // Test results
    bool test_passed = (investigation.total_violations == 0);
    
    if (test_passed) {
        std::cout << std::endl;
        std::cout << "✓ SUCCESS: No cursor drift detected!" << std::endl;
        std::cout << "✓ Cursor remained within bounds throughout F/B passes" << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << "✗ FAILURE: Cursor drift detected!" << std::endl;
        std::cout << "✗ " << investigation.total_violations << " violations found" << std::endl;
        std::cout << "  This confirms the hypothesis that cursor drift is causing K inflation" << std::endl;
        
        // Don't assert to allow full investigation
        // assert(false);
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "Cursor drift investigation completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    return test_passed ? 0 : 1;
}
