ifndef SHAPEIT5_COMMON_INCLUDED
SHAPEIT5_COMMON_INCLUDED := 1
# (common makefile already included above)
endif

CXXFLAGS += -I../common/src -I../phase_common/src -DSHAPEIT5_TEST_TRACE
LOCAL_INC := -I../common/src -I../phase_common/src

include ../common/makefile_common.mk

objdir := obj
bindir := bin

# Phase_common objects required by some tests (define BEFORE using it)
PHASE_COMMON_OBJS := \
	$(objdir)/phase_variant.o \
	$(objdir)/phase_variant_map.o \
	$(objdir)/phase_bitmatrix.o \
	$(objdir)/phase_bitvector.o \
	$(objdir)/phase_genotype_set.o \
	$(objdir)/phase_haplotype_set.o \
	$(objdir)/phase_conditioning_set_managment.o \
	$(objdir)/phase_conditioning_set_selection.o \
	$(objdir)/phase_super_site_builder.o \
	$(objdir)/phase_ibd2_tracks.o \
	$(objdir)/haplotype_segment_single.o \
	$(objdir)/haplotype_segment_double.o \
	$(objdir)/hmm_parameters.o \
	$(objdir)/compute_job.o \
	$(objdir)/window_set.o \
	$(objdir)/genotype_build.o \
	$(objdir)/genotype_managment_from_pc.o \
	$(objdir)/genotype_prune.o \
	$(objdir)/genotype_sweep.o \
	$(objdir)/genotype_mendel.o \
	$(objdir)/test_toolbox.o

TEST_BINARIES := \
    $(bindir)/test_supersite_emissions \
    $(bindir)/test_supersite_emissions_real \
    $(bindir)/test_supersite_accessor \
    $(bindir)/test_supersite_unpack \
    $(bindir)/test_supersite_builder \
	$(bindir)/test_supersite_anchor_encoding \
    $(bindir)/test_supersite_hmm \
    $(bindir)/test_supersite_hmm_states \
    $(bindir)/test_missing_multiallelic_multinomial \
    $(bindir)/test_supersite_vs_biallelic_simple \
    $(bindir)/test_forward_backward_comprehensive \
    $(bindir)/test_packing_format_diagnostic \
    $(bindir)/test_supersite_anchor_gating \
    $(bindir)/test_supersite_backward_sc \
    $(bindir)/test_supersite_mutual_exclusivity_make \
    $(bindir)/test_supersite_backward_projection \
    $(bindir)/test_supersite_boundary_conditions \
    $(bindir)/test_supersite_float_double_parity \
    $(bindir)/test_supersite_class_microcases \
    $(bindir)/test_supersite_no_double_counting \
    $(bindir)/test_supersite_backward_parity \
    $(bindir)/test_supersite_perf_smoke \
    $(bindir)/test_supersite_combine \
	$(bindir)/test_supersite_outer_product \
    $(bindir)/test_supersite_representation_parity \
    $(bindir)/test_supersite_expansion_parity \
    $(bindir)/test_supersite_expansion_epochs \
    $(bindir)/test_supersite_expansion_epochs_w_missing \
    $(bindir)/test_segment_boundary_multiallelic \
    $(bindir)/test_segment_boundary_invariants \
    $(bindir)/test_supersite_transition_underflow \
    $(bindir)/test_empty_segment \
    $(bindir)/test_empty_segment_underflow \
	$(bindir)/test_hmm_indexing_supersite \
	$(bindir)/test_supersite_cursor_parity \
	$(bindir)/test_k_inflation_measurement \
	$(bindir)/test_hmm_state_inflation \
	$(bindir)/test_emission_pattern_validation \
	$(bindir)/test_phase_decision_accuracy

# Thread Safety Tests
THREADING_TEST_BINARIES := \
    $(bindir)/test_phaser_threading_safety \
    $(bindir)/test_conditioning_set_threading \
    $(bindir)/test_supersite_threading_safety \
    $(bindir)/test_ibd2_threading_safety \
    $(bindir)/test_supersite_parity_threading

# Default goal: only build test binaries (do not build bin/tests)
.DEFAULT_GOAL := test-binaries
.PHONY: test-binaries threading-tests all-tests
test-binaries: $(TEST_BINARIES)
threading-tests: $(THREADING_TEST_BINARIES)
all-tests: $(TEST_BINARIES) $(THREADING_TEST_BINARIES)

# Link rules for each test
$(bindir)/test_supersite_emissions: $(objdir)/test_supersite_emissions.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_emissions_real: $(objdir)/test_supersite_emissions_real.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_accessor: $(objdir)/test_supersite_accessor.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_unpack: $(objdir)/test_supersite_unpack.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_builder: $(objdir)/test_supersite_builder.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_anchor_encoding: $(objdir)/test_supersite_anchor_encoding.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) --coverage -lgcov

$(bindir)/test_supersite_hmm: $(objdir)/test_supersite_hmm.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_hmm_states: $(objdir)/test_supersite_hmm_states.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_anchor_gating: $(objdir)/test_supersite_anchor_gating.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_backward_sc: $(objdir)/test_supersite_backward_sc.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_mutual_exclusivity_make: $(objdir)/test_supersite_mutual_exclusivity_make.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_backward_projection: $(objdir)/test_supersite_backward_projection.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_boundary_conditions: $(objdir)/test_supersite_boundary_conditions.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_float_double_parity: $(objdir)/test_supersite_float_double_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_class_microcases: $(objdir)/test_supersite_class_microcases.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_no_double_counting: $(objdir)/test_supersite_no_double_counting.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_backward_parity: $(objdir)/test_supersite_backward_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_perf_smoke: $(objdir)/test_supersite_perf_smoke.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_combine: $(objdir)/test_supersite_combine.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_outer_product: $(objdir)/test_supersite_outer_product.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)
$(bindir)/test_missing_multiallelic_multinomial: $(objdir)/test_missing_multiallelic_multinomial.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_vs_biallelic_simple: $(objdir)/test_supersite_vs_biallelic_simple.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_representation_parity: $(objdir)/test_supersite_representation_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_expansion_parity: $(objdir)/test_supersite_expansion_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_expansion_epochs: $(objdir)/test_supersite_expansion_epochs.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_expansion_epochs_w_missing: $(objdir)/test_supersite_expansion_epochs_w_missing.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_forward_backward_comprehensive: $(objdir)/test_forward_backward_comprehensive.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_packing_format_diagnostic: $(objdir)/test_packing_format_diagnostic.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_segment_boundary_multiallelic: $(objdir)/test_segment_boundary_multiallelic.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_segment_boundary_invariants: $(objdir)/test_segment_boundary_invariants.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_transition_underflow: $(objdir)/test_supersite_transition_underflow.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_empty_segment: $(objdir)/test_empty_segment.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_empty_segment_underflow: $(objdir)/test_empty_segment_underflow.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_cursor_drift_investigation: $(objdir)/test_cursor_drift_investigation.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_k_inflation_measurement: $(objdir)/test_k_inflation_measurement.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_hmm_state_inflation: $(objdir)/test_hmm_state_inflation.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_emission_pattern_validation: $(objdir)/test_emission_pattern_validation.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_phase_decision_accuracy: $(objdir)/test_phase_decision_accuracy.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_emission_logic_differences: $(objdir)/test_emission_logic_differences.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_pbwt_conditioning_parity: $(objdir)/test_pbwt_conditioning_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_hmm_indexing_supersite: $(objdir)/test_hmm_indexing_supersite.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_cursor_parity: $(objdir)/test_supersite_cursor_parity.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_expected_codes: $(objdir)/test_supersite_expected_codes.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

# Thread Safety Test Link Rules
$(bindir)/test_phaser_threading_safety: $(objdir)/test_phaser_threading_safety.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) -lpthread

$(bindir)/test_conditioning_set_threading: $(objdir)/test_conditioning_set_threading.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) -lpthread

$(bindir)/test_supersite_threading_safety: $(objdir)/test_supersite_threading_safety.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) -lpthread

$(bindir)/test_ibd2_threading_safety: $(objdir)/test_ibd2_threading_safety.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) -lpthread

$(bindir)/test_supersite_parity_threading: $(objdir)/test_supersite_parity_threading.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS) -lpthread

# Explicit compile rules for external sources (not under tests/src)
$(objdir)/genotype_managment_from_pc.o: ../phase_common/src/objects/genotype/genotype_managment.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/genotype_build.o: ../phase_common/src/objects/genotype/genotype_build.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/genotype_prune.o: ../phase_common/src/objects/genotype/genotype_prune.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_variant.o: ../phase_common/src/objects/variant.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_variant_map.o: ../phase_common/src/containers/variant_map.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_bitmatrix.o: ../phase_common/src/containers/bitmatrix.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_bitvector.o: ../phase_common/src/containers/bitvector.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_haplotype_set.o: ../phase_common/src/containers/haplotype_set.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)


$(objdir)/phase_conditioning_set_managment.o: ../phase_common/src/containers/conditioning_set/conditioning_set_managment.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_conditioning_set_selection.o: ../phase_common/src/containers/conditioning_set/conditioning_set_selection.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_super_site_builder.o: ../phase_common/src/objects/super_site_builder.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_ibd2_tracks.o: ../phase_common/src/containers/ibd2_tracks.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/haplotype_segment_single.o: ../phase_common/src/models/haplotype_segment_single.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/haplotype_segment_double.o: ../phase_common/src/models/haplotype_segment_double.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/hmm_parameters.o: ../phase_common/src/objects/hmm_parameters.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/compute_job.o: ../phase_common/src/objects/compute_job.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/window_set.o: ../phase_common/src/containers/window_set.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

# Compile rule for genotype_set
$(objdir)/phase_genotype_set.o: ../phase_common/src/containers/genotype_set.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/genotype_sweep.o: ../phase_common/src/objects/genotype/genotype_sweep.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/genotype_mendel.o: ../phase_common/src/objects/genotype/genotype_mendel.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

# Provide toolbox globals exactly once
$(objdir)/test_toolbox.o: src/test_toolbox.cpp
	$(CXX) $(CXXFLAG) -D_DECLARE_TOOLBOX_HERE -c $< -o $@ -Isrc $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

# Pattern rule for all other test source files
$(objdir)/test_%.o: src/test_%.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ -Isrc $(LOCAL_INC) -I$(HTSLIB_INC) -I$(BOOST_INC)

# Pattern rule above handles all test_*.cpp compilation automatically

.PHONY: test-run test-threading-run
test-run: test-binaries
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_emissions
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_emissions_real
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_accessor
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_unpack
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_builder
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_missing_multiallelic_multinomial
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_vs_biallelic_simple
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_forward_backward_comprehensive
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_packing_format_diagnostic

# Regression check: ensure backward IMPUTE guard is present (prevents AlphaMissing OOB regressions)
	@echo "Running regression check: backward IMPUTE guard script"
	./scripts/check_backward_impute_guard.sh || (echo "ERROR: backward IMPUTE guard regression failed" && exit 1)
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_expansion_parity
	# Strict unit test: verify ambiguous-cursor parity PASS diagnostics for both precisions
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_cursor_parity || (echo "ERROR: test_supersite_cursor_parity failed" && exit 1)

# Run threading safety tests
test-threading-run: threading-tests
	@echo "=== Running Thread Safety Tests ==="
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_phaser_threading_safety
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_conditioning_set_threading
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_threading_safety
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_ibd2_threading_safety
	@echo "=== Thread Safety Tests Complete ==="

# Add extra cleaning without conflicting with common 'clean'
.PHONY: clean-extra
clean: clean-extra

clean-extra:
	$(RM) $(TEST_BINARIES) $(THREADING_TEST_BINARIES) $(PHASE_COMMON_OBJS)
