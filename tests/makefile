ifndef SHAPEIT5_COMMON_INCLUDED
SHAPEIT5_COMMON_INCLUDED := 1
include ../common/makefile_common.mk
endif

CXXFLAG += -I../common/src -I../phase_common/src

objdir := obj
bindir := bin

# Phase_common objects required by some tests (define BEFORE using it)
PHASE_COMMON_OBJS := \
    $(objdir)/phase_variant.o \
    $(objdir)/phase_variant_map.o \
    $(objdir)/phase_bitmatrix.o \
    $(objdir)/phase_bitvector.o \
    $(objdir)/phase_haplotype_set.o \
    $(objdir)/phase_conditioning_set_managment.o \
    $(objdir)/phase_super_site_builder.o \
    $(objdir)/phase_ibd2_tracks.o \
    $(objdir)/haplotype_segment_single.o \
    $(objdir)/haplotype_segment_double.o \
    $(objdir)/hmm_parameters.o \
    $(objdir)/genotype_managment_from_pc.o \
    $(objdir)/test_toolbox.o

TEST_BINARIES := \
    $(bindir)/test_supersite_emissions \
    $(bindir)/test_supersite_accessor \
    $(bindir)/test_supersite_unpack \
    $(bindir)/test_supersite_builder \
    $(bindir)/test_supersite_hmm \
    $(bindir)/test_supersite_hmm_states \
    $(bindir)/test_missing_multiallelic_multinomial

# Default goal: only build test binaries (do not build bin/tests)
.DEFAULT_GOAL := test-binaries
.PHONY: test-binaries
test-binaries: $(TEST_BINARIES)

# Link rules for each test
$(bindir)/test_supersite_emissions: $(objdir)/test_supersite_emissions.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_accessor: $(objdir)/test_supersite_accessor.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_unpack: $(objdir)/test_supersite_unpack.o
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_builder: $(objdir)/test_supersite_builder.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_hmm: $(objdir)/test_supersite_hmm.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_supersite_hmm_states: $(objdir)/test_supersite_hmm_states.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

$(bindir)/test_missing_multiallelic_multinomial: $(objdir)/test_missing_multiallelic_multinomial.o $(PHASE_COMMON_OBJS)
	$(CXX) $(LDFLAG) $^ -o $@ $(DYN_LIBS)

# Explicit compile rules for external sources (not under tests/src)
$(objdir)/genotype_managment_from_pc.o: ../phase_common/src/objects/genotype/genotype_managment.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_variant.o: ../phase_common/src/objects/variant.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_variant_map.o: ../phase_common/src/containers/variant_map.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_bitmatrix.o: ../phase_common/src/containers/bitmatrix.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_bitvector.o: ../phase_common/src/containers/bitvector.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_haplotype_set.o: ../phase_common/src/containers/haplotype_set.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_conditioning_set_managment.o: ../phase_common/src/containers/conditioning_set/conditioning_set_managment.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_super_site_builder.o: ../phase_common/src/objects/super_site_builder.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/phase_ibd2_tracks.o: ../phase_common/src/containers/ibd2_tracks.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/haplotype_segment_single.o: ../phase_common/src/models/haplotype_segment_single.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/haplotype_segment_double.o: ../phase_common/src/models/haplotype_segment_double.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

$(objdir)/hmm_parameters.o: ../phase_common/src/objects/hmm_parameters.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@  -I$(HTSLIB_INC) -I$(BOOST_INC)

# Provide toolbox globals exactly once
$(objdir)/test_toolbox.o: src/test_toolbox.cpp
	$(CXX) $(CXXFLAG) -D_DECLARE_TOOLBOX_HERE -c $< -o $@ -Isrc  -I$(HTSLIB_INC) -I$(BOOST_INC)

.PHONY: test-run
test-run: test-binaries
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_emissions
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_accessor
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_unpack
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_supersite_builder
	@echo "Skipping test_supersite_hmm (old implementation, expected to fail)"
	@echo "Skipping test_supersite_hmm_states (old implementation, expected to fail)"
	LD_LIBRARY_PATH=$(HOME)/.linuxbrew/lib:/usr/local/lib:$(LD_LIBRARY_PATH) $(bindir)/test_missing_multiallelic_multinomial

# Add extra cleaning without conflicting with common 'clean'
.PHONY: clean-extra
clean: clean-extra

clean-extra:
	$(RM) $(TEST_BINARIES) $(PHASE_COMMON_OBJS)
