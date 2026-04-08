#
# Makefile for heliolinx
#
# To build all programs (in parallel), run:      make -j4
# To install to ../bin, run:                     sudo make install
# To install to (say) /usr/local, run:           PREFIX=/usr/local make install
#

PREFIX ?= ..
CXXFLAGS ?= -O3

###

PROGRAMS_PATH = $(PREFIX)/bin
PROGRAMS = make_tracklets make_outim heliolinc heliolinc_omp heliolinc_lowmem heliolinc_lowmem_omp heliovane link_purify link_planarity link_planarity_omp link_purify_chisq parse_clust2det_MPC80 parse_clust2det modsplit_hlfile merge_tracklet_files merge_det_catalogs make_trailed_tracklets parse_trk2det calc_heliohypmat label_hldet label_hldet_mpc helio_highgrade helio_highgrade_omp analyze_linkage01a make_hypmask prefilter_tracklets split_tracklets_by_time

LIB = libheliolinx.a
LIB_SOURCES = solarsyst_dyn_geo01.cpp

SOURCES = $(LIB_SOURCES) $(PROGRAMS:%=%.cpp)

.PHONY: all
all: $(PROGRAMS)

# Shell scripts installed alongside compiled binaries
SCRIPTS = heliolinc_zonesplit

.PHONY: install
install: $(PROGRAMS)
	mkdir -p $(PREFIX)/bin
	cp -p $(PROGRAMS) $(PREFIX)/bin
	cp -p $(SCRIPTS) $(PREFIX)/bin
	chmod +x $(addprefix $(PREFIX)/bin/, $(SCRIPTS))

$(PROGRAMS): %: %.o $(LIB)
	$(CXX) $< -L. -lheliolinx -fopenmp $(OUTPUT_OPTION)

$(LIB): $(LIB_SOURCES:%.cpp=%.o)
	$(AR) -rc $@ $^

clean:
	rm -f $(PROGRAMS) $(LIB) $(SOURCES:%.cpp=%.o)
	rm -rf .deps

###

DEPDIR := .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d && touch $@

%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -std=c++11 -fopenmp -I../include -c $(OUTPUT_OPTION) $<
	@ $(POSTCOMPILE)

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SOURCES:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))
