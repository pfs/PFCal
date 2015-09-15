# name of the library
LIBNAME = PFCalEEAnalysis

#Necessary to use shell built-in commands
SHELL=bash

# Define include paths
USERINCLUDES += -I$(ROOTSYS)/include/
USERINCLUDES += -I$(ROOFITSYS)/include/
USERINCLUDES += -isystem $(BOOSTSYS)/include

USERINCLUDES += -Iinclude/ -I../userlib/include/ -I$(FASTJET_INSTALL)/include/
USERINCLUDES += -I$(HEPMC_DIR)/include/

# Define libraries to link
USERLIBS += $(shell root-config --glibs) -lRooFit -lGenVector #-lTreePlayer -lTMVA
USERLIBS += -L../userlib/lib -lPFCalEEuserlib
USERLIBS += -Wl,-rpath,$(FASTJET_INSTALL)/lib -lm  -L$(FASTJET_INSTALL)/lib -lfastjettools -lfastjet -lfastjetplugins -lsiscone_spherical -lsiscone
USERLIBS += -L$(BOOSTSYS)/lib -lboost_regex -lboost_program_options -lboost_filesystem

#CXXFLAGS = -Wall -W -Wno-unused-function -Wno-parentheses -Wno-char-subscripts -Wno-unused-parameter -O2 
CXXFLAGS = -g -Wall -W -O2 -std=c++0x -fPIC# -std=c++11 
#CXXFLAGS = -Wall -W -O2 -std=c++11 
LDFLAGS = -shared -Wall -W


# If possible we'll use the clang compiler, it's faster and gives more helpful error messages
# If it's not available, then fallback to gcc.
CXX=g++
LD=g++
# CLANGPATH := $(shell type -p clang++)
# ifeq ($(CLANGPATH),)
# $(warning clang++ not found, reverting to g++!)
# CXX=g++
# LD=g++
# endif

CXXFLAGS += $(USERINCLUDES)
LIBS += $(USERLIBS)

# A list of directories
BASEDIR = $(shell pwd)
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

# Build a list of srcs and bins to build
SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))
#BINS=$(subst $(TESTDIR), $(EXEDIR),$(subst .$(TEST_EXT),,$(EXES)))

#BINS=$(EXEDIR)/MixPUSignal $(EXEDIR)/egammaResolution $(EXEDIR)/higgsResolution $(EXEDIR)/timeResolution $(EXEDIR)/plotMomentum $(EXEDIR)/plotTime

#BINS=$(EXEDIR)/plotParticleDensity $(EXEDIR)/validation $(EXEDIR)/timeResolution $(EXEDIR)/hepmcHggChargedTracks

#BINS=$(EXEDIR)/mipStudy $(EXEDIR)/mipProfile $(EXEDIR)/mipSelection $(EXEDIR)/ootpuOnHiggsPhotons $(EXEDIR)/higgsResolution $(EXEDIR)/plotNabove100fC

#BINS=$(EXEDIR)/vbfResolution $(EXEDIR)/calibration

#BINS=$(EXEDIR)/maxChargeStudy $(EXEDIR)/plotMomentum $(EXEDIR)/egammaResolution $(EXEDIR)/higgsResolution $(EXEDIR)/ootpuOnHiggsPhotons $(EXEDIR)/getAbsorberWeight

#BINS=$(EXEDIR)/egammaResolution $(EXEDIR)/higgsResolution
BINS=$(EXEDIR)/egammaResoWithTruth $(EXEDIR)/egammaResolution $(EXEDIR)/hadronResolution

SUBDIRS = ../userlib 

.PHONY: $(SUBDIRS) all
all: $(SUBDIRS) lib $(BINS)

docs: all
	doxygen Doxyfile

$(SUBDIRS):
	$(MAKE) -C $@

#$(EXEDIR)/%:  $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
#	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/plotMomentum:  $(TESTDIR)/plotMomentum.cpp $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(HEPMC_DIR)lib -lHepMC #-lHepMCfio

$(EXEDIR)/hepmcHggChargedTracks:  $(TESTDIR)/hepmcHggChargedTracks.cpp $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(HEPMC_DIR)lib -lHepMC #-lHepMCfio

$(EXEDIR)/maxChargeStudy:  $(TESTDIR)/maxChargeStudy.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/plotNabove100fC:  $(TESTDIR)/plotNabove100fC.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/plotTime:  $(TESTDIR)/plotTime.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/plotParticleDensity:  $(TESTDIR)/plotParticleDensity.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/ootpuOnHiggsPhotons:  $(TESTDIR)/ootpuOnHiggsPhotons.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/MixPUSignal:  $(TESTDIR)/MixPUSignal.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/egammaResolution:  $(TESTDIR)/egammaResolution.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/hadronResolution:  $(TESTDIR)/hadronResolution.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/egammaResoWithTruth:  $(TESTDIR)/egammaResoWithTruth.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/vbfResolution:  $(TESTDIR)/vbfResolution.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/calibration:  $(TESTDIR)/calibration.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/higgsResolution:  $(TESTDIR)/higgsResolution.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/timeResolution:  $(TESTDIR)/timeResolution.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/getAbsorberWeight:  $(TESTDIR)/getAbsorberWeight.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/mipStudy:  $(TESTDIR)/mipStudy.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/mipProfile:  $(TESTDIR)/mipProfile.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/mipSelection:  $(TESTDIR)/mipSelection.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(EXEDIR)/validation:  $(TESTDIR)/validation.cpp $(LIBDIR)/lib$(LIBNAME).so $(wildcard $(BASEDIR)/include/*.h*)
	$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME)

$(OBJDIR)/%.$(OBJ_EXT):  $(SRCDIR)/%.cc $(BASEDIR)/include/%.h*
	$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so:  $(OBJS) $(LIBDEPENDS) 
	$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so $(OBJS) $(LIBS)

lib: $(LIBDIR)/lib$(LIBNAME).so

# info:
# 	@echo "LIBS: " $(LIBS)
# 	@echo "CXXFLAGS: " $(CXXFLAGS)
# 	@echo "Source files: " $(SRCS) 
# 	@echo "Object files: " $(OBJS)
# 	@echo "Executables:  " $(TARGETS)

clean:
	rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)

userclean:
	$(MAKE) -C userlib/ clean

