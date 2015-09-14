# $Id: GNUmakefile,v 1.2 2000-10-19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := PFCalEE
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

CPPFLAGS += -Iuserlib/include/ $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --glibs) -Luserlib/lib -lPFCalEEuserlib

.PHONY: $(SUBDIRS) all
all: $(SUBDIRS) lib bin

#use this line to also use pythia6
#all : $(SUBDIRS) pythia6 lib bin

include $(G4INSTALL)/config/binmake.gmk

 # if you want to use Pythia library, uncomment out the next lines.
#
#  G4LIB_USE_PYTHIA := 1
#ifdef G4LIB_USE_PYTHIA
#    CPPFLAGS += -DG4LIB_USE_PYTHIA
#endif

INCFLAGS += -I$(HEPMC_DIR)/include/

ifdef G4LIB_USE_PYTHIA
      LDLIBS1 += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio -L$(G4TMPDIR) -lPythia6 -lg2c
else
      LDLIBS1 += -L$(HEPMC_DIR)/lib -lHepMC -lHepMCfio #$(G4TMPDIR)/HEPEvtcom.o
endif

SUBDIRS = userlib

#pythia stuff
#FCFLAGS += -c
#pythia6: $(G4TMPDIR)/libPythia6.so
#$(G4TMPDIR)/libPythia6.so: $(G4TMPDIR)/pythia6.o
#        $(FC) -shared -Wl,-soname,libPythia6.so -o $(G4TMPDIR)/libPythia6.so $(G4TMPDIR)/pythia6.o
#$(G4TMPDIR)/pythia6.o:
#        $(FC) $(FCFLAGS) $(PYTHIA6)/pythia-$(PYTHIA6_VERSION).f -o $(G4TMPDIR)/pythia6.o


$(SUBDIRS):
	$(MAKE) -C $@

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

userclean:
	$(MAKE) -C userlib/ clean
