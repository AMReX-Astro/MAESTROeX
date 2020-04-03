PRECISION = DOUBLE
PROFILE = FALSE
DEBUG = TRUE
DIM = 3

COMP = gnu

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE

USE_ACC = FALSE

# programs to be compiled
ALL: radial_$(DIM)d.ex

Bpack   := ./Make.package
Blocs   := .
# EXTERN_SEARCH = .

MAESTROEX_HOME := ../..

EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MICROPHYSICS_HOME)/networks/general_null/simple.net

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/amrdata
include $(AMREX_HOME)/Src/Extern/amrdata/Make.package
vpathdir += $(AMREX_HOME)/Src/Extern/amrdata

include $(MAESTROEX_HOME)/Exec/Make.Maestro

radial_$(DIM)d.ex: $(objForExecs)
	@echo Linking $@ ...
	$(SILENT) $(PRELINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(libraries)
