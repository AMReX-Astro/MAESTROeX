NDEBUG := t
MPI    :=
OMP    :=
MKVERBOSE :=t 
COMP := gfortran

FBOXLIB_HOME := $(FBOXLIB_HOME)

# some routines need an eos/network (i.e. to compute thermodynamic 
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := t

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# include the main Makefile stuff
include $(FBOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# Microphysics core directories
MICROPHYS_CORE := 

# core BoxLib directories
FBOXLIB_CORE := Src/BaseLib

# other packages needed for data_processing
Fmdirs := 

# directories containing files that are 'include'-d via Fortran
Fmincludes := 

ifdef NEED_EOS_NETWORK
  EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS

  NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks

  EOS_DIR := helmholtz

  NETWORK_DIR := ignition_reaclib/URCA-simple

  MICROPHYS_CORE += $(MAESTRO_HOME)/Microphysics/EOS \
                    $(MAESTRO_HOME)/Microphysics/networks \
                    $(EOS_TOP_DIR) \
                    $(EOS_TOP_DIR)/$(EOS_DIR) \
                    $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \
                    $(MICROPHYSICS_HOME)/screening \
                    $(MICROPHYSICS_HOME)/neutrinos

  f90sources += probin.f90

  include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES
  include $(NETWORK_TOP_DIR)/GNetwork.mak

  Fmdirs += Microphysics/networks
  Fmdirs += Microphysics/EOS

  ifdef NEED_VODE
    Fmdirs += Util/VODE
  endif

  ifdef NEED_BLAS
    Fmdirs += Util/BLAS
  endif

endif

# Add EXTRA_THERMO to FPP_DEFINES to calculate dpdA, etc.
EXTRA_THERMO := t
FPP_DEFINES += -DEXTRA_THERMO=$(EXTRA_THERMO)

# we are not using the CUDA stuff
FPP_DEFINES += -DAMREX_DEVICE=""

# the helmeos has a table
ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
  EOS_PATH := $(EOS_TOP_DIR)/helmholtz
  ALL: table
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi

# URCA network has tables
ifeq ($(findstring URCA-simple, $(NETWORK_DIR)), URCA-simple)
  ALL: urcatables
endif

urcatables:
	@if [ ! -f 23Ne-23Na_betadecay.dat ]; then echo ${bold}Linking 23Ne-23Na_betadecay.dat${normal}; ln -s $(NETWORK_TOP_DIR)/$(NETWORK_DIR)/23Ne-23Na_betadecay.dat .;  fi
	@if [ ! -f 23Na-23Ne_electroncapture.dat ]; then echo ${bold}Linking 23Na-23Ne_electroncapture.dat${normal}; ln -s $(NETWORK_TOP_DIR)/$(NETWORK_DIR)/23Na-23Ne_electroncapture.dat .;  fi

# any MAESTRO stuff
Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))

# Microphysics
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# BoxLib stuff
Fmpack += $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(FBOXLIB_CORE), $(FBOXLIB_HOME)/$(dir))

# any include files
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

#programs := fad_excess
#programs := fconv_slopes
programs := fconv_radial

f90sources += probin.f90

PROBIN_TEMPLATE := $(MICROPHYSICS_HOME)/unit_test/dummy.probin.template
PROBIN_PARAMETER_DIRS += $(MICROPHYSICS_HOME)/unit_test/  
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)


PROBIN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(FBOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "

ALL: $(pnames)

# include the fParallel Makefile rules
include $(FBOXLIB_HOME)/Tools/F_mk/GMakerules.mak

#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(FBOXLIB_HOME)/Tools/F_scripts/makebuildinfo.py \
           --modules "$(Fmdirs) $(MICROPHYS_CORE) $(UNIT_DIR)" \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --fboxlib_home "$(FBOXLIB_HOME)" \
           --source_home "$(MICROPHYSICS_HOME)" \
           --network "$(NETWORK_DIR)" \
           --integrator "$(INTEGRATOR_DIR)" \
           --eos "$(EOS_DIR)"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90



#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean ::
	$(RM) probin.f90
	$(RM) build_info.f90


realclean ::
	$(RM) helm_table.dat
	$(RM) 23Na-23Ne_electroncapture.dat
	$(RM) 23Ne-23Na_betadecay.dat


