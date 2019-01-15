DEBUG      = FALSE
DIM        = 3
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_REACT  = TRUE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..

# if not already defined, point to Microphysics
MICROPHYSICS_HOME ?= ../../../../Microphysics

# define the physics packages to build this problem
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar

#NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
#NETWORK_DIR := rprox
# NETWORK_DIR := triple_alpha
#NETWORK_DIR := ignition_chamulak
NETWORK_DIR := aprox13
# NETWORK_DIR := ignition_simple

# a special case for the test suite
ifeq ($(NETWORK_DIR), ignition_chamulak)
  NETWORK_TOP_DIR := $(MAESTROEX_HOME)/Microphysics/networks
endif

Bpack   := ./Make.package
Blocs   := .

PROBIN_PARAMETER_DIRS := .

# include the MAESTRO build stuff
include $(MAESTROEX_HOME)/Exec/Make.Maestro
