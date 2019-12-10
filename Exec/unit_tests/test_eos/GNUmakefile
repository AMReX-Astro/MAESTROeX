DEBUG      = FALSE
DIM        = 3
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = TRUE
USE_REACT  = TRUE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..

# if not already defined, point to Microphysics
MICROPHYSICS_HOME ?= ../../../../Microphysics

# define the physics packages to build this problem
# EOS_DIR := gamma_law_general
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTROEX_HOME)/Microphysics/networks/$(NETWORK_DIR)/H_core.net

Bpack   := ./Make.package
Blocs   := .

PROBIN_PARAMETER_DIRS := .

# include the MAESTRO build stuff
include $(MAESTROEX_HOME)/Exec/Make.Maestro
