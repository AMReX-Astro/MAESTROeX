DEBUG      = FALSE
DIM        = 1
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = TRUE
USE_REACT  = TRUE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..

# if not already defined, point to Microphysics
MICROPHYSICS_HOME ?= ../../../../Microphysics

# Set the EOS, conductivity, and network directories
# We first check if these exist in $(MAESTROEX_HOME)/Microphysics/(EOS/conductivity/networks)
# If not we use the version in $(MICROPHYSICS_HOME)/Microphysics/(EOS/conductivity/networks)
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := constant
NETWORK_DIR := general_null
NETWORK_FILE := triple_alpha_plus_o.net
GENERAL_NET_INPUTS := $(MAESTROEX_HOME)/Microphysics/networks/$(NETWORK_DIR)/$(NETWORK_FILE)

Bpack   := ./Make.package
Blocs   := .

PROBIN_PARAMETER_DIRS := .

# include the MAESTRO build stuff
include $(MAESTROEX_HOME)/Exec/Make.Maestro
