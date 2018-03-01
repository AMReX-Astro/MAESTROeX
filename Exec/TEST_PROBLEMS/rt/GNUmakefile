PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..

# Set the EOS, conductivity, and network directories
# We first check if these exist in $(MAESTROEX_HOME)/Microphysics/(EOS/conductivity/networks)
# If not we use the version in $(MICROPHYSICS_HOME)/Microphysics/(EOS/conductivity/networks)
EOS_DIR          := gamma_law
CONDUCTIVITY_DIR := constant
NETWORK_DIR      := general_null

# for the general_null network, we read in the species names, aion, and zion from this file
GENERAL_NET_INPUTS := $(MAESTROEX_HOME)/Microphysics/networks/$(NETWORK_DIR)/simple.net

Blocs   := . 

include $(MAESTROEX_HOME)/Exec/Make.Maestro
