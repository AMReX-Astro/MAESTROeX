PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

# define the location of the MAESTRO_C home directory
MAESTRO_C_HOME  := ../../..

# Set the EOS, conductivity, and network directories
# We first check if these exist in $(MAESTRO_C_HOME)/Microphysics/(EOS/conductivity/networks)
# If not we use the version in $(MICROPHYSICS_HOME)/Microphysics/(EOS/conductivity/networks)
EOS_DIR          := gamma_law
CONDUCTIVITY_DIR := constant
NETWORK_DIR      := general_null

# for the general_null network, we read in the species names, aion, and zion from this file
GENERAL_NET_INPUTS := $(MAESTRO_C_HOME)/Microphysics/networks/$(NETWORK_DIR)/simple.net

Blocs   := . 

include $(MAESTRO_C_HOME)/Exec/Make.Maestro
