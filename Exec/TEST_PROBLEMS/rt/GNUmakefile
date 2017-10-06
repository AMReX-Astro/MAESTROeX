PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE

# define the location of the MAESTRO_C top directory
MAESTRO_C_HOME  := ../../..

# This sets the EOS directory in $(MAESTRO_C_HOME)/Microphysics/EOS
EOS_DIR     := gamma_law

# This sets the conductivity directory in $(MAESTRO_C_HOME)/Microphysics/conductivity
CONDUCTIVITY_DIR := constant

# This sets the network directory in $(MAESTRO_C_HOME)/Microphysics/Networks
NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTRO_C_HOME)/Microphysics/networks/$(NETWORK_DIR)/simple.net

Bpack   := ./Make.package
Blocs   := . 

include $(MAESTRO_C_HOME)/Exec/Make.Maestro
