PRECISION  = DOUBLE
PROFILE    = FALSE
DEBUG      = TRUE
DIM        = 2
COMP	   = gnu
USE_MPI    = TRUE
USE_OMP    = FALSE
USE_REACT  = TRUE

# define the location of the MAESTRO_C home directory
MAESTRO_C_HOME  := ../../..

# Set the EOS, conductivity, and network directories
# We first check if these exist in $(MAESTRO_C_HOME)/Microphysics/(EOS/conductivity/networks)
# If not we use the version in $(MICROPHYSICS_HOME)/Microphysics/(EOS/conductivity/networks)
EOS_DIR          := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_DIR      := ignition_simple

Blocs   := . 

include $(MAESTRO_C_HOME)/Exec/Make.Maestro
