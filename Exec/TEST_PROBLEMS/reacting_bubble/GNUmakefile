DEBUG      = FALSE
DIM        = 3
#COMP	   = gnu
#USE_MPI    = TRUE
USE_OMP    = FALSE
USE_REACT  = FALSE
# for GPU compiling
COMP     = PGI
USE_AMREX_CUDA = TRUE
USE_GPU_PRAGMA = TRUE
CUDA_ARCH=60
COMPILE_CUDA_PATH=/usr/local/cuda-9.2
USE_MPI  = FALSE

TINY_PROFILE = FALSE
PROFILE      = FALSE # TRUE overrides TINY_PROFILE


# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..

# if not already defined, point to Microphysics
MICROPHYSICS_HOME ?= ../../../../Microphysics

# Set the EOS, conductivity, and network directories
# We first check if these exist in $(MAESTROEX_HOME)/Microphysics/(EOS/conductivity/networks)
# If not we use the version in $(MICROPHYSICS_HOME)/Microphysics/(EOS/conductivity/networks)
EOS_DIR          := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_DIR      := ignition_simple

Bpack   := ./Make.package
Blocs   := . 

include $(MAESTROEX_HOME)/Exec/Make.Maestro
