DEBUG      = FALSE
DIM        = 2
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_REACT  = TRUE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..


# define the physics packages to build this problem
EOS_DIR := gamma_law
CONDUCTIVITY_DIR := constant

NETWORK_DIR := general_null
NETWORK_INPUTS := triple_alpha_plus_o.net

Bpack   := ./Make.package
Blocs   := .

PROBIN_PARAMETER_DIRS := .

# include the MAESTRO build stuff
include $(MAESTROEX_HOME)/Exec/Make.Maestro
