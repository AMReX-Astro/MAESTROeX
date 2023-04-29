DEBUG      = FALSE
DIM        = 3
COMP	   = gnu
USE_MPI    = FALSE
USE_OMP    = FALSE
USE_REACT  = TRUE

# define the location of the MAESTROEX home directory
MAESTROEX_HOME  := ../../..


# define the physics packages to build this problem
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_DIR := aprox13

Bpack   := ./Make.package
Blocs   := .

PROBIN_PARAMETER_DIRS := .

# include the MAESTRO build stuff
include $(MAESTROEX_HOME)/Exec/Make.Maestro
