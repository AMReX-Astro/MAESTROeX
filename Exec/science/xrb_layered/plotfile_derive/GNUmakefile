PRECISION = DOUBLE
PROFILE = FALSE

DEBUG = FALSE

DIM = 2

COMP = g++
FCOMP = gfortran

USE_MPI = FALSE
USE_OMP = FALSE

USE_REACT = TRUE
USE_CXX_EOS = TRUE

MAX_ZONES := 16384

DEFINES += -DNPTS_MODEL=$(MAX_ZONES)

# programs to be compiled
EBASE := derive

# EOS and network
EOS_DIR := helmholtz
NETWORK_DIR := CNO_extras

# Conductivity
CONDUCTIVITY_HOME  = $(MICROPHYSICS_HOME)/conductivity
USE_CONDUCTIVITY := TRUE
CONDUCTIVITY_DIR := stellar

Bpack := ./Make.package
Blocs := . ..

EXTERN_SEARCH += . ..

USE_AMR_CORE = TRUE

include $(MICROPHYSICS_HOME)/Make.Microphysics
