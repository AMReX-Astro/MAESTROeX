#ifndef _Maestro_F_H_
#define _Maestro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

// these are interfaces to fortran subroutines
void runtime_pretty_print(int* jobinfo_file_name,
                          const int* jobinfo_file_length);
void runtime_init();

//////////////////////
// in maestro_init.f90

void maestro_network_init();
void maestro_eos_init();
//////////////////////

//////////////////////
// in meth_params.F90

// read in shared C++/Fortran parameters from inputs file
void read_method_params();

//////////////////////

#ifdef __cplusplus
}
#endif

#endif
