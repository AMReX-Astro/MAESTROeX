#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void get_min_timestep(amrex::Real* min_timestep);


#ifdef __cplusplus
}
#endif
#endif
