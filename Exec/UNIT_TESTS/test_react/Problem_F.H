#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void get_min_timestep(amrex::Real* min_timestep);
void get_react_its(int* react_its);

void initdata_thermal(const int* lo, const int* hi,
              const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
              const int* domlo, const int* domhi,
              const amrex::Real* s0_init,
              const amrex::Real* p0_init,
              const amrex::Real* dx);


#ifdef __cplusplus
}
#endif
#endif
