#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void make_analytic_solution(const int* lo, const int* hi,
                            amrex::Real* solution, const int* s_lo, const int* s_hi,
                            const amrex::Real* dx, const amrex::Real time);

void initscaldata(const int* lo, const int* hi,
              const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
              const amrex::Real* s0_init,
              amrex::Real* p0_init,
              const amrex::Real* dx);


#ifdef __cplusplus
}
#endif
#endif
