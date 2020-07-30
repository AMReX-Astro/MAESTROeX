#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void do_tests(const int* lo, const int* hi, amrex::Real* scal, const int* s_lo,
              const int* s_hi, amrex::Real* error, const int* e_lo,
              const int* e_hi, const int* nc_e, const int* dom_lo,
              const int* dom_hi);

#ifdef __cplusplus
}
#endif
#endif
