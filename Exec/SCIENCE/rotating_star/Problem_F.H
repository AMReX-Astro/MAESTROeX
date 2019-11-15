#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

void set_enuc_omegadot(const amrex::Real* rho0, 
    amrex::Real* rho_omegadot, amrex::Real* rho_Hnuc);

void put_1d_spec_array_on_cart(const int* lo, const int* hi, const int lev,
                              amrex::Real* s0_cart,
                              const int* s0_cart_lo, const int* s0_cart_hi, 
                              const amrex::Real* s0, const int is_input_edge_centered);

void put_1d_spec_array_on_cart_sphr(const int* lo, const int* hi,
                amrex::Real* s0_cart,
                const int* s0_cart_lo, const int* s0_cart_hi, 
                const amrex::Real* s0, const amrex::Real* dx,
                const int is_input_edge_centered,
                const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
                const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);


#ifdef __cplusplus
}
#endif

#endif