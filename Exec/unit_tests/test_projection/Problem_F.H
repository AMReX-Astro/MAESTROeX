#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void init_vel(const int* lo, const int* hi, amrex::Real* vel, const int* v_lo,
              const int* v_hi, const amrex::Real* dx);

void init_mac_vel(const int* lo, const int* hi, amrex::Real* umac,
                  const int* u_lo, const int* u_hi, amrex::Real* vmac,
                  const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
                  amrex::Real* wmac, const int* w_lo, const int* w_hi,
#endif
                  const amrex ::Real* dx);

void get_project_type(int* project_type);

void add_grad_scalar(const int* lo, const int* hi, amrex::Real* gphi,
                     const int* g_lo, const int* g_hi, const int* nc_g,
                     amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const int* nc_v, const int* domain_phys_bc,
                     const amrex::Real* dx);

void add_grad_scalar_mac(
    const int* lo, const int* hi, amrex::Real* gphi_x, const int* gx_lo,
    const int* gx_hi, amrex::Real* gphi_y, const int* gy_lo, const int* gy_hi,
#if (AMREX_SPACEDIM == 3)
    amrex::Real* gphi_z, const int* gz_lo, const int* gz_hi,
#endif
    amrex::Real* u, const int* u_lo, const int* u_hi, amrex::Real* v,
    const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
    amrex::Real* w, const int* w_lo, const int* w_hi,
#endif
    const int* domain_phys_bc, const int* box_phys_bc, const amrex::Real* dx);

void convert_MAC_to_cc(const int* lo, const int* hi, const amrex::Real* umac,
                       const int* m_lo, const int* m_hi,
                       const amrex::Real* vmac, const int* v_lo,
                       const int* v_hi,
#if (AMREX_SPACEDIM == 3)
                       const amrex::Real* wmac, const int* w_lo,
                       const int* w_hi,
#endif
                       amrex::Real* u, const int* u_lo, const int* u_hi);

#ifdef __cplusplus
}
#endif
#endif
