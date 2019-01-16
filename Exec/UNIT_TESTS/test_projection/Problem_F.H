#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void init_vel(const int* lo, const int* hi,
              amrex::Real* vel, const int* v_lo, const int* v_hi,
              const int* nc_v, const amrex::Real* dx);

void get_project_type(int* project_type);

void add_grad_scalar(const int* lo, const int* hi,
                     amrex::Real* gphi, const int* g_lo, const int* g_hi,
                     const int* nc_g,
                     amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const int* nc_v,
                     const int* domain_phys_bc, const amrex::Real* dx);

void add_grad_scalar_mac(const int* lo, const int* hi,
                         amrex::Real* gphi, const int* g_lo, const int* g_hi,
                         const int* nc_g,
                         amrex::Real* vel, const int* v_lo, const int* v_hi,
                         const int* nc_v,
                         const int* domain_phys_bc, const int* box_phys_bc,
                         const amrex::Real* dx);

void convert_MAC_to_cc(const int* lo, const int* hi,
                       const amrex::Real* umac, const int* m_lo, const int* m_hi,
                       const int* nc_m,
                       amrex::Real* u, const int* u_lo, const int* u_hi,
                       const int* nc_u);

#ifdef __cplusplus
}
#endif
#endif
