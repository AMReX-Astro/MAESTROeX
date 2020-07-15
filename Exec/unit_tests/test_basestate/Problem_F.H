#ifndef _Problem_F_H_
#define _Problem_F_H_
#include <AMReX_BLFort.H>

/* problem-specific stuff goes here */

#ifdef __cplusplus
extern "C" {
#endif

void get_heating(amrex::Real* Hbar, const amrex::Real* rho0,
                 const amrex::Real* tempbar, const amrex::Real* rhoX0,
                 const amrex::Real time, const amrex::Real dt,
                 const amrex::Real* r_cc_loc);

void make_Sbar(amrex::Real* Sbar, const amrex::Real* rho0,
               const amrex::Real* tempbar, const amrex::Real* rhoX0,
               const amrex::Real* Hbar);

void compute_gamma1bar(amrex::Real* gamma1bar, const amrex::Real* rho0,
                       const amrex::Real* tempbar, const amrex::Real* rhoX0,
                       const amrex::Real* p0);

void update_temp(const amrex::Real* rho0, amrex::Real* tempbar,
                 const amrex::Real* rhoX0, const amrex::Real* p0);

void update_species(const amrex::Real* rho0,
                    const amrex::Real* rho0_predicted_edge,
                    const amrex::Real* rhoX0_old, amrex::Real* rhoX0_new,
                    const amrex::Real* w0, const amrex::Real* r_cc_loc,
                    const amrex::Real* r_edge_loc, const amrex::Real dt);

void set_initial_dt(amrex::Real* dt);

void make_rhoX0(const amrex::Real* s0, amrex::Real* rhoX0);

void check_hseness(const amrex::Real* rho0, const amrex::Real* p0,
                   amrex::Real* max_hse_error);

#ifdef __cplusplus
}
#endif
#endif
