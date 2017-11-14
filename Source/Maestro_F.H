#ifndef _Maestro_F_H_
#define _Maestro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    //////////////////////
    // in meth_params.F90

    // read in shared C++/Fortran parameters from inputs file
    void read_method_params();

    // deallocate strings
    void finalize_meth_params();

    //////////////////////

    //////////////////////
    // in maestro.f90

    void maestro_network_init();
    void get_num_spec(int* nspec);
    void get_spec_names(int* spec_names, int* ispec, int* len);

    // set fortran-specific parameters in meth_params.F90
    void set_method_params(const int& Density,
                           const int& Enthalpy,
                           const int& FirstSpec,
                           const int& Temperature,
                           const int& Pressure,
                           const int& Nscalars);

    void init_base_state (const amrex::Real* s0_init, 
                          const amrex::Real* p0_init,
                          const int& nlevs,
                          const amrex::Real* prob_lo);

    //////////////////////


    //////////////////////
    // in maestro.f90

    void init_base_state_geometry(const int& max_radial_level_in,
                                  const int& nr_fine_in,
                                  const amrex::Real& dr_fine_in,
                                  const amrex::Real* r_cc_loc,
                                  const amrex::Real* r_edge_loc,
                                  const amrex::Real* dx_fine,
                                  const int* domhi_fine,
                                  const amrex::Real* prob_lo,
                                  const amrex::Real* prob_hi);


    void destroy_base_state_geometry();

    //////////////////////


    //////////////////////
    // in make_S.f90
    void make_S_cc(const int* lo, const int* hi, 
                   amrex::Real* S_cc,  const int* s_lo, const int* s_hi, const int& nc_s,
                   const amrex::Real* scal,  const int* c_lo, const int* c_hi, const int& nc_c,
                   const amrex::Real* rodot, const int* r_lo, const int* r_hi, const int& nc_r,
                   const amrex::Real* rHnuc, const int* n_lo, const int* n_hi, const int& nc_n,
                   const amrex::Real* rHext, const int* e_lo, const int* e_hi, const int& nc_e,
                   const amrex::Real* therm, const int* t_lo, const int* t_hi, const int& nc_t);
    //////////////////////

    void initdata(const int& lev,
                  const amrex::Real& time, 
                  const int* lo, const int* hi,
                  BL_FORT_FAB_ARG_3D(scal),
                  BL_FORT_FAB_ARG_3D(vel),
                  const amrex::Real* s0_init, 
                  const amrex::Real* p0_init,
                  const amrex::Real* dx,
                  const amrex::Real* prob_lo);

    void get_face_velocity(const int& level, const amrex::Real& time, 
                           AMREX_D_DECL(BL_FORT_FAB_ARG(xvel),
                                        BL_FORT_FAB_ARG(yvel),
                                        BL_FORT_FAB_ARG(zvel)),
                           const amrex::Real* dx, const amrex::Real* problo);

    //////////////////////
    // in tagging.f90

    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
                     const BL_FORT_FAB_ARG_3D(state),
                     const int* tagval, const int* clearval,
                     const int* lo, const int* hi,
                     const amrex::Real* dx, const amrex::Real* problo,
                     const amrex::Real* time, const amrex::Real* phierr);

    //////////////////////

    void advect(const int* lo, const int*hi,
                const BL_FORT_FAB_ARG_3D(statein),
                BL_FORT_FAB_ARG_3D(stateout),
                AMREX_D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
                             const BL_FORT_FAB_ARG_3D(yvel),
                             const BL_FORT_FAB_ARG_3D(zvel)),
                AMREX_D_DECL(BL_FORT_FAB_ARG_3D(fx),
                             BL_FORT_FAB_ARG_3D(fy),
                             BL_FORT_FAB_ARG_3D(fz)),
                const amrex::Real* dx, const amrex::Real& dt, const int& ncomp);

    //////////////////////
    // in bc_fill_nd.F90

    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const amrex::Real* dx, const amrex::Real* glo, 
		  const amrex::Real* time, const int* bc);

    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
