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
                           const int& Nscalars,
                           const amrex::Real* problo,
                           const amrex::Real* probhi);

    void init_base_state (const amrex::Real* s0_init, 
                          const amrex::Real* p0_init,
                          const int& nlevs);

    //////////////////////


    //////////////////////
    // in maestro.f90

    void init_base_state_geometry(const int& nlevs,
                                  const amrex::Real* dx_fine,
                                  const int* ncell_fine);


    void destroy_base_state_geometry();

    //////////////////////


    //////////////////////
    // in make_S.f90
    void make_S_cc(const int* lo, const int* hi, 
                   amrex::Real* S_cc, const int* s_lo, const int* s_hi,
                   const int& ng, const int& ncomp);
    //////////////////////

    void initdata(const int& level, const int& max_levs,
                  const amrex::Real& time, 
                  const int* lo, const int* hi,
                  BL_FORT_FAB_ARG_3D(scal),
                  BL_FORT_FAB_ARG_3D(vel),
                  const amrex::Real* s0_init, 
                  const amrex::Real* p0_init,
                  const amrex::Real* dx);

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
