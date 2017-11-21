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
                          const amrex::Real* rho0,
                          const amrex::Real* rhoh0,
                          const amrex::Real* p0,
                          const amrex::Real* tempbar,
                          const int& nlevs,
                          const amrex::Real* prob_lo);

    //////////////////////


    //////////////////////
    // in base_state_geometry.f90

    void init_base_state_geometry(const int& max_radial_level_in,
                                  const int& nr_fine_in,
                                  const amrex::Real& dr_fine_in,
                                  const amrex::Real* r_cc_loc,
                                  const amrex::Real* r_edge_loc,
                                  const amrex::Real* dx_fine,
                                  const int* domhi_fine,
                                  const amrex::Real* prob_lo,
                                  const amrex::Real* prob_hi);

    void compute_cutoff_coords(const amrex::Real* rho0);

    void init_multilevel(const int& finest_radial_level_in);

    void restrict_base(amrex::Real* s0, const int& is_cell_centered);

    void fill_ghost_base(amrex::Real* s0, const int& is_cell_centered);

    void destroy_base_state_geometry();

    //////////////////////

    //////////////////////
    // in make_gamma.F90
    void make_gamma(const int& lev, const int* lo, const int* hi, 
                    amrex::Real* gamma,  const int* g_lo, const int* g_hi, const int& nc_g,
                    const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int& nc_s,
                    const amrex::Real* p0);
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

    //////////////////////
    // in burner_loop.f90
    void burner_loop(const int& lev, const int* lo, const int* hi, 
                     const amrex::Real* s_in,  const int* i_lo, const int* i_hi, const int& nc_i,
                     amrex::Real* s_out,  const int* o_lo, const int* o_hi, const int& nc_o,
                     const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi, const int& nc_e,
                     amrex::Real* rho_odot, const int* r_lo, const int* r_hi, const int& nc_r,
                     amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi, const int& nc_n,
                     const amrex::Real* tempbar_init_in, const amrex::Real& dt_in);
    //////////////////////

    //////////////////////
    // in average.f90
    void average(const int& lev, const int* lo, const int* hi,
                 const amrex::Real* phi, const int* p_lo, const int* p_hi, 
                 amrex::Real* phisum);

    void divide_phisum_by_ncell(amrex::Real* phisum, int* ncell);
    //////////////////////

    //////////////////////
    // in make_grav.f90
    void make_grav_cell(amrex::Real* grav_cell,
                        const amrex::Real* rho0,
                        const amrex::Real* r_cc_loc,
                        const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in make_div_coeff.f90
    void make_div_coeff(amrex::Real* div_coeff,
                        const amrex::Real* rho0,
                        const amrex::Real* p0,
                        const amrex::Real* gamma1bar_in,
                        const amrex::Real* grav_cell_in);
    //////////////////////


    //////////////////////
    // in initdata.f90
    void initdata(const int& lev,
                  const amrex::Real& time, 
                  const int* lo, const int* hi,
                  BL_FORT_FAB_ARG_3D(scal),
                  BL_FORT_FAB_ARG_3D(vel),
                  const amrex::Real* s0_init, 
                  const amrex::Real* p0_init,
                  const amrex::Real* dx,
                  const amrex::Real* prob_lo);
    //////////////////////

    //////////////////////
    // in tagging.f90
    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
                     const BL_FORT_FAB_ARG_3D(state),
                     const int* tagval, const int* clearval,
                     const int* lo, const int* hi,
                     const amrex::Real* dx, const amrex::Real* problo,
                     const amrex::Real* time, const amrex::Real* temperr);
    //////////////////////

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
