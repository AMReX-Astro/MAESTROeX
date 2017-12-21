#ifndef _Maestro_F_H_
#define _Maestro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" 
{
#endif

    // these are interfaces to fortran subroutines

    //////////////////////
    // in average.f90
    void average(const int& lev, const int* lo, const int* hi,
                 const amrex::Real* phi, const int* p_lo, const int* p_hi, 
                 amrex::Real* phisum);

    void divide_phisum_by_ncell(amrex::Real* phisum, int* ncell);
    //////////////////////

    //////////////////////
    // in base_state_geometry.f90

    void init_base_state_geometry(const int& max_radial_level_in,
                                  const int& nr_fine_in,
                                  const amrex::Real& dr_fine_in,
                                  const amrex::Real* r_cc_loc,
                                  const amrex::Real* r_edge_loc,
                                  const amrex::Real* dx_fine,
                                  const int* domhi_fine);

    void compute_cutoff_coords(const amrex::Real* rho0);

    void init_multilevel(const int& finest_radial_level_in);

    void restrict_base(amrex::Real* s0, const int& is_cell_centered);

    void fill_ghost_base(amrex::Real* s0, const int& is_cell_centered);

    void destroy_base_state_geometry();

    //////////////////////

    //////////////////////
    // in bc_fill_nd.F90

    void phifill (amrex::Real* phi, const int* phi_lo, const int* phi_hi,
		  const int* domlo, const int* domhi,
		  const amrex::Real* dx, const amrex::Real* gridlo, 
		  const amrex::Real* time, const int* bc);

    //////////////////////

    //////////////////////
    // in burner.F90
    void burner_init();
    //////////////////////

    //////////////////////
    // in burner_loop.f90
    void burner_loop(const int& lev, const int* lo, const int* hi, 
                     const amrex::Real* s_in,     const int* i_lo, const int* i_hi, const int* nc_i,
                           amrex::Real* s_out,    const int* o_lo, const int* o_hi, const int* nc_o,
                     const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi,
                           amrex::Real* rho_odot, const int* r_lo, const int* r_hi, const int* nc_r,
                           amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi,
                     const amrex::Real* tempbar_init_in,
                     const amrex::Real& dt_in);
    //////////////////////

    //////////////////////
    // in compute_dt.F90
    void firstdt(amrex::Real& dt, amrex::Real& umax, const int* lo, const int* hi, 
                 const amrex::Real* dx,
                 const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
                 const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u);
    //////////////////////

    //////////////////////
    // in enforce_HSE.f90
    void enforce_HSE(const amrex::Real* rho0,
                           amrex::Real* p0,
                     const amrex::Real* grav_cell,
                     const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in fill_umac_ghost.f90
    void fill_umac_ghost(const int* domlo, const int* domhi,
                         const int* lo, const int* hi,
                         amrex::Real* umac, const int* umac_lo, const int* umac_hi,
                         amrex::Real* vmac, const int* vmac_lo, const int* vmac_hi,
#if (AMREX_SPACEDIM == 3)
                         amrex::Real* wmac, const int* wmac_lo, const int* wmac_hi,
#endif
                         const int* lo_bc, const int* hi_bc);
    //////////////////////

    //////////////////////
    // in fill_3d_data.f90
    void put_1d_array_on_cart(const int& lev, const int* lo, const int* hi,
                              amrex::Real* s0_cart,
                              const int* s0_cart_lo, const int* s0_cart_hi, const int* nc_s,
                              const amrex::Real* s0, const int& is_input_edge_centered,
                              const int& is_output_a_vector);
    //////////////////////

    //////////////////////
    // in initdata.f90
    void initdata(const int& lev,
                  const amrex::Real& time, 
                  const int* lo, const int* hi,
                  const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
                  const amrex::Real* vel,  const int* vel_lo,  const int* vel_hi,  const int* nc_v,
                  const amrex::Real* s0_init, 
                  const amrex::Real* p0_init,
                  const amrex::Real* dx);
    //////////////////////

    //////////////////////
    // in maestro.f90

    void maestro_network_init();
    void get_num_spec(int* nspec);
    void get_spec_names(int* spec_names, int* ispec, int* len);

    // set fortran-specific parameters in meth_params.F90
    void set_method_params (const int& Density,
                            const int& Enthalpy,
                            const int& FirstSpec,
                            const int& Temperature,
                            const int& Pressure,
                            const int& Nscalars,
                            const amrex::Real* prob_lo_in,
                            const amrex::Real* prob_hi_in);

    void set_rel_eps (const amrex::Real& rel_eps_in);

    void init_base_state (const amrex::Real* s0_init, 
                          const amrex::Real* p0_init,
                          const amrex::Real* rho0,
                          const amrex::Real* rhoh0,
                          const amrex::Real* p0,
                          const amrex::Real* tempbar,
                          const int& nlevs);

    //////////////////////

    //////////////////////
    // in make_beta0.f90
    void make_beta0(      amrex::Real* beta0,
                    const amrex::Real* rho0,
                    const amrex::Real* p0,
                    const amrex::Real* gamma1bar,
                    const amrex::Real* grav_cell);
    //////////////////////

    //////////////////////
    // in make_gamma.F90
    void make_gamma(const int& lev, const int* lo, const int* hi, 
                          amrex::Real* gamma, const int* g_lo, const int* g_hi,
                    const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
                    const amrex::Real* p0);
    //////////////////////

    //////////////////////
    // in make_grav.f90
    void make_grav_cell(      amrex::Real* grav_cell,
                        const amrex::Real* rho0,
                        const amrex::Real* r_cc_loc,
                        const amrex::Real* r_edge_loc);

    void make_grav_edge(      amrex::Real* grav_edge,
                        const amrex::Real* rho0,
                        const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in make_S.f90
    void make_S_cc(const int* lo, const int* hi, 
                         amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
                   const amrex::Real* scal,  const int* c_lo, const int* c_hi, const int* nc_c,
                   const amrex::Real* rodot, const int* r_lo, const int* r_hi, const int* nc_r,
                   const amrex::Real* rHnuc, const int* n_lo, const int* n_hi,
                   const amrex::Real* rHext, const int* e_lo, const int* e_hi,
                   const amrex::Real* therm, const int* t_lo, const int* t_hi);


    void make_ccrhs(const int& lev, const int* lo, const int* hi, 
                          amrex::Real* ccrhs, const int* c_lo, const int* c_hi,
                    const amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
                    const amrex::Real* Sbar,
                    const amrex::Real* beta0);


    void make_nodalrhs(const int& lev, const int* lo, const int* hi, 
                       amrex::Real* nodalrhs,  const int* n_lo, const int* n_hi,
                       const amrex::Real* ccrhs,  const int* c_lo, const int* c_hi);

    //////////////////////

    //////////////////////
    // in make_vel_force.f90
    void make_vel_force(const int& lev, const int* lo, const int* hi,
                        amrex::Real* vel_force, const int* f_lo, const int* f_hi, const int* nc_f,
                        const amrex::Real* gpi, const int* g_lo, const int* g_hi, const int* nc_g,
                        const amrex::Real* rho, const int* r_lo, const int* r_hi,
                        const amrex::Real* uedge, const int* u_lo, const int* u_hi,
                        const amrex::Real* vedge, const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
                        const amrex::Real* wedge, const int* w_lo, const int* w_hi,
#endif
                        const amrex::Real* w0,
                        const amrex::Real* w0_force,
                        const amrex::Real* rho0,
                        const amrex::Real* grav,
                        const int& is_final_update,
                        const int& do_add_utilde_force);
    //////////////////////

    //////////////////////
    // in make_w0.f90
    void make_w0(amrex::Real* w0, const amrex::Real* w0_old,
                 const amrex::Real* w0_force, const amrex::Real* Sbar_in, 
                 const amrex::Real* rho0_old, const amrex::Real* rho0_new,
                 const amrex::Real* p0_old, const amrex::Real* p0_new,
                 const amrex::Real* gamma1bar_old, const amrex::Real* gamma1bar_new,
                 const amrex::Real* p0_minus_peosbar, const amrex::Real* psi,
                 const amrex::Real* etarho_ec, const amrex::Real* etarho_cc,
                 const amrex::Real* delta_chi_w0, const amrex::Real* r_cc_loc,
                 const amrex::Real* r_edge_loc, const amrex::Real& dt,
                 const amrex::Real& dtold, const int& is_predictor);
    //////////////////////

    //////////////////////
    // in meth_params.F90

    // read in shared C++/Fortran parameters from inputs file
    void read_method_params();

    // deallocate strings
    void finalize_meth_params();

    //////////////////////

    //////////////////////
    // in rhoh_vs_t.F90
    void makeTfromRhoH (const int& lev, const int* lo, const int* hi, 
                        amrex::Real* state,  const int* s_lo, const int* s_hi, const int* nc_s,
                        const amrex::Real* p0);

    void makeTfromRhoP (const int& lev, const int* lo, const int* hi, 
                        amrex::Real* state,  const int* s_lo, const int* s_hi, const int* nc_s,
                        const amrex::Real* p0, const int& updateRhoH);
    //////////////////////

    //////////////////////
    // in sponge.F90
    void init_sponge(const amrex::Real* rho0);

    void mk_sponge (const int* lo, const int* hi, 
                    amrex::Real* sponge,  const int* s_lo, const int* s_hi,
                    const amrex::Real* dx, const amrex::Real& dt);
    //////////////////////


    //////////////////////
    // in tagging.f90
    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
                     const amrex::Real* state, const int* state_lo, const int* state_hi,
                     const int* tagval, const int* clearval,
                     const int* lo, const int* hi,
                     const amrex::Real* dx,
                     const amrex::Real* time, const amrex::Real* temperr);
    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
