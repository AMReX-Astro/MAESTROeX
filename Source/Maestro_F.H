#ifndef _Maestro_F_H_
#define _Maestro_F_H_
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    // these are interfaces to fortran subroutines

    void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

    //////////////////////
    // in advect_base.F90
    void advect_base_dens(const amrex::Real* w0,
			  const amrex::Real* rho0_old,
			        amrex::Real* rho0_new,
			        amrex::Real* rho0_predicted_edge,
			  const amrex::Real dt,
			  const amrex::Real* r_cc_loc,
			  const amrex::Real* r_edge_loc);

    void advect_base_enthalpy(const amrex::Real* w0,
			      const amrex::Real* rho0_old,
			      const amrex::Real* rhoh0_old,
			            amrex::Real* rhoh0_new,
			      const amrex::Real* rho0_predicted_edge,
			      const amrex::Real* psi,
			      const amrex::Real dt,
			      const amrex::Real* r_cc_loc,
			      const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in average.f90
    void average(const int* lev, const int* lo, const int* hi,
                 const amrex::Real* phi, const int* p_lo, const int* p_hi,
                 amrex::Real* phisum);

    void divide_phisum_by_ncell(amrex::Real* phisum, int* ncell);

    void average_sphr_irreg(const int* lev, const int* lo, const int* hi,
			    const amrex::Real* phi, const int* p_lo, const int* p_hi,
			    amrex::Real* phisum, int* ncell,
			    const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);

    void divide_phisum_by_ncell_irreg(amrex::Real* phisum, int* ncell);

    void average_sphr(amrex::Real* phisum, amrex::Real* phibar,
		      int* ncell, amrex::Real* radii, const int* finest_level);

    void sum_phi_3d_sphr(const int* lev, const int* lo, const int* hi,
			 const amrex::Real* phi, const int* p_lo, const int* p_hi,
			 amrex::Real* phisum,
			 const amrex::Real* radii, const int* finest_level,
			 const amrex::Real* dx, int* ncell,
                         const int* mask, const int* m_lo, const int* m_hi,
                         const int* use_mask);

    void compute_radii_sphr(const int* lev, amrex::Real* radii, const int* finest_level,
			    const amrex::Real* dx);
    //////////////////////

    //////////////////////
    // in base_state_geometry.f90

    void init_base_state_geometry(const int* max_radial_level_in,
                                  const int* nr_fine_in,
                                  const amrex::Real* dr_fine_in,
                                  amrex::Real* r_cc_loc,
                                  amrex::Real* r_edge_loc,
                                  const amrex::Real* dx_fine,
				  int* nr_irreg_in);

    void init_base_state_map_sphr(amrex::Real* cc_to_r, const int* lo, const int* hi,
				  const amrex::Real* dx_fine,
				  const amrex::Real* dx_lev);

    void compute_cutoff_coords(const amrex::Real* rho0);

    void init_multilevel(const int* tag_array,
			 const int* finest_radial_level_in);

    void destroy_base_state_geometry();

    //////////////////////

    //////////////////////
    // in bc_fill_nd.F90

    void scalarfill (amrex::Real* phi, const int* phi_lo, const int* phi_hi,
		  const int* domlo, const int* domhi,
		  const amrex::Real* dx, const amrex::Real* gridlo,
		  const amrex::Real* time, const int* bc, const int icomp);

    void velfill (amrex::Real* vel, const int* vel_lo, const int* vel_hi,
    	  const int* domlo, const int* domhi,
    	  const amrex::Real* dx, const amrex::Real* gridlo,
    	  const amrex::Real* time, const int* bc, const int icomp);

    //////////////////////

    //////////////////////
    // in burner.F90
    void burner_init();
    //////////////////////

    //////////////////////
    // in burner_loop.f90
    void burner_loop(const int* lo, const int* hi,
                     const int lev,
                     const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                           amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                     const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi,
                           amrex::Real* rho_odot, const int* r_lo, const int* r_hi,
                           amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi,
                     const amrex::Real* tempbar_init_in,
                     const amrex::Real dt_in,
                     const int* mask, const int* m_lo, const int* m_hi,
                     const int use_mask);

    void burner_loop_sphr(const int* lo, const int* hi,
                          const amrex::Real* s_in,     const int* i_lo, const int* i_hi,
                                amrex::Real* s_out,    const int* o_lo, const int* o_hi,
                          const amrex::Real* rho_Hext, const int* e_lo, const int* e_hi,
                                amrex::Real* rho_odot, const int* r_lo, const int* r_hi,
                                amrex::Real* rho_Hnuc, const int* n_lo, const int* n_hi,
                          const amrex::Real* tempbar_init_cart, const int* t_lo, const int* t_hi,
                          const amrex::Real dt_in,
                          const int* mask, const int* m_lo, const int* m_hi,
                          const int use_mask);

    //////////////////////

    //////////////////////
    // in cell_to_edge.f90
    void cell_to_edge(const amrex::Real* s0_cell, amrex::Real* s0_edge);

    //////////////////////

    //////////////////////
    // in compute_dt.F90
    void estdt(const int lev, amrex::Real* dt, amrex::Real* umax,
               const int* lo, const int* hi,
               const amrex::Real* dx,
               const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
               const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
               const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
               const amrex::Real* divu,  const int* d_lo, const int* d_hi,
               const amrex::Real* dSdt,  const int* t_lo, const int* t_hi,
               const amrex::Real* w0, const amrex::Real* p0, const amrex::Real* gamma1bar);

    void estdt_sphr(amrex::Real* dt, amrex::Real* umax,
		    const int* lo, const int* hi, const amrex::Real* dx,
		    const amrex::Real* scal, const int* s_lo, const int* s_hi, const int* nc_s,
		    const amrex::Real* u, const int* u_lo, const int* u_hi, const int* nc_u,
		    const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
		    const amrex::Real* divu, const int* d_lo, const int* d_hi,
		    const amrex::Real* dSdt, const int* t_lo, const int* t_hi,
		    const amrex::Real* w0,
		    const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
		    const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
		    const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
		    const amrex::Real* gp0_cart, const int* g_lo, const int* g_hi);

    void firstdt(const int lev, amrex::Real* dt, amrex::Real* umax,
                 const int* lo, const int* hi,
                 const amrex::Real* dx,
                 const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
                 const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
                 const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
                 const amrex::Real* divu,  const int* d_lo, const int* d_hi,
                 const amrex::Real* p0, const amrex::Real* gamma1bar);

    void firstdt_sphr(amrex::Real* dt, amrex::Real* umax,
		      const int* lo, const int* hi,
		      const amrex::Real* dx,
		      const amrex::Real* scal,  const int* s_lo, const int* s_hi, const int* nc_s,
		      const amrex::Real* u,     const int* u_lo, const int* u_hi, const int* nc_u,
		      const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f,
		      const amrex::Real* divu,  const int* d_lo, const int* d_hi,
		      const amrex::Real* gp0_cart, const int* g_lo, const int* g_hi);

    void estdt_divu(amrex::Real* gp0,
		      const amrex::Real* p0, const amrex::Real* gamma1bar,
		      const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in compute_grad_phi.F90
    void compute_grad_phi(const int* lo, const int* hi,
                          const amrex::Real* phi,  const int* p_lo, const int* p_hi,
                                amrex::Real* gphi, const int* g_lo, const int* g_hi,
                          const amrex::Real* dx);

    void compute_grad_phi_rad(const amrex::Real* phi,
			      amrex::Real* gphi_rad);
    //////////////////////

    //////////////////////
    // in debug.F90
    void print_base_cc(const amrex::Real* base);
    void print_base_edge(const amrex::Real* base);

    void print_mf(const int* lev, const int* lo, const int* hi,
                  const amrex::Real* mf,
                  const int* m_lo, const int* m_hi, const int* nc_m);

    void print_edge(const int* lev, const int* lo, const int* hi,
                    const amrex::Real* mf,
                    const int* m_lo, const int* m_hi, const int* nc_m);
    //////////////////////

    //////////////////////
    // in diag.F90
    void diag(const int* lev, const int* lo, const int* hi,
	      const amrex::Real* scal, const int* s_lo, const int* s_hi, const int* nc_s,
	      const amrex::Real* rho_Hnuc, const int* hn_lo, const int* hn_hi,
	      const amrex::Real* rho_Hext, const int* he_lo, const int* he_hi,
	      const amrex::Real* rho0, const amrex::Real* p0,
	      const amrex::Real* u, const int* u_lo, const int* u_hi,
	      const amrex::Real* w0, const amrex::Real* dx,
	      amrex::Real* Mach_max, amrex::Real* temp_max,
	      amrex::Real* enuc_max, amrex::Real* Hext_max,
              const int* mask, const int* m_lo, const int* m_hi,
              const int* use_mask);

    void diag_sphr(const int* lev, const int* lo, const int* hi,
		   const amrex::Real* scal, const int* s_lo, const int* s_hi, const int* nc_s,
		   const amrex::Real* rho_Hnuc, const int* hn_lo, const int* hn_hi,
		   const amrex::Real* rho_Hext, const int* he_lo, const int* he_hi,
		   const amrex::Real* u, const int* u_lo, const int* u_hi,
		   const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
		   const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
		   const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
		   const amrex::Real* w0r, const int* wr_lo, const int* wr_hi,
		   const amrex::Real* dx,
		   const amrex::Real* normal, const int* n_lo, const int* n_hi,
		   amrex::Real* T_max, amrex::Real* coord_Tmax, amrex::Real* vel_Tmax,
		   amrex::Real* enuc_max, amrex::Real* coord_enucmax, amrex::Real* vel_enucmax,
		   amrex::Real* kin_ener, amrex::Real* int_ener, amrex::Real* nuc_ener,
		   amrex::Real* U_max, amrex::Real* Mach_max,
		   int* ncenter, amrex::Real* T_center, amrex::Real* vel_center,
                   const int* mask, const int* m_lo, const int* m_hi,
                   const int* use_mask);

    void diag_grav_energy(amrex::Real* grav_ener,
			  const amrex::Real* rho0,
                          const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in enforce_HSE.f90
    void enforce_HSE(const amrex::Real* rho0,
                           amrex::Real* p0,
                     const amrex::Real* grav_cell,
		     const amrex::Real* r_cc_loc,
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
                         const int* phys_bc);

    void PC_EDGE_INTERP(const int* flo, const int* fhi,
                        const int* nc, const int* ratio, const int* dir,
                        amrex::Real* crse, const int* c_lo, const int* c_hi, const int* nc_c,
                        amrex::Real* fine, const int* f_lo, const int* f_hi, const int* nc_f);

    void EDGE_INTERP(const int* flo, const int* fhi,
                     const int* nc, const int* ratio, const int* dir,
                     amrex::Real* fine, const int* f_lo, const int* f_hi, const int* nc_f);
    //////////////////////

    //////////////////////
    // in fill_3d_data.F90
    void put_1d_array_on_cart(const int* lo, const int* hi, const int lev,
                              amrex::Real* s0_cart,
                              const int* s0_cart_lo, const int* s0_cart_hi, const int nc_s,
                              const amrex::Real* s0, const int is_input_edge_centered,
                              const int is_output_a_vector);

    void put_1d_array_on_cart_sphr(const int* lo, const int* hi,
				   amrex::Real* s0_cart,
				   const int* s0_cart_lo, const int* s0_cart_hi, const int nc_s,
				   const amrex::Real* s0, const amrex::Real* dx,
				   const int is_input_edge_centered,
				   const int is_output_a_vector,
				   const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
				   const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);

    void addw0(const int* lev, const int* lo, const int* hi,
               amrex::Real* uedge, const int* u_lo, const int* u_hi,
               amrex::Real* vedge, const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
               amrex::Real* wedge, const int* w_lo, const int* w_hi,
#endif
               const amrex::Real* w0, const amrex::Real* mult);

    void addw0_sphr(const int* lo, const int* hi,
		    amrex::Real* umac, const int* u_lo, const int* u_hi,
		    amrex::Real* vmac, const int* v_lo, const int* v_hi,
		    amrex::Real* wmac, const int* w_lo, const int* w_hi,
		    const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
		    const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
		    const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
		    const amrex::Real* mult);

    void make_w0mac_sphr(const int* lo, const int* hi,
			 const amrex::Real* w0,
			 amrex::Real* w0macx, const int* x_lo, const int* x_hi,
			 amrex::Real* w0macy, const int* y_lo, const int* y_hi,
			 amrex::Real* w0macz, const int* z_lo, const int* z_hi,
			 amrex::Real* w0_cart, const int* w0_lo, const int* w0_hi, const int* nc_w0,
			 const amrex::Real* dx,
			 const amrex::Real* r_edge_loc);

    void make_s0mac_sphr(const int* lo, const int* hi,
			 const amrex::Real* s0,
			 amrex::Real* s0macx, const int* x_lo, const int* x_hi,
			 amrex::Real* s0macy, const int* y_lo, const int* y_hi,
			 amrex::Real* s0macz, const int* z_lo, const int* z_hi,
			 amrex::Real* s0_cart, const int* s0_lo, const int* s0_hi,
			 const amrex::Real* dx,
			 const amrex::Real* r_cc_loc);

    void make_s0mac_sphr_irreg(const int* lo, const int* hi,
			       const amrex::Real* s0,
			       amrex::Real* s0macx, const int* x_lo, const int* x_hi,
			       amrex::Real* s0macy, const int* y_lo, const int* y_hi,
			       amrex::Real* s0macz, const int* z_lo, const int* z_hi,
			       amrex::Real* s0_cart, const int* s0_lo, const int* s0_hi,
			       const amrex::Real* dx,
			       const amrex::Real* r_cc_loc);

    void make_normal(amrex::Real* normal, const int* n_lo, const int* n_hi,
		     const amrex::Real* dx);

    void put_data_on_faces(const int* lo, const int* hi,
			   const amrex::Real* scc, const int* cc_lo, const int* cc_hi,
			   amrex::Real* facex, const int* x_lo, const int* x_hi,
			   amrex::Real* facey, const int* y_lo, const int* y_hi,
#if (AMREX_SPACEDIM == 3)
			   amrex::Real* facez, const int* z_lo, const int* z_hi,
#endif
			   const int* harmonic_avg);
    //////////////////////

    //////////////////////
    // in initdata.f90
    void initdata(const int* lev,
                  const amrex::Real* time,
                  const int* lo, const int* hi,
                  const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
                  const amrex::Real* vel,  const int* vel_lo,  const int* vel_hi,  const int* nc_v,
                  const amrex::Real* s0_init,
                  const amrex::Real* p0_init,
                  const amrex::Real* dx);

    void initdata_sphr(const amrex::Real* time,
		       const int* lo, const int* hi,
		       const amrex::Real* scal, const int* scal_lo, const int* scal_hi, const int* nc_s,
		       const amrex::Real* vel,  const int* vel_lo,  const int* vel_hi,  const int* nc_v,
		       const amrex::Real* s0_init,
		       const amrex::Real* p0_init,
		       const amrex::Real* dx,
		       const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
		       const amrex::Real* cell_cc_to_r, const int* ccr_lo, const int* ccr_hi);
    //////////////////////

    //////////////////////
    // in maestro_init.f90

    void maestro_network_init();
    void maestro_extern_init();
    void maestro_conductivity_init();
    void get_num_spec(int* nspec);
    void get_spec_names(int* spec_names, int* ispec, int* len);
    void get_spec_az(const int* ispec, amrex::Real* A, amrex::Real* Z);

    // set fortran-specific parameters in meth_params.F90
    void set_method_params (const int* Density,
                            const int* Enthalpy,
                            const int* FirstSpec,
                            const int* Temperature,
                            const int* Pressure,
                            const int* Nscalars,
                            const amrex::Real* prob_lo_in,
                            const amrex::Real* prob_hi_in);

    void set_rel_eps (const amrex::Real* rel_eps_in);
    void get_rel_eps (const amrex::Real* rel_eps_in);

    void init_base_state (const amrex::Real* s0_init,
                          const amrex::Real* p0_init,
                          const amrex::Real* rho0,
                          const amrex::Real* rhoh0,
                          const amrex::Real* p0,
                          const amrex::Real* tempbar,
			  const amrex::Real* tempbar_init);

    void init_base_state_irreg (const amrex::Real* s0_init,
				const amrex::Real* p0_init,
				const amrex::Real* rho0,
				const amrex::Real* rhoh0,
				const amrex::Real* p0,
				const amrex::Real* tempbar,
				const amrex::Real* tempbar_init,
				const amrex::Real* r_cc_loc,
				const amrex::Real* r_edge_loc);
    //////////////////////

    //////////////////////
    // in make_beta0.f90
    void make_beta0(      amrex::Real* beta0,
                    const amrex::Real* rho0,
                    const amrex::Real* p0,
                    const amrex::Real* gamma1bar,
                    const amrex::Real* grav_cell);

    void make_beta0_irreg(      amrex::Real* beta0,
			   const amrex::Real* rho0,
			   const amrex::Real* p0,
			   const amrex::Real* gamma1bar,
			   const amrex::Real* grav_cell,
		      	   const amrex::Real* r_cc_loc,
			   const amrex::Real* r_edge_loc);
    //////////////////////


    //////////////////////
    // in make_edge_scal.F90
#if (AMREX_SPACEDIM == 2)
    void make_edge_scal_2d(
#elif (AMREX_SPACEDIM == 3)
    void make_edge_scal_3d(
#endif
                      const int* domlo, const int* domhi,
                      const int* lo, const int* hi,
                      const amrex::Real* s,      const int* s_lo, const int* s_hi, const int* nc_s, int ng_s,
                            amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
                            amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
#if (AMREX_SPACEDIM == 3)
                            amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
#endif
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
                      const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
#endif
                      int ng_um,
                      const amrex::Real* force,  const int* f_lo, const int* f_hi, const int* nc_f,
                      const amrex::Real* dx, const amrex::Real* dt, const int* is_vel,
                      const int* adv_bc, const int* nbccomp, const int* comp, const int* bccomp,
                      const int* is_conservative);

    //////////////////////

    //////////////////////
    // in make_eta.F90
    void make_etarho_planar(amrex::Real* etarho_ec, amrex::Real* etarho_cc,
			    const amrex::Real* etarhosum, const amrex::Real* ncell);

    void sum_etarho(const int* lev, const int* domlo, const int* domhi,
		    const int* lo,  const int* hi,
		    const amrex::Real* etarhoflux, const int* x_lo, const int* x_hi,
                    amrex::Real* etarhosum);

    void construct_eta_cart(const int* lo, const int* hi,
			    const amrex::Real* rho_old, const int* ro_lo, const int* ro_hi,
			    const amrex::Real* rho_new, const int* rn_lo, const int* rn_hi,
			    const amrex::Real* umac,    const int* u_lo,  const int* u_hi,
			    const amrex::Real* vmac,    const int* v_lo,  const int* v_hi,
			    const amrex::Real* wmac,    const int* w_lo,  const int* w_hi,
			    const amrex::Real* w0macx,  const int* x_lo,  const int* x_hi,
			    const amrex::Real* w0macy,  const int* y_lo,  const int* y_hi,
			    const amrex::Real* w0macz,  const int* z_lo,  const int* z_hi,
			    const amrex::Real* normal,  const int* n_lo,  const int* n_hi,
			    amrex::Real* eta_cart, const int* e_lo,  const int* e_hi,
			    const amrex::Real* rho0_old, const amrex::Real* rho0_new,
			    const amrex::Real* dx,
			    const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc,
			    const amrex::Real* cc_to_r, const int* ccr_lo, const int* ccr_hi);
    //////////////////////

    //////////////////////
    // in make_explicit_thermal.F90
    void make_thermal_coeffs(const int* lo,  const int* hi,
			     const amrex::Real* scal, const int* s_lo, const int* s_hi,
			     amrex::Real* Tcoeff, const int* t_lo, const int* t_hi,
			     amrex::Real* hcoeff, const int* h_lo, const int* h_hi,
			     amrex::Real* Xkcoeff, const int* xk_lo, const int* xk_hi,
			     amrex::Real* pcoeff, const int* p_lo, const int* p_hi);
    //////////////////////

    //////////////////////
    // in make_flux.F90
#if (AMREX_SPACEDIM == 2)
    void make_rhoX_flux_2d(
#elif (AMREX_SPACEDIM == 3)
    void make_rhoX_flux_3d(
#endif
			   const int* lev, const int* lo, const int* hi,
			   amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
			   amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
#if (AMREX_SPACEDIM == 3)
			   amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
#endif
			   amrex::Real* etarhoflux, const int* eta_lo, const int* eta_hi,
			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
			   amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
#if (AMREX_SPACEDIM == 3)
			   amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
#endif
			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
			   const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
			   const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
#endif
			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
			   const amrex::Real* rho0_predicted_edge, const amrex::Real* w0,
			   const int* startcomp, const int* endcomp);

#if (AMREX_SPACEDIM == 3)
    void make_rhoX_flux_3d_sphr(const int* lo, const int* hi,
				amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
				amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
				amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
				amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
				amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
				amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
				const amrex::Real* umac,   const int* u_lo, const int* u_hi,
				const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
				const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
				const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
				const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
				const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
				const amrex::Real* rho0_edgex, const int* rx_lo, const int* rx_hi,
				const amrex::Real* rho0_edgey, const int* ry_lo, const int* ry_hi,
				const amrex::Real* rho0_edgez, const int* rz_lo, const int* rz_hi,
				const int* startcomp, const int* endcomp);
#endif

#if (AMREX_SPACEDIM == 2)
    void make_rhoh_flux_2d(
#elif (AMREX_SPACEDIM == 3)
    void make_rhoh_flux_3d(
#endif
			   const int* lev, const int* lo, const int* hi,
			   amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
			   amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
#if (AMREX_SPACEDIM == 3)
			   amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
#endif
			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
			   amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
#if (AMREX_SPACEDIM == 3)
			   amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
#endif
			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
			   const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
			   const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
#endif
			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
			   const amrex::Real* rhoh0_old, const amrex::Real* rhoh0_edge_old,
			   const amrex::Real* rhoh0_new, const amrex::Real* rhoh0_edge_new,
			   const amrex::Real* w0);

#if (AMREX_SPACEDIM == 3)
    void make_rhoh_flux_3d_sphr(const int* lo, const int* hi,
				amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
				amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
				amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
				amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
				amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
				amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
				const amrex::Real* umac,   const int* u_lo, const int* u_hi,
				const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
				const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
				const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
				const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
				const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
				const amrex::Real* rho0_edgex, const int* rx_lo, const int* rx_hi,
				const amrex::Real* rho0_edgey, const int* ry_lo, const int* ry_hi,
				const amrex::Real* rho0_edgez, const int* rz_lo, const int* rz_hi,
				const amrex::Real* h0_edgex,   const int* hx_lo, const int* hx_hi,
				const amrex::Real* h0_edgey,   const int* hy_lo, const int* hy_hi,
				const amrex::Real* h0_edgez,   const int* hz_lo, const int* hz_hi);

    void make_rhoh_flux_3d_sphr_irreg(const int* lo, const int* hi,
				amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
				amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
				amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
				amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
				amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
				amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
				const amrex::Real* umac,   const int* u_lo, const int* u_hi,
				const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
				const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
				const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
				const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
				const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
				const amrex::Real* rhoh0_edgex,   const int* hx_lo, const int* hx_hi,
				const amrex::Real* rhoh0_edgey,   const int* hy_lo, const int* hy_hi,
				const amrex::Real* rhoh0_edgez,   const int* hz_lo, const int* hz_hi);
#endif
    //////////////////////

    //////////////////////
    // in make_gamma.F90
    void make_gamma(const int* lo, const int* hi, const int lev,
                          amrex::Real* gamma, const int* g_lo, const int* g_hi,
                    const amrex::Real* scal,  const int* s_lo, const int* s_hi,
                    const amrex::Real* p0);

    void make_gamma_sphr(const int* lo, const int* hi,
			 amrex::Real* gamma, const int* g_lo, const int* g_hi,
			 const amrex::Real* scal, const int* s_lo, const int* s_hi,
			 const amrex::Real* p0_cart, const int* p0_lo, const int* p0_hi);
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
    // in make_heating.f90
    void make_heating(const int* lo, const int* hi,
                            amrex::Real* rho_Hext, const int* r_lo, const int* r_hi,
                      const amrex::Real* scal,     const int* s_lo, const int* s_hi, const int* nc_s,
                      const amrex::Real* dx, const amrex::Real* time);
    //////////////////////

    //////////////////////
    // in make_pi_cc.F90
    void make_pi_cc(const int* lo, const int* hi,
                    const amrex::Real* pi,         const int* p_lo, const int* p_hi,
                          amrex::Real* pi_cc,      const int* c_lo, const int* c_hi,
                    const amrex::Real* beta0_cart, const int* b_lo, const int* b_hi);
    //////////////////////

    //////////////////////
    // in make_plot_variables.F90
    void make_magvel(const int* lo, const int* hi, const int lev,
                     const amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const amrex::Real* w0,
                     amrex::Real* magvel, const int* m_lo, const int* m_hi);

    void make_magvel_sphr(const int* lo, const int* hi,
                          const amrex::Real* vel, const int* v_lo, const int* v_hi,
                          amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                          amrex::Real* w0macy, const int* y_lo, const int* y_hi,
                          amrex::Real* w0macz, const int* z_lo, const int* z_hi,
                          amrex::Real* magvel, const int* m_lo, const int* m_hi);

    void make_velrc(const int* lo, const int* hi,
                        const amrex::Real* vel, const int* v_lo, const int* v_hi,
                        const amrex::Real* w0rcart, const int* w_lo, const int* w_hi,
                        const amrex::Real* normal, const int* n_lo, const int* n_hi,
                        amrex::Real* rad_vel, const int* r_lo, const int* r_hi,
                        amrex::Real* circ_vel, const int* c_lo, const int* c_hi);

    void make_ad_excess(const int* lo, const int* hi,
                        const amrex::Real* state, const int* s_lo, const int* s_hi,
                        const int nc_s,
                        amrex::Real* ad_excess, const int* a_lo, const int* a_hi);

    void make_ad_excess_sphr(const int* lo, const int* hi,
                        const amrex::Real* state, const int* s_lo, const int* s_hi,
                        const int nc_s,
                        const amrex::Real* normal, const int* n_lo, const int* n_hi,
                        amrex::Real* ad_excess, const int* a_lo, const int* a_hi);

    void make_vorticity(const int* lo, const int* hi,
                        const amrex::Real* vel, const int* v_lo, const int* v_hi,
                        const amrex::Real* dx,
                        amrex::Real* vorticity, const int* d_lo, const int* d_hi,
                        const int* bc);

    void make_deltagamma(const int* lo, const int* hi, const int lev,
                        const amrex::Real* state, const int* s_lo, const int* s_hi,
                        const int nc_s,
                        const amrex::Real* p0, const amrex::Real* gamma1bar,
                        amrex::Real* deltagamma, const int* d_lo, const int* d_hi);

    void make_deltagamma_sphr(const int* lo, const int* hi,
                            const amrex::Real* state, const int* s_lo, const int* s_hi,
                            const int nc_s,
                            const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
                            const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi,
                            amrex::Real* deltagamma, const int* d_lo, const int* d_hi);

    void make_entropy(const int* lo, const int* hi, const int lev,
                     const amrex::Real* state, const int* s_lo, const int* s_hi,
                     const int nc_s,
                     amrex::Real* entropy, const int* e_lo, const int* e_hi);

    void make_divw0(const int* lo, const int* hi, const int lev,
                    const amrex::Real* w0, const amrex::Real* dx,
                    amrex::Real* divw0, const int* d_lo, const int* d_hi);

    void make_divw0_sphr(const int* lo, const int* hi,
                         const amrex::Real* w0macx, const int* x_lo, const int* x_hi,
                         const amrex::Real* w0macy, const int* y_lo, const int* y_hi,
                         const amrex::Real* w0macz, const int* z_lo, const int* z_hi,
                         const amrex::Real* dx,
                         amrex::Real* divw0, const int* d_lo, const int* d_hi);

     void make_pidivu(const int* lo, const int* hi,
                     const amrex::Real* vel, const int* v_lo, const int* v_hi,
                     const amrex::Real* dx,
                     const amrex::Real* pi_cc, const int* p_lo, const int* p_hi, const int nc,
                     amrex::Real* pidivuw0, const int* d_lo, const int* d_hi);

    void make_abar(const int* lo, const int* hi,
                 const amrex::Real* state, const int* s_lo, const int* s_hi, const int nc_s,
                 amrex::Real* abar, const int* a_lo, const int* a_hi);
/////////////////////

    //////////////////////
    // in make_psi.F90
    void make_psi_planar(const amrex::Real* etarho_cc,
			       amrex::Real* psi);

    void make_psi_spherical(amrex::Real* psi,
			    const amrex::Real* w0,
			    const amrex::Real* gamma1bar,
			    const amrex::Real* p0_avg,
			    const amrex::Real* Sbar_in,
			    const amrex::Real* r_cc_loc,
			    const amrex::Real* r_edge_loc);

    void make_psi_irreg(const amrex::Real* etarho_cc,
			const amrex::Real* grav_cell,
			      amrex::Real* psi);

    // void make_psi_irreg(amrex::Real* psi,
    // 			const amrex::Real* p0_old,
    // 			const amrex::Real* p0_new,
    // 			const amrex::Real* dt);
    //////////////////////

    //////////////////////
    // in make_S.f90

    // compute S from constraint div(u) = S
    void make_S_cc(const int* lo, const int* hi, const int lev,
                         amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
            amrex::Real* delta_gamma1_term,  const int* dg_lo, const int* dg_hi,
                 amrex::Real* delta_gamma1,  const int* df_lo, const int* df_hi,
                   const amrex::Real* scal,  const int* c_lo, const int* c_hi,
                   const amrex::Real* u,     const int* u_lo, const int* u_hi,
                   const amrex::Real* rodot, const int* r_lo, const int* r_hi,
                   const amrex::Real* rHnuc, const int* n_lo, const int* n_hi,
                   const amrex::Real* rHext, const int* e_lo, const int* e_hi,
                   const amrex::Real* therm, const int* t_lo, const int* t_hi,
                   const amrex::Real* p0, const amrex::Real* gamma1bar, const amrex::Real* dx);

     void make_S_cc_sphr(const int* lo, const int* hi,
                    amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
                    amrex::Real* delta_gamma1_term,const int* dg_lo, const int* dg_hi,
                    amrex::Real* delta_gamma1,const int* df_lo, const int* df_hi,
                    const amrex::Real* scal,  const int* c_lo, const int* c_hi,
                    const amrex::Real* u,     const int* u_lo, const int* u_hi,
                    const amrex::Real* rodot, const int* r_lo, const int* r_hi,
                    const amrex::Real* rHnuc, const int* n_lo, const int* n_hi,
                    const amrex::Real* rHext, const int* e_lo, const int* e_hi,
                    const amrex::Real* therm, const int* t_lo, const int* t_hi,
                    const amrex::Real* p0_cart, const int* p0_lo, const int* p0_hi,
                    const amrex::Real* gradp0_cart, const int* gp0_lo, const int* gp0_hi,
                    const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi,
                    const amrex::Real* normal, const int* no_lo, const int* no_hi);

    // rhcc = beta0 * (S - Sbar)
    void make_rhcc_for_nodalproj(const int* lo, const int* hi, const int lev,
                                 amrex::Real* rhcc,  const int* c_lo, const int* c_hi,
                                 const amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
                                 const amrex::Real* Sbar,
                                 const amrex::Real* beta0,
                                 const amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi);

    void make_rhcc_for_nodalproj_sphr(const int* lo, const int* hi,
				      amrex::Real* rhcc,  const int* c_lo, const int* c_hi,
				      const amrex::Real* S_cc, const int* s_lo, const int* s_hi,
				      const amrex::Real* Sbar_cart, const int* sb_lo, const int* sb_hi,
				      const amrex::Real* beta0_cart, const int* b_lo, const int* b_hi,
              const amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi);

    void create_correction_cc(const int* lo, const int* hi, const int lev,
			      amrex::Real* correction_cc, const int* c_lo, const int* c_hi,
			      const amrex::Real* delta_p_term, const int* dp_lo, const int* dp_hi,
			      const amrex::Real* beta0,
			      const amrex::Real* gamma1bar,
			      const amrex::Real* p0,
			      const amrex::Real dt);

    void create_correction_cc_sphr(const int* lo, const int* hi,
				   amrex::Real* correction_cc, const int* c_lo, const int* c_hi,
				   const amrex::Real* delta_p_term, const int* dp_lo, const int* dp_hi,
				   const amrex::Real* beta0_cart, const int* b_lo, const int* b_hi,
				   const amrex::Real* gamma1bar_cart, const int* g_lo, const int* g_hi,
				   const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
				   const amrex::Real* rho0_cart, const int* r_lo, const int* r_hi,
				   const amrex::Real dt);

    // rhcc = beta0 * (S - Sbar) + beta0 * delta_chi
    void make_rhcc_for_macproj(const int* lo, const int* hi, const int lev,
                               amrex::Real* rhcc,  const int* c_lo, const int* c_hi,
                               const amrex::Real* S_cc,  const int* s_lo, const int* s_hi,
                               const amrex::Real* Sbar,
                               const amrex::Real* beta0,
                               const amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi,
                               const amrex::Real* gamma1bar,
			                         const amrex::Real* p0,
                               const amrex::Real* delta_p_term, const int* dp_lo, const int* dp_hi,
			             amrex::Real* delta_chi, const int* dc_lo, const int* dc_hi,
			       const amrex::Real dt, const int is_predictor);

    void make_rhcc_for_macproj_sphr(const int* lo, const int* hi,
				    amrex::Real* rhcc,  const int* c_lo, const int* c_hi,
				    const amrex::Real* S_cc, const int* s_lo, const int* s_hi,
				    const amrex::Real* Sbar_cart, const int* sb_lo, const int* sb_hi,
				    const amrex::Real* beta0_cart, const int* b_lo, const int* b_hi,
				    const amrex::Real* rho0_cart, const int* r_lo, const int* r_hi,
                    const amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi,
				    const amrex::Real* gamma1bar_cart, const int* g1_lo, const int* g1_hi,
				    const amrex::Real* p0_cart, const int* p0_lo, const int* p0_hi,
				    const amrex::Real* delta_p_term, const int* dp_lo, const int* dp_hi,
				    amrex::Real* delta_chi, const int* dc_lo, const int* dc_hi,
				    const amrex::Real dt, const int is_predictor);

    void create_correction_delta_gamma1_term(const int* lo, const int* hi, const int lev,
				    amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi,
				    const amrex::Real* delta_gamma1, const int* df_lo, const int* df_hi,
			      const amrex::Real* gamma1bar,
			      const amrex::Real* psi,
			      const amrex::Real* p0);

    void create_correction_delta_gamma1_term_sphr(const int* lo, const int* hi,
				   amrex::Real* delta_gamma1_term, const int* dg_lo, const int* dg_hi,
                   const amrex::Real* delta_gamma1, const int* df_lo, const int* df_hi,
				   const amrex::Real* gamma1bar_cart, const int* g1_lo, const int* g1_hi,
				   const amrex::Real* psi_cart, const int* ps_lo, const int* ps_hi,
				   const amrex::Real* p0_cart, const int* p0_lo, const int* p0_hi);

    //////////////////////


    //////////////////////
    // in make_scal_force.F90
    void modify_scal_force(const int* lo, const int* hi, const int lev,
                           amrex::Real* force, const int* f_lo, const int* f_hi,
                           const amrex::Real* scal,  const int* s_lo, const int* s_hi,
                           const amrex::Real* umac,  const int* u_lo, const int* u_hi,
                           const amrex::Real* vmac,  const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
                           const amrex::Real* wmac,  const int* w_lo, const int* w_hi,
#endif
                           const amrex::Real* s0, const amrex::Real* s0_edge, const amrex::Real* w0,
                           const amrex::Real* dx, const int do_fullform);

    void modify_scal_force_sphr(const int* lo, const int* hi,
                const int* domlo, const int* domhi,
				amrex::Real* force, const int* f_lo, const int* f_hi,
				const amrex::Real* scal,  const int* s_lo, const int* s_hi,
				const amrex::Real* umac,  const int* u_lo, const int* u_hi,
				const amrex::Real* vmac,  const int* v_lo, const int* v_hi,
				const amrex::Real* wmac,  const int* w_lo, const int* w_hi,
				const amrex::Real* s0,  const int* s0_lo, const int* s0_hi,
				const amrex::Real* w0,
				const amrex::Real* dx, const int do_fullform,
				const amrex::Real* divu_cart,  const int* d_lo, const int* d_hi);

    void mkrhohforce(const int* lo, const int* hi, const int lev,
                     amrex::Real* rhoh_force, const int* f_lo, const int* f_hi,
#if (AMREX_SPACEDIM == 2)
                     const amrex::Real* vmac,  const int* v_lo, const int* v_hi,
#elif (AMREX_SPACEDIM == 3)
                     const amrex::Real* wmac,  const int* w_lo, const int* w_hi,
#endif
                     const amrex::Real* thermal, const int* t_lo, const int* t_hi,
                     const amrex::Real* p0, const amrex::Real* rho0, const amrex::Real* grav,
                     const amrex::Real* psi, const int is_prediction, const int add_thermal);

    void mkrhohforce_sphr(const int* lo, const int* hi,
			  amrex::Real* rhoh_force, const int* f_lo, const int* f_hi,
			  const amrex::Real* umac,  const int* u_lo, const int* u_hi,
			  const amrex::Real* vmac,  const int* v_lo, const int* v_hi,
			  const amrex::Real* wmac,  const int* w_lo, const int* w_hi,
			  const amrex::Real* thermal, const int* t_lo, const int* t_hi,
			  const amrex::Real* p0_cart, const int* p_lo, const int* p_hi,
			  const amrex::Real* p0macx,  const int* x_lo, const int* x_hi,
			  const amrex::Real* p0macy,  const int* y_lo, const int* y_hi,
			  const amrex::Real* p0macz,  const int* z_lo, const int* z_hi,
			  const amrex::Real* psi_cart, const int* ps_lo, const int* ps_hi,
              const amrex::Real* dx,
			  const int is_prediction, const int add_thermal);
    //////////////////////

    //////////////////////
    // in make_vel_force.f90
    void make_vel_force(const int* lo, const int* hi, const int lev,
                        amrex::Real* vel_force, const int* f_lo, const int* f_hi,
                        const amrex::Real* gpi, const int* g_lo, const int* g_hi,
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
                        const int do_add_utilde_force);

    void make_vel_force_sphr(const int* lo, const int* hi,
			     amrex::Real* vel_force, const int* f_lo, const int* f_hi,
			     const amrex::Real* gpi, const int* g_lo, const int* g_hi,
			     const amrex::Real* rho, const int* r_lo, const int* r_hi,
			     const amrex::Real* uedge, const int* u_lo, const int* u_hi,
			     const amrex::Real* vedge, const int* v_lo, const int* v_hi,
			     const amrex::Real* wedge, const int* w_lo, const int* w_hi,
			     const amrex::Real* normal, const int* n_lo, const int* n_hi,
			     const amrex::Real* gradw0_cart, const int* gw_lo, const int* gw_hi,
			     const amrex::Real* w0_force_cart, const int* wf_lo, const int* wf_hi,
			     const amrex::Real* rho0_cart, const int* r0_lo, const int* r0_hi,
                 const amrex::Real* grav_cart, const int* gr_lo, const int* gr_hi,
			     const amrex::Real* dx,
			     const int do_add_utilde_force);
    //////////////////////

    //////////////////////
    // in make_w0.f90
    void make_w0(amrex::Real* w0, const amrex::Real* w0_old,
                 const amrex::Real* w0_force, const amrex::Real* Sbar_in,
                 const amrex::Real* rho0_old, const amrex::Real* rho0_new,
                 const amrex::Real* p0_old, const amrex::Real* p0_new,
                 const amrex::Real* gamma1bar_old, const amrex::Real* gamma1bar_new,
                 const amrex::Real* p0_minus_peosbar,
                 const amrex::Real* etarho_ec, const amrex::Real* etarho_cc,
                 const amrex::Real* delta_chi_w0, const amrex::Real* r_cc_loc,
                 const amrex::Real* r_edge_loc, const amrex::Real* dt,
                 const amrex::Real* dtold, const int* is_predictor);
    //////////////////////

    //////////////////////
    // in meth_params.F90

    // read in shared C++/Fortran parameters from inputs file
    void read_method_params();

    // deallocate strings
    void finalize_meth_params();

    //////////////////////

    //////////////////////
    // in mkutrans.F90
#if (AMREX_SPACEDIM == 2)
    void mkutrans_2d(
#elif (AMREX_SPACEDIM == 3)
    void mkutrans_3d(
#endif
                      const int* lev, const int* domlo, const int* domhi,
                      const int* lo, const int* hi,
                      const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int* nc_ut, int ng_ut,
                      const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int* nc_uf, int ng_uf,
                      amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                      amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
#if (AMREX_SPACEDIM == 3)
                      amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
		      const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
		      const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
		      const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
#endif
                      const amrex::Real* w0, const amrex::Real* dx, const amrex::Real* dt,
                      const int* adv_bc, const int* phys_bc);

    //////////////////////

    //////////////////////
    // in mac_solver.F90
    void mac_solver_rhs(const int* lo, const int* hi, const int lev,
			amrex::Real* newrhs, const int* nrhs_lo, const int* nrhs_hi,
			const amrex::Real* oldrhs, const int* orhs_lo, const int* orhs_hi,
                        const amrex::Real* uedge, const int* u_lo, const int* u_hi,
                        const amrex::Real* vedge, const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
			const amrex::Real* wedge, const int* w_lo, const int* w_hi,
#endif
			const amrex::Real* dx);

    void mult_beta0(const int* lev, const int* lo, const int* hi,
               amrex::Real* uedge, const int* u_lo, const int* u_hi,
               amrex::Real* vedge, const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
               amrex::Real* wedge, const int* w_lo, const int* w_hi,
#endif
	       const amrex::Real* beta0_edge,
               const amrex::Real* beta0, const int* mult_or_div);

    void mac_bcoef_face(const int* lev, const int* lo, const int* hi,
			amrex::Real* xface, const int* x_lo, const int* x_hi,
			amrex::Real* yface, const int* y_lo, const int* y_hi,
#if (AMREX_SPACEDIM == 3)
			amrex::Real* zface, const int* z_lo, const int* z_hi,
#endif
			const amrex::Real* rhocc, const int* r_lo, const int* r_hi);
    //////////////////////

    //////////////////////
    // in regrid.f90
    void regrid_base_state_cc(amrex::Real* state_cc);
    void regrid_base_state_edge(amrex::Real* state_ec);
    //////////////////////

    //////////////////////
    // in rhoh_vs_t.F90
    void makeTfromRhoH (const int* lo, const int* hi, const int lev,
                        amrex::Real* state,  const int* s_lo, const int* s_hi,
                        const amrex::Real* p0);

    void makeTfromRhoH_sphr (const int* lo, const int* hi,
			     amrex::Real* state, const int* s_lo, const int* s_hi,
			     const amrex::Real* p0, const int* p0_lo, const int* p0_hi);

    void makeTfromRhoP (const int* lo, const int* hi, const int lev,
                        amrex::Real* state, const int* s_lo, const int* s_hi,
                        const amrex::Real* p0, const int updateRhoH);

    void makeTfromRhoP_sphr (const int* lo, const int* hi,
    			     amrex::Real* state, const int* s_lo, const int* s_hi,
    			     const amrex::Real* p0, const int* p0_lo, const int* p0_hi,
                     const int updateRhoH);

    void makePfromRhoH(const int* lo, const int* hi,
		       const amrex::Real* state, const int* s_lo, const int* s_hi,
		       const amrex::Real* temp_old, const int* t_lo, const int* t_hi,
		       amrex::Real* peos, const int* p_lo, const int* p_hi);

    void makeMachfromRhoH (const int* lo, const int* hi, const int lev,
                           const amrex::Real* state, const int* s_lo, const int* s_hi,
                           const amrex::Real* u, const int* u_lo, const int* u_hi,
                           const amrex::Real* p0, const amrex::Real* w0,
                           amrex::Real* mach, const int* m_lo, const int* m_hi);

    void makeMachfromRhoH_sphr (const int* lo, const int* hi,
                                const amrex::Real* state,  const int* s_lo, const int* s_hi,
                                const amrex::Real* u, const int* u_lo, const int* u_hi,
                                const amrex::Real* p0, const int* p0_lo, const int* p0_hi,
                                const amrex::Real* w0, const int* w_lo, const int* w_hi,
                                amrex::Real* mach, const int* m_lo, const int* m_hi);


    void makeCsfromRhoH (const int* lo, const int* hi, const int lev,
                           const amrex::Real* state, const int* s_lo, const int* s_hi,
                           const amrex::Real* p0,
                           amrex::Real* cs, const int* c_lo, const int* c_hi);

    void makeCsfromRhoH_sphr (const int* lo, const int* hi,
                          const amrex::Real* state, const int* s_lo, const int* s_hi,
                          const amrex::Real* p0cart, const int* p_lo, const int* p_hi,
                          amrex::Real* cs, const int* c_lo, const int* c_hi);
    //////////////////////

    //////////////////////
    // in sponge.F90
    void init_sponge(const amrex::Real* rho0);

    void init_sponge_irreg(const amrex::Real* rho0,
			   const amrex::Real* r_cc_loc, const amrex::Real* r_edge_loc);

    void mk_sponge (const int* lo, const int* hi,
                    amrex::Real* sponge,  const int* s_lo, const int* s_hi,
                    const amrex::Real* dx, const amrex::Real dt);
    //////////////////////


    //////////////////////
    // in tagging.f90
    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
                     const amrex::Real* state, const int* state_lo, const int* state_hi,
                     const int* tagval, const int* clearval,
                     const int* lo, const int* hi,
                     const amrex::Real* dx,
                     const amrex::Real* time,
		     const int* lev, int* tag_array);

    void tag_boxes(int* tag, const int* tag_lo, const int* tag_hi,
		   const int* tagval, const int* clearval,
		   const int* lo, const int* hi,
		   const amrex::Real* dx,
		   const amrex::Real* time,
		   const int* lev, const int* tag_array);

    void retag_array(const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const int* lev, const int* tag_array);
    //////////////////////

    //////////////////////
    // in update_scal.F90
    void update_rhoX(const int* lo, const int* hi,
			const amrex::Real* sold, const int* so_lo, const int* so_hi,
			      amrex::Real* snew, const int* sn_lo, const int* sn_hi,
			const amrex::Real* sfluxx, const int* x_lo, const int* x_hi,
			const amrex::Real* sfluxy, const int* y_lo, const int* y_hi,
#if (AMREX_SPACEDIM == 3)
			const amrex::Real* sfluxz, const int* z_lo, const int* z_hi,
#endif
			const amrex::Real* force,  const int* f_lo, const int* f_hi,
			const amrex::Real* dx, const amrex::Real dt,
			const int startcomp, const int endcomp);


    void update_rhoh(const int* lo, const int* hi, const int lev,
			const amrex::Real* sold, const int* so_lo, const int* so_hi,
			      amrex::Real* snew, const int* sn_lo, const int* sn_hi,
			const amrex::Real* sfluxx, const int* x_lo, const int* x_hi,
			const amrex::Real* sfluxy, const int* y_lo, const int* y_hi,
#if (AMREX_SPACEDIM == 3)
			const amrex::Real* sfluxz, const int* z_lo, const int* z_hi,
#endif
			const amrex::Real* force,  const int* f_lo, const int* f_hi,
			const amrex::Real* p0_new,
			const amrex::Real* dx, const amrex::Real dt,
			const int nspec);

#if (AMREX_SPACEDIM == 3)
    void update_rhoh_3d_sphr(const int* lo, const int* hi,
			 const amrex::Real* sold, const int* so_lo, const int* so_hi,
			       amrex::Real* snew, const int* sn_lo, const int* sn_hi,
			 const amrex::Real* sfluxx, const int* x_lo, const int* x_hi,
			 const amrex::Real* sfluxy, const int* y_lo, const int* y_hi,
			 const amrex::Real* sfluxz, const int* z_lo, const int* z_hi,
			 const amrex::Real* force,  const int* f_lo, const int* f_hi,
			 const amrex::Real* p0_new_cart, const int* p_lo, const int* p_hi,
			 const amrex::Real* dx, const amrex::Real dt,
			 const int nspec);
#endif
    //////////////////////

    //////////////////////
    // in update_velocity.F90
    void update_velocity(const int* lo, const int* hi, const int lev,
			 const amrex::Real* uold, const int* uo_lo, const int* uo_hi,
			       amrex::Real* unew, const int* un_lo, const int* un_hi,
			 const amrex::Real* umac, const int* u_lo, const int* u_hi,
			 const amrex::Real* vmac, const int* v_lo, const int* v_hi,
#if (AMREX_SPACEDIM == 3)
			 const amrex::Real* wmac, const int* w_lo, const int* w_hi,
#endif
			 const amrex::Real* uedgex, const int* x_lo, const int* x_hi,
			 const amrex::Real* uedgey, const int* y_lo, const int* y_hi,
#if (AMREX_SPACEDIM == 3)
			 const amrex::Real* uedgez, const int* z_lo, const int* z_hi,
#endif
			 const amrex::Real* force, const int* f_lo, const int* f_hi,
			 const amrex::Real* sponge, const int* s_lo, const int* s_hi,
			 const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt);

    void update_velocity_sphr(const int* lo, const int* hi,
			      const amrex::Real* uold, const int* uo_lo, const int* uo_hi,
			       amrex::Real* unew, const int* un_lo, const int* un_hi,
			      const amrex::Real* umac, const int* u_lo, const int* u_hi,
			      const amrex::Real* vmac, const int* v_lo, const int* v_hi,
			      const amrex::Real* wmac, const int* w_lo, const int* w_hi,
			      const amrex::Real* uedgex, const int* x_lo, const int* x_hi,
			      const amrex::Real* uedgey, const int* y_lo, const int* y_hi,
			      const amrex::Real* uedgez, const int* z_lo, const int* z_hi,
			      const amrex::Real* force, const int* f_lo, const int* f_hi,
			      const amrex::Real* sponge, const int* s_lo, const int* s_hi,
			      const amrex::Real* w0,
			      const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
			      const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
			      const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
			      const amrex::Real* dx, const amrex::Real dt);
    //////////////////////

    //////////////////////
    // in velpred.F90
#if (AMREX_SPACEDIM == 2)
    void velpred_2d(
#elif (AMREX_SPACEDIM == 3)
    void velpred_3d(
#endif
                      const int* lev, const int* domlo, const int* domhi,
                      const int* lo, const int* hi,
                      const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int* nc_ut, int ng_ut,
                      const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int* nc_uf, int ng_uf,
                      const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                      const amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
#if (AMREX_SPACEDIM == 3)
                      const amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
#endif
                      amrex::Real* umac, const int* mu_lo, const int* mu_hi,
                      amrex::Real* vmac, const int* mv_lo, const int* mv_hi,
#if (AMREX_SPACEDIM == 3)
                      amrex::Real* wmac, const int* mw_lo, const int* mw_hi,
		      const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
		      const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
		      const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
#endif
                      const amrex::Real* force, const int* f_lo, const int* f_hi, const int* nc_f, int ng_f,
                      const amrex::Real* w0, const amrex::Real* dx, const amrex::Real* dt,
                      const int* adv_bc, const int* phys_bc);

    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
