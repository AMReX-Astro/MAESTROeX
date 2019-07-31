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
    // in make_edge_scal.F90
#if (AMREX_SPACEDIM == 1)
    void make_edge_scal_1d(
                      const int* domlo, const int* domhi,
                      const int* lo, const int* hi,
                      const amrex::Real* s,      const int* s_lo, const int* s_hi, const int nc_s, int ng_s,
                            amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int nc_x,
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const int ng_um,
                      const amrex::Real* force,  const int* f_lo, const int* f_hi, const int nc_f,
                      const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                      const int* adv_bc, const int nbccomp, const int comp, const int bccomp,
                      const int is_conservative);

#elif (AMREX_SPACEDIM == 2)
  void make_edge_scal_predictor_2d(const int* lo, const int* hi, const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* s, const int* s_lo, const int* s_hi, const int nc_s, int ng_s,
                    const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                    const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                    const amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                    const amrex::Real* Im, const int* im_lo, const int* im_hi,
                    amrex::Real* sl, const int* sl_lo, const int* sl_hi,
                    amrex::Real* sr, const int* sr_lo, const int* sr_hi,
                    amrex::Real* simh, const int* si_lo, const int* si_hi,
                    const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                    const int* adv_bc, const int nbccomp, const int comp, const int bccomps);

    void make_edge_scal_2d(const int* lo, const int* hi, const int idir,
                      const int* domlo, const int* domhi,
                    const amrex::Real* s,      const int* s_lo, const int* s_hi, const int nc_s, int ng_s,
                          amrex::Real* sedge, const int* x_lo, const int* x_hi, const int nc_x,
                    const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                    const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                    const amrex::Real* Ipf, const int* ipf_lo, const int* ipf_hi,
                    const amrex::Real* Imf, const int* imf_lo, const int* imf_hi,
                    const amrex::Real* sl, const int* sl_lo, const int* sl_hi,
                    const amrex::Real* sr, const int* sr_lo, const int* sr_hi,
                    const amrex::Real* simh, const int* si_lo, const int* si_hi,
                    const amrex::Real* force,  const int* f_lo, const int* f_hi, const int nc_f,
                    const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                    const int* adv_bc, const int nbccomp, const int comp, const int bccomp,
                    const int is_conservative);

#elif (AMREX_SPACEDIM == 3)

    void make_divu(const int* lo, const int* hi,
                      amrex::Real* divu, const int* d_lo, const int* d_hi,
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                      const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
                      const amrex::Real* dx, const int is_conservative);

    void make_edge_scal_predictor_3d(const int* lo, const int* hi,
                      const int idir,
                      const int* domlo, const int* domhi,
                      const amrex::Real* s,      const int* s_lo, const int* s_hi, const int nc_s,
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                      const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
                      const amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                      const amrex::Real* Im, const int* im_lo, const int* im_hi,
                      const amrex::Real* slopez, const int* slo_lo, const int* slo_hi,
                      amrex::Real* sl, const int* sl_lo, const int* sl_hi,
                      amrex::Real* sr, const int* sr_lo, const int* sr_hi,
                      amrex::Real* simh, const int* si_lo, const int* si_hi,
                      const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                      const int* adv_bc, const int nbccomp, const int comp, const int bccomp);

    void make_edge_scal_transverse_3d(const int* lo, const int* hi,
                      const int norm_dir, const int trans_dir,
                      const int* domlo, const int* domhi,
                      const amrex::Real* s,      const int* s_lo, const int* s_hi, const int nc_s,
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                      const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
                      const amrex::Real* divu, const int* d_lo, const int* d_hi,
                      const amrex::Real* slx, const int* slx_lo, const int* slx_hi,
                      const amrex::Real* srx, const int* srx_lo, const int* srx_hi,
                      const amrex::Real* simhx, const int* six_lo, const int* six_hi,
                      const amrex::Real* sly, const int* sly_lo, const int* sly_hi,
                      const amrex::Real* sry, const int* sry_lo, const int* sry_hi,
                      const amrex::Real* simhy, const int* siy_lo, const int* siy_hi,
                      const amrex::Real* slz, const int* slz_lo, const int* slz_hi,
                      const amrex::Real* srz, const int* srz_lo, const int* srz_hi,
                      const amrex::Real* simhz, const int* siz_lo, const int* siz_hi,
                      amrex::Real* simh_trans, const int* sit_lo, const int* sit_hi,
                      const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                      const int* adv_bc, const int nbccomp, const int comp, const int bccomp,
                      const int is_conservative);

    void make_edge_scal_3d(const int* lo, const int* hi,
                      const int idir,
                      const int* domlo, const int* domhi,
                      const amrex::Real* s,      const int* s_lo, const int* s_hi, const int nc_s,
                            amrex::Real* sedge, const int* x_lo, const int* x_hi, const int nc_x,
                      const amrex::Real* umac,   const int* u_lo, const int* u_hi,
                      const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
                      const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
                      const amrex::Real* Ipf, const int* ipf_lo, const int* ipf_hi,
                      const amrex::Real* Imf, const int* imf_lo, const int* imf_hi,
                      const amrex::Real* sl, const int* sl_lo, const int* sl_hi,
                      const amrex::Real* sr, const int* slr_lo, const int* sr_hi,
                      const amrex::Real* simhxy, const int* xy_lo, const int* xy_hi,
                      const amrex::Real* simhxz, const int* xz_lo, const int* xz_hi,
                      const amrex::Real* simhyx, const int* yx_lo, const int* yx_hi,
                      const amrex::Real* simhyz, const int* yz_lo, const int* yz_hi,
                      const amrex::Real* simhzx, const int* zx_lo, const int* zx_hi,
                      const amrex::Real* simhzy, const int* zy_lo, const int* zy_hi,
                      const amrex::Real* force,  const int* f_lo, const int* f_hi, const int nc_f,
                      const amrex::Real* dx, const amrex::Real dt, const int is_vel,
                      const int* adv_bc, const int nbccomp, const int comp, const int bccomp,
                      const int is_conservative);

#endif
    //////////////////////

    //////////////////////
    // in make_flux.F90
#if (AMREX_SPACEDIM == 1)
    void make_rhoX_flux_1d(
			   const int* lev, const int* lo, const int* hi,
			   amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
			   amrex::Real* etarhoflux, const int* eta_lo, const int* eta_hi,
			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
			   const amrex::Real* rho0_predicted_edge, const amrex::Real* w0,
			   const int* startcomp, const int* endcomp);

#elif (AMREX_SPACEDIM == 2)
   void make_rhoX_flux_2d(const int* lo, const int* hi,
               const int lev, const int idir,
			   amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int nc_fx,
			   amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int nc_fy,
			   amrex::Real* etarhoflux, const int* eta_lo, const int* eta_hi,
			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int nc_x,
			   amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int nc_y,
			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
			   const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
			   const amrex::Real* rho0_predicted_edge, const amrex::Real* w0,
			   const int startcomp, const int endcomp);

#elif (AMREX_SPACEDIM == 3)
   void make_rhoX_flux_3d(
			   const int* lo, const int* hi, const int lev, const int idir,
			   amrex::Real* sflux, const int* f_lo, const int* f_hi, const int nc_f,
			   amrex::Real* etarhoflux, const int* eta_lo, const int* eta_hi,
			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int nc_x,
			   amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int nc_y,
			   amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int nc_z,
			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
			   const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
			   const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
			   const amrex::Real* rho0_predicted_edge, const amrex::Real* w0,
			   const int startcomp, const int endcomp);
#endif

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

// #if (AMREX_SPACEDIM == 1)
//     void make_rhoh_flux_1d(
// #elif (AMREX_SPACEDIM == 2)
//     void make_rhoh_flux_2d(
// #elif (AMREX_SPACEDIM == 3)
//     void make_rhoh_flux_3d(
// #endif
// 			   const int* lev, const int* lo, const int* hi,
// 			   amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
// #if (AMREX_SPACEDIM >= 2)
// 			   amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
// #if (AMREX_SPACEDIM == 3)
// 			   amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
// #endif
// #endif
// 			   amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
// #if (AMREX_SPACEDIM >= 2)
// 			   amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
// #if (AMREX_SPACEDIM == 3)
// 			   amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
// #endif
// #endif
// 			   const amrex::Real* umac,   const int* u_lo, const int* u_hi,
// #if (AMREX_SPACEDIM >= 2)
// 			   const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
// #if (AMREX_SPACEDIM == 3)
// 			   const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
// #endif
// #endif
// 			   const amrex::Real* rho0_old, const amrex::Real* rho0_edge_old,
// 			   const amrex::Real* rho0_new, const amrex::Real* rho0_edge_new,
// 			   const amrex::Real* rhoh0_old, const amrex::Real* rhoh0_edge_old,
// 			   const amrex::Real* rhoh0_new, const amrex::Real* rhoh0_edge_new,
// 			   const amrex::Real* w0);
//
// #if (AMREX_SPACEDIM == 3)
//     void make_rhoh_flux_3d_sphr(const int* lo, const int* hi,
// 				amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
// 				amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
// 				amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
// 				amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
// 				amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
// 				amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
// 				const amrex::Real* umac,   const int* u_lo, const int* u_hi,
// 				const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
// 				const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
// 				const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
// 				const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
// 				const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
// 				const amrex::Real* rho0_edgex, const int* rx_lo, const int* rx_hi,
// 				const amrex::Real* rho0_edgey, const int* ry_lo, const int* ry_hi,
// 				const amrex::Real* rho0_edgez, const int* rz_lo, const int* rz_hi,
// 				const amrex::Real* h0_edgex,   const int* hx_lo, const int* hx_hi,
// 				const amrex::Real* h0_edgey,   const int* hy_lo, const int* hy_hi,
// 				const amrex::Real* h0_edgez,   const int* hz_lo, const int* hz_hi);
//
//     void make_rhoh_flux_3d_sphr_irreg(const int* lo, const int* hi,
// 				amrex::Real* sfluxx, const int* fx_lo, const int* fx_hi, const int* nc_fx,
// 				amrex::Real* sfluxy, const int* fy_lo, const int* fy_hi, const int* nc_fy,
// 				amrex::Real* sfluxz, const int* fz_lo, const int* fz_hi, const int* nc_fz,
// 				amrex::Real* sedgex, const int* x_lo, const int* x_hi, const int* nc_x,
// 				amrex::Real* sedgey, const int* y_lo, const int* y_hi, const int* nc_y,
// 				amrex::Real* sedgez, const int* z_lo, const int* z_hi, const int* nc_z,
// 				const amrex::Real* umac,   const int* u_lo, const int* u_hi,
// 				const amrex::Real* vmac,   const int* v_lo, const int* v_hi,
// 				const amrex::Real* wmac,   const int* w_lo, const int* w_hi,
// 				const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
// 				const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
// 				const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
// 				const amrex::Real* rhoh0_edgex,   const int* hx_lo, const int* hx_hi,
// 				const amrex::Real* rhoh0_edgey,   const int* hy_lo, const int* hy_hi,
// 				const amrex::Real* rhoh0_edgez,   const int* hz_lo, const int* hz_hi);
// #endif
    //////////////////////

    //////////////////////
    // in mkutrans.F90
#if (AMREX_SPACEDIM == 1)
    void mkutrans_1d(const int* lev, const int* domlo, const int* domhi,
                      const int* lo, const int* hi,
                      const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int* nc_ut, int ng_ut,
                      const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int* nc_uf, int ng_uf,
                      amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                      const amrex::Real* w0, const amrex::Real* dx, const amrex::Real* dt,
                      const int* adv_bc, const int* phys_bc);

#elif (AMREX_SPACEDIM == 2)
  void mkutrans_2d(const int* lo, const int* hi,
                    const int lev, const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf,
                    amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    const amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                    const amrex::Real* Im, const int* im_lo, const int* im_hi,
                    const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt,
                    const int* adv_bc, const int* phys_bc);

#elif (AMREX_SPACEDIM == 3)
  void mkutrans_3d(const int* lo, const int* hi,
                    const int lev, const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf,
                    amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
                    const amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                    const amrex::Real* Im, const int* im_lo, const int* im_hi,
    		        const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
    		        const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
    		        const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
                    const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt,
                    const int* adv_bc, const int* phys_bc);
#endif

    //////////////////////

    //////////////////////
    // in ppm.F90
    void ppm_2d(const int* lo, const int* hi,
        const amrex::Real* s,  const int* s_lo, const int* s_hi,
        const int nc_s,
        const amrex::Real* u,  const int* u_lo, const int* u_hi,
        const amrex::Real* v,  const int* v_lo, const int* v_hi,
        amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
        amrex::Real* Im, const int* im_lo, const int* im_hi,
        const int* domlo, const int* domhi,
        const int* adv_bc,
        const amrex::Real* dx, const amrex::Real dt,
        const int is_umac, const int comp, const int bccomp);

    void ppm_3d(const int* lo, const int* hi,
        const amrex::Real* s,  const int* s_lo, const int* s_hi,
        const int nc_s,
        const amrex::Real* u,  const int* u_lo, const int* u_hi,
        const amrex::Real* v,  const int* v_lo, const int* v_hi,
        const amrex::Real* w,  const int* w_lo, const int* w_hi,
        amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
        amrex::Real* Im, const int* im_lo, const int* im_hi,
        const int* domlo, const int* domhi,
        const int* adv_bc,
        const amrex::Real* dx, const amrex::Real dt,
        const int is_umac, const int comp, const int bccomp);

    //////////////////////

    //////////////////////
    // in slope.F90

    void slopex_2d(const int* lo, const int* hi,
                   const amrex::Real* state, const int* s_lo, const int* s_hi, const int nc_s,
                   const amrex::Real* slopex, const int* sl_lo, const int* sl_hi,
                   const int nc_sl,
                   const int* domlo, const int* domhi,
                   const int nvar, const int* adv_bc, const int nbccomp, const int bccomp);

   void slopey_2d(const int* lo, const int* hi,
                  const amrex::Real* state, const int* s_lo, const int* s_hi, const int nc_s,
                  const amrex::Real* slopey, const int* sl_lo, const int* sl_hi,
                  const int nc_sl,
                  const int* domlo, const int* domhi,
                  const int nvar, const int* adv_bc, const int nbccomp, const int bccomp);

    void slopez_3d(const int* lo, const int* hi,
                  const amrex::Real* state, const int* s_lo, const int* s_hi, const int nc_s,
                  const amrex::Real* slopez, const int* sl_lo, const int* sl_hi,
                  const int nc_sl,
                  const int* domlo, const int* domhi,
                  const int nvar, const int* adv_bc, const int nbccomp, const int bccomp);

    //////////////////////

    //////////////////////
    // in velpred.F90
#if (AMREX_SPACEDIM == 1)
    void velpred_1d(const int* lo, const int* hi,
                      const int lev, const int* domlo, const int* domhi,
                      const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                      const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf, int ng_uf,
                      const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                      amrex::Real* umac, const int* mu_lo, const int* mu_hi,
                      const amrex::Real* force, const int* f_lo, const int* f_hi, const int nc_f, int ng_f,
                      const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt,
                      const int* adv_bc, const int* phys_bc);

#elif (AMREX_SPACEDIM == 2)
  void velpred_interface_2d(const int* lo, const int* hi,
                    const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf, int ng_uf,
                    const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    const amrex::Real* Imu, const int* imu_lo, const int* imu_hi,
                    const amrex::Real* Ipu, const int* ipu_lo, const int* ipu_hi,
                    const amrex::Real* Imv, const int* imv_lo, const int* imv_hi,
                    const amrex::Real* Ipv, const int* ipv_lo, const int* ipv_hi,
                    amrex::Real* ul, const int* ul_lo, const int* ul_hi,
                    amrex::Real* ur, const int* ur_lo, const int* ur_hi,
                    amrex::Real* uimh, const int* ui_lo, const int* ui_hi,
                    const amrex::Real* dx, const amrex::Real dt,
                    const int* adv_bc, const int* phys_bc);

    void velpred_2d(const int* lo, const int* hi,
                      const int lev, const int idir,
                      const int* domlo, const int* domhi,
                      const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                      const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                      const amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                      amrex::Real* umac, const int* mu_lo, const int* mu_hi,
                      amrex::Real* vmac, const int* mv_lo, const int* mv_hi,
                      const amrex::Real* Imf, const int* imf_lo, const int* imf_hi,
                      const amrex::Real* Ipf, const int* ipf_lo, const int* ipf_hi,
                      const amrex::Real* ulx, const int* ulx_lo, const int* ulx_hi,
                      const amrex::Real* urx, const int* urx_lo, const int* urx_hi,
                      const amrex::Real* uimhx, const int* uix_lo, const int* uix_hi,
                      const amrex::Real* uly, const int* uly_lo, const int* uly_hi,
                      const amrex::Real* ury, const int* ury_lo, const int* ury_hi,
                      const amrex::Real* uimhy, const int* uiy_lo, const int* uiy_hi,
                      const amrex::Real* force, const int* f_lo, const int* f_hi, const int nc_f, int ng_f,
                      const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt,
                      const int* adv_bc, const int* phys_bc);

#elif (AMREX_SPACEDIM == 3)
    void velpred_interface_3d(const int* lo, const int* hi,
                    const int idir, const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf, int ng_uf,
                    const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    const amrex::Real* Imu, const int* imu_lo, const int* imu_hi,
                    const amrex::Real* Ipu, const int* ipu_lo, const int* ipu_hi,
                    const amrex::Real* Imv, const int* imv_lo, const int* imv_hi,
                    const amrex::Real* Ipv, const int* ipv_lo, const int* ipv_hi,
                    const amrex::Real* Imw, const int* imw_lo, const int* imw_hi,
                    const amrex::Real* Ipw, const int* ipw_lo, const int* ipw_hi,
                    amrex::Real* ul, const int* ul_lo, const int* ul_hi,
                    amrex::Real* ur, const int* ur_lo, const int* ur_hi,
                    amrex::Real* uimh, const int* ui_lo, const int* ui_hi,
                    const amrex::Real* dx, const amrex::Real dt,
                    const int* phys_bc);

    void velpred_transverse_3d(const int* lo, const int* hi,
                    const int base_dir, const int norm_dir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                    const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    const amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    const amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
                    const amrex::Real* ulx, const int* ulx_lo, const int* ulx_hi,
                    const amrex::Real* urx, const int* urx_lo, const int* urx_hi,
                    const amrex::Real* uimhx, const int* uix_lo, const int* uix_hi,
                    const amrex::Real* uly, const int* uly_lo, const int* uly_hi,
                    const amrex::Real* ury, const int* ury_lo, const int* ury_hi,
                    const amrex::Real* uimhy, const int* uiy_lo, const int* uiy_hi,
                    const amrex::Real* ulz, const int* ulz_lo, const int* ulz_hi,
                    const amrex::Real* urz, const int* urz_lo, const int* urz_hi,
                    const amrex::Real* uimhz, const int* uiz_lo, const int* uiz_hi,
                    const amrex::Real* uimh_trans, const int* uit_lo, const int* uit_hi,
		            const amrex::Real* dx, const amrex::Real dt,
                    const int* phys_bc);

    void velpred_3d(const int* lo, const int* hi,
                    const int lev, const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, int ng_ut,
                    const amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    const amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    const amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
                    amrex::Real* umac, const int* mu_lo, const int* mu_hi,
                    amrex::Real* vmac, const int* mv_lo, const int* mv_hi,
                    amrex::Real* wmac, const int* mw_lo, const int* mw_hi,
		      const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
		      const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
		      const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
              const amrex::Real* Imf, const int* imf_lo, const int* imf_hi,
              const amrex::Real* Ipf, const int* ipf_lo, const int* ipf_hi,
              const amrex::Real* ul, const int* ul_lo, const int* ul_hi,
              const amrex::Real* ur, const int* ur_lo, const int* ur_hi,
              const amrex::Real* uimhyz, const int* uyz_lo, const int* uyz_hi,
              const amrex::Real* uimhzy, const int* uzy_lo, const int* uzy_hi,
              const amrex::Real* vimhxz, const int* vxz_lo, const int* vxz_hi,
              const amrex::Real* vimhzx, const int* vzx_lo, const int* vzx_hi,
              const amrex::Real* wimhxy, const int* wxy_lo, const int* wxy_hi,
              const amrex::Real* wimhyx, const int* wyx_lo, const int* wyx_hi,
                    const amrex::Real* force, const int* f_lo, const int* f_hi, const int nc_f, int ng_f,
                    const amrex::Real* w0, const amrex::Real* dx, const amrex::Real dt,
                    const int* phys_bc);
#endif

    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
