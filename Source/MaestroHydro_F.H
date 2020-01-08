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
#if (AMREX_SPACEDIM == 2)
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
    // in mkutrans.F90
#if (AMREX_SPACEDIM == 2)
  void mkutrans_2d(const int* lo, const int* hi,
                    const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, const int ng_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf, int const ng_uf,
                    amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                    amrex::Real* Im, const int* im_lo, const int* im_hi,
                    const amrex::Real* w0_cart, const int* w_lo, const int* w_hi, 
                    const amrex::Real* dx, const amrex::Real dt,
                    const int* adv_bc, const int* phys_bc);

#elif (AMREX_SPACEDIM == 3)
  void mkutrans_3d(const int* lo, const int* hi,
                    const int idir,
                    const int* domlo, const int* domhi,
                    const amrex::Real* utilde,  const int* ut_lo, const int* ut_hi, const int nc_ut, const int ng_ut,
                    const amrex::Real* ufull,   const int* uf_lo, const int* uf_hi, const int nc_uf, int const ng_uf,
                    amrex::Real* utrans, const int* uu_lo, const int* uu_hi,
                    amrex::Real* vtrans, const int* uv_lo, const int* uv_hi,
                    amrex::Real* wtrans, const int* uw_lo, const int* uw_hi,
                    amrex::Real* Ip, const int* ip_lo, const int* ip_hi,
                    amrex::Real* Im, const int* im_lo, const int* im_hi,
    		        const amrex::Real* w0macx, const int* wx_lo, const int* wx_hi,
    		        const amrex::Real* w0macy, const int* wy_lo, const int* wy_hi,
    		        const amrex::Real* w0macz, const int* wz_lo, const int* wz_hi,
                    const amrex::Real* w0_cart, const int* w_lo, const int* w_hi, 
                    const amrex::Real* dx, const amrex::Real dt,
                    const int* adv_bc, const int* phys_bc);
#endif

    //////////////////////

    //////////////////////
    // in ppm.F90
#if (AMREX_SPACEDIM == 2)
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
        const int is_umac, const int comp, const int bccomp, 
        const int nbccomp);
#elif (AMREX_SPACEDIM == 3)
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
        const int is_umac, const int comp, const int bccomp, 
        const int nbccomp);
#endif

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
#if (AMREX_SPACEDIM == 2)
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
                    const amrex::Real* w0_cart, const int* w_lo, const int* w_hi,
                    const amrex::Real* dx, const amrex::Real dt,
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
                    const amrex::Real* w0_cart, const int* w_lo, const int* w_hi, 
                    const amrex::Real* dx, const amrex::Real dt,
                    const int* phys_bc);
#endif

    void get_rel_eps (const amrex::Real* rel_eps_in);

    //////////////////////

#ifdef __cplusplus
}
#endif

#endif
