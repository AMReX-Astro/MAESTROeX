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

    void get_rel_eps (const amrex::Real* rel_eps_in);

#ifdef __cplusplus
}
#endif

#endif
