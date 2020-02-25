#ifndef _eos_F_H_
#define _eos_F_H
#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void eos_get_small_temp(amrex::Real* small_temp_out);

  void eos_get_small_dens(amrex::Real* small_dens_out);

#ifdef __cplusplus
}
#endif

#endif
